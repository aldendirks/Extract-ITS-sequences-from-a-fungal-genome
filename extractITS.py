#!/usr/bin/env python

import argparse
import os
import subprocess
import sys

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


ALLOWED_REGIONS = ["SSU", "ITS1", "5.8S", "ITS2", "LSU", "all", "none"]


def get_params(argv):
    parser = argparse.ArgumentParser(description="Utility to extract rDNA sequences from fungal genomes")
    parser.add_argument("-i", "--input", help="input genome file", required=True)
    parser.add_argument("-o", "--output", help="output directory", required=True)
    parser.add_argument("--which", help="specify which rDNA sequences to extract (SSU|ITS1|5.8S|ITS2|LSU|all|none); default=all", default="all", required=False)
    parser.add_argument("-t", "--threads", help=" Number of threads/cores to use", default="1", required=False,)
    parser.add_argument("-p", "--prefix", help="prefix for output filenames", required=False)
    return parser.parse_args(argv)


def terminal(cmd):
    process = subprocess.run(cmd, capture_output=True, text=True)
    return process


def touch_file(file_path):
    with open(file_path, "a"):
        os.utime(file_path, None)


if __name__ == "__main__":
    args = get_params(sys.argv[1:])
    input_filename = os.path.basename(args.input)

    args.output = os.path.abspath(args.output)

    os.makedirs(args.output, exist_ok=True)

    if not args.prefix:
        output_prefix = os.path.splitext(input_filename)[0]
    else:
        output_prefix = args.prefix

    if args.which not in ALLOWED_REGIONS:
        print(
            "ERROR, argument --which not recognized. Please choose between "
            "'SSU', 'ITS1', '5.8S', 'ITS2', 'LSU', 'all', or 'none'"
        )
        sys.exit(1)

    barrnap_tmp = os.path.join(args.output, "barrnap.tmp")
    with open(barrnap_tmp, "w") as barrnap_output:
        barrnap_result = subprocess.run(
            ["barrnap", "--kingdom", "euk", "--threads", args.threads, args.input],
            stdout=barrnap_output,
            stderr=subprocess.PIPE,
            text=True,
        )
    if barrnap_result.returncode != 0:
        print("ERROR running barrnap:")
        print(barrnap_result.stderr.strip())
        sys.exit(barrnap_result.returncode)

    gff = pd.read_csv(barrnap_tmp, header=None, comment="#", sep="\t")
    gff = (
        gff[gff[8].str.contains("18S") | gff[8].str.contains("28S")]
        .sort_values(by=[0, 6, 3, 4])
        .reset_index(drop=True)
    )
    regions = []
    print(gff)

    for i in gff.index:
        try:
            if "partial" in gff.loc[i, 8]:
                pass
            elif gff.loc[i, 6] == "+":
                if (
                    "18S_rRNA" in gff.loc[i, 8]
                    and "28S_rRNA" in gff.loc[i + 1, 8]
                    and gff.loc[i, 0] == gff.loc[i + 1, 0]
                    and gff.loc[i, 6] == gff.loc[i + 1, 6]
                ):
                    regions.append(
                        (
                            gff.loc[i, 0],
                            gff.loc[i, 6],
                            gff.loc[i, 3] - 1,
                            gff.loc[i + 1, 4] - 1,
                        )
                    )
            elif gff.loc[i, 6] == "-":
                if (
                    "28S_rRNA" in gff.loc[i, 8]
                    and "18S_rRNA" in gff.loc[i + 1, 8]
                    and gff.loc[i, 0] == gff.loc[i + 1, 0]
                    and gff.loc[i, 6] == gff.loc[i + 1, 6]
                ):
                    regions.append(
                        (
                            gff.loc[i, 0],
                            gff.loc[i, 6],
                            gff.loc[i, 3] - 1,
                            gff.loc[i + 1, 4] - 1,
                        )
                    )
            else:
                print(f"Error, impossible to read {gff.loc[i, 6]}")
        except (KeyError, IndexError):
            pass

    print(regions)

    with open(args.input, "r") as handle:
        fasta = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    extracted_regions = []
    for region in regions:
        print(region)
        seq = fasta[region[0]][region[2] : region[3]]
        seq.id = f"{region[0]}_{region[1]}_{region[2]}_{region[3]}"
        seq.name = seq.id
        seq.description = seq.id
        extracted_regions.append(seq)
        print(seq)

    regions_fasta_path = os.path.join(args.output, "regions.fasta")
    with open(regions_fasta_path, "w") as output_handle:
        SeqIO.write(extracted_regions, output_handle, "fasta")

    if extracted_regions:
        itsx_result = terminal(
            [
                "ITSx",
                "-t",
                "F",
                "-i",
                regions_fasta_path,
                "-o",
                os.path.join(args.output, output_prefix),
                "--cpu",
                args.threads,
                "--save_regions",
                args.which,
            ]
        )
        if itsx_result.returncode != 0:
            print("ERROR running ITSx:")
            print(itsx_result.stderr.strip())
            sys.exit(itsx_result.returncode)

        ribcistron = ["SSU", "ITS1", "5_8S", "ITS2", "LSU"]
        for cist in ribcistron:
            itsx_output = os.path.join(args.output, f"{output_prefix}.{cist}.fasta")
            if not os.path.isfile(itsx_output):
                continue

            print("removing duplicates from", cist)

            with open(itsx_output, "r") as handle:
                its1s = list(SeqIO.parse(handle, "fasta"))

            seqs = [str(record.seq) for record in its1s]
            diff = list(set(seqs))
            unique = []

            for n in range(len(diff)):
                if len(diff) == 1:
                    record_name = output_prefix
                else:
                    record_name = output_prefix + "_#" + str(n + 1)

                unique.append(
                    SeqRecord(
                        Seq(diff[n]),
                        record_name,
                        record_name,
                        str(seqs.count(diff[n]))
                        + " copies of this sequence found in "
                        + input_filename,
                    )
                )

            filtered_output = os.path.join(
                args.output, f"{output_prefix}.{cist}_filtered.fasta"
            )
            with open(filtered_output, "w") as output_handle:
                SeqIO.write(unique, output_handle, "fasta")

            if len(diff) == 1:
                copy_count = seqs.count(diff[0])
                print(
                    "\nDONE. A single "
                    + cist
                    + " sequence was found in "
                    + str(copy_count)
                    + " copies, in file "
                    + input_filename
                    + "."
                )
            elif len(diff) > 1:
                print(
                    "\nDONE. "
                    + str(len(diff))
                    + " different "
                    + cist
                    + " sequences were found in file "
                    + input_filename
                    + ". See output files for more info."
                )
            elif len(diff) == 0:
                print(
                    "\nDONE. No "
                    + cist
                    + " sequence found in file "
                    + input_filename
                )
    else:
        print(
            "\nDONE. No " + args.which + " sequence found in file " + input_filename + "."
        )
        touch_file(
            os.path.join(args.output, f"{output_prefix}.{args.which}_filtered.fasta")
        )

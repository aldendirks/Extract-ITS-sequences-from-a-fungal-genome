#!/usr/bin/env python

import atexit
import argparse
import os
import shlex
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
    parser.add_argument("-t", "--threads", help="number of threads/cores to use", default="8", required=False,)
    parser.add_argument("-p", "--prefix", help="prefix for output filenames", required=False)
    return parser.parse_args(argv)


def touch_file(file_path):
    with open(file_path, "a"):
        os.utime(file_path, None)


class TeeStream:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
        return len(data)

    def flush(self):
        for stream in self.streams:
            stream.flush()

    def isatty(self):
        return any(getattr(stream, "isatty", lambda: False)() for stream in self.streams)


if __name__ == "__main__":
    args = get_params(sys.argv[1:])
    input_filename = os.path.basename(args.input)
    
    args.output = os.path.abspath(args.output)

    os.makedirs(args.output, exist_ok=True)

    if not args.prefix:
        output_prefix = os.path.splitext(input_filename)[0]
    else:
        output_prefix = args.prefix

    log_path = os.path.join(args.output, f"{output_prefix}.log")
    log_handle = open(log_path, "w", encoding="utf-8")
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = TeeStream(original_stdout, log_handle)
    sys.stderr = TeeStream(original_stderr, log_handle)

    def close_log_file():
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        log_handle.flush()
        log_handle.close()

    atexit.register(close_log_file)
    print(f"Logging output to: {log_path}")

    if args.which not in ALLOWED_REGIONS:
        print(
            "\nERROR, argument --which not recognized. Please choose between "
            "'SSU', 'ITS1', '5.8S', 'ITS2', 'LSU', 'all', or 'none'."
        )
        sys.exit(1)

    print("\nRun summary:")
    print(f"  Input genome: {os.path.abspath(args.input)}")
    print(f"  Output directory: {args.output}")
    print(f"  Output prefix: {output_prefix}")
    print(f"  Regions requested: {args.which}")
    print(f"  Threads: {args.threads}")

    barrnap_tmp = os.path.join(args.output, "barrnap.tmp")
    barrnap_cmd = ["barrnap", "--kingdom", "euk", "--threads", args.threads, args.input]
    print("\nRunning barrnap to identify rDNA regions in the genome assembly...")
    print(f"\t{shlex.join([str(part) for part in barrnap_cmd])}")
    with open(barrnap_tmp, "w") as barrnap_output:
        barrnap_result = subprocess.run(
            barrnap_cmd,
            stdout=barrnap_output,
            stderr=subprocess.PIPE,
            text=True,
        )
    if barrnap_result.returncode != 0:
        print("\nERROR running barrnap:")
        print(barrnap_result.stderr.strip())
        sys.exit(barrnap_result.returncode)

    gff = pd.read_csv(barrnap_tmp, header=None, comment="#", sep="\t")
    print("\nBarrnap results:")
    print(gff)
    gff = (
        gff[gff[8].str.contains("18S") | gff[8].str.contains("28S")]
        .sort_values(by=[0, 6, 3, 4])
        .reset_index(drop=True)
    )
    regions = []
    print("\nBarrnap results filtered and sorted:")
    print(gff)

    print("\nIdentifying rDNA regions based on barrnap results...")
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

    with open(args.input, "r") as handle:
        fasta = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    extracted_regions = []
    for region in regions:
        seq = fasta[region[0]][region[2] : region[3]]
        seq.id = f"{region[0]}_{region[1]}_{region[2]}_{region[3]}"
        seq.name = seq.id
        seq.description = seq.id
        extracted_regions.append(seq)
        print(f"\t{seq.description}")

    regions_fasta_path = os.path.join(args.output, "regions.fasta")
    with open(regions_fasta_path, "w") as output_handle:
        SeqIO.write(extracted_regions, output_handle, "fasta")

    print("\nExtracting rDNA regions from the genome assembly with ITSx...")
    if extracted_regions:
        itsx_cmd = [
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
        print(f"\t{shlex.join([str(part) for part in itsx_cmd])}")
        itsx_result = subprocess.run(itsx_cmd, capture_output=True, text=True)
        if itsx_result.returncode != 0:
            print("\nERROR running ITSx:")
            print(itsx_result.stderr.strip())
            sys.exit(itsx_result.returncode)

        ribcistron = ["SSU", "ITS1", "5_8S", "ITS2", "LSU"]
        for cist in ribcistron:
            itsx_output = os.path.join(args.output, f"{output_prefix}.{cist}.fasta")
            if not os.path.isfile(itsx_output):
                continue

            print(f"\nRemoving duplicates from {cist}")

            with open(itsx_output, "r") as handle:
                its1s = list(SeqIO.parse(handle, "fasta"))

            seqs = [str(record.seq) for record in its1s]
            diff = list(set(seqs))
            unique = []

            for n in range(len(diff)):
                if len(diff) == 1:
                    record_name = output_prefix
                else:
                    record_name = f"{output_prefix}_#{n + 1}"

                unique.append(
                    SeqRecord(
                        Seq(diff[n]),
                        record_name,
                        record_name,
                        f"{seqs.count(diff[n])} copies of this sequence found in {input_filename}",
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
                    f"A single {cist} sequence was found in {copy_count} "
                    f"copies, in file {input_filename}."
                )
            elif len(diff) > 1:
                print(
                    f"{len(diff)} different {cist} sequences were found in "
                    f"file {input_filename}. See output files for more info."
                )
            elif len(diff) == 0:
                print(
                    f"No {cist} sequence found in file {input_filename}"
                )
    else:
        print(
            f"No {args.which} sequence found in file {input_filename}."
        )
        touch_file(
            os.path.join(args.output, f"{output_prefix}.{args.which}_filtered.fasta")
        )
    
    print("\nDone!")

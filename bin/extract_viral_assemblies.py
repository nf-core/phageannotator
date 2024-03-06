#!/usr/bin/env python

import argparse
import pandas as pd
from Bio import SeqIO
import sys
import gzip


def parse_args(args=None):
    Description = "Extract genomes identified as contained via mash screen, and append to the sample's FASTA assembly."
    Epilog = "Example usage: python append_screen_hits.py <FILES_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-v",
        "--virus_fasta",
        help="Path to virus FASTA file (gzipped) containing virus sequences.",
    )
    parser.add_argument(
        "-e",
        "--exclude",
        help="Exclude entries containing this substring from virus assemblies file",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output TSV file containing the names of all viral assemblies.",
    )
    return parser.parse_args(args)


def extract_viral_assemblies(virus_fasta, exclude, output):

    viral_assemblies = []
    with gzip.open(virus_fasta, "rt") as reference_fasta_gunzip:
        for record in SeqIO.parse(reference_fasta_gunzip, "fasta"):
            if exclude in record.id:
                continue
            else:
                viral_assemblies.append(record.id)

    out = open(output, "w")
    for sequence in viral_assemblies:
        out.write(sequence + "\n")


def main(args=None):
    args = parse_args(args)
    extract_viral_assemblies(args.virus_fasta, args.exclude, args.output)


if __name__ == "__main__":
    sys.exit(main())

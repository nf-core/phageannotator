#!/usr/bin/env python

import sys
import argparse
import pandas as pd
from Bio import SeqIO

def parse_args(args=None):
    Description = "Create a scaffold-to-bin file for input to inStrain."
    Epilog = "Example usage: python create_instrain_stb.py -f fasta.fasta -o fasta.stb"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-f",
        "--fasta",
        help="Path to FASTA file that will be input into inStrain.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to the TSV file, that will be scaffold-to-bin file for inStrain",
    )
    return parser.parse_args(args)


def create_instrain_stb(fasta, output):
    # create list to store contig names in
    contig_names = []

    # read in the fasta
    for record in SeqIO.parse(fasta, "fasta"):
        contig_names.append(record.id)

    stb_df = pd.DataFrame()
    stb_df['scaffold'] = contig_names
    stb_df['bin'] = stb_df['scaffold']
    stb_df.to_csv(output, sep='\t', index=False, header=False)


def main(args=None):
    args = parse_args(args)
    create_instrain_stb(args.fasta, args.output)


if __name__ == "__main__":
    sys.exit(main())

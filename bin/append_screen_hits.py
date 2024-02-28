#!/usr/bin/env python

import argparse
from Bio import SeqIO
import sys
import gzip


def parse_args(args=None):
    Description = "Extract viral assembly names from an virus fasta input."
    Epilog = "Example usage: python extract_viral_assemblies.py <FILES_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-r",
        "--reference_fasta",
        help="Path to FASTA file (gzipped) that was sketched with mash sketch, and searched for containment with mash screen.",
    )
    parser.add_argument(
        "-s",
        "--mash_screen_results",
        help="Path to the TSV file output by running mash screen.",
    )
    parser.add_argument(
        "-f",
        "--assembly_fasta",
        help="Path to the FASTA file (gzipped) containing assemblies from the screened sample.",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="The sample ID for identifying sample of origin when combining across assemblies.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output FASTA file containing assemblies and appended mash screen hits.",
    )
    return parser.parse_args(args)


def append_screen_hits(reference_fasta, mash_screen_results, assembly_fasta, prefix, output):
    mash_screen_results_df = pd.read_csv(
        mash_screen_results,
        sep="\t",
        header=None,
        index_col=False,
        names=["identity", "shared-hashes", "median-multiplicity", "p-value", "query-id", "query-comment"],
    )
    reference_hits = set(mash_screen_results_df["query-id"])
    contained_genomes = []
    tested_genomes = set()
    with gzip.open(reference_fasta, "rt") as reference_fasta_gunzip:
        for record in SeqIO.parse(reference_fasta_gunzip, "fasta"):
            if record.id in reference_hits:
                if record.id in tested_genomes:
                    continue
                else:
                    record.id = "mash_screen|" + record.id
                    contained_genomes.append(record)
                    tested_genomes.add(record.id)
    with gzip.open(assembly_fasta, "rt") as assembly_fasta_gunzip:
        for record in SeqIO.parse(assembly_fasta_gunzip, "fasta"):
            record.id = prefix + "|" + record.id
            contained_genomes.append(record)
    SeqIO.write(contained_genomes, output, "fasta")


def main(args=None):
    args = parse_args(args)
    append_screen_hits(args.reference_fasta, args.mash_screen_results, args.assembly_fasta, args.prefix, args.output)


if __name__ == "__main__":
    sys.exit(main())

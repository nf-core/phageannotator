#!/usr/bin/env python

from Bio import SeqIO
import argparse
import sys
import gzip

def parse_args(args=None):
    Description = "Extract cluster representatives into a single FASTA file."
    Epilog = "Example usage: python extract_cluster_representatives.py <FILES_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-f",
        "--fasta",
        help="Path to FASTA file (gzipped) that contains phages used in clustering.",
    )
    parser.add_argument(
        "-c",
        "--clusters",
        help="Path to the TSV file containing cluster representatives and member sequences.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Path to the where cluster representative FASTA file should be output."
    )
    return parser.parse_args(args)

def extract_cluster_representatives(fasta, clusters, output):

    # open clustering results
    clusters = open(clusters, 'r')

    cluster_reps= []
    for line in clusters:
        stripped = line.strip()
        rep, nodes = stripped.split('\t')
        cluster_reps.append(rep)

    cluster_reps_set = set(cluster_reps)

    # extract representative sequences from fasta file
    cluster_rep_sequences = []
    already_added = set()

    fasta_gunzipped = gzip.open(fasta, 'rt')
    for record in SeqIO.parse(fasta_gunzipped, "fasta"):
        if record.id in already_added:
            continue
        if record.id in cluster_reps_set:
            cluster_rep_sequences.append(record)
            already_added.add(record.id)

    # save all sequences to specified file
    SeqIO.write(cluster_rep_sequences, output, "fasta")


def main(args=None):
    args = parse_args(args)
    extract_cluster_representatives(args.fasta, args.clusters, args.output)


if __name__ == "__main__":
    sys.exit(main())

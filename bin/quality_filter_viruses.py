#!/usr/bin/env python

import argparse
import pandas as pd
from Bio import SeqIO
import sys
import gzip
import os

def parse_args(args=None):
    Description = "Filter virus sequence based on CheckV output"
    Epilog = "Example usage: python quality_filter_viruses.py -v viruses.fna.gz -p proviruses.fna.gz -q quality_summary.tsv -o filtered_viruses.fna"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-v",
        "--viruses",
        help="Path to viruses FASTA (gzipped) file output by CheckV.",
    )
    parser.add_argument(
        "-p",
        "--proviruses",
        help="Path to proviruses FASTA (gzipped) file output by CheckV.",
    )
    parser.add_argument(
        "-q",
        "--quality_summary",
        help="Path to quality summary TSV file output by CheckV.",
    )
    parser.add_argument(
        "-l",
        "--min_length",
        help="Minimum length of viral sequence to pass filtering (unless min_completeness threshold is met).",
    )
    parser.add_argument(
        "-c",
        "--min_completeness",
        help="Minimum completeness of viral sequence to pass filtering (unless min_length threshold is met).",
    )
    parser.add_argument(
        "-r",
        "--remove_proviruses",
        help="Remove sequences that have been labeled as proviruses.",
        action=argparse.BooleanOptionalAction,
        default=False
    )
    parser.add_argument(
        "-w",
        "--remove_warnings",
        help="Remove sequences with k-mer or length warnings from CheckV.",
        action=argparse.BooleanOptionalAction,
        default=False
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output FASTA file containing filtered viral sequences.",
    )
    return parser.parse_args(args)

def quality_filter_viruses( viruses, proviruses, quality_summary, min_length, min_completeness, remove_proviruses, remove_warnings, output):
    if os.stat(quality_summary).st_size != 0:
        quality_summary_df = pd.read_csv(quality_summary, sep="\t")

        quality_summary_df['length'] = quality_summary_df.apply(lambda x: x.proviral_length if x.proviral_length > 0 else x.contig_length, axis=1)

        quality_summary_filtered = quality_summary_df[(quality_summary_df["length"] >= float(min_length))|(quality_summary_df["completeness"] >= float(min_completeness))]

        if remove_proviruses:
            quality_summary_filtered = quality_summary_filtered[(quality_summary_filtered["provirus"] != 'Yes')]

        # remove genomes > 1.5x longer than expected
        if len(quality_summary_filtered[quality_summary_filtered['warnings'].notnull()]) > 1 and remove_warnings:
            quality_summary_filtered = quality_summary_filtered[(quality_summary_filtered['warnings'].str.contains('contig >1.5x longer than expected genome length') == False) | (quality_summary_filtered['warnings'].isnull())]

        filtered_viruses = set(quality_summary_filtered["contig_id"])
        filtered_virus_seqs = []

        # parse through and combine virus sequences for each sample
        with gzip.open(viruses, "rt") as viruses_gunzip:
            for record in SeqIO.parse(viruses_gunzip, "fasta"):
                if record.id in filtered_viruses:
                    filtered_virus_seqs.append(record)

        with gzip.open(proviruses, "rt") as proviruses_gunzip:
            for record in SeqIO.parse(proviruses_gunzip, "fasta"):
                if record.id.rpartition('_')[0] in filtered_viruses:
                    if '|provirus' in record.id:
                        genomad_provirus = record.id.split('|provirus')[1]
                        genomad_start = genomad_provirus.split('_')[1]
                        checkv_provirus = record.description.split(' ')[1]
                        checkv_provirus_coords = checkv_provirus.split('/')[0]
                        checkv_start, checkv_stop = checkv_provirus_coords.split('-')
                        checkv_start_total = int(checkv_start) + int(genomad_start) -1
                        checkv_stop_total = int(checkv_stop) + int(genomad_start) - 1
                        record.id = record.id + "|checkv_provirus_" + str(checkv_start_total) + "_" + str(checkv_stop_total)
                    else:
                        checkv_provirus = record.description.split(' ')[1]
                        checkv_provirus_coords = checkv_provirus.split('/')[0]
                        checkv_start, checkv_stop = checkv_provirus_coords.split('-')
                        checkv_start_total = int(checkv_start)
                        checkv_stop_total = int(checkv_stop)
                        record.id = record.id + "|checkv_provirus_" + str(checkv_start_total) + "_" + str(checkv_stop_total)
                    filtered_virus_seqs.append(record)

        # save all sequences to specified file
        SeqIO.write(filtered_virus_seqs, output, "fasta")

    else:
        output = open(output, 'x')
        output.close()

def main(args=None):
    args = parse_args(args)
    quality_filter_viruses(args.viruses, args.proviruses, args.quality_summary, args.min_length, args.min_completeness, args.remove_proviruses, args.remove_warnings, args.output)


if __name__ == "__main__":
    sys.exit(main())

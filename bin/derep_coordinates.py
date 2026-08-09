#!/usr/bin/env python

from Bio import SeqIO
import argparse
import sys
import gzip
import pandas as pd
import math


def parse_args(args=None):
    Description = "Dereplicate provirus coordinates to match dereplicated provirus assemblies"
    Epilog = "Example usage: python derep_coordinates.py --coordinates <COORDS_TSV> --clusters <CLUSTERS_TSV> --output <DEREP_COORDS_TSV>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-r",
        "--coordinates",
        help="Path to TSV file containing provirus coordinates for each assembly.",
    )
    parser.add_argument(
        "-c",
        "--clusters",
        help="Path to the TSV file containing cluster information from dereplicating assemblies.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output TSV file containing provirus coordinates for dereplicated assembly.",
    )
    return parser.parse_args(args)

def derep_coordinates(coords_tsv, clusters_tsv, output):
    # open coordinates file
    coords = pd.read_csv(coords_tsv, sep='\t')
    # open cluster results
    clusters = pd.read_csv(clusters_tsv, sep='\t', header=None)
    # identify coords contained in derep clusters
    derep_coords = coords[coords['scaffold'].isin(set(clusters[0]))]

    # save coords file
    derep_coords.to_csv(output, sep='\t', index=False)

def main(args=None):
    args = parse_args(args)
    derep_coordinates(args.coordinates, args.clusters, args.output )


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python

from Bio import SeqIO
import argparse
import sys
import gzip
import pandas as pd
import math


def parse_args(args=None):
    Description = "Extract proviruses present in assemblies based on geNomad and CheckV."
    Epilog = "Example usage: python extractproviruses.py --fasta <FASTA> --genomad <virus_summary> --checkv <quality_summary> --output_fasta <proviruses.fasta.gz> --output_meta <provirus_metadata.tsv>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-f",
        "--fasta",
        help="Path to FASTA file (gzipped) that contains contigs.",
    )
    parser.add_argument(
        "-g",
        "--genomad",
        help="Path to the TSV file containing geNomad's virus summary.",
    )
    parser.add_argument(
        "-c",
        "--checkv",
        help="Path to the TSV file containing CheckV's contamination summary.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output FASTA file containing provirus scaffolds.",
    )
    parser.add_argument(
        "-t",
        "--tsv",
        help="Output TSV file containing provirus coordinates.",
    )
    return parser.parse_args(args)

def extract_proviruses(fasta, genomad, checkv, out_fasta, out_tsv):
    # open genomad results
    genomad = pd.read_csv(genomad, sep='\t')

    # identify genomad proviruses and their coordinates
    genomad_proviruses = genomad[genomad['topology'] == 'Provirus']

    # open checkv results
    checkv = pd.read_csv(checkv, sep='\t')

    # filter to checkv proviruses
    checkv_proviruses = checkv[checkv['provirus'] == 'Yes']

    # identify checkv provirus coordinates
    checkv_coords = pd.DataFrame()
    for index, row in checkv_proviruses.iterrows():
        contig_id = row['contig_id']
        region_types = row['region_types'].split(',')
        provirus_count = 0
        # parse though regions for each contig
        for i in range(len(region_types)):
            if region_types[i] == 'viral':
                # if a region is viral, extract contig name, assign a provirus id, and add checkv start/end coords
                provirus_info = pd.DataFrame()
                provirus_count += 1
                provirus_info['seq_name'] = [contig_id]
                provirus_info['provirus_id'] = [contig_id + '|checkv_provirus_' + str(provirus_count)]
                provirus_info[['provirus_start', 'provirus_stop']] = [row['region_coords_bp'].split(',')[i].split('-')]
                checkv_coords = pd.concat([checkv_coords, provirus_info], axis=0)

    # merge genomad and checkv
    genomad_checkv = genomad_proviruses.merge(checkv_coords, on='seq_name', how='outer')

    # if checkv provirus exists, add it to df and skip associated genomad provirus. Else, use geNomad provirus
    provirus_combined_coords = pd.DataFrame()
    already_added = {}
    for index, row in genomad_checkv.iterrows():
        if str(row['seq_name']) + '_' + str(row['provirus_id']) in already_added:
            continue
        else:
            provirus_coords = pd.DataFrame()
            # remove prefix added above or by geNomad
            provirus_coords['scaffold'] = [row['seq_name'].split('|provirus')[0]]
            # if checkv provirus only; use checkv coordinates
            if float(row['provirus_start']) > 0 and '-' not in str(row['coordinates']):
                provirus_coords['start'] = [row['provirus_start']]
                provirus_coords['stop'] = [row['provirus_stop'] ]
            # if checkv and genomad provirus; add genomad start and checkv start (-1 to set checkv start to 0)
            # find stop by adding length (checkv_stop - checkv_start) to start coord found above
            elif float(row['provirus_start']) > 0 and '-' in row['coordinates']:
                provirus_coords['start'] = [int(row['coordinates'].split('-')[0]) + (int(row['provirus_start']) - 1)]
                provirus_coords['stop'] = [provirus_coords['start'][0] + int(row['provirus_stop']) - int(row['provirus_start'])]
            # if genomad provirus only; use genomad coordinates
            elif math.isnan(row['provirus_start']) and '-' in row['coordinates']:
                provirus_coords['start'] = [int(row['coordinates'].split('-')[0])]
                provirus_coords['stop'] = [int(row['coordinates'].split('-')[1])]
            # rename provirus fragment based on final start/end
            if '|provirus_' in row['seq_name']:
                provirus_coords['fragment'] = row['seq_name'].split('|provirus_')[0] + '|provirus_' + str(provirus_coords['start'][0]) + '_' + str(provirus_coords['stop'][0])
            else:
                provirus_coords['fragment'] = row['seq_name'] + '|provirus_' + str(provirus_coords['start'][0]) + '_' + str(provirus_coords['stop'][0])
            # combine coordinates with previous row coordinates
            provirus_combined_coords = pd.concat([provirus_combined_coords, provirus_coords], axis=0)

    # reorder columns to match propagate input
    provirus_combined_coords_reorg = provirus_combined_coords[['scaffold', 'fragment', 'start', 'stop']]
    # save coords file
    provirus_combined_coords_reorg.to_csv(out_tsv, sep='\t', index=False)
    # identify scaffolds to extract
    provirus_scaffolds = set(provirus_combined_coords_reorg['scaffold'])

    # extract provirus scaffolds from fasta
    exracted_scaffolds = []
    fasta_gunzipped = gzip.open(fasta, "rt")
    for record in SeqIO.parse(fasta_gunzipped, "fasta"):
        if record.id in provirus_scaffolds:
            record.description = record.id
            exracted_scaffolds.append(record)
    # save all extracted provirus scaffolds to specified file
    SeqIO.write(exracted_scaffolds, out_fasta, "fasta")


def main(args=None):
    args = parse_args(args)
    extract_proviruses(args.fasta, args.genomad, args.checkv, args.output, args.tsv)


if __name__ == "__main__":
    sys.exit(main())

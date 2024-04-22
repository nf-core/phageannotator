#! /usr/bin/env python3
# Author: Kristopher Kieft, UW-Madison
# kieft@wisc.edu

# PropagAtE: Prophage Activity Estimator
# Version: v1.1.0
# Release date: January 2022

import subprocess
import os
import sys
import pysam
import numpy as np
from numba import jit


def prophages_check(vibe):
    with open(vibe, 'r') as phage_list:
        header = phage_list.readline().strip('\n').split('\t')
        check_line = phage_list.readline().strip('\n')
        if not check_line:
            return 'no prophages'

        try:
            if header[0] == "scaffold" and header[1] == "fragment" and header[5] == "nucleotide start" and header[6] == "nucleotide stop":
                return 'vibrant header'
        except IndexError:
            pass
        try:
            if header[0] == "scaffold" and (header[1] == "fragment" or header[1] == "prophage") and header[2] == "start" and header[3] == "stop":
                return 'custom header'
        except IndexError:
            pass

        sys.stderr.write("\nError: -v coordinates file is formatted incorrectly. See README for details. Exiting.\n")
        return ''

def does_exist(file, descript):
    if os.path.exists(file):
        sys.stderr.write(f"\nError: the {descript} already exists. Exiting.\n\n")
        return True
    return False

def not_exist(file, descript):
    if not os.path.exists(file):
        sys.stderr.write(f"\nError: the {descript} file does not exist. Exiting.\n\n")
        return True
    return False

@jit(nopython=True)
def quick_stats(depth, start, stop, min_cov, all_phages, l, host_depth):
    host_depth[start:stop] = np.nan # do not count any prophages in host coverage
    all_phages += l
    if l >= 1000: # minimum length allowed
        cov = depth[start:stop]
        a = np.nanmean(cov)
        m = np.nanmedian(cov)
        s = np.nanstd(cov)
        d = np.count_nonzero(cov >= min_cov) / l

    return a,m,s,d,all_phages,host_depth

@jit(nopython=True)
def host_stats(depth):
    eff_host = depth.size
    if eff_host > 0:
        a = np.nanmean(depth)
        m = np.nanmedian(depth)
        s = np.nanstd(depth)
        return a,m,s,eff_host
    else:
        return 0,0,0,0

def coverage_stats(depth, genome, prophage_dict, prophage_lengths, min_cov, length):
    '''
    Calculate average, median and standard deviation
    '''
    prophages = prophage_dict[genome]
    phage_covs = {}
    all_phages = 0 # length of all phages on this scaffold
    host_depth = depth.copy() # in case of prophage overlap
    for p in prophages:
        phage,start,stop = p
        l = prophage_lengths[phage]
        a,m,s,d,all_phages,host_depth = quick_stats(depth, start, stop, min_cov, all_phages, l, host_depth)
        phage_covs[phage] = (a,m,s,l,d)

    host_depth_filter = host_depth[~np.isnan(host_depth)]
    a,m,s,eff_host = host_stats(host_depth_filter) # host
    return phage_covs, a, m, s, eff_host

@jit(nopython=True)
def add_depth(depth, start, end, ed, rl, read_id):
    '''
    Add depth with read alignment filter
    '''
    if ed/rl <= read_id:
        for i in range(start,end):
            depth[i] += 1
    return depth

@jit(nopython=True)
def add_depth_no_ed(depth, start, end):
    '''
    Add depth, ignore read alignment identity
    '''
    for i in range(start,end):
        depth[i] += 1
    return depth

def fasta_parse(infile):
    '''
    Standard fasta parser.
    '''
    with open(infile, 'r') as fasta:
        for line in fasta:
            if line[0] == '>':
                try:
                    yield header, seq
                except NameError:
                    pass # first line
                seq = ''
                header = line[1:].strip("\n")
            else:
                seq += line.strip("\n")

        # last one
        yield header, seq

def bowtie2_build(fasta, base, outpath, spaces, u):
    '''
    Build bowtie2 index.
    '''
    build = False
    if not os.path.exists(fasta.rsplit('.',1)[0] + ".bowtie2.index.1.bt2"):
        if spaces:
            with open(f"{outpath}{base}.{u}.fasta", 'w') as fasta_out:
                for name,seq in fasta_parse(fasta):
                    name = name.replace(" ", "~!!~")
                    fasta_out.write(f">{name}\n{seq}\n")
            fasta = f"{outpath}{base}.{u}.fasta"

        subprocess.run(f"bowtie2-build {fasta} {outpath}{base}.bowtie2.index > /dev/null 2> /dev/null", shell=True)
        build = True

    return build

def run_bowtie2_paired(base, outpath, threads, forward, reverse):
    subprocess.run(f"bowtie2 -x {outpath}{base}.bowtie2.index -1 {forward} -2 {reverse} -S {outpath}{base}.sam -q -p {threads} --no-unal --no-discordant > /dev/null 2> /dev/null", shell=True)

def run_bowtie2_interleaved(base, outpath, threads, interleaved):
    subprocess.run(f"bowtie2 -x {outpath}{base}.bowtie2.index --interleaved {interleaved} -S {outpath}{base}.sam -q -p {threads} --no-unal --no-discordant > /dev/null 2> /dev/null", shell=True)

def run_bowtie2_unpaired(base, outpath, threads, unpaired):
    subprocess.run(f"bowtie2 -x {outpath}{base}.bowtie2.index -U {unpaired} -S {outpath}{base}.sam -q -p {threads} --no-unal > /dev/null 2> /dev/null", shell=True)

def post_bowtie2(outpath, base, clean, u, build, threads):
    '''
    Clean up bowtie2 and SAM -> BAM
    '''
    sam = f'{outpath}{base}.sam'
    bam = f'{outpath}{base}.bam'
    subprocess.run(f"samtools view -@ {threads} -h -b {sam} > {bam} 2> /dev/null", shell=True)

    if clean:
        subprocess.run(f"rm {outpath}{base}.bowtie2.index.* 2> /dev/null", shell=True)
        subprocess.run(f"rm {sam} 2> /dev/null", shell=True)
        if build:
            subprocess.run(f"rm {outpath}{base}.{u}.fasta 2> /dev/null", shell=True)

    return bam

def sam_bam(sam, outpath, threads):
    '''
    SAM -> BAM
    '''
    try:
        temp = sam.rsplit("/",1)[1]
        base = temp.rsplit(".",1)[0]
    except Exception:
        base = sam.rsplit(".",1)[0]

    bam = f'{outpath}{base}.bam'
    subprocess.run(f"samtools view -@ {threads} -h -b {sam} > {bam} 2> /dev/null", shell=True)

    return bam


def sort_bam(bam, outpath, threads, clean, clean_bam):
    '''
    Sort BAM file.
    '''
    if int(os.stat(bam).st_size) == 0:
        return False, None
    # check if anything aligned
    check_aligned = subprocess.check_output(f"samtools view -@ {threads} {bam} | head -n 1", shell=True)
    if len(check_aligned) == 0:
        return False, None

    try:
        temp = bam.rsplit("/",1)[1]
        base = temp.rsplit(".",1)[0]
    except Exception:
        base = bam.rsplit(".",1)[0]

    sort_check = False
    if clean_bam: # only need to check if bamfile was the input
        try:
            check = subprocess.check_output(f"samtools view -@ {threads} -H {bam} | grep '@HD'", shell=True)
            if "coordinate" in str(check):
                sort_check = True
        except Exception:
            # no @HD line, retain sort_check = False
            pass
    if not sort_check:
        subprocess.run(f"samtools sort -@ {threads} -o {outpath}{base}.sorted.bam {bam} 2> /dev/null", shell=True)

        if clean and not clean_bam:
            subprocess.run(f'rm {bam}', shell=True)

        bam = f'{outpath}{base}.sorted.bam'

    return True, bam

@jit(nopython=True)
def get_len(start,stop):
    '''
    Length of prophages.
    Python is zero-based and end-exclusive in range(),
    so even though start = start-1, leave this math
    as stop-start to include the start and stop base precisely
    '''
    return stop-start


@jit(nopython=True)
def reverse_and_zero_based(start, stop):
    '''
    Reverse start and stop if the prophage is in reverse coordinates.
    -1 to start for Python zero-based format.
    Leave stop for Python range exclusive format.
    '''
    if start > stop:
        start, stop = stop, start
    start -= 1
    return start,stop


def process_coordinates_file(vibe, spaces, vibe_header):
    '''
    Convert coordinates file into usable dictionary and list of names.
    '''
    prophage_dict = {}
    prophage_lengths = {}
    prophage_dict_frags = {}
    if vibe_header == 'vibrant header':
        x = 5
        y = 6
    elif vibe_header == 'custom header':
        x = 2
        y = 3

    with open(vibe, 'r') as phage_list:
        next(phage_list)
        for line in phage_list:
            line = line.strip('\n').split('\t')
            if not line: continue
            name = line[0]
            frag = line[1]

            if spaces:
                name = name.replace(" ", "~!!~")
                frag = frag.replace(" ", "~!!~")

            start = int(line[x])
            stop = int(line[y])
            start, stop = reverse_and_zero_based(start, stop)
            prophage_dict.setdefault(name, []).append((frag, start, stop))
            prophage_lengths[frag] = get_len(start,stop)
            prophage_dict_frags[frag] = name

        genomes_full = set(list(prophage_dict.keys()))
        prophages = len(prophage_dict_frags.keys())

        return True, prophage_dict, prophage_dict_frags, genomes_full, prophages, prophage_lengths

def get_lengths(fasta, spaces, genomes_full):
    '''
    Get genome lengths for initialization of np.zeros(length)
    '''
    lengths = {}
    if spaces:
        for name,seq in fasta_parse(fasta):
            name = name.replace(" ", "~!!~")
            if name in genomes_full:
                lengths[name] = len(seq)
    else:
        for name,seq in fasta_parse(fasta):
            if name in genomes_full:
                lengths[name] = len(seq)

    return lengths


def extract_coverage(bam, read_id, lengths, mask, outfile, prophage_dict, effect, ratio_cutoff, prophage_dict_frags, prophage_lengths, min_breadth, min_cov, clean, spaces):
    '''
    Code mainly from vRhyme (same author).
    Extract coverage information from BAM file.
    '''
    with open(outfile, 'w') as output:
        output.write("prophage\thost\tactive\tCohenD\tprophage-host_ratio\tmean_difference\tprophage_len\tprophage_mean_cov\tprophage_median_cov\tprophage_sd_cov\tprophage_cov_breadth\thost_len\thost_mean_cov\thost_median_cov\thost_sd_cov\n")
    written = []
    total = 0

    bai = False
    if not os.path.exists(bam + '.bai'):
        subprocess.run(f'samtools index {bam}', shell=True)
        bai = True

    if mask != 0:
        mask_f = mask
        mask_r = -mask
    else:
        mask_f = False
        mask_r = False

    bamfile = pysam.AlignmentFile(bam, "rb")

    if read_id != 0:
        for x in bamfile.fetch(until_eof=True):
            genome = x.reference_name
            length = lengths.get(genome,False)
            if length: # not in keep
                try:
                    if genome != prev:
                        if mask_f:
                            depth[:mask_f] = np.nan
                            depth[mask_r:] = np.nan
                        prophage_covs, avg, med, sd, eff_host = coverage_stats(depth, prev, prophage_dict, prophage_lengths, min_cov, length)
                        written, total = write_coverages(outfile, prophage_covs, prev, avg, med, sd, eff_host, effect, ratio_cutoff, written, total, min_breadth, min_cov, spaces)
                        depth = np.zeros(length)
                except NameError:
                    # prev not defined
                    depth = np.zeros(length)

                prev = genome

                ed = 0
                for t in x.tags:
                    if t[0] == 'NM':
                        ed = t[1]
                        break
                rl = x.query_length

                start = x.reference_start # 0-based
                end = x.reference_end
                if end:
                    depth = add_depth(depth, start, end, ed, rl, read_id)

        if length: # not in keep
            # last one
            if mask_f:
                depth[:mask_f] = np.nan
                depth[mask_r:] = np.nan
            prophage_covs, avg, med, sd, eff_host = coverage_stats(depth, prev, prophage_dict, prophage_lengths, min_cov, length)
            written, total = write_coverages(outfile, prophage_covs, prev, avg, med, sd, eff_host, effect, ratio_cutoff, written, total, min_breadth, min_cov, spaces)
            depth = None
    else: # no mismatches
        for x in bamfile.fetch(until_eof=True):
            genome = x.reference_name
            length = lengths.get(genome,False)
            if length: # not in keep
                try:
                    if genome != prev:
                        if mask_f:
                            depth[:mask_f] = np.nan
                            depth[mask_r:] = np.nan
                        prophage_covs, avg, med, sd, eff_host = coverage_stats(depth, prev, prophage_dict, prophage_lengths, min_cov, length)
                        written, total = write_coverages(outfile, prophage_covs, prev, avg, med, sd, eff_host, effect, ratio_cutoff, written, total, min_breadth, min_cov, spaces)
                        depth = np.zeros(length)
                except NameError:
                    depth = np.zeros(length)

                prev = genome
                start = x.reference_start # 0-based
                end = x.reference_end
                if end:
                    depth = add_depth_no_ed(depth, start, end)

        if length: # not in keep
            # last one
            if mask_f:
                depth[:mask_f] = np.nan
                depth[mask_r:] = np.nan
            prophage_covs, avg, med, sd, eff_host = coverage_stats(depth, prev, prophage_dict, prophage_lengths, min_cov, length)
            written, total = write_coverages(outfile, prophage_covs, prev, avg, med, sd, eff_host, effect, ratio_cutoff, written, total, min_breadth, min_cov, spaces)
            depth = None

    bamfile.close()
    if clean and bai:
        subprocess.run(f'rm {bam}.bai', shell=True)

    include_zeros(outfile, prophage_lengths, written, prophage_dict_frags, lengths, spaces)

    return total

@jit(nopython=True)
def cohenD(phage_mean, phage_sd, host_mean, host_sd):
    """
    Cohen's d equation
    """
    try:
        pool = ((phage_sd**2+host_sd**2)/2)**0.5
        d = abs((host_mean-phage_mean)/pool)
    except Exception:
        # host has 0 coverage
        d = 0
    return d

@jit(nopython=True)
def activity(phage_mean,cov_depth,avg,d,effect,ratio_cutoff,min_breadth,min_cov,total):
    '''
    Determine activity based on cutoffs
    '''
    try:
        ratio = phage_mean/avg
    except Exception:
        ratio = phage_mean
    diff = phage_mean-avg
    active = 'dormant'
    if d >= effect and ratio >= ratio_cutoff:
        if cov_depth >= min_breadth and phage_mean >= min_cov:
            active = 'active'
            total += 1
        else:
            active = 'ambiguous'
    return active,total,diff,ratio

def write_coverages(outfile, prophage_covs, host, avg, med, sd, length, effect, ratio_cutoff, written, total, min_breadth, min_cov, spaces):
    '''
    Within the BAM loop.
    Perform statistical analyses and write out final results.
    '''
    with open(outfile, 'a') as output:
        for key,p in prophage_covs.items():
            phage_mean,phage_med,phage_sd,l,cov_depth = p
            d = cohenD(phage_mean, phage_sd, avg, sd)
            active,total,diff,ratio = activity(phage_mean,cov_depth,avg,d,effect,ratio_cutoff,min_breadth,min_cov,total)
            if spaces:
                key = key.replace("~!!~", " ")
                host = host.replace("~!!~", " ")
            output.write(f'{key}\t{host}\t{active}\t{d}\t{ratio}\t{diff}\t{l}\t{phage_mean}\t{phage_med}\t{phage_sd}\t{cov_depth}\t{length}\t{avg}\t{med}\t{sd}\n')
            written.append(key)
    return written, total

def include_zeros(outfile, prophage_lengths, written, prophage_dict_frags, lengths, spaces):
    '''
    Write out all the prophages not identified within the BAM file
    '''
    written = set(written)
    with open(outfile, 'a') as output:
        for key,host in prophage_dict_frags.items():
            l = prophage_lengths[key]
            length = lengths[host]
            if spaces:
                key = key.replace("~!!~", " ")
                host = host.replace("~!!~", " ")
            if key in written: continue
            output.write(f'{key}\t{host}\tnot present\tNA\tNA\tNA\t{l}\tNA\tNA\tNA\tNA\t{length}\tNA\tNA\tNA\n')

try:
    import warnings
    warnings.filterwarnings("ignore")
    import argparse
    import subprocess
    import sys
    import time
    import datetime
    import logging
    import os

except Exception as e:
    sys.stderr.write("\nError: please verify dependancy imports are installed and up to date:\n\n")
    sys.stderr.write(str(e) + "\n\n")
    exit(1)

# Set up variables
start = time.time()

descript = '''
    PropagAtE: Prophage Activity Estimator (v1.1.0)

    Using a prophage coordinates file, fasta file and coverage information,
    calculate if prophages were active in the given sample.
    Prophages should be connected (integrated) to the host scaffold/genome.

    Example: input paired reads and run Bowtie2
    Propagate -f scaffolds.fasta -r forward.fastq reverse.fastq -o output_folder -v prophage_coordinates.tsv -t threads

    Example: input BAM alignment file
    Propagate -f scaffolds.fasta -b alignment.bam -o output_folder -v prophage_coordinates.tsv

'''

propagate = argparse.ArgumentParser(description=descript, formatter_class=argparse.RawTextHelpFormatter, usage=argparse.SUPPRESS)
propagate.add_argument('--version', action='version', version='PropagAtE v1.1.0')
required = propagate.add_argument_group('REQUIRED')
coverage = propagate.add_argument_group('PICK ONE')
common = propagate.add_argument_group('COMMON')
edit = propagate.add_argument_group('EDIT METHODS')

# Input / Output
required.add_argument('-f', metavar='', type=str, nargs=1, required = True, help='input genomes/scaffolds (can have extra sequences in file)')
required.add_argument('-v', metavar='', type=str, nargs=1, required = True, help='VIBRANT "integrated_prophage_coordinates" file or custom file (see README)')
#
coverage.add_argument('-b', metavar='', type=str, nargs=1, default = [''], help='input BAM sequence alignment file')
coverage.add_argument('-s', metavar='', type=str, nargs=1, default = [''], help='input SAM sequence alignment file')
coverage.add_argument('-r', metavar='', type=str, nargs=2, default = ['',''], help='input paired read files separated by a space (forward reverse)')
coverage.add_argument('-i', metavar='', type=str, nargs=1, default = [''], help='input interleaved paired read file')
coverage.add_argument('-u', metavar='', type=str, nargs=1, default = [''], help='input unpaired read file')
#
common.add_argument('-o', metavar='', type=str, nargs=1, default = [''], help='name of output folder [default = Propagate_results_(-v)]')
common.add_argument('-t', metavar='', type=str, nargs=1, default = ['1'], help='threads [1]')
#
edit.add_argument('-p', metavar='', type=str, nargs=1, default = ['0.97'], help='minimum percent identity per aligned read for calculating coverage [0.97]')
edit.add_argument('-e', metavar='', type=str, nargs=1, default = ['0.70'], help="minimum effect size for significance by Cohen's d test [default=0.70, minimum=0.60]")
edit.add_argument('-c', metavar='', type=str, nargs=1, default = ['2.0'], help="minimum prophage:host coverage ratio for significance [default=2.0, minimum=1.50]")
edit.add_argument('--mask', metavar='', type=str, nargs=1, default = ['150'], help="mask coverage values <int> bases on each end of a scaffold [150]")
edit.add_argument('--min', metavar='', type=str, nargs=1, default = ['1.0'], help="minimum average coverage to consider prophage present and for --breadth [1.0]")
edit.add_argument('--breadth', metavar='', type=str, nargs=1, default = ['0.5'], help="minimum breadth of coverage as fraction of bases >= minimum coverage (--min) [0.5]")
edit.add_argument('--clean', action='store_true', help='remove generated SAM, unsorted BAM, Bowtie2 index. Retain user input files and sorted BAM [off]')

# Parse arguments
args = propagate.parse_args()
samfile = str(args.s[0])
bamfile = str(args.b[0])
clean_bam = False
if bamfile:
    clean_bam = True
fasta = str(args.f[0])
forward = str(args.r[0])
reverse = str(args.r[1])
if forward and not reverse:
    sys.stderr.write("\nA reverse read set must be provided with the forward. Exiting.\n")
    exit(1)
interleaved = str(args.i[0])
unpaired = str(args.u[0])
threads = str(args.t[0])
vibe = str(args.v[0])
try:
    temp = vibe.rsplit("/",1)[1]
    u = temp.rsplit(".",1)[0]
except Exception:
    u = vibe.rsplit(".",1)[0]
mask = int(args.mask[0])
min_cov = float(args.min[0])
min_breadth = float(args.breadth[0])
effect = float(args.e[0])
read_id = float(args.p[0])
ratio_cutoff = float(args.c[0])
outpath = str(args.o[0])
if not outpath:
    outpath = f'PropagAtE_results_{u}/'
if outpath[-1] != '/':
    outpath += '/'

# make sure inputs are correct
exist = [
    does_exist(outpath, 'output folder'),
    not_exist(fasta, 'fasta'),
    not_exist(vibe, 'coordinates')
]

if any(exist):
    exit(1)

vibe_header = prophages_check(vibe)
if not vibe_header:
    exit(1)

# verify inputs
check = [samfile, bamfile, forward, interleaved, unpaired]
check = [c for c in check if c != '']
if len(check) > 1 or not check:
    sys.stderr.write(f"\nOnly one input file (-s, -b, -r, -i, -u) is allowed. {len(check)} provided. Exiting.\n")
    exit(1)

if forward and reverse:
    if not forward.endswith('.fastq') and not forward.endswith('.fastq.gz'):
        sys.stderr.write("\nError: Provided paired reads files must both have the extension .fastq or .fastq.gz. Exiting.\n")
        sys.stderr.write(f"{forward}\n")
        exit(1)
    if not reverse.endswith('.fastq') and not reverse.endswith('.fastq.gz'):
        sys.stderr.write("\nError: Provided paired reads files must both have the extension .fastq or .fastq.gz. Exiting.\n")
        sys.stderr.write(f"{reverse}\n\n")
        exit(1)
if interleaved:
    if not interleaved.endswith('.fastq') and not interleaved.endswith('.fastq.gz'):
        sys.stderr.write("\nError: Provided interleaved reads file must have the extension .fastq or .fastq.gz. Exiting.\n")
        sys.stderr.write(f"{interleaved}\n\n")
        exit(1)
if unpaired:
    if not unpaired.endswith('.fastq') and not unpaired.endswith('.fastq.gz'):
        sys.stderr.write("\nError: Provided unpaired reads file must have the extension .fastq or .fastq.gz. Exiting.\n")
        sys.stderr.write(f"{unpaired}\n\n")
        exit(1)

if samfile:
    if not_exist(samfile, 'sam file'):
        exit(1)
    if not bamfile.endswith('.sam'):
        sys.stderr.write("\nError: Provided sam file must have the extension .sam. Exiting.\n")
        exit(1)
if bamfile:
    if not_exist(bamfile, 'bam file'):
        exit(1)
    if not bamfile.endswith('.bam'):
        sys.stderr.write("\nError: Provided bam file must have the extension .bam. Exiting.\n")
        exit(1)

if effect < 0.6:
    sys.stderr.write("\nError: Cohen's d effect size (-e) should not be set below 0.6. Exiting.\n")
    exit(1)
if min_breadth > 1:
    sys.stderr.write("\nError: breadth (--breadth) should be a decimal value <= 1. Exiting.\n")
    exit(1)
if ratio_cutoff < 1.5:
    sys.stderr.write("\nError: ratio cutoff (-c) should not be set below 1.5. Exiting.\n")
    exit(1)
if read_id > 1:
    sys.stderr.write("\nError: percent identity (-p) should be a decimal value <= 1. Exiting.\n")
    exit(1)
read_id = 1.0 - read_id


# set up folder and log
subprocess.run(f'mkdir {outpath}', shell=True)
if outpath.count('/') > 1:
    base = outpath.rsplit("/",2)[1]
else:
    base = outpath[:-1]
outfile = f'{outpath}{base}.tsv'
logfilename = f'{outpath}{base}.log'
logging.basicConfig(filename=logfilename, level=logging.INFO, format='%(message)s')

##### ----------------------------------------------------------------------------------------------------------------------- #####
logging.info("Command:    %s" % ' '.join(sys.argv))
logging.info("")
logging.info("Date:       %s" % str(datetime.date.today()))
logging.info("Time:       %s" % str(datetime.datetime.now().time()).rsplit(".",1)[0])
logging.info("Program:    PropagAtE v1.1.0\n")

logging.info("Time (min) |  Log                                                   ")
logging.info("--------------------------------------------------------------------")

if vibe_header == 'no prophages':
    logging.info("%s         No prophages were found in the input coordinates file" % str(round((time.time() - float(start))/60,1)))
    logging.info("%s         Analysis finished" % str(round((time.time() - float(start))/60,1)))
    logging.info("")
    logging.info("")
    logging.info("Results file:      %s" % outfile.replace(outpath,''))
    logging.info("Active prophages:   0")
    logging.info("")
    os.mkdir(outpath)
    open(outpath + "/" + outpath + ".tsv").close()

# check for spaces in fasta
spaces = False
try:
    check = subprocess.check_output(f'grep -c " " {fasta}', shell=True)
    check = int(spaces.strip("'").strip("b"))
    if check > 0:
        spaces = True
except Exception:
    spaces = True

# If input is reads/fasta run Bowtie2
if forward or interleaved or unpaired:
    logging.info("%s         Reads input identified, using %s threads to run Bowtie2." % (str(round((time.time() - float(start))/60,1)),threads))

    try:
        subprocess.check_output("which bowtie2", shell=True)
    except Exception:
        sys.stderr.write("\nError: Bowtie2 does not appear to be installed or is not in the system's PATH. Exiting.\n")
        logging.info("\nError: Bowtie2 does not appear to be installed or is not in the system's PATH. Exiting.\n")
        exit(1)
    try:
        temp = fasta.rsplit("/",1)[1]
        base = temp.rsplit(".",1)[0]
    except Exception:
        base = fasta.rsplit(".",1)[0]

    build = bowtie2_build(fasta, base, outpath, spaces, u)
    if forward:
        if not os.path.exists(forward) or not os.path.exists(reverse):
            sys.stderr.write("\nError: the forward and/or reverse reads files do not exist. Exiting.\n\n")
            exit(1)
        run_bowtie2_paired(base, outpath, threads, forward, reverse)
    elif interleaved:
        if not os.path.exists(interleaved):
            sys.stderr.write("\nError: the interleaved reads file does not exist. Exiting.\n\n")
            exit(1)
        run_bowtie2_interleaved(base, outpath, threads, interleaved)
    elif unpaired:
        if not os.path.exists(unpaired):
            sys.stderr.write("\nError: the unpaired reads file does not exist. Exiting.\n\n")
            exit(1)
        run_bowtie2_unpaired(base, outpath, threads, unpaired)
    bamfile = post_bowtie2(outpath, base, args.clean, u, build, threads)

if samfile:
    logging.info("%s         Converting SAM file to BAM format" % str(round((time.time() - float(start))/60,1)))
    bamfile = sam_bam(samfile, outpath, threads)

if bamfile:
    logging.info("%s         Checking if BAM file needs to be sorted" % str(round((time.time() - float(start))/60,1)))
    check_bam, bamfile = sort_bam(bamfile, outpath, threads, args.clean, clean_bam)
    if not check_bam:
        sys.stderr.write("\nError: The SAM/BAM file appears to be empty, was not converted properly, or no reads aligned. Exiting.\n")
        logging.info("\nError: The SAM/BAM file appears to be empty, was not converted properly, or no reads aligned. Exiting.\n")
        exit(1)

# 'bamfile' is now from read alignment, sam, unsorted bam, or direct input

# Read in prophage coordinate data
logging.info("%s         Generating a list of all prophage regions" % str(round((time.time() - float(start))/60,1)))
check, prophage_dict, prophage_dict_frags, genomes_full, prophages, prophage_lengths = process_coordinates_file(vibe, spaces, vibe_header)

number_hosts = len(genomes_full)
logging.info("%s         Number of prophage regions identified: %s" % (str(round((time.time() - float(start))/60,1)),prophages))
logging.info("%s         Number of unique host regions identified: %s" % (str(round((time.time() - float(start))/60,1)),number_hosts))

# process sam/bam files
logging.info("%s         Extracting coverage and performing statistical analyses" % str(round((time.time() - float(start))/60,1)))
lengths = get_lengths(fasta, spaces, genomes_full)
total = extract_coverage(bamfile, read_id, lengths, mask, outfile, prophage_dict, effect, ratio_cutoff, prophage_dict_frags, prophage_lengths, min_breadth, min_cov, args.clean, spaces)

logging.info("%s         Analysis finished" % str(round((time.time() - float(start))/60,1)))
logging.info("")
logging.info("")
logging.info("Results file:      %s" % outfile.replace(outpath,''))
logging.info("Active prophages:  %s" % total)
logging.info("")
logging.info('                                                               ##')
logging.info('                                                             ##  ##')
logging.info('                                                           ##      ##')
logging.info('######   ##  ##     ##     #######   ######    #####       ##      ##')
logging.info('##  ##   ##  ##   ##  ##   ##        ##       ##             ##  ##')
logging.info('######   ######   ######   ##  ###   ######    ###             ##')
logging.info('##       ##  ##   ##  ##   ##   ##   ##           ##           ##')
logging.info('##       ##  ##   ##  ##   #######   ######   #####            ##')
logging.info('                                                            #  ##  #')
logging.info('                                                           # # ## # #')
logging.info('                                                          #   #  #   #')
logging.info('                                                         #            #')
logging.info("")


#
#
#

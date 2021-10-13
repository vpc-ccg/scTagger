#!/usr/bin/env python3
import sys
import argparse
from collections import Counter
from itertools import combinations,groupby
from operator import itemgetter

import edlib
from tqdm import tqdm
from matplotlib import pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract long-read barcodes using short-read adapter alignment")
    parser.add_argument("-r",
                        "--reads",
                        nargs="+",
                        type=str,
                        required=True,
                        help="Space separated paths to reads in FASTQ")
    parser.add_argument("-t",
                        "--threads",
                        default=1,
                        type=int,
                        help="Number of threads. Default: 1")
    parser.add_argument("-sa",
                        "--short-read-adapter",
                        type=str,
                        default='CTACACGACGCTCTTCCGATCT',
                        help="Short-read adapter. Default: CTACACGACGCTCTTCCGATCT")
    parser.add_argument("-o",
                        "--outfile",
                        type=str,
                        default=None,
                        help="Path to output file. Default: stdout")
    args = parser.parse_args()
    assert args.threads > 0
    return args


rev_compl_l = [chr(i) for i in range(128)]
rev_compl_l[ord('A')] = 'T'
rev_compl_l[ord('C')] = 'G'
rev_compl_l[ord('G')] = 'C'
rev_compl_l[ord('T')] = 'A'

def rev_compl(s):
    return ''.join(rev_compl_l[ord(c)] for c in reversed(s))

def get_lr_bc_matches(reads, a1):
    a2 = rev_compl(a1)
    lr_bc_matches = list()
    for fastq in reads:
        for idx,l in tqdm(enumerate(open(fastq))):
            if idx % 0 == 0:
                rname = l.split()[0][1:]
                b = ''
                d = '-1'
                loc = 'NA'
            if idx % 4 != 1:
                continue
            l = l.rstrip()
            aln_1 = edlib.align(a1, l,'HW', 'locations')
            aln_2 = edlib.align(a2, l,'HW', 'locations')
            if aln_1['editDistance'] == aln_2['editDistance']:
                pass
            elif aln_1['editDistance'] < aln_2['editDistance']:
                locs = [x[1] for x in aln_1['locations']]
                if 35 <= min(locs) <= max(locs) <=50 or 120 <= min(locs) <= max(locs) <= 160:
                    b = l[loc-2:max(locs)+20]
                    d = aln_1['editDistance']
                    loc = min(locs)
            else:
                locs = [x[0]-len(l)-1 for x in aln_2['locations']]
                if -40 >= max(locs) >= min(locs) >= -60 or -80 >= max(locs) >= min(locs) >= -120:
                    b = l[min(locs)-20:loc+2]
                    d = aln_2['editDistance']
                    loc = max(locs)
            lr_bc_matches.append((
                rname,
                b,
                d,
                l,
            ))
    return lr_bc_matches

def output_lr_bc_matches(lr_bc_matches, outfile):
    for rname,b,d,l in tqdm(lr_bc_matches):
        print(rname,d,l,b, sep='\t', file=outfile)

def main():
    args = parse_args()

    lr_bc_matches = get_lr_bc_matches(args.reads, args.short_read_adapter)
    if args.outfile:
        args.outfile = open(args.outfile, 'w+')
    else:
        args.outfile = sys.stdout
    output_lr_bc_matches(lr_bc_matches, args.outfile)
    args.outfile.close()

if __name__ == "__main__":
    main()
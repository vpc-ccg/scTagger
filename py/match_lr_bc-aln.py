#!/usr/bin/env python3
import sys
import argparse
from multiprocessing import Pool
import gzip
# from collections import Counter
# from itertools import combinations,groupby
# from operator import itemgetter

import edlib
from tqdm import tqdm
from matplotlib import pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract long-read barcodes using short-read adapter alignment")
    parser.add_argument("-lr",
                        "--long-read-segments",
                        type=str,
                        required=True,
                        help="Long-read segments TSV file")
    parser.add_argument("-sr",
                        "--short-read-barcodes",
                        type=str,
                        required=True,
                        help="Short-read barcode list TSV file")
    parser.add_argument("-t",
                        "--threads",
                        default=1,
                        type=int,
                        help="Number of threads. Default: 1")
    parser.add_argument("-o",
                        "--outfile",
                        type=str,
                        default=None,
                        help="Path to output file. Output file is gzipped. STDOUT is in normal text. Default: stdout")
    args = parser.parse_args()
    assert args.threads > 0
    return args

sr_barcodes = set()

rev_compl_l = [chr(i) for i in range(128)]
rev_compl_l[ord('A')] = 'T'
rev_compl_l[ord('C')] = 'G'
rev_compl_l[ord('G')] = 'C'
rev_compl_l[ord('T')] = 'A'

def rev_compl(s):
    return ''.join(rev_compl_l[ord(c)] for c in reversed(s))

def get_matches(lb):
    matches = list()
    d = float('inf')
    if lb == '':
        return tuple(matches), d
    for idx, sb in enumerate(sr_barcodes):
        aln = edlib.align(sb, lb,'HW', 'distance')
        if aln['editDistance'] < d:
            matches = list()
            d = aln['editDistance']
        if aln['editDistance'] == d:
            matches.append(idx)
    return tuple(matches), d


def get_lr_segments(long_read_segments_tsv):
    gz = long_read_segments_tsv.endswith('gz')
    lr_names = list()
    lr_segments = list()
    if gz:
        f = gzip.open(long_read_segments_tsv, 'rt')
    else:
        f = open(long_read_segments_tsv, 'r')
    for l in tqdm(f):
        # if len(lr_names) >= 10000: break
        l = l.rstrip('\n').split('\t')
        lr_names.append(l[0])
        lr_segments.append(l[3])
    f.close()
    return lr_names,lr_segments



def run_get_matches(lr_segments, threads):
    lr_matches = list()
    lr_dists = list()
    with Pool(threads) as p:
        for matches, d in tqdm(p.imap(get_matches, lr_segments, chunksize=100), total=len(lr_segments)):
            lr_matches.append(matches)
            lr_dists.append(d)
    return lr_matches, lr_dists

def get_sr_barcodes(barcodes_tsv):
    gz = barcodes_tsv.endswith('gz')
    global sr_barcodes
    sr_barcodes = set()
    if gz:
        f = gzip.open(barcodes_tsv)
    else:
        f = open(barcodes_tsv)
    for l in f:
        if gz: l = l.decode()
        l = l.rstrip().split('\t')
        b = l[0]
        sr_barcodes.add(b)
        sr_barcodes.add(rev_compl(b))
    sr_barcodes = list(sr_barcodes)

def output_lr_bc_matches(lr_names, lr_segments, lr_matches, lr_dists, outfile):
    for n,s,m,d in tqdm(zip(lr_names, lr_segments, lr_matches, lr_dists), total=len(lr_names)):
        c = len(m)
        m = ','.join([sr_barcodes[bid] for bid in m])
        outfile.write(f'{n}\t{d}\t{c}\t{s}\t{m}\n')

def main():
    args = parse_args()
    print(args)

    get_sr_barcodes(args.short_read_barcodes)
    lr_names,lr_segments = get_lr_segments(args.long_read_segments)
    lr_matches, lr_dists = run_get_matches(lr_segments, args.threads)
 
    if args.outfile:
        args.outfile = gzip.open(args.outfile, 'wt')
    else:
        args.outfile = sys.stdout
    output_lr_bc_matches(
        lr_names=lr_names, 
        lr_segments=lr_segments,
        lr_matches=lr_matches, 
        lr_dists=lr_dists, 
        outfile=args.outfile,
    )
    args.outfile.close()

if __name__ == "__main__":
    main()
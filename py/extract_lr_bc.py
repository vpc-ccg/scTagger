#!/usr/bin/env python3
import enum
import sys
import argparse
from multiprocessing import Pool
import gzip
from collections import deque

import numpy as np
import pandas as pd
import matplotlib.ticker as mtick
import matplotlib.pyplot as plt
from tqdm import tqdm
import edlib


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract long-read barcodes using short-read adapter alignment")
    parser.add_argument("-r",
                        "--reads",
                        nargs="+",
                        type=str,
                        required=True,
                        help="Space separated paths to reads in FASTQ")
    parser.add_argument("-g",
                        "--ranges",
                        nargs="+",
                        type=str,
                        default=list(),
                        help="""Space separated of the ranges of where SR adapter should be found on the LR's.
                                Intervals should have start with "f" or "r" to indicated which side of the read the interval is suppused to be.
                                E.g. f20:40 r1:30 means SR adapter needs to be either on the range 20-40bp (20 and 40 both included) 
                                and on the forward strand OR the last 30bp from the end of the read and on the reverse strand.
                                If adapters are found on multiple ranges on the read, then all of them are considered invalid.
                                Ranges are 1-indexed and inclusive.
                                Default: Detect from data.""")
    parser.add_argument("-z",
                        "--gzipped",
                        dest='gzipped',
                        action='store_true',
                        help="Indicate input is gzipped. Default: Assume input is gzipped if it ends with \".gz\".")
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
                        help="Path to output file. Output file is gzipped. STDOUT is in normal text. Default: stdout")
    parser.add_argument("-p",
                        "--plotfile",
                        type=str,
                        default=None,
                        help="Path to plot file. Default: no plotting")

    args = parser.parse_args()
    ranges = [list(),list()]
    ranges_dict = [dict(), dict()]
    for r in args.ranges:
        assert r[0] in 'fr', r
        strand = r[0]
        r = r[1:].split(':')
        assert len(r) == 2,r
        s = int(r[0])
        e = int(r[1])
        assert 0 < s <= e, (s,e) 
        if strand == 'f':
            ranges_idx = 0
            s, e = s - 1, e
        elif strand == 'r':
            ranges_idx = 1
            s, e = -e, -s+1
        else:
            assert False, (strand)
        # print(strand, r, (s,e))
        for i in np.arange(s,e):
            assert not i in ranges_dict[ranges_idx], (ranges_idx, i, ranges_dict[ranges_idx])
            # print(i, end=', ')
            ranges_dict[ranges_idx][i] = len(ranges[ranges_idx])
        ranges[ranges_idx].append((s,e))
        # print()

    args.ranges = ranges        
    args.ranges_dict = ranges_dict        
    assert args.threads > 0
    return args


rev_compl_l = [chr(i) for i in range(128)]
rev_compl_l[ord('A')] = 'T'
rev_compl_l[ord('C')] = 'G'
rev_compl_l[ord('G')] = 'C'
rev_compl_l[ord('T')] = 'A'

def rev_compl(s):
    return ''.join(rev_compl_l[ord(c)] for c in reversed(s))

def read_fastqs(fastqs, gzipped):
    rnames = list()
    seqs = list()
    for fastq in fastqs:
        print(f'Reading {fastq}', file=sys.stderr)
        if gzipped or fastq.endswith('.gz'):
            f = gzip.open(fastq, 'rt')
        else:
            f = open(fastq, 'r')
        for idx,l in tqdm(enumerate(f)):
            if idx % 4 == 0:
                rnames.append(l.split()[0][1:])
            if idx % 4 == 1:
                seqs.append(l.rstrip())
    return rnames,seqs

def get_alns(seq):
    d = -1
    s = 'NA'
    locs = ['NA']
    aln_1 = edlib.align(a1, seq,'HW', 'locations')
    aln_2 = edlib.align(a2, seq,'HW', 'locations')
    if aln_1['editDistance'] == aln_2['editDistance']:
        pass
    elif aln_1['editDistance'] < aln_2['editDistance']:
        locs = [x[1] for x in aln_1['locations']]
        d = aln_1['editDistance']
        s = '+'
    else:
        locs = [x[0]-len(seq)-1 for x in aln_2['locations']]
        d = aln_2['editDistance']
        s = '-'
    return (
        s,
        d,
        locs,
    )

def get_ranges(data):
    ranges = list()
    min_l = min(data)
    max_l = max(data)
    L = np.arange(min_l, max_l+1)
    F = np.zeros(max_l-min_l+1)
    for l in data:
        F[l-min_l]+=1
    T = np.sum(F)
    while True:
        PEAK = np.argmax(F)
        neighborhood_sum = sum(F[max(0,PEAK-20):min(PEAK+20,len(F))])
        print(f'--> {neighborhood_sum/T: 5.2%} of strend reads fall around {L[PEAK]}', file=sys.stderr)
        if sum(F[max(0,PEAK-20):min(PEAK+20,len(F))]) < 0.01*T:
            break
        Q = deque()
        Q.append(PEAK)
        first = PEAK
        last = PEAK
        while len(Q) > 0:
            i = Q.popleft()
            F[i] = 0
            if i<= PEAK and i-1>0 and F[i-1] > T*0.001:
                Q.append(i-1)
                first=i-1
            if i>= PEAK and i+1<len(F) and F[i+1] > T*0.001:
                Q.append(i+1)
                last=i+1
        for i in range(max(0,first-20),min(last+20,len(F))):
            F[i]=0
        ranges.append((L[first],L[last]))
    return ranges

def get_possible_ranges(alns):
    data_f = list()
    data_r = list()
    for s,d,locs in alns:
        if d < 0 or d > 5:
            continue
        if s == '+': 
            data = data_f
        else:
            data = data_r
        for l in locs:
            data.append(l)
    ranges_f = get_ranges(data_f)
    print(f'Found these ranges on + strand:\t{ranges_f}', file=sys.stderr)
    ranges_r = get_ranges(data_r)
    print(f'Found these ranges on - strand:\t{ranges_r}', file=sys.stderr)
    return ranges_f,ranges_r

def get_range_dicts(ranges):
    ranges_dicts = list()
    for R in ranges:
        ranges_dict = dict()
        for idx,(s,e) in enumerate(R):
            for i in np.arange(s,e):
                assert not i in ranges_dict, (i,ranges_dict)
                ranges_dict[i] = idx
        ranges_dicts.append(ranges_dict)
    return ranges_dicts

def get_lr_bc_alns(seqs, threads):
    print(f'Aligning {a1} to {len(seqs)} reads on {threads} threads', file=sys.stderr)
    alns = list()
    with Pool(threads) as p:
        for s,d,locs in p.imap(get_alns, seqs, chunksize=10000):
            alns.append((s,d,locs))
    return alns


def output_matches(rnames, seqs, alns, ranges_dicts, outfile, num_bp_after=20):
    for rname,seq,(strand,dist,locs) in tqdm(zip(rnames, seqs, alns), total=len(rnames)):
        # print(locs)
        range_idx = int(strand == '-')
        range_hits = {ranges_dicts[range_idx].get(l, -1) for l in locs}
        # print(range_idx, range_hits)
        if -1 in range_hits or len(range_hits) != 1:
            outfile.write(f'{rname}\t{-1}\t{"NA"}\t{""}\n')
        else:
            if strand == '+':
                s = max(0,        min(locs)-2)
                e = min(len(seq), max(locs)+num_bp_after)
                l = s
            else:
                s = max(-len(seq), min(locs)-num_bp_after)
                e = min(0, max(locs)+2)
                l = e
            outfile.write(f'{rname}\t{dist}\t{l}\t{seq[s:e or None]}\n')

def show_plots(lr_bc_matches, outfile):
    names = ["read_id", "distance", "strand", "segment"]
    data = pd.DataFrame(lr_bc_matches, columns=names)

    new_data = data.groupby("distance").count().reset_index()
    new_data = new_data.rename(index={1:"0", 2: "1", 3: "2", 4: "3", 5: "4", 6: "5", 7: "6", 8: "7", 9: "8",
                                      10: "9", 11: "10", 0: 'NA'})

    target_row = 0
    idx =  [i for i in range(len(new_data)) if i != target_row] + [target_row]
    new_data = new_data.iloc[idx]
    new_data["read_id"].cumsum(axis=0)
    new_data["cumSum"] = new_data["read_id"].cumsum(axis=0)
    new_data["cumSumPer"] = new_data["read_id"].cumsum(axis=0) / len(data) * 100
    fig = plt.figure()

    ax = fig.add_subplot(111)
    ax2 = ax.twinx()

    width = 0.2

    new_data.read_id.plot(kind='bar', color='red', ax=ax, width=width, position=1)
    new_data.cumSum.plot(kind='bar', color='blue', ax=ax, width=width, position=0)
    new_data.cumSumPer.plot(kind='bar', color='blue', ax=ax2, width=width, position=0)

    ax.set_ylabel('Number of Long-reads')
    ax.set_xlabel("Edit distance")
    ax2.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax2.set_ylabel('Percentage of Long-reads')

    plt.savefig(outfile)

def main():
    args = parse_args()
    print(args, file=sys.stderr)
    global a1,a2
    a1 = args.short_read_adapter
    a2 = rev_compl(a1)

    rnames,seqs = read_fastqs(args.reads, args.gzipped)
    alns = get_lr_bc_alns(seqs, args.threads)
    if len(args.ranges[0]) + len(args.ranges[1]) == 0:
        print('No ranges for SR adapters have been preset. Detecting directly from data...', file=sys.stderr)
        args.ranges = get_possible_ranges(alns)
        args.ranges_dict = get_range_dicts(args.ranges)
    # print(args.ranges)
    # print(args.ranges_dict)
    if args.outfile:
        args.outfile = gzip.open(args.outfile, 'wt+')
    else:
        args.outfile = sys.stdout
    output_matches(rnames, seqs, alns, args.ranges_dict, args.outfile)
    args.outfile.close()
    # output_lr_bc_matches(lr_bc_matches, args.outfile)
    # args.outfile.close()
    # if args.plotfile != None:
    #     show_plots(lr_bc_matches, args.plotfile)

if __name__ == "__main__":
    main()

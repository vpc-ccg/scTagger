#!/usr/bin/env python3
import sys
import argparse
from multiprocessing import Pool
import gzip
import edlib
import pandas as pd
import matplotlib.ticker as mtick

from tqdm import tqdm
import matplotlib.pyplot as plt


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
                        help="Path to output file. Output file is gzipped. STDOUT is in normal text. Default: stdout")
    parser.add_argument("-p",
                        "--plotfile",
                        type=str,
                        default=None,
                        help="Path to plot file.")

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

def get_aln(seq):
    b = ''
    d = '-1'
    loc = 'NA'
    aln_1 = edlib.align(a1, seq,'HW', 'locations')
    aln_2 = edlib.align(a2, seq,'HW', 'locations')
    if aln_1['editDistance'] == aln_2['editDistance']:
        pass
    elif aln_1['editDistance'] < aln_2['editDistance']:
        locs = [x[1] for x in aln_1['locations']]
        if 35 <= min(locs) <= max(locs) <=50 or 120 <= min(locs) <= max(locs) <= 160:
            loc = min(locs)
            b = seq[loc-2:max(locs)+20]
            d = aln_1['editDistance']
    else:
        locs = [x[0]-len(seq)-1 for x in aln_2['locations']]
        if -40 >= max(locs) >= min(locs) >= -60 or -80 >= max(locs) >= min(locs) >= -120:
            loc = max(locs)
            b = seq[min(locs)-20:loc+2]
            d = aln_2['editDistance']
    return (
        b,
        d,
        loc,
    )


def get_lr_bc_matches(fastqs, threads):
    lr_bc_matches = list()
    rnames = list()
    seqs = list()
    for fastq in fastqs:
        print(f'Reading {fastq}')
        for idx,l in tqdm(enumerate(open(fastq))):
            if idx % 4 == 0:
                rnames.append(l.split()[0][1:])
            if idx % 4 == 1:
                seqs.append(l.rstrip())
    print(f'Aligning {a1} to {len(seqs)} reads on {threads} threads')
    with Pool(threads) as p:
        for n,(b,d,loc) in tqdm(zip(rnames,p.imap(get_aln, seqs, chunksize=100)), total=len(seqs)):
            lr_bc_matches.append((n,b,d,loc))
    return lr_bc_matches

def output_lr_bc_matches(lr_bc_matches, outfile):
    for rname,b,d,l in tqdm(lr_bc_matches):
        outfile.write(f'{rname}\t{d}\t{l}\t{b}\n')

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
    print(args)
    global a1,a2
    a1 = args.short_read_adapter
    a2 = rev_compl(a1)


    lr_bc_matches = get_lr_bc_matches(args.reads, threads=args.threads)
    if args.outfile:
        args.outfile = gzip.open(args.outfile, 'wt')
    else:
        args.outfile = sys.stdout
    output_lr_bc_matches(lr_bc_matches, args.outfile)
    args.outfile.close()
    show_plots(lr_bc_matches, args.plotfile)

if __name__ == "__main__":
    main()

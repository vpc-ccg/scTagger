#!/usr/bin/env python3
from math import ceil
import sys
import argparse
import gzip
from collections import Counter

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract selected short-read barcodes")
    parser.add_argument("-i",
                        "--input",
                        type=str,
                        required=True,
                        help="Input file")
    parser.add_argument("-o",
                        "--outfile",
                        type=str,
                        default=None,
                        help="Path to output file. Default: STDOUT")
    parser.add_argument("-p",
                        "--plotfile",
                        type=str,
                        default=None,
                        help="Path to plot file")
    parser.add_argument("--thresh",
                        type=float,
                        default=0.005,
                        help="Percentage theshold required per step to continue adding read barcodes. Default: 0.005")
    parser.add_argument("--step-size",
                        type=int,
                        default=1000,
                        help="Number of barcodes processed at a time and whose sum is used to check against the theshold. Default: 1000")
    parser.add_argument("--max-barcode-cnt",
                        type=int,
                        default=25_000,
                        help="Max number of barcodes to keep. Default: 25000")
    args = parser.parse_args()
    assert 0 <= args.thresh <= 1
    assert 0 < args.step_size 
    assert 0 < args.max_barcode_cnt 
    return args

def read_short_reads(path):
    if path.endswith('.gz'):
        f = gzip.open(path, 'rt')
    else:
        f = open(path, 'r')
    barcodes = list()
    total = 0
    for line in tqdm(f):
        total+=1
        _,cb,_ = line.rstrip('\n').split('\t')
        if cb == 'NA':
            continue
        barcodes.append(cb)
    return total,barcodes

def get_barcode_hist(barcode_cnts, total, step_size):
    remaining = total
    distribution = {}
    for idx, x in enumerate(barcode_cnts, start=1):
        if idx % step_size == 0:
            distribution[idx] = 1 - remaining / total
        remaining -= x[1]
    if not idx % step_size == 0:
        distribution[idx] = 1 - remaining / total
    return distribution

def show_plot(distribution, step_size, last_idx, outfile):
    x = sorted(distribution.keys())
    y1 = [distribution[idx]*100 for idx in x]
    y2 = [distribution[idx]*100 for idx in x]
    for idx in range(1, len(y2)):
        y2[idx] = y1[idx] - y1[idx-1] 
    fig = plt.figure(figsize=(10,5))
    fig.suptitle(f'SR coverage with each additional {step_size} unique barcodes')
    ax1 = fig.add_subplot(111)
    plt.xticks(
        range(step_size, max(x), step_size*ceil(max(x)/step_size/18)),
        rotation=45,
    )
    ax2 = ax1.twinx()

    plot_lines = list()
    plot_lines.extend(
        ax1.plot(x, y1, color='#1b9e77', label='Cumulative % coverage (left y-axis)')
    )
    plot_lines.extend(
        ax2.plot(x, y2, color='#d95f02', label=f'Coverage (right y-axis)')
    )
    ax2.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax1.yaxis.set_major_formatter(mtick.PercentFormatter())
    plot_lines.extend(
        ax2.plot([last_idx,last_idx], [min(y2), max(y2)], color='#7570b3', label='Selected barcodes',ls='dashed')
    )
    plt.legend(plot_lines, [l.get_label() for l in plot_lines], loc='center right')
    plt.savefig(outfile)

def main():
    args = parse_args()
    total, barcodes = read_short_reads(args.input)
    barcode_cnts = Counter(barcodes)
    barcode_cnts_tuple = sorted(barcode_cnts.items(), key=lambda x: x[1], reverse=True)[:args.max_barcode_cnt]
    barcodes,barcode_cnts = tuple(zip(*barcode_cnts_tuple))
    barcode_hist = get_barcode_hist(
        barcode_cnts=barcode_cnts_tuple,
        step_size=args.step_size,
        total=total,
    )

    last_idx = len(barcodes)
    last_f = 0
    for idx,f in sorted(barcode_hist.items()):
        if idx == 0:
            continue
        if f-last_f < args.thresh:
            last_idx = min(
                idx,
                len(barcodes)
            )
            break
        last_f = f
    if args.plotfile != None:
        show_plot(
            distribution=barcode_hist,
            step_size=args.step_size,
            last_idx=last_idx,
            outfile=args.plotfile,
        )
    if args.outfile:
        outfile = gzip.open(args.outfile, 'wt+')
    else:
        outfile = sys.stdout
    for b, c in zip(barcodes[:last_idx],barcode_cnts[:last_idx]):
        outfile.write(f'{b}\t{c}\n')
    outfile.close()

if __name__ == "__main__":
    main()

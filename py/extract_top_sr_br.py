#!/usr/bin/env python3
import sys
import argparse
import gzip
import pandas as pd

from tqdm import tqdm

from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


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
                        default=50_000,
                        help="Max number of barcodes to keep. Default: 50,000")
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

def get_barcode_hist(barcode_cnts, step_size, max_barcode_cnt):
    distribution = list()
    for idx in range(0, min(max_barcode_cnt,len(barcode_cnts)), step_size):
        cur_barcode_cnts = barcode_cnts[idx:idx+step_size]
        distribution.append(
            sum(cur_barcode_cnts)
        )        
    return distribution


def show_plot(short_read_barcode, distribution, outfile):
    x = []
    y1 = []
    y2 = []
    start = 0

    total = len(short_read_barcode)
    for key in distribution:
        x.append(key)
        y1.append(total - ((1 - distribution[key]) * total))
        y2.append(distribution[key] - start)
        start = distribution[key]

    df = pd.DataFrame({'cumulative': y1,
                       'Change': y2}, index=x)
    df["cumulative"] = df["cumulative"] / total * 100
    df["Change"] = df["Change"] * 100
    df = df.drop(labels=0, axis=0)


    fig = plt.figure()

    ax = fig.add_subplot(111)
    ax2 = ax.twinx()

    df.cumulative.plot(kind='line', color='red', ax=ax)
    df.Change.plot(kind='line', color='blue', ax=ax2)

    ax.set_ylabel('Percentage of cumulative reads cover')
    ax.set_xlabel("Selected barcode")
    ax2.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax2.set_ylabel('Percentage of cover change in reads')

    plt.savefig(outfile)

def main():
    args = parse_args()
    total, barcodes = read_short_reads(args.input)
    barcode_cnts = Counter(barcodes)
    barcode_cnts = sorted(barcode_cnts.items(), key=lambda x: x[1], reverse=True)
    barcodes,barcode_cnts = tuple(zip(*barcode_cnts))
    barcode_hist = get_barcode_hist(
        barcode_cnts=barcode_cnts,
        step_size=args.step_size,
        max_barcode_cnt=args.max_barcode_cnt,
    )
    if args.plotfile != None:
        show_plot(
            total=total,
            distribution=distribution,
            outfile=args.plotfile,
        )

    last_idx = args.max_barcode_cnt
    for idx,c in enumerate(barcode_hist):
        if c / total < args.thresh:
            last_idx = min(
                (idx+1) * args.step_size,
                args.max_barcode_cnt,
                len(barcodes)
            )
            break
    if args.outfile:
        outfile = gzip.open(args.outfile, 'wt+')
    else:
        outfile = sys.stdout
    for b, c in zip(barcodes[:last_idx],barcode_cnts[:last_idx]):
        outfile.write(f'{b}\t{c}\n')
    outfile.close()

if __name__ == "__main__":
    main()

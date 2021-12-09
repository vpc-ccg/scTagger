#!/usr/bin/env python3
import sys
import argparse
import gzip
import edlib
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
    parser.add_argument("-t",
                        "--threads",
                        default=1,
                        type=int,
                        help="Number of threads. Default: 1")
    parser.add_argument("-o",
                        "--outfile",
                        type=str,
                        default=None,
                        help="Path to output file")
    parser.add_argument("-p",
                        "--plotfile",
                        type=str,
                        default=None,
                        help="Path to plot file")

    args = parser.parse_args()
    assert args.threads > 0
    return args

def rev(x):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(x))

def read_short_reads(path):
    file1 = open(path, 'r')
    Lines = file1.readlines()

    barcodes_list = []
    distribution = dict()
    count = 0
    for line in tqdm(Lines):
        b = line.replace("\n", "").split("\t")
        if len(set(b[2]) - set("ATCG")) == 0:
            count += (b[1] == b[2])
            if b[1] != b[2]:
                distance = edlib.align(b[2], b[1], mode="HW")["editDistance"]
                distance_rev = edlib.align(b[2], rev(b[1]), mode="HW")["editDistance"]
                barcodes_list.append((b[2], min(distance, distance_rev)))
            else:
                barcodes_list.append((b[2], 0))

    short_read_barcode = [x[0] for x in barcodes_list]
    return short_read_barcode


def select_barcode(inputPath, thresh, steps, endpoint):
    sr_barcodes_list = read_short_reads(inputPath)
    diff = 1
    start = 0
    short_reads_count = Counter(sr_barcodes_list)
    sort_orders = sorted(short_reads_count.items(), key=lambda x: x[1], reverse=True)
    total = sum([x[1] for x in sort_orders])
    remaining = total
    selected_barcode = []
    distribution = {}
    for idx, x in enumerate(sort_orders):
        if diff > thresh:
            selected_barcode.append(x)

        if idx % steps == 0:
            distribution[idx] = 1 - remaining / total
            diff = (1 - remaining / total) - start
            start = (1 - remaining / total)
        remaining -= x[1]

        if idx > endpoint:
            break

    return selected_barcode, distribution, sort_orders


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
    sr_barcodes_list = read_short_reads(args.input)
    short_read_barcode = [x[0] for x in sr_barcodes_list[0]]
    selected_barcode, distribution, sort_orders = select_barcode(short_read_barcode, 0.005, 1000, 50000)
    show_plot(short_read_barcode, distribution, args.plotfile)

    outfile = gzip.open(args.outfile, 'wt')
    for b, c in sort_orders:
        outfile.write(f'{b}\t{c}\n')
    outfile.close()

if __name__ == "__main__":
    main()

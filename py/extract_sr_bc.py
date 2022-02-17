#!/usr/bin/env python3
import sys
import argparse
import gzip

import pysam
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract long-read barcodes using short-read adapter alignment")
    parser.add_argument("-i",
                        "--input",
                        type=str,
                        required=True,
                        help="BAM file of short reads alignmed by 10x Genomics Cell Ranger")
    parser.add_argument("-t",
                        "--threads",
                        default=1,
                        type=int,
                        help="Number of threads. Default: 1")
    parser.add_argument("-o",
                        "--outfile",
                        type=str,
                        required=True,
                        help="Path to output file. Output file is gzipped.")
    args = parser.parse_args()
    assert args.threads > 0
    return args

def output_sr_bc(bam, outpath, threads):
    outfile = gzip.open(outpath, 'wt+')
    for aln in tqdm(pysam.AlignmentFile(bam, 'rb', threads=threads)):
        if aln.flag > 256:
            continue
        tags = dict(aln.tags)
        C = tags.get('CB', 'NA').split('-')[0]
        U = tags.get('UB', 'NA').split('-')[0]
        qname = aln.query_name
        outfile.write(f'{qname}\t{C}\t{U}\n')
    outfile.close()


def main():
    args = parse_args()
    print(args, file=sys.stderr)
    output_sr_bc(args.bam, args.output, args.threads)
    
if __name__ == "__main__":
    main()

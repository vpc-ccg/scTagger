#!/usr/bin/env python3
import sys
import argparse
from multiprocessing import Pool, RawArray
import gzip
from collections import deque, Counter
from math import ceil

from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib.ticker as mtick
import matplotlib.pyplot as plt

import edlib
import pysam
import ahocorasick


def parse_args():
    parser = argparse.ArgumentParser(description="scTagger pipeline!")
    subparsers = parser.add_subparsers(dest='subcommand')
    subparsers.required = True

    parser_lr_bc = subparsers.add_parser('extract_lr_bc')
    # add a required argument
    parser_lr_bc.add_argument("-r", "--reads", nargs="+", type=str, required=True,
                              help="Space separated paths to reads in FASTQ")
    parser_lr_bc.add_argument("-g", "--ranges", nargs="+", type=str, default=list(),
                              help="""Space separated of the ranges of where SR adapter should be found on the LR's.
                                Intervals should have start with "f" or "r" to indicated which side of the read the interval is suppused to be.
                                E.g. f20:40 r1:30 means SR adapter needs to be either on the range 20-40bp (20 and 40 both included) 
                                and on the forward strand OR the last 30bp from the end of the read and on the reverse strand.
                                If adapters are found on multiple ranges on the read, then all of them are considered invalid.
                                Ranges are 1-indexed and inclusive.
                                Default: Detect from data.""")
    parser_lr_bc.add_argument("-z", "--gzipped", dest='gzipped', action='store_true',
                              help="Indicate input is gzipped. Default: Assume input is gzipped if it ends with \".gz\".")
    parser_lr_bc.add_argument("-t", "--threads", default=1, type=int,
                              help="Number of threads. Default: 1")
    parser_lr_bc.add_argument("-sa", "--short-read-adapter", type=str,
                              default='CTACACGACGCTCTTCCGATCT',
                              help="Short-read adapter. Default: CTACACGACGCTCTTCCGATCT")
    parser_lr_bc.add_argument("-o", "--outfile", type=str, default=None,
                              help="Path to output file. Output file is gzipped. STDOUT is in normal text. Default: stdout")
    parser_lr_bc.add_argument("-p", "--plotfile", type=str, default=None,
                              help="Path to plot file. Default: no plotting")
    parser_lr_bc.add_argument("--num-bp-after", type=int, default=20,
                              help="Number of bases after the end of the SR adapter alignment to generate. Default: 20")

    parser_top_sr = subparsers.add_parser('extract_sr_bc')
    parser_top_sr.add_argument("-i", "--input", type=str, required=True,
                               help="Input BAM file")
    parser_top_sr.add_argument("-o", "--outfile", type=str, default=None,
                               help="Path to output file. Default: STDOUT")
    parser_top_sr.add_argument("-p", "--plotfile", type=str, default=None,
                               help="Path to plot file")
    parser_top_sr.add_argument("-t", "--threads", default=1, type=int,
                               help="Number of threads. Default: 1")
    parser_top_sr.add_argument("--thresh", type=float, default=0.005,
                               help="Percentage theshold required per step to continue adding read barcodes. Default: 0.005")
    parser_top_sr.add_argument("--step-size", type=int, default=1000,
                               help="Number of barcodes processed at a time and whose sum is used to check against the theshold. Default: 1000")
    parser_top_sr.add_argument("--max-barcode-cnt", type=int, default=25_000,
                               help="Max number of barcodes to keep. Default: 25000")

    parser_sr_bc_from_lr = subparsers.add_parser('extract_sr_bc_from_lr')
    parser_sr_bc_from_lr.add_argument("-i", "--input", type=str, required=True,
                               help="Input TSV file generated from extract_lr_bc step")
    parser_sr_bc_from_lr.add_argument("-o", "--outfile", type=str, default=None,
                               help="Path to output file. Default: STDOUT")
    parser_sr_bc_from_lr.add_argument("-wl", "--barcode-whitelist", type=str, required=True,
                               help="Path to TXT barcode whiltelist such as 10x Chromimum v3 3M-february-2018.txt.gz barcode file")
    parser_sr_bc_from_lr.add_argument("--thresh", type=float, default=0.005,
                               help="Percentage theshold required per step to continue adding read barcodes. Default: 0.005")
    parser_sr_bc_from_lr.add_argument("--step-size", type=int, default=1000,
                               help="Number of barcodes processed at a time and whose sum is used to check against the theshold. Default: 1000")
    parser_sr_bc_from_lr.add_argument("--max-barcode-cnt", type=int, default=25_000,
                               help="Max number of barcodes to keep. Default: 25000")

    parser_match_trie = subparsers.add_parser('match_trie')
    parser_match_trie.add_argument("-lr", "--long-read-segments", type=str,
                                   required=True, help="Long-read segments TSV file")
    parser_match_trie.add_argument("-sr", "--short-read-barcodes", type=str,
                                   required=True, help="Short-read barcode list TSV file")
    parser_match_trie.add_argument("-mr", "--max-error", default=2, type=int,
                                   help="Maximum number of errors allowed for barcode matching. Default: 2")
    parser_match_trie.add_argument("-m", "--mem", default=16.0, type=float,
                                   help="Maximum number of GB of RAM to be used. Default: 16.0")
    parser_match_trie.add_argument("-bl", "--barcode-length", default=16, type=int,
                                   help="Length of barcodes. Default: 16")
    parser_match_trie.add_argument("-t", "--threads", default=16, type=int,
                                   help="Number of threads to use for searching. Default: 16")
    parser_match_trie.add_argument("-p", "--plotfile", default=None, type=str,
                                   help="Path of plot file. Default: no plotting")
    parser_match_trie.add_argument("-o", "--outfile", type=str, default=None,
                                   help="Path to output file. Output file is gzipped. STDOUT is in normal text. Default: stdout")
    args = parser.parse_args()

    if args.subcommand == 'extract_lr_bc':
        assert 0 < args.num_bp_after
        ranges = [list(), list()]
        global ranges_dicts, num_bp_after
        num_bp_after = args.num_bp_after
        ranges_dicts = [dict(), dict()]
        for r in args.ranges:
            assert r[0] in 'fr', r
            strand = r[0]
            r = r[1:].split(':')
            assert len(r) == 2, r
            s = int(r[0])
            e = int(r[1])
            assert 0 < s <= e, (s, e)
            if strand == 'f':
                ranges_idx = 0
                s, e = s - 1, e
            elif strand == 'r':
                ranges_idx = 1
                s, e = -e, -s + 1
            else:
                assert False, (strand)
            for i in np.arange(s, e):
                assert not i in ranges_dicts[ranges_idx], (
                    ranges_idx, i, ranges_dicts[ranges_idx])
                ranges_dicts[ranges_idx][i] = len(ranges[ranges_idx])
            ranges[ranges_idx].append((s, e))
        args.ranges = ranges
        assert args.threads > 0

    if args.subcommand == 'extract_sr_bc':
        assert 0 <= args.thresh <= 1
        assert 0 < args.step_size
        assert 0 < args.max_barcode_cnt

    if args.subcommand == 'extract_sr_bc_from_lr':
        assert 0 <= args.thresh <= 1
        assert 0 < args.step_size
        assert 0 < args.max_barcode_cnt

    if args.subcommand == 'match_trie':
        assert args.mem > 0
        assert args.barcode_length > 0
        assert args.barcode_length > args.max_error >= 0

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
        for idx, l in tqdm(enumerate(f)):
            if idx % 4 == 0:
                rnames.append(l.split()[0][1:])
            if idx % 4 == 1:
                seqs.append(l.rstrip())
    return rnames, seqs


def get_alns(seq):
    d = -1
    s = 'NA'
    locs = ['NA']
    aln_1 = edlib.align(a1, seq, 'HW', 'locations')
    aln_2 = edlib.align(a2, seq, 'HW', 'locations')
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
        F[l-min_l] += 1
    T = np.sum(F)
    while True:
        PEAK = np.argmax(F)
        neighborhood_sum = sum(F[max(0, PEAK-20):min(PEAK+20, len(F))])
        print(
            f'--> {neighborhood_sum/T: 5.2%} of strend reads fall around {L[PEAK]}', file=sys.stderr)
        if sum(F[max(0, PEAK-20):min(PEAK+20, len(F))]) < 0.01*T:
            break
        Q = deque()
        Q.append(PEAK)
        first = PEAK
        last = PEAK
        while len(Q) > 0:
            i = Q.popleft()
            F[i] = 0
            if i <= PEAK and i-1 > 0 and F[i-1] > T*0.001:
                Q.append(i-1)
                first = i-1
            if i >= PEAK and i+1 < len(F) and F[i+1] > T*0.001:
                Q.append(i+1)
                last = i+1
        for i in range(max(0, first-20), min(last+20, len(F))):
            F[i] = 0
        ranges.append((L[first], L[last]))
    return ranges


def get_possible_ranges(alns):
    data_f = list()
    data_r = list()
    for s, d, locs in alns:
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
    return ranges_f, ranges_r


def set_global_range_dicts(ranges):
    global ranges_dicts
    ranges_dicts = list()
    for R in ranges:
        ranges_dict = dict()
        for idx, (s, e) in enumerate(R):
            for i in np.arange(s, e):
                assert not i in ranges_dict, (i, ranges_dict)
                ranges_dict[i] = idx
        ranges_dicts.append(ranges_dict)


def get_lr_bc_alns(seqs, threads):
    print(
        f'Aligning {a1} to {len(seqs)} reads on {threads} threads', file=sys.stderr)
    alns = list()
    with Pool(threads) as p:
        for s, d, locs in p.imap(get_alns, seqs, chunksize=10000):
            alns.append((s, d, locs))
    return alns


def get_filter_alns(alns, threads):
    print(
        f'Filtering alignments using ranges on {threads} threads', file=sys.stderr)
    fin_alns = list()
    with Pool(threads) as p:
        for dist, loc, s, e in tqdm(p.imap(filter_aln, alns, chunksize=10000), total=len(alns)):
            fin_alns.append((dist, loc, s, e))
    return fin_alns


def filter_aln(aln):
    strand, dist, locs = aln
    range_idx = int(strand == '-')
    range_hits = {ranges_dicts[range_idx].get(l, -1) for l in locs}
    if -1 in range_hits or len(range_hits) != 1:
        dist = -1
        loc = 'NA'
        s = e = -1
    else:
        if strand == '+':
            s = max(0, min(locs)-2)
            e = max(locs)+num_bp_after
            loc = s
        else:
            s = min(locs)-num_bp_after
            e = min(0, max(locs)+2)
            loc = e
    return dist, loc, s, e


def get_lr_bc_alns(seqs, threads):
    print(
        f'Aligning {a1} to {len(seqs)} reads on {threads} threads', file=sys.stderr)
    alns = list()
    with Pool(threads) as p:
        for s, d, locs in tqdm(p.imap(get_alns, seqs, chunksize=10000), total=len(seqs)):
            alns.append((s, d, locs))
    return alns


def output_matches(rnames, seqs, alns, outfile):
    print(f'Writng to {outfile}', file=sys.stderr)
    for rname, seq, (dist, loc, s, e) in tqdm(zip(rnames, seqs, alns), total=len(rnames)):
        outfile.write(f'{rname}\t{dist}\t{loc}\t{seq[s:e or None]}\n')


def show_plots_extract_lr_bc(rnames, alns, outfile):
    names = ["read_id", "distance"]
    data = pd.DataFrame([
        (r, d) for r, (d, _, _, _) in zip(rnames, alns)
    ], columns=names)

    new_data = data.groupby("distance").count().reset_index()
    new_data = new_data.rename(index={1: "0", 2: "1", 3: "2", 4: "3", 5: "4", 6: "5", 7: "6", 8: "7", 9: "8",
                                      10: "9", 11: "10", 0: 'NA'})

    target_row = 0
    idx = [i for i in range(len(new_data)) if i != target_row] + [target_row]
    new_data = new_data.iloc[idx]
    new_data["read_id"].cumsum(axis=0)
    new_data["cumSum"] = new_data["read_id"].cumsum(axis=0)
    new_data["cumSumPer"] = new_data["read_id"].cumsum(
        axis=0) / len(data) * 100
    fig = plt.figure(figsize=(10, 5))

    ax = fig.add_subplot(111)
    ax2 = ax.twinx()

    width = 0.2

    new_data.read_id.plot(kind='bar', color='red',
                          ax=ax, width=width, position=1)
    new_data.cumSum.plot(kind='bar', color='blue',
                         ax=ax, width=width, position=0)
    new_data.cumSumPer.plot(kind='bar', color='blue',
                            ax=ax2, width=width, position=0)

    ax.set_ylabel('Number of Long-reads')
    ax.set_xlabel("Edit distance")
    ax2.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax2.set_ylabel('Percentage of Long-reads')

    plt.savefig(outfile)


def extract_lr_bc(args):
    global a1, a2
    a1 = args.short_read_adapter
    a2 = rev_compl(a1)

    rnames, seqs = read_fastqs(args.reads, args.gzipped)
    alns = get_lr_bc_alns(seqs, args.threads)
    if len(args.ranges[0]) + len(args.ranges[1]) == 0:
        print('No ranges for SR adapters have been preset. Detecting directly from data...', file=sys.stderr)
        args.ranges = get_possible_ranges(alns)
        set_global_range_dicts(args.ranges)
    alns = get_filter_alns(alns, args.threads)
    # print(args.ranges)
    # print(args.ranges_dict)
    if args.outfile:
        args.outfile = gzip.open(args.outfile, 'wt+')
    else:
        args.outfile = sys.stdout
    output_matches(rnames, seqs, alns, args.outfile)
    args.outfile.close()
    if args.plotfile != None:
        show_plots_extract_lr_bc(rnames, alns, args.plotfile)


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


def plot_sr_bc_coverage(distribution, step_size, last_idx, outfile):
    x = sorted(distribution.keys())
    y1 = [distribution[idx]*100 for idx in x]
    y2 = [distribution[idx]*100 for idx in x]
    for idx in range(1, len(y2)):
        y2[idx] = y1[idx] - y1[idx-1]
    fig = plt.figure(figsize=(10, 5))
    fig.suptitle(
        f'SR coverage with each additional {step_size} unique barcodes')
    ax1 = fig.add_subplot(111)
    plt.xticks(
        range(step_size, max(x), step_size*ceil(max(x)/step_size/18)),
        rotation=45,
    )
    ax2 = ax1.twinx()

    plot_lines = list()
    plot_lines.extend(
        ax1.plot(x, y1, color='#1b9e77',
                 label='Cumulative % coverage (left y-axis)')
    )
    plot_lines.extend(
        ax2.plot(x, y2, color='#d95f02', label=f'Coverage (right y-axis)')
    )
    ax2.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax1.yaxis.set_major_formatter(mtick.PercentFormatter())
    plot_lines.extend(
        ax2.plot([last_idx, last_idx], [min(y2), max(y2)],
                 color='#7570b3', label='Selected barcodes', ls='dashed')
    )
    plt.legend(plot_lines, [l.get_label()
               for l in plot_lines], loc='center right')
    plt.savefig(outfile)


def extract_sr_barcode(bamfile, threads):
    print(f'\n====\nExtracting SR barcodes from {bamfile}:')
    barcodes = list()
    total = 0
    alns_per_contig = {x.contig: x.total for x in pysam.AlignmentFile(
        bamfile, 'rb').get_index_statistics()}
    imap_args = [
        (bamfile, x['SN'])
        for x in pysam.AlignmentFile(bamfile, 'rb').header['SQ']
    ]
    p = Pool(threads)
    with tqdm(total=sum(alns_per_contig.values())) as tqdm_bar:
        for contig, c_total, c_barcodes in tqdm(p.imap_unordered(read_bam_contig, imap_args, chunksize=1), total=len(imap_args)):
            tqdm_bar.update(alns_per_contig[contig])
            barcodes.extend(c_barcodes)
            total += c_total
    p.close()
    return total, barcodes


def read_bam_contig(args):
    bamfile, contig = args
    barcodes = list()
    total = 0
    for aln in pysam.AlignmentFile(bamfile, 'rb').fetch(contig=contig):
        if aln.flag > 256:
            continue
        tags = dict(aln.tags)
        C = tags.get('CB', 'NA').split('-')[0]
        total += 1
        if C == 'NA':
            continue
        barcodes.append(C)
    return contig, total, barcodes


def extract_sr_bc(args):
    total, barcodes = extract_sr_barcode(args.input, args.threads)

    print("\n=====\nCounting and sorting barcodes")
    barcode_cnts = Counter(barcodes)
    barcode_cnts_tuple = sorted(barcode_cnts.items(
    ), key=lambda x: x[1], reverse=True)[:args.max_barcode_cnt]
    barcodes, barcode_cnts = tuple(zip(*barcode_cnts_tuple))
    barcode_hist = get_barcode_hist(
        barcode_cnts=barcode_cnts_tuple,
        step_size=args.step_size,
        total=total,
    )

    last_idx = len(barcodes)
    last_f = 0
    for idx, f in sorted(barcode_hist.items()):
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
        plot_sr_bc_coverage(
            distribution=barcode_hist,
            step_size=args.step_size,
            last_idx=last_idx,
            outfile=args.plotfile,
        )
    print(f"\n=====\nWriting the top {last_idx} barcodes")
    if args.outfile:
        outfile = gzip.open(args.outfile, 'wt+')
    else:
        outfile = sys.stdout
    for b, c in tqdm(zip(barcodes[:last_idx], barcode_cnts[:last_idx]), total=last_idx):
        outfile.write(f'{b}\t{c}\n')
    outfile.close()


map_char = [0 for _ in range(128)]
map_char[ord('A')] = 0
map_char[ord('C')] = 1
map_char[ord('G')] = 2
map_char[ord('T')] = 3

alphabet_size = 4
barcodes = None
barcodes_rc = None
lr_segs = None
lr_segs_idxs = None
lr_names = None


class Trie(object):
    def __init__(self, max_error, barcode_length, capacity=100_000_000):
        self.nid_to_rids = [list()]
        self.nid_to_child = np.zeros(
            (alphabet_size, capacity),
            dtype=np.uint32,
        )
        self.size = 1
        self.max_error = max_error
        self.barcode_length = barcode_length

    def __increase_capacity__(self):
        self.nid_to_child = np.concatenate(
            (
                self.nid_to_child,
                np.zeros(
                    self.nid_to_child.shape,
                    dtype=np.uint32,
                )
            ),
            axis=1
        )
        print(f'Doubled capacity. New shape: {self.nid_to_child.shape}')

    def insert(self, word, read_id):
        nid = 0
        for index, char in enumerate(word):
            cid = map_char[ord(char)]
            # Needs to create node
            if self.nid_to_child[cid][nid] == 0:
                new_nid = self.size
                self.size += 1
                if self.size >= self.nid_to_child.shape[1]:
                    self.__increase_capacity__()
                self.nid_to_child[cid][nid] = new_nid
                self.nid_to_rids.append(list())
            nid = self.nid_to_child[cid][nid]
            if index >= self.barcode_length - self.max_error - 1:
                self.nid_to_rids[nid].append(read_id)

    def dfs(self, nid, barcode, error, index, result):
        if index == len(barcode):
            edit_distance = self.max_error - error
            result[edit_distance].extend(
                self.nid_to_rids[nid]
            )
            return
        # delete
        if error > 0:
            self.dfs(nid, barcode, error - 1, index + 1, result)
        for cid in range(alphabet_size):
            child_id = self.nid_to_child[cid][nid]
            if child_id == 0:
                continue
            # dont change
            if cid == map_char[ord(barcode[index])]:
                self.dfs(child_id, barcode, error, index + 1, result)
            # mutation
            if cid != map_char[ord(barcode[index])] and error > 0:
                self.dfs(child_id, barcode, error - 1, index + 1, result)
            # insert
            if error > 0:
                self.dfs(child_id, barcode, error - 1, index, result)

    def query(self, barcode):
        result = [list() for _ in range(self.max_error + 1)]
        self.dfs(
            nid=0,
            barcode=barcode,
            error=self.max_error,
            index=0,
            result=result,
        )
        return result


def read_long_reads(path):
    global lr_segs_idxs, lr_segs, lr_names
    print('Reading long reads barcodes')
    if path.endswith('.gz'):
        f = gzip.open(path, 'rt')
    else:
        f = open(path, 'r')
    lr_segs_tmp = list()
    lr_names = list()
    for l in f:
        l = l.rstrip('\n').split('\t')
        lr_names.append(l[0])
        lr_segs_tmp.append(l[3])

    N = sum(len(s) for s in lr_segs_tmp)
    lr_segs = RawArray('c', N)
    lr_segs_idxs = RawArray('i', len(lr_segs_tmp) + 1)
    a = 0
    b = 0
    for s in lr_segs_tmp:
        lr_segs_idxs[b] = a
        for c in s:
            lr_segs[a] = ord(c)
            a += 1
        b += 1
    lr_segs_idxs[b] = a
    print(f'There are {len(lr_names):,} LRs')


def read_sr_barcodes(path):
    print('Reading short reads barcodes')
    global barcodes, barcodes_rc
    if path.endswith('.gz'):
        f = gzip.open(path, 'rt')
    else:
        f = open(path)
    barcodes = [l.rstrip('\n').split('\t')[0] for l in f]
    barcodes_rc = [rev_compl(b) for b in barcodes]
    print(f'There are {len(barcodes):,} SR barcodes')


def trie_search(trie):
    result = dict()
    for bid in range(len(barcodes)):
        for e, rids in enumerate(trie.query(barcodes[bid])):
            for rid in rids:
                if not rid in result:
                    result[rid] = (
                        len(barcodes[0]) + 1,  # Edit dist
                        set(),  # Tuples of (bid,strand)
                    )
                if e < result[rid][0]:
                    result[rid] = (e, set())
                if e == result[rid][0]:
                    result[rid][1].add((bid, True))
        for e, rids in enumerate(trie.query(barcodes_rc[bid])):
            for rid in rids:
                if not rid in result:
                    result[rid] = (
                        len(barcodes[0]) + 1,  # Edit dist
                        set(),  # Tuples of (bid,strand)
                    )
                if e < result[rid][0]:
                    result[rid] = (e, set())
                if e == result[rid][0]:
                    result[rid][1].add((bid, False))
    return result


def get_matches_prefix(args):
    max_error, prefix = args
    barcode_length = len(barcodes[0])
    trie = Trie(max_error, barcode_length)
    trie.insert(prefix.decode(), -1)
    result = {}
    for rid, (idx_1, idx_2) in enumerate(zip(lr_segs_idxs[:-1], lr_segs_idxs[1:])):
        segment = lr_segs[idx_1:idx_2]
        for i in range(len(segment)):
            sub_segment = segment[i:i + barcode_length + max_error]
            if sub_segment.startswith(prefix) and len(sub_segment) >= barcode_length - max_error:
                trie.insert(sub_segment.decode(), rid)
    result = trie_search(trie)
    del trie
    return result


def run_get_matches_prefix(max_error, threads):
    prefixes = [b'']
    while threads > len(prefixes):  # *.5:
        new_prefixes = list()
        for p in prefixes:
            for c in 'ACGT':
                new_prefixes.append(
                    (p.decode() + c).encode()
                )
            prefixes = new_prefixes
    args = [(max_error, p) for p in prefixes]
    print(args)
    result = dict()
    print(
        f'Running get_matches_prefix on {threads} threads and {len(args)} tasks')
    with Pool(threads) as p:
        for proc_result in tqdm(p.imap(get_matches_prefix, args, chunksize=1), total=len(args)):
            for rid, (e, X) in proc_result.items():
                if not rid in result:
                    result[rid] = (
                        len(barcodes[0]) + 1,  # Edit dist
                        set(),  # Tuples of (bid,strand)
                    )
                if e < result[rid][0]:
                    result[rid] = (e, set())
                if e == result[rid][0]:
                    for bid, strand in X:
                        result[rid][1].add((bid, strand))
    return result


def show_plot_match_trie(full_data, plotfile, max_error):
    read_id = []
    barcodes = []
    distance = []

    for key in full_data:

        read_id.append(key)
        find_dict = full_data[key]
        tmp_barcodes = []
        not_find = True
        for index in range(max_error + 1):
            if len(find_dict[index]) > 0:
                tmp_barcodes = list(find_dict[index])
                not_find = False
                distance.append(index)
                break

        if not_find:
            distance.append(-1)

        barcodes.append(tmp_barcodes)

    trie_dataframe = pd.DataFrame(
        {"read_id": read_id, "barcodes": barcodes, "distance": distance})
    new_data = trie_dataframe.groupby("distance").count()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    width = 0.2

    new_data.read_id.plot(kind='bar', color='red',
                          ax=ax, width=width, position=1)

    ax.set_ylabel('Number of long-reads')
    ax.set_xlabel("Edit distance")

    plt.savefig(plotfile)


def match_trie(args):
    read_sr_barcodes(args.short_read_barcodes)
    barcode_lens = {len(b) for b in barcodes}
    assert barcode_lens == {args.barcode_length}, barcode_lens
    del barcode_lens

    read_long_reads(args.long_read_segments)
    result = run_get_matches_prefix(
        max_error=args.max_error,
        threads=args.threads,
    )

    if args.outfile:
        if args.outfile.endswith('gz'):
            outfile = gzip.open(args.outfile, 'wt+')
        else:
            outfile = open(args.outfile, 'w+')
    else:
        outfile = sys.stdout
    for rid, (e, bids) in tqdm(sorted(result.items()), total=len(result)):
        name = lr_names[rid]
        idx_1 = lr_segs_idxs[rid]
        idx_2 = lr_segs_idxs[rid + 1]
        seg = lr_segs[idx_1:idx_2].decode()
        outfile.write(f'{name}\t')
        if len(bids) == 0:
            e = 'inf'
        outfile.write(f'{e}\t')
        outfile.write(f'{len(bids)}\t')
        outfile.write(f'{seg}\t')
        outfile.write(
            f'{",".join([barcodes[b] if s else rev_compl(barcodes[b]) for b, s in sorted(bids)])}\n')
    outfile.close()

def extract_sr_bc_from_lr(args):
    if args.barcode_whitelist.endswith('.gz'):
        infile = gzip.open(args.barcode_whitelist, 'rt')
    else:
        infile = open(args.barcode_whitelist)
    print(f'Reading whiltelist barcodes from: {args.barcode_whitelist}')
    barcodes = [l[:-1] for l in tqdm(infile)]
    infile.close()

    A = ahocorasick.Automaton()
    print(f'\n=====\nAdding forward barcodes...')
    for idx,bc in tqdm(enumerate(barcodes), total=len(barcodes)):
        A.add_word(bc,idx)
    print(f'\n=====\nAdding reverse compliment barcodes...')
    for idx,bc in tqdm(enumerate(barcodes), total=len(barcodes)):
        A.add_word(rev_compl(bc),-idx)
    print(f'\n=====\nBuilding Aho-Corasick automaton...')
    A.make_automaton()    


    if args.input.endswith('.gz'):
        infile = gzip.open(args.input, 'rt')
    else:
        infile = open(args.input)
    print(f'\n=====\nMatching exact barcodes on long-reads: {args.input}')
    C = Counter()
    for l in tqdm(infile):
        _,_,p,seg = l.rstrip('\n').split('\t')
        if p=='NA':
            continue
        hits = tuple(A.iter(seg))
        if len(hits) > 1:
            continue
        for _,bc in hits:
            C[abs(bc)]+=1
    print(f'\n=====\nFound {len(C):,} unique barcodes on long-reads') 
    sorted_bc = sorted(C.items(), reverse=True, key=lambda x: x[1])[:args.max_barcode_cnt]
    
    total = sum(c for _,c in sorted_bc)
    for last_idx in range(0, len(sorted_bc), args.step_size):
        percent = sum(c for _,c in sorted_bc[last_idx:last_idx+args.step_size])/total
        if percent < args.thresh:
            break
    sorted_bc = sorted_bc[:last_idx+args.step_size]

    print(f"\n=====\nWriting the top {len(sorted_bc)} barcodes")
    if args.outfile:
        outfile = gzip.open(args.outfile, 'wt+')
    else:
        outfile = sys.stdout
    for bc,c in tqdm(sorted_bc):
        outfile.write(f'{barcodes[bc]}\t{c}\n')
    outfile.close()


def main():
    args = parse_args()
    print(args)

    if args.subcommand == 'extract_lr_bc':
        extract_lr_bc(args)

    if args.subcommand == 'extract_sr_bc':
        extract_sr_bc(args)

    if args.subcommand == 'extract_sr_bc_from_lr':
        extract_sr_bc_from_lr(args)

    if args.subcommand == 'match_trie':
        match_trie(args)


if __name__ == "__main__":
    main()

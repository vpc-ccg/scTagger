#!/usr/bin/env python3
import sys
import argparse
import gzip
import math

import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import psutil
from multiprocessing import Pool

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
    parser.add_argument("-mr",
                        "--max-error",
                        default=2,
                        type=int,
                        help="Maximum number of errors allowed for barcode matching. Default: 2")
    parser.add_argument("-m",
                        "--mem",
                        default=16.0,
                        type=float,
                        help="Maximum number of GB of RAM to be used. Default: 16.0")
    parser.add_argument("-bl",
                        "--barcode-length",
                        default=16,
                        type=int,
                        help="Length of barcodes. Default: 16")
    parser.add_argument("-t",
                        "--threads",
                        default=16,
                        type=int,
                        help="Number of threads to use for searching. Default: 16")
    parser.add_argument("-p",
                        "--plotfile",
                        default=None,
                        type=str,
                        help="Path of plot file. Default: no plotting")
    parser.add_argument("-o",
                        "--outfile",
                        type=str,
                        default=None,
                        help="Path to output file. Output file is gzipped. STDOUT is in normal text. Default: stdout")
    args = parser.parse_args()
    assert args.mem > 0
    assert args.barcode_length > 0
    assert args.barcode_length > args.max_error >= 0
    
    return args

map_char = [None for i in range(128)]
map_char[ord('A')] = 0
map_char[ord('C')] = 1
map_char[ord('G')] = 2
map_char[ord('T')] = 3
map_char[ord('N')] = 4

rev_compl_l = [chr(i) for i in range(128)]
rev_compl_l[ord('A')] = 'T'
rev_compl_l[ord('C')] = 'G'
rev_compl_l[ord('G')] = 'C'
rev_compl_l[ord('T')] = 'A'

trie = None
barcodes = None
def rev_compl(s):
    return ''.join(rev_compl_l[ord(c)] for c in reversed(s))

class Trie(object):

    def __init__(self, max_error, barcode_length, capacity=100_000_000):
        self.nid_to_rids = [list()]
        self.nid_to_child = np.zeros(
            (5,capacity),
            dtype=np.uint32,
        )
        self.size = 1
        self.max_error = max_error
        self.barcode_length = barcode_length

    def __increase_capacity__(self):
        new_capacity = self.nid_to_child.shape[1]*2
        print(f'Increasing trie capacity to {new_capacity}')
        self.nid_to_child.resize(
            (5,new_capacity),
        )

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
        for cid in range(5):
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
        result = [list() for _ in range(self.max_error+1)]
        self.dfs(
            nid=0, 
            barcode = barcode, 
            error = self.max_error,
            index = 0,
            result = result,
        )
        return result

def read_long_reads(path):

    names=["read_id", "distance", "strand", "segment"]
    print('Reading long reads segments')
    long_segment = pd.read_csv(path, delimiter = "\t", names=names)

    return long_segment

def read_short_reads(path):
    print('Reading short reads barcodes')
    global barcodes
    if path.endswith('.gz'):
        f = gzip.open(path, 'rt')
    else:
        f = open(path)
    barcodes = [l.rstrip().split('\t')[0] for l in f]
    print(f'There {len(barcodes)} SR barcodes')

def trie_run_query(bid):
    result = list()
    for e,rids in enumerate(trie.query(barcodes[bid])):
        result.append((e,True,rids))
    for e,rids in enumerate(trie.query(rev_compl(barcodes[bid]))):
        result.append((e,False,rids))
    return bid,result

def trie_search(read_total, threads):
    result = [(100,set()) for _ in range(read_total)]
    with Pool(threads) as p:
        for bid,X in tqdm(p.imap(trie_run_query, range(len(barcodes))), total=len(barcodes)):
            for (e,strand,rids) in X:
                for rid in rids:
                    if result[rid][0] > e:
                        result[rid] = (e,set())
                    if result[rid][0] == e:
                        result[rid][1].add((bid,strand))
    return result

def run_get_matches_memory(long_reads, max_error, barcode_length, memory_GB, threads):
    global trie
    trie = Trie(max_error, barcode_length)
    result = {}
    print('Building trie')
    for index, row in tqdm(long_reads.iterrows(), total=long_reads.shape[0]):
        segment = row["segment"]
        if type(segment)!=str and math.isnan(segment):
            continue
        for i in range(len(segment)):
            if len(segment[i:i + barcode_length + max_error]) >= barcode_length - max_error:
                trie.insert(segment[i:i + barcode_length + max_error], index)
    print(f'Memory use: {psutil.Process().memory_info().rss / (1024**3):.2f}GB')
    result = trie_search(long_reads.shape[0], threads)
    del trie
    print(f'Memory use: {psutil.Process().memory_info().rss / (1024**3):.2f}GB')
    return result


def show_plot(full_data, plotfile, max_error):
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

    trie_dataframe = pd.DataFrame({"read_id": read_id, "barcodes": barcodes, "distance": distance})
    new_data = trie_dataframe.groupby("distance").count()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    width = 0.2

    new_data.read_id.plot(kind='bar', color='red', ax=ax, width=width, position=1)

    ax.set_ylabel('Number of long-reads')
    ax.set_xlabel("Edit distance")

    plt.savefig(plotfile)


def main():
    args = parse_args()
    print(args)
    global barcodes
    read_short_reads(args.short_read_barcodes)
    long_reads = read_long_reads(args.long_read_segments)
    result = run_get_matches_memory(
        long_reads=long_reads,
        max_error=args.max_error,
        barcode_length=args.barcode_length,
        memory_GB=args.mem,
        threads=args.threads,
    )
    memory_useage_GB = psutil.Process().memory_info().rss / (1024**3)
    print(f"Memory usage: {memory_useage_GB:.2f}GB")
    # if args.plotfile != None:
    #     show_plot(result, args.plotfile, args.max_error)
    if args.outfile:
        outfile = gzip.open(args.outfile, 'wt')
    else:
        outfile = sys.stdout
    for rid,(e,bids) in enumerate(result):
        segment = long_reads.iloc[rid]["segment"]
        if type(segment)!=str and math.isnan(segment):
            continue
        outfile.write(f'{long_reads.iloc[rid]["read_id"]}\t')
        if len(bids) == 0:
            e = 'inf'
        outfile.write(f'{e}\t')
        outfile.write(f'{len(bids)}\t')
        outfile.write(f'{segment}\t')
        outfile.write(f'{",".join([barcodes[b] if s else rev_compl(barcodes[b]) for b,s in sorted(bids)])}\n')
    outfile.close()

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import sys
import argparse
import gzip
import math
import psutil
from multiprocessing import Pool,RawArray

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

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

map_char = [0 for _ in range(128)]
map_char[ord('A')] = 0
map_char[ord('C')] = 1
map_char[ord('G')] = 2
map_char[ord('T')] = 3

rev_compl_l = [chr(i) for i in range(128)]
rev_compl_l[ord('A')] = 'T'
rev_compl_l[ord('C')] = 'G'
rev_compl_l[ord('G')] = 'C'
rev_compl_l[ord('T')] = 'A'

alphabet_size = 4 
barcodes = None
barcodes_rc = None
lr_segs = None
lr_segs_idxs = None
lr_names = None

def rev_compl(s):
    return ''.join(rev_compl_l[ord(c)] for c in reversed(s))

def mem_use_gb():
    return psutil.Process().memory_info().rss/(1024**3)

class Trie(object):
    def __init__(self, max_error, barcode_length, capacity=100_000_000):
        self.nid_to_rids = [list()]
        self.nid_to_child = np.zeros(
            (alphabet_size,capacity),
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
            axis = 1
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
    lr_segs_idxs = RawArray('i', len(lr_segs_tmp)+1)
    a = 0
    b = 0
    for s in lr_segs_tmp:
        lr_segs_idxs[b] = a
        for c in s:
            lr_segs[a]=ord(c)
            a+=1
        b+=1
    lr_segs_idxs[b] = a
    print(f'There are {len(lr_names):,} LRs')

def read_sr_barcodes(path):
    print('Reading short reads barcodes')
    global barcodes,barcodes_rc
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
        for e,rids in enumerate(trie.query(barcodes[bid])):
            for rid in rids:
                if not rid in result:
                    result[rid] = (
                        len(barcodes[0])+1, # Edit dist
                        set(),              # Tuples of (bid,strand)
                    )
                if  e < result[rid][0]:
                    result[rid] = (e,set())
                if  e == result[rid][0]:
                    result[rid][1].add((bid,True))
        for e,rids in enumerate(trie.query(barcodes_rc[bid])):
            for rid in rids:
                if not rid in result:
                    result[rid] = (
                        len(barcodes[0])+1, # Edit dist
                        set(),              # Tuples of (bid,strand)
                    )
                if  e < result[rid][0]:
                    result[rid] = (e,set())
                if  e == result[rid][0]:
                    result[rid][1].add((bid,False))
    return result

def get_matches_prefix(args):
    max_error, prefix = args
    barcode_length = len(barcodes[0])
    trie = Trie(max_error, barcode_length)
    trie.insert(prefix.decode(), -1)
    result = {}
    for rid, (idx_1,idx_2) in enumerate(zip(lr_segs_idxs[:-1],lr_segs_idxs[1:])):
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
    while threads > len(prefixes):#*.5:
        new_prefixes = list()
        for p in prefixes:
            for c in 'ACGT':
                new_prefixes.append(
                    (p.decode()+c).encode()
                )
            prefixes = new_prefixes
    args = [(max_error,p) for p in prefixes]
    print(args)
    result = dict()
    print(f'Running get_matches_prefix on {threads} threads and {len(args)} tasks')
    with Pool(threads) as p:
        for proc_result in tqdm(p.imap(get_matches_prefix, args, chunksize=1), total=len(args)):
            for rid,(e,X) in proc_result.items():
                if not rid in result:
                    result[rid] = (
                        len(barcodes[0])+1, # Edit dist
                        set(),              # Tuples of (bid,strand)
                    )
                if  e < result[rid][0]:
                    result[rid] = (e,set())
                if  e == result[rid][0]:
                    for bid,strand in X:
                        result[rid][1].add((bid,strand))
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
    read_sr_barcodes(args.short_read_barcodes)
    barcode_lens = {len(b) for b in barcodes}
    assert barcode_lens == {args.barcode_length}, barcode_lens
    del barcode_lens

    long_reads = read_long_reads(args.long_read_segments)
    result = run_get_matches_prefix(
        max_error=args.max_error,
        threads=args.threads,
    )
    memory_useage_GB = mem_use_gb()
    print(f"Memory usage: {memory_useage_GB:.2f}GB")
    # if args.plotfile != None:
    #     show_plot(result, args.plotfile, args.max_error)
    if args.outfile:
        if args.outfile.endswith('gz'):
            outfile = gzip.open(args.outfile, 'wt+')
        else:
            outfile = open(args.outfile, 'w+')
    else:
        outfile = sys.stdout
    for rid,(e,bids) in tqdm(sorted(result.items()),total=len(result)):
        name = lr_names[rid]
        idx_1 = lr_segs_idxs[rid]
        idx_2 = lr_segs_idxs[rid+1]
        seg = lr_segs[idx_1:idx_2].decode() 
        outfile.write(f'{name}\t')
        if len(bids) == 0:
            e = 'inf'
        outfile.write(f'{e}\t')
        outfile.write(f'{len(bids)}\t')
        outfile.write(f'{seg}\t')
        outfile.write(f'{",".join([barcodes[b] if s else rev_compl(barcodes[b]) for b,s in sorted(bids)])}\n')
    outfile.close()

if __name__ == "__main__":
    main()

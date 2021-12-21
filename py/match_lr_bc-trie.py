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
                        help="Number of maximum error for barcode matching. Default: 22")
    parser.add_argument("-m",
                        "--max-memory",
                        default=1e10,
                        type=int,
                        help="Number of maximum memory that can use.")
    parser.add_argument("-bl",
                        "--barcode-length",
                        default=16,
                        type=int,
                        help="Length of barcodes. Default: 16")
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
    return args

class TrieNode:
    """A node in the trie structure"""

    def __init__(self):
        self.reads = set()
        self.children = {}


class Trie(object):

    def __init__(self, max_error, barcode_length):
        self.root = TrieNode()
        self.max_error = max_error
        self.barcode_length = barcode_length

    def insert(self, word, read_id):
        node = self.root

        for index, char in enumerate(word):
            if char in node.children:
                node = node.children[char]
            else:
                new_node = TrieNode()
                node.children[char] = new_node
                node = new_node
            if index >= self.barcode_length - self.max_error - 1:
                node.reads.add(read_id)

    def dfs(self, node, barcode, error, index):
        if index == len(barcode):
            edit_distance = self.max_error - error
            self.output[edit_distance].update(node.reads)
            return
        # delete
        if error > 0:
            self.dfs(node, barcode, error - 1, index + 1)
        for char, child in node.children.items():

            # dont change
            if char == barcode[index]:
                self.dfs(child, barcode, error, index + 1)

            # mutation
            if char != barcode[index] and error > 0:
                self.dfs(child, barcode, error - 1, index + 1)

            # insert
            if error > 0:
                self.dfs(child, barcode, error - 1, index)

    def query(self, x):

        node = self.root

        self.output = {}
        for i in range(self.max_error + 1):
            self.output[i] = set()

        self.dfs(node, x, self.max_error, 0)
        for i in range(0, self.max_error):
            for j in range(i + 1, self.max_error + 1):
                self.output[j] -= self.output[i]

        return self.output

def read_long_reads(path):

    names=["read_id", "distance", "strand", "segment"]
    print('Reading long reads segments')
    long_segment = pd.read_csv(path, delimiter = "\t", names=names)

    return long_segment

def read_short_reads(path):
    print('Reading short reads barcodes')
    names=["barcode", "frequency"]
    selected_barcode = pd.read_csv(path, delimiter = "\t", names=names)
    print(f'There {selected_barcode.shape[0]} SR barcodes')
    return selected_barcode

rev_compl_l = [chr(i) for i in range(128)]
rev_compl_l[ord('A')] = 'T'
rev_compl_l[ord('C')] = 'G'
rev_compl_l[ord('G')] = 'C'
rev_compl_l[ord('T')] = 'A'

def rev_compl(s):
    return ''.join(rev_compl_l[ord(c)] for c in reversed(s))

def trie_search(barcodes, trie, max_error):
    all_barcode = set()

    for index, row in barcodes.iterrows():
        b = row['barcode']
        all_barcode.add(b)
        all_barcode.add(rev_compl(b))

    print('all_barcode')
    result = {}

    for b in tqdm(all_barcode):
        result[b] = trie.query(b)

    full_result = {}
    print('result')
    for key in tqdm(result):
        for index in result[key]:
            for read_id in result[key][index]:
                if read_id not in full_result:
                    full_result[read_id] = [set() for _ in range(max_error+1)]
                full_result[read_id][int(index)].add(key)
    return full_result

def run_get_matches(selected_barcode, long_reads, max_error, barcode_length):
    trie = Trie(max_error, barcode_length)

    print('Building trie')
    for index, row in tqdm(long_reads.iterrows(), total=long_reads.shape[0]):
        segment = row["segment"]
        if type(segment)!=str and math.isnan(segment):
            continue
        for i in range(len(segment)):
            if len(segment[i:i + barcode_length + max_error]) >= barcode_length - max_error:
                trie.insert(segment[i:i + barcode_length + max_error], index)

    return trie_search(selected_barcode, trie, max_error)


def run_get_matches_memory(selected_barcode, long_reads, max_error, barcode_length, memory):

    trie = Trie(max_error, barcode_length)
    result = {}
    print('Building trie')
    for index, row in tqdm(long_reads.iterrows(), total=long_reads.shape[0]):
        segment = row["segment"]

        if type(segment)!=str and math.isnan(segment):
            continue

        memory_useage = psutil.Process().memory_info().rss / (1024 * 1024)
        if memory < memory_useage * 1.1:
            print("Try to remove memory")

            result.update(trie_search(selected_barcode, trie, max_error))
            del trie
            trie = Trie(max_error, barcode_length)

        for i in range(len(segment)):
            if len(segment[i:i + barcode_length + max_error]) >= barcode_length - max_error:
                trie.insert(segment[i:i + barcode_length + max_error], index)
    
    result.update(trie_search(selected_barcode, trie, max_error))
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

    selected_barcode = read_short_reads(args.short_read_barcodes)
    long_reads = read_long_reads(args.long_read_segments)
    result = run_get_matches_memory(selected_barcode, long_reads,
                                    args.max_error, args.barcode_length,
                                    args.max_memory)

    memory_useage = psutil.Process().memory_info().rss / (1024 * 1024)
    print("Memory usage: ", memory_useage)
    if args.plotfile != None:
        show_plot(result, args.plotfile, args.max_error)
    if args.outfile:
        outfile = gzip.open(args.outfile, 'wt')
    else:
        outfile = sys.stdout
    for rid in sorted(result):
        outfile.write(f'{long_reads.iloc[rid]["read_id"]}\t')
        for e in range(args.max_error+1):
            if len(result[rid][e]) > 0:
                # s = ','.join([str(x) for x in result[rid][e]])
                break   
        if len(result[rid][e]) == 0:
            e = 'inf'
        outfile.write(f'{e}\t')
        outfile.write(f'{len(result[rid][e])}\t')
        outfile.write(f'{long_reads.iloc[rid]["segment"]}\t')
        outfile.write(f'{",".join(str(x) for x in sorted(result[rid][e]))}\n')
    outfile.close()
        

if __name__ == "__main__":
    main()

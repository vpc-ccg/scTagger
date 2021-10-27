#!/usr/bin/env python3
import sys
import argparse

from tqdm import tqdm
import pandas as pd



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
                        default=0,
                        type=int,
                        help="Number of maximum error for barcode matching. Default: 0")
    parser.add_argument("-bl",
                        "--barcode-length",
                        default=16,
                        type=int,
                        help="Length of barcodes. Default: 16")
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

        return self.output

def read_long_reads(path):

    names=["read_id", "distance", "stand","segment"]
    long_segment = pd.read_csv(path, delimiter = "\t", names=names)

    return long_segment

def read_short_reads(path):

    names=["barcode", "trash"]
    selected_barcode = pd.read_csv(path, delimiter = "\t", names=names)

    return selected_barcode

def rev(x):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(x))


def run_get_matches(selected_barcode, long_reads, max_error, barcode_length):
    trie = Trie(max_error, barcode_length)

    for index, row in tqdm(long_reads.iterrows()):
        segment = row["segment"]
        try:
            for i in range(len(segment)):
                if len(segment[i:i + barcode_length + max_error]) >= barcode_length - max_error:
                    trie.insert(segment[i:i + barcode_length + max_error], row["read_id"])

        except:
            continue

    result = {}
    all_barcode = set()

    for b in selected_barcode:
        all_barcode.add(b)
        all_barcode.add(rev(b))

    for b in tqdm(all_barcode):
        result[b] = trie.query(b)

    full_result = {}
    for key in tqdm(result):
        for index in result[key]:
            for read_id in result[key][index]:
                if read_id not in full_result:
                    full_result[read_id] = {"0": set(), "1": set(), "2": set()}

                full_result[read_id][str(index)].add(key)

     return full_result

def main():
    args = parse_args()
    print(args)

    selected_barcode = read_short_reads(args.short_read_barcodes)
    long_reads = read_long_reads(args.long_read_segments)

    result = run_get_matches(selected_barcode, long_reads, args.max_error, args.barcode_length)


if __name__ == "__main__":
    main()
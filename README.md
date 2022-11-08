# scTagger paper branch
scTagger matches cellular barcodes from 10x Genomics single cell experiments with hybrid short- and long-read RNA sequencing. This branch reproduces the results presented in the paper manuscript of scTagger.

# Installation

To install the needed dependencies, use conda (or mamba for a faster environment resolution) to create and activate the needed environment:

```bash
mamba env create -f envs/run.yml
conda activate run_scTagger
```

# Scripts

scTagger has three main scripts:

- `extract_sr_br.py`: Extracts the cellular barcodes of short-reads from the 10x Genomics Cell Ranger BAM alignment file of short-reads.
- `extract_top_sr_br.py`: Extracts the top most frequent short-read barcodes 
- `extract_lr_br.py`: Extracts the cellular barcode segments from the long-reads using the alignment of the short-read sequencing adapter.
- `match_lr_bc-trie.py`: Matches the short-read barcodes to the long-read segments using a trie data structure

Additionally, scTagger has an alternative, brute-force implementation of the matching in:

- `match_lr_bc-aln.py`



# Reproducing scTagger results

## Available datasets

The real and simulated datasets used in the paper are accessible at our Figshare project page:

https://figshare.com/projects/scTagger_data/138886

While the simulated datasets are provided as-is, the long-read sequences of the real datasets had to be *trimmed* to remove any sequences after the short-read adapter. This allows users to run the matching stages of scTagger (and of FLAMES, the other tool we compare against) without exposing any private patient data. 

In total, there are 6 samples available to download:

- `N_trim`, `NOA1_trim`, and `NOA2_trim`: Trimmed versions of the `N`, `NOA1` and `NOA2` real datasets.
- `S`, `M`, and `L`: the simulated small, medium and large datasets used in the manuscript.

Note, that in addition to these datasets, the user can generate new simulated datasets using a real dataset if they have their one. 

## Running the benchmarking

You can run the benchmark using Snakemake:

`snakemake -j <threads>`

However, in order to run the simulation datasets, Snakemake need to run Cell Ranger v3 which is not installed in the conda environment. Please ensure that Cell Ranger is installed and `cellranger` is in the environment `PATH` variable. Additionally, edit `cellranger_ref` variable in the `config.yaml` file to point to the Cell Ranger genome reference index. 

Note, that you can avoid running the simulation datasets by commenting their lines in the `config.yaml` file.

### Benchmarks

- `analysis/benchmark-out/<sample>/gtime.tsv`: GNU Time captured CPU/real time and max memory use for different processes
- `analysis/evaluation/<real sample>.all_vs_all.stats`: Comparison of FLAMES and scTagger matching results vs. the brute-force method
- `analysis/evaluation/<sim sample>.<tool>.stats`: The accuracy statistics of the tool based the ground truth of the simulated dataset.

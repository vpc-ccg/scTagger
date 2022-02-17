# scTagger
scTagger matches cellular barcodes from 10x Genomics single cell experiments with hybrid short- and long-read RNA sequencing.

# Installation

To install the dependency, use conda (or mamba for faster environment resolution) to create and activate the needed environment:

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



# Using scTagger

## TLDR

To run scTagger scripts, you need two main input file:

- `lr.fq`: FASTQ file of the untrimmed long-reads. This can be gzipped.
- `sr.bam`: A BAM file output by 10x Cell Ranger. You should be able to find this under the output directory of 10x Cell Ranger with the path: `outs/possorted_genome_bam.bam`

With these two files at hand, you can run the following scripts:

```
py/extract_sr_br.py \
    -i sr.bam \
    -o sr_bc.tsv.gz \
    -t <threads>

py/extract_top_sr_br.py \
    -i sr_bc.tsv.gz      
    -o sr_bc.TOP.tsv.gz \

py/extract_lr_br.py \
    -r lr.fq           
    -o lr_bc.tsv.gz \

py/match_lr_bc-trie.py \
    -lr lr_bc.tsv.gz \
    -sr sr_bc.TOP.tsv.gz \
    -o lr_bc_matches.tsv.gz \
    -t <threads>
```

## Snakemake

The advantage of running the Snakemake is that you can use it also run Cell Ranger as well all scTagger's steps. To use Snakemake, you need to first modify the `config.yaml` file. Use one of the existing samples as a template. Note that Cell Ranger is very finicky about the file path of the short-reads. 

Once all paths are updated in `config.yaml`, you can run the Snakemake:

```bash
snakemake --use-conda -j <threads>
```

You can also use the Snakemake to generate simulation data based on your real data:

```bash
snakemake -s Snakefile-simulate --use-conda -j <threads>
```

Make sure to uncomment the `sim_samples` samples list. This list lists the real samples that are going to be used as basis for gene isoform expression in Minnow single cell simulator. For each sample of the `sim_samples`, one simulated dataset will be generated using the simulation parameters listed by the `gene_cell_seed_sr_lr` list. 

Note that the Snakemake will also run the other tool we compare against, `FLAMES`.

## Script parameters

To view the list of parameters and their descriptions, run any of the scripts followed by `--help`.
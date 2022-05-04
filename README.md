# scTagger
scTagger matches barcodes of short- and long-reads of single-cell RNA-seq experiments to achieve the information of both datasets. 

## Installation

### Running with Snakemake
Snakemake provides access to the entire scTagger pipeline. You can check the pipeline at the Snakefile. Make sure to add your file paths to the config.yaml file. You can also add a new sample to run scTagger in various samples. After editing config.yaml, you can run Snakemake with your specific settings.

### Running manually
scTagger has a single python file in the py directory containing different functions to match long-reads and short-reads barcodes. 

The hole pipeline contains three steps that you can run each part separately:

#### Extract long-reads segment
The first step of the scTagger pipeline is to extract a segment where the probability of seeing a barcode is more than in other places.  To run this step, you can use the following command. 

```
py/scTagger.py extract_lr_bc -r "path/to/long/read/fastq" -o "path/to/output/file" -p "path/to/output/plots"
```

**Augments**

* `-r`: Space separated paths to reads in FASTQ
* `-g`: Space separated of the ranges of where SR adapter should be found on the LR's (Optional, Default: Detect from data)
* `-z`: Indicate input is gzipped (Optional, Default: Assume input is gzipped if it ends with \".gz\")
* `-t`: Number of threads (Optional, Default: 1)
* `-sa`: Short-read adapter (Optional, Default: "CTACACGACGCTCTTCCGATCT")
* `--num-bp-afte`: Number of bases after the end of the SR adapter alignment to generate (Optional, Default: 20)
* `-o`: Path to output file
* `-p`: Path to plot file (Optional, Default: No plotting)

**Inputs**
* A list of fastQ files of long reads

**Outputs**
* A Tsv file: 
  * First column is read-id 
  * Second column is the best edit distance with the short-read adapter
  * Third column is the starting point of long-read that matches with the adapter
  * Fourth column is the long-read segment that find. 
* A plot of optimal alignment locations of the short read adapter to the long reads. 

#### Extract short-reads barcodes

The second step is to extract the top short-reads barcodes that cover most of the reads.

```
py/scTagger.py extract_top_sr_bc -i "path/to/bam/file" -o "path/to/output/file" -p "path/to/output/plot"
```

**Arguments**
* `-i`: Input file
* `-o`: Path to output file.
* `-p`: Path to plot file (Optional, Default: No plotting)
* `--thresh`: Percentage theshold required per step to continue adding read barcodes (Optional, Default: 0.005)
* `--step-size`: Number of barcodes processed at a time and whose sum is used to check against the theshold (Optional, Default: 1000)
* `--max-barcode-cnt`: Max number of barcodes to keep (Optional, Default: 25000)

**Input**
* A bam file of short reads data

**Output**
* A TSV file
  * First column is barcodes
  * Second column is the number of appearances of the barcode
* A cumulative plot of SR coverage with batches of 1,000 barcodes 

#### Match long-reads segment with short-reads barcode
The last step is to match long read segments with selected barcodes from short reads
```
py/scTagger.py match_trie -lr "path/to/output/extract/long-read/segment" -sr "path/to/output/extract/top/short-read" -o "path/to/output/file" -t "number of threads"
```

**Arguments**
* `-lr`: Long-read segments TSV file
* `-sr`: Short-read barcode list TSV file
* `-mr`: Maximum number of errors allowed for barcode matching (Optional, Default: 2)
* `-m`: Maximum number of GB of RAM to be used (Optional, Default: 16.0)
* `-bl`: Length of barcodes (Optional, Default: 16)
* `-t`: Number of threads to use for searching (Optional, Default: 16)
* `-p`: Path of plot file
* `-o`: Path to output file. Output file is gzipped


**Inputs**
* Use the output of extracting long read segment and selecting top barcodes part as the inputs of this section 

**Outputs**
* A TSV file
  *  First column is the read id
  *  Second column is the minimum edit distance
  *  Third column is the number of short reads barcodes that match with the long read
  *  Fourth column is the long-read segment, and the Fifth column is a list of all short-read barcodes with minimum edit distance 
* A bar plot that shows the number of long reads by the minimum edit distance of their match barcode



Note: This is an active development branch.
Please check the [paper](https://github.com/vpc-ccg/scTagger/tree/paper) branch of this repository for the archived paper experiements and implementation. 


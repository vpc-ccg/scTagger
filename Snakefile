from collections import Counter
import gzip
import os

configfile: 'config.yaml'

outpath = config['outpath'].rstrip('/')
clrg_d = f'{os.path.abspath(outpath)}/cellranger-out'

rule all:
    input:
        expand(f'{clrg_d}/{{sample}}/{{sample}}/outs/possorted_genome_bam.bam', sample=config['samples']),
        expand(f'{outpath}/{{sample}}/{{sample}}.sr_bc.tsv.gz', sample=config['samples']),
        expand(f'{outpath}/{{sample}}/{{sample}}.lr_bc.tsv.gz', sample=config['samples']),
        expand(f'{outpath}/{{sample}}/{{sample}}.lr_bc_matches.tsv.gz', sample=config['samples']),

rule cellranger_count:
    input:
        ref = lambda wildcards: os.path.abspath(config['references'][config['samples'][wildcards.sample]['ref']]['cellranger_ref']),
        i1 = lambda wildcards: config['samples'][wildcards.sample]['sr']['I1'],
        r1 = lambda wildcards: config['samples'][wildcards.sample]['sr']['R1'],
        r2 = lambda wildcards: config['samples'][wildcards.sample]['sr']['R2'],
    output:
        bam = f'{clrg_d}/{{sample}}/{{sample}}/outs/possorted_genome_bam.bam',
        matrix = directory(f'{clrg_d}/{{sample}}/{{sample}}/outs/raw_feature_bc_matrix'),
    params:
        outdir = f'{clrg_d}/{{sample}}',
        sr_dir = lambda wildcards: os.path.abspath(config['samples'][wildcards.sample]['sr']['dir']),
        sr_prefix = lambda wildcards: config['samples'][wildcards.sample]['sr']['prefix'],
        mem_gb = 512
    threads:
        32
    resources:
        mem  = '512G',
        time = 60*12-1,
    shell:
        'rm -r {params.outdir} && mkdir -p {params.outdir} && cd {params.outdir} && '
        'cellranger count '
        ' --id={wildcards.sample}'
        ' --chemistry=SC3Pv3'
        ' --transcriptome={input.ref}'
        ' --fastq={params.sr_dir}'
        ' --sample={params.sr_prefix}'
        ' --localcores {threads}'        
        ' --localmem {params.mem_gb}'

rule extract_sr_bc:
    input:
        bam = f'{clrg_d}/{{sample}}/{{sample}}/outs/possorted_genome_bam.bam',
    output:
        tsv = f'{outpath}/{{sample}}/{{sample}}.sr_bc.tsv.gz',
    params:
        script = config['exec']['scTagger'],
    resources:
        mem  = "1G",
        time = 59,
    shell:
        '{params.script} extract_sr_bc -i {input.bam} -o {output.tsv}'

rule extract_lr_bc:
    input:
        reads = lambda wildcards: config['samples'][wildcards.sample]['lr_fastqs'],
    wildcard_constraints:
        sample = '|'.join([re.escape(s) for s in config['samples']]+['^$'])
    output:
        tsv = f'{outpath}/{{sample}}/{{sample}}.lr_bc.tsv.gz',
    params:
        script = config['exec']['scTagger'],
    threads:
        32
    resources:
        mem  = "256G",
        time = 59,
    shell:
        '{params.script} extract_lr_bc -r {input.reads} -o {output.tsv} -t {threads}'

rule match_trie:
    input:
        lr_tsv = f'{outpath}/{{sample}}/{{sample}}.lr_bc.tsv.gz',
        sr_tsv = f'{outpath}/{{sample}}/{{sample}}.sr_bc.tsv.gz',
    output:
        lr_tsv = f'{outpath}/{{sample}}/{{sample}}.lr_bc_matches.tsv.gz',
    params:
        script = config['exec']['scTagger'],
    threads:
        32
    resources:
        mem  = "64G",
        time = 60*5-1,
    shell:
        '{params.script} match_trie -lr {input.lr_tsv} -sr {input.sr_tsv} -o {output.lr_tsv} -t {threads}'
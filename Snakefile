from collections import Counter
import gzip
import os

import pysam
from tqdm import tqdm

configfile: 'config.yaml'

for k in config.keys():
    if not k.startswith('override_'):
        continue
    keys = k[len('override_'):].split('_')
    top_dict = eval('config{}'.format(''.join(['["{}"]'.format(x) for x in keys[:-1]])))
    assert keys[-1] in top_dict
    top_dict[keys[-1]]=config[k]

outpath = config['outpath'].rstrip('/')
bnch_d = f'{outpath}/benchmark-out'
clrg_d = f'{os.path.abspath(outpath)}/cellranger-out'
extr_d = f'{outpath}/extract-out'
fred_d = f'{outpath}/freddie-out'
seurat_d = f'{outpath}/seurat-out'
flames_d = f'{outpath}/flames-out'
simu_d = f'{outpath}/simulation'

for sample in list(config['samples'].keys()):
    for gcsSL in config['gene_cell_seed_sr_lr']:
        k = f'sim_{sample}_{gcsSL}'
        config['samples'][k] = dict(
            ref = config['samples'][sample]['ref'],
            reads = os.path.abspath(f'{simu_d}/reads/{sample}/{gcsSL}/long-reads.fastq'),
            sr_fastq_dir = os.path.abspath(f'{simu_d}/reads/{sample}/{gcsSL}'),
            sr_fastq_prefix = 'hg_100',
            sr_reads = dict(
                R1 = os.path.abspath(f'{simu_d}/reads/{sample}/{gcsSL}/hg_100_S1_L001_R1_001.fastq.gz'),
                R2 = os.path.abspath(f'{simu_d}/reads/{sample}/{gcsSL}/hg_100_S1_L001_R2_001.fastq.gz'),
            ),
            cell_count = int(gcsSL.split('-')[1]),
            seurat_config = config['samples'][sample]['seurat_config'],
        )
    del config['samples'][sample]


rule all:
    input:
        expand('{}/{{sample}}/{{sample}}/outs/possorted_genome_bam.bam'.format(clrg_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.lr_bc_matches.TRIE.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.lr_bc_matches.Validation.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.lr_bc_matches.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.sr_bc.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.lr_bc.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.sorted.bam'.format(fred_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.match.csv'.format(flames_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.match.fastq.gz'.format(flames_d), sample=config['samples']),
        # expand('{}/{{sample}}/freddie.split'.format(fred_d),         sample=config['samples']),
        # expand('{}/{{sample}}/freddie.segment'.format(fred_d),       sample=config['samples']),
        # expand('{}/{{sample}}/freddie.cluster'.format(fred_d),       sample=config['samples']),
        # expand('{}/{{sample}}/freddie.isoforms.gtf'.format(fred_d),  sample=config['samples']),
        # expand('{}/{{sample}}/'.format(seurat_d), sample=config['samples'])

rule extract_lr_br:
    input:
        script = config['exec']['split'],
        reads = lambda wildcards: config['samples'][wildcards.sample]['reads'],
    benchmark:
        '{}/extract_lr_br/{{sample}}.txt'.format(bnch_d)
    wildcard_constraints:
        sample = '|'.join([re.escape(s) for s in config['samples']]+['^$'])
    output:
        tsv = protected('{}/{{sample}}/{{sample}}.lr_bc.tsv.gz'.format(extr_d)),
        plot = protected('{}/{{sample}}/{{sample}}.lr_sa_distance.jpg'.format(extr_d)),
    threads:
        16
    resources:
        mem  = "256G",
        time = 59,
    shell:
        '{input.script} -r {input.reads} -o {output.tsv} -t {threads} -p {output.plot}'

rule minimap2:
    input:
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        genome = lambda wildcards: config['references'][config['samples'][wildcards.sample]['ref']]['genome'],
    output:
        bam=protected('{}/{{sample}}/{{sample}}.sorted.bam'.format(fred_d)),
        bai=protected('{}/{{sample}}/{{sample}}.sorted.bam.bai'.format(fred_d)),
    conda:
        'extern/freddie/envs/minimap2.yml'
    threads:
        32
    resources:
        mem  = "128G",
        time = 1439,
    wildcard_constraints:
        sample = '|'.join([re.escape(s) for s in config['samples']]+['^$'])
    shell:
        'minimap2 -a -x splice -t {threads} {input.genome} {input.reads} | '
        '  samtools sort -T {output.bam}.tmp -m 2G -@ {threads} -O bam - > {output.bam} && '
        '  samtools index {output.bam} '

rule freddie_split:
    input:
        script = config['exec']['freddie']['split'],
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        bam = '{}/{{sample}}/{{sample}}.sorted.bam'.format(fred_d),
    output:
        split = directory('{}/{{sample}}/freddie.split'.format(fred_d)),
    conda:
        'extern/freddie/envs/freddie.yml'
    threads:
        32
    resources:
        mem  = "16G",
        time = 359,
    shell:
        '{input.script} -b {input.bam} -r {input.reads} -o {output.split} -t {threads}'

rule freddie_segment:
    input:
        script = config['exec']['freddie']['segment'],
        split  = '{}/{{sample}}/freddie.split'.format(fred_d),
    output:
        segment = directory('{}/{{sample}}/freddie.segment'.format(fred_d)),
    conda:
        'extern/freddie/envs/freddie.yml'
    threads:
        32
    resources:
        mem  = "32G",
        time = 599,
    shell:
        '{input.script} -s {input.split} -o {output.segment} -t {threads}'

rule freddie_cluster:
    input:
        license = config['gurobi']['license'],
        script  = config['exec']['freddie']['cluster'],
        segment = '{}/{{sample}}/freddie.segment'.format(fred_d),
    output:
        cluster = directory('{}/{{sample}}/freddie.cluster'.format(fred_d)),
    conda:
        'extern/freddie/envs/freddie.yml'
    params:
        logs    = directory('{}/{{sample}}/freddie.cluster_logs'.format(fred_d)),
        log     = '{}/{{sample}}/freddie.cluster.log'.format(fred_d),
        timeout = config['gurobi']['timeout'],
    conda:
        'extern/freddie/envs/freddie.yml'
    threads:
        32
    resources:
        mem  = "32G",
        time = 999,
    shell:
        'export GRB_LICENSE_FILE={input.license}; '
        '{input.script} -s {input.segment} -o {output.cluster} -l {output.logs} -t {threads} -to {params.timeout} > {output.log}'

rule freddie_isoforms:
    input:
        script  = config['exec']['freddie']['isoforms'],
        split   = '{}/{{sample}}/freddie.split'.format(fred_d),
        cluster = '{}/{{sample}}/freddie.cluster'.format(fred_d),
    output:
        isoforms = protected('{}/{{sample}}/freddie.isoforms.gtf'.format(fred_d)),
    conda:
        'extern/freddie/envs/freddie.yml'
    threads:
        8
    resources:
        mem  = "16G",
        time = 359,
    shell:
        '{input.script} -s {input.split} -c {input.cluster} -o {output.isoforms} -t {threads}'


rule extract_sr_br:
    input:
        bam = '{}/{{sample}}/{{sample}}/outs/possorted_genome_bam.bam'.format(clrg_d),
    benchmark:
        '{}/extract_sr_br/{{sample}}.txt'.format(bnch_d)
    output:
        tsv = protected('{}/{{sample}}/{{sample}}.sr_bc.tsv.gz'.format(extr_d)),
    threads:
        8
    resources:
        mem  = "1G",
        time = 59,
    run:
        outfile = gzip.open(output.tsv, 'wt')
        for aln in tqdm(pysam.AlignmentFile(input.bam, 'rb', threads=threads)):
            if aln.flag > 256:
                continue
            tags = dict(aln.tags)
            C = tags.get('CB', 'NA').split('-')[0]
            U = tags.get('UB', 'NA').split('-')[0]
            qname = aln.query_name
            outfile.write(f'{qname}\t{C}\t{U}\n')
        outfile.close()

rule get_top_sr_br:
    input:
        script = config['exec']['select'],
        tsv = '{}/{{sample}}/{{sample}}.sr_bc.tsv.gz'.format(extr_d),
    output:
        tsv = protected('{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d)),
        # plot = protected('{}/{{sample}}/{{sample}}.sr_select_barcode.jpg'.format(extr_d))
    benchmark:
        '{}/get_top_sr_br/{{sample}}.txt'.format(bnch_d)
    resources:
        mem  = "8G",
        time = 59,
    shell:
        '{input.script} -i {input.tsv} -o {output.tsv}'

rule match_aln:
    input:
        script  = config['exec']['match_aln'],
        lr_tsv = '{}/{{sample}}/{{sample}}.lr_bc.tsv.gz'.format(extr_d),
        sr_tsv = '{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d),
    benchmark:
        '{}/match_aln/{{sample}}.txt'.format(bnch_d)
    output:
        lr_tsv = protected('{}/{{sample}}/{{sample}}.lr_bc_matches.tsv.gz'.format(extr_d)),
    threads:
        32
    resources:
        mem  = "250G",
        time = 60*6-1,
    shell:
        '{input.script} -lr {input.lr_tsv} -sr {input.sr_tsv} -t {threads} -o {output.lr_tsv}'

rule match_trie:
    input:
        script  = config['exec']['match_trie'],
        lr_tsv = '{}/{{sample}}/{{sample}}.lr_bc.tsv.gz'.format(extr_d),
        sr_tsv = '{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d),
    benchmark:
        '{}/match_trie/{{sample}}.txt'.format(bnch_d)
    output:
        lr_tsv = protected('{}/{{sample}}/{{sample}}.lr_bc_matches.TRIE.tsv.gz'.format(extr_d)),
        # plot = protected('{}/{{sample}}/{{sample}}.lr_bc_match_distance.jpg'.format(extr_d)),
    threads:
        1
    resources:
        mem  = "128G",
        time = 60*5-1,
    shell:
        '{input.script} -lr {input.lr_tsv} -sr {input.sr_tsv} -o {output.lr_tsv} '#-p {output.plot}'


rule validate_trie:
    input:
        lr_aln_tsv = '{}/{{sample}}/{{sample}}.lr_bc_matches.tsv.gz'.format(extr_d),
        lr_trie_tsv = '{}/{{sample}}/{{sample}}.lr_bc_matches.TRIE.tsv.gz'.format(extr_d),
    output:
        check = '{}/{{sample}}/{{sample}}.lr_bc_matches.Validation.tsv.gz'.format(extr_d),
    threads:
        1
    resources:
        mem  = "250G",
        time = 60*6-1,
    run:
        lr = dict()
        for l in tqdm(gzip.open(input.lr_aln_tsv, 'rt')):
            l = l.rstrip('\n').split('\t')
            rid = l[0]
            aln_matches = tuple(sorted(l[4].split(',')))
            assert not rid in lr, rid
            lr[rid] = [aln_matches]
        for l in tqdm(gzip.open(input.lr_trie_tsv, 'rt')):
            l = l.rstrip('\n').split('\t')
            rid = l[0]
            for trie_matches in l[1:]:
                if trie_matches=='.':
                    continue
                assert len(lr[rid])==1, rid
                trie_matches = tuple(sorted(trie_matches.split(',')))
                lr[rid].append(trie_matches)
                break             
        for am,tm in lr.values():
            assert am==tm
        
rule cell_ranger_count:
    input:
        sr_fastq_dir = lambda wildcards: config['samples'][wildcards.sample]['sr_fastq_dir'],
        ref = lambda wildcards: config['references'][config['samples'][wildcards.sample]['ref']]['cellranger_ref'],
    output:
        bam = '{}/{{sample}}/{{sample}}/outs/possorted_genome_bam.bam'.format(clrg_d),
        matrix = directory('{}/{{sample}}/{{sample}}/outs/raw_feature_bc_matrix'.format(clrg_d)),
    params:
        outdir = '{}/{{sample}}'.format(clrg_d),
        sample_prefix = lambda wildcards: config['samples'][wildcards.sample]['sr_fastq_prefix'],
        mem_gb = 512
    threads:
        32
    resources:
        mem  = '512G',
        time = 60*12-1,
    benchmark:
        '{}/cell_ranger_count/{{sample}}.txt'.format(bnch_d)
    shell:
        'rm -r {params.outdir} && mkdir -p {params.outdir} && cd {params.outdir} && '
        '  cellranger count '
        ' --id={wildcards.sample}'
        ' --chemistry=SC3Pv3'
        ' --transcriptome={input.ref}'
        ' --fastq={input.sr_fastq_dir}'
        ' --sample={params.sample_prefix}'
        ' --localcores {threads}'        
        ' --localmem {params.mem_gb}'

rule seurat:
    input:
        config_r = lambda wildcards: config['samples'][wildcards.sample]['seurat_config'],
        script = config['exec']['seurat']['preprocessing'],
        input_path = '{}/{{sample}}/{{sample}}/outs/raw_feature_bc_matrix'.format(clrg_d),
    output:
        out_path = protected('{}/{{sample}}/'.format(seurat_d)),
    benchmark:
        '{}/seurat/{{sample}}.txt'.format(bnch_d)
    shell:
        'Rscript {input.script} {input.config_r} {input.input_path} {output.out_path}/'

rule flames:
    input:
        script = 'extern/FLAMES/src/bin/match_cell_barcode',
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        whitelist = '{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d),
    output:
        csv   = protected('{}/{{sample}}/{{sample}}.match.csv'.format(flames_d)),
        fastq = protected('{}/{{sample}}/{{sample}}.match.fastq.gz'.format(flames_d)),
    params:
        tmp_in_dir = '{}/{{sample}}_fastq'.format(flames_d),
        edit_dist = 4,
    conda:
        'envs/flames.yml'
    shell:
        'rm -rf {params.tmp_in_dir}  && \n'
        'mkdir -p {params.tmp_in_dir}  && \n'
        'ln -s {input.reads} {params.tmp_in_dir}/  && \n'
        '{input.script} {params.tmp_in_dir}/ {output.csv} {output.fastq} <(less {input.whitelist}) {params.edit_dist}'

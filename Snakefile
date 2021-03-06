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
bnch_d = f'{os.path.abspath(outpath)}/benchmark-out'
clrg_d = f'{os.path.abspath(outpath)}/cellranger-out'
extr_d = f'{outpath}/extract-out'
fred_d = f'{outpath}/freddie-out'
seurat_d = f'{outpath}/seurat-out'
flames_d = f'{outpath}/flames-out'
simu_d = f'{outpath}/simulation'
eval_d = f'{outpath}/evaluation'

# for sample in list(config['sim_samples']):
#     for gcsSL in config['gene_cell_seed_sr_lr']:
#         k = f'sim_{sample}_{gcsSL}'
#         config['samples'][k] = dict(
#             ref = config['samples'][sample]['ref'],
#             reads = os.path.abspath(f'{simu_d}/reads/{sample}/{gcsSL}/lr.fq'),
#             sr_fastq_dir = os.path.abspath(f'{simu_d}/reads/{sample}/{gcsSL}'),
#             sr_fastq_prefix = 'hg_100',
#             sr_reads = dict(
#                 R1 = os.path.abspath(f'{simu_d}/reads/{sample}/{gcsSL}/hg_100_S1_L001_R1_001.fastq.gz'),
#                 R2 = os.path.abspath(f'{simu_d}/reads/{sample}/{gcsSL}/hg_100_S1_L001_R2_001.fastq.gz'),
#             ),
#             cell_count = int(gcsSL.split('-')[1]),
#             seurat_config = config['samples'][sample]['seurat_config'],
#         )

rev_compl_l = [chr(i) for i in range(128)]
rev_compl_l[ord('A')] = 'T'
rev_compl_l[ord('C')] = 'G'
rev_compl_l[ord('G')] = 'C'
rev_compl_l[ord('T')] = 'A'

def rev_compl(s):
    return ''.join(rev_compl_l[ord(c)] for c in reversed(s))

rule all:
    input:
        expand('{}/{{sample}}/{{sample}}/outs/possorted_genome_bam.bam'.format(clrg_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.lr_bc_matches.TRIE.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.lr_bc_matches.Validation.txt'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.lr_bc_matches.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.sr_bc.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.lr_bc.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.sorted.bam'.format(fred_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.match.csv'.format(flames_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.match.fastq.gz'.format(flames_d), sample=config['samples']),
        expand('{}/{{sample}}.{{tool}}.stats'.format(eval_d), sample=[s for s in config['samples'] if s.startswith('sim_')], tool=[
            'flames',
            'alns',
            'trie',
        ]),
        expand('{}/{{sample}}.{{tool}}.tsv.gz'.format(eval_d), sample=[s for s in config['samples'] if s.startswith('sim_')], tool=[
            'flames',
            'alns',
            'trie',
        ]),
        expand('{}/{{sample}}/freddie.split'.format(fred_d),         sample=config['samples']),
        expand('{}/{{sample}}/freddie.segment'.format(fred_d),       sample=config['samples']),
        expand('{}/{{sample}}/freddie.cluster'.format(fred_d),       sample=config['samples']),
        expand('{}/{{sample}}/freddie.isoforms.gtf'.format(fred_d),  sample=config['samples']),
        # expand('{}/{{sample}}/'.format(seurat_d), sample=config['samples'])

rule make_time:
    output:
        time = '{}/{{sample}}/gtime.tsv'.format(bnch_d)
    run:
        outfile = open(output.time, 'w+')
        record = list()
        record.append('method')
        record.append('real')
        record.append('user')
        record.append('memory')
        record.append('threads')
        outfile.write('\t'.join(record))
        outfile.write('\n')
        outfile.close()

rule extract_lr_br:
    input:
        reads = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        time = ancient('{}/{{sample}}/gtime.tsv'.format(bnch_d))
    wildcard_constraints:
        sample = '|'.join([re.escape(s) for s in config['samples']]+['^$'])
    output:
        tsv = protected('{}/{{sample}}/{{sample}}.lr_bc.tsv.gz'.format(extr_d)),
        plot = protected('{}/{{sample}}/{{sample}}.lr_sa_distance.pdf'.format(extr_d)),
    threads:
        32
    resources:
        mem  = "256G",
        time = 59,
    params:
        script = config['exec']['split'],
    shell:
        'GNU_TIME=$(which time) && $GNU_TIME -f "{rule}\\t%e\\t%U\\t%M\\t{threads}" -a -o {input.time} '
        '{params.script} -r {input.reads} -o {output.tsv} -t {threads} -p {output.plot}'

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
        '{input.script} -s {input.segment} -o {output.cluster} -l {params.logs} -t {threads} -to {params.timeout} > {params.log}'

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

rule filter_bc:
    input:
        script = config['exec']['select'],
        tsv = '{}/{{sample}}/{{sample}}.sr_bc.tsv.gz'.format(extr_d),
        time = ancient('{}/{{sample}}/gtime.tsv'.format(bnch_d))
    output:
        tsv = protected('{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d)),
        plot = protected('{}/{{sample}}/{{sample}}.sr_bc.TOP.pdf'.format(extr_d)),
    resources:
        mem  = "8G",
        time = 59,
    shell:
        'GNU_TIME=$(which time) && $GNU_TIME -f "{rule}\\t%e\\t%U\\t%M\\t{threads}" -a -o {input.time} '
        '{input.script} -i {input.tsv} -o {output.tsv} -p {output.plot}'

rule match_aln:
    input:
        script  = config['exec']['match_aln'],
        lr_tsv = '{}/{{sample}}/{{sample}}.lr_bc.tsv.gz'.format(extr_d),
        sr_tsv = '{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d),
        time = ancient('{}/{{sample}}/gtime.tsv'.format(bnch_d))
    output:
        lr_tsv = protected('{}/{{sample}}/{{sample}}.lr_bc_matches.tsv.gz'.format(extr_d)),
    threads:
        32
    resources:
        mem  = "250G",
        time = 60*6-1,
    shell:
        'GNU_TIME=$(which time) && $GNU_TIME -f "{rule}\\t%e\\t%U\\t%M\\t{threads}" -a -o {input.time} '
        '{input.script} -lr {input.lr_tsv} -sr {input.sr_tsv} -t {threads} -o {output.lr_tsv}'

rule match_trie:
    input:
        script  = config['exec']['match_trie'],
        lr_tsv = '{}/{{sample}}/{{sample}}.lr_bc.tsv.gz'.format(extr_d),
        sr_tsv = '{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d),
        time = ancient('{}/{{sample}}/gtime.tsv'.format(bnch_d))
    output:
        lr_tsv = protected('{}/{{sample}}/{{sample}}.lr_bc_matches.TRIE.tsv.gz'.format(extr_d)),
        # plot = protected('{}/{{sample}}/{{sample}}.lr_bc_matches.TRIE.pdf'.format(extr_d)),
    threads:
        32
    resources:
        mem  = "64G",
        time = 60*5-1,
    shell:
        'GNU_TIME=$(which time) && $GNU_TIME -f "{rule}\\t%e\\t%U\\t%M\\t{threads}" -a -o {input.time} '
        '{input.script} -lr {input.lr_tsv} -sr {input.sr_tsv} -o {output.lr_tsv} -t {threads}'#-p {output.plot} -m {params.mem}'


rule validate_trie:
    input:
        lr_aln_tsv = '{}/{{sample}}/{{sample}}.lr_bc_matches.tsv.gz'.format(extr_d),
        lr_trie_tsv = '{}/{{sample}}/{{sample}}.lr_bc_matches.TRIE.tsv.gz'.format(extr_d),
    output:
        check = '{}/{{sample}}/{{sample}}.lr_bc_matches.Validation.txt'.format(extr_d),
    threads:
        1
    resources:
        mem  = "250G",
        time = 60*6-1,
    run:
        lr = dict()
        # stats = dict(
        #     matches = 0,
        #     mismatch = 0,
        #     error_too_low = 0,
        #     skipped = 0,
        # )
        for l in tqdm(gzip.open(input.lr_aln_tsv, 'rt')):
            l = l.rstrip('\n').split('\t')
            rid = l[0]
            if l[4] !='':
                matches = tuple(sorted(l[4].split(',')))
            else:
                matches = tuple()
            assert len(matches) == int(l[2])
            assert not rid in lr, (rid,l)
            lr[rid] = dict(
                aln=((l[1],matches)),
                trie=(('inf'),tuple()),
            )
        for l in tqdm(gzip.open(input.lr_trie_tsv, 'rt')):
            l = l.rstrip('\n').split('\t')
            rid = l[0]
            if l[4] !='':
                matches = tuple(sorted(l[4].split(',')))
            else:
                matches = tuple()
            assert len(matches) == int(l[2])
            lr[rid]['trie'] = (l[1],matches)
        C = Counter(
            (X['aln']==X['trie'], X['aln'][0], len(X['aln'][1])<2, X['trie'][0], len(X['trie'][1])<2,) for X in lr.values()
        )
        outfile = open(output.check, 'w+')
        outfile.write('match\taln_e\taln_u\ttrie_e\ttrie_u\t#\t%\n')
        for k,v in sorted(C.items(), key=lambda x:x[1], reverse=True):
            k = '\t'.join(str(x) for x in k)
            outfile.write(f'{k}\t{v}\t{v/len(lr):.2%}\n')
        outfile.close()
        
rule cellranger_count:
    input:
        sr_fastq_dir = lambda wildcards: config['samples'][wildcards.sample]['sr_fastq_dir'],
        ref = lambda wildcards: config['references'][config['samples'][wildcards.sample]['ref']]['cellranger_ref'],
        time = ancient('{}/{{sample}}/gtime.tsv'.format(bnch_d))
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
    shell:
        'INPUT_TIME=$(readlink -f {input.time}) && '
        'rm -r {params.outdir} && mkdir -p {params.outdir} && cd {params.outdir} && '
        'GNU_TIME=$(which time) && $GNU_TIME -f "{rule}\\t%e\\t%U\\t%M\\t{threads}" -a -o $INPUT_TIME '
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
        time = ancient('{}/{{sample}}/gtime.tsv'.format(bnch_d))
    output:
        out_path = protected('{}/{{sample}}/'.format(seurat_d)),
    shell:
        'GNU_TIME=$(which time) && $GNU_TIME -f "{rule}\\t%e\\t%U\\t%M\\t{threads}" -a -o {input.time} '
        'Rscript {input.script} {input.config_r} {input.input_path} {output.out_path}/'

rule flames:
    input:
        script = 'extern/FLAMES/src/bin/match_cell_barcode',
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        whitelist = '{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d),
        time = ancient('{}/{{sample}}/gtime.tsv'.format(bnch_d))
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
        'GNU_TIME=$(which time) && $GNU_TIME -f "{rule}\\t%e\\t%U\\t%M\\t{threads}" -a -o {input.time} '
        '{input.script} {params.tmp_in_dir}/ {output.csv} {output.fastq} <(less {input.whitelist}) {params.edit_dist}'

rule convert_flames:
    input:
        fastq = '{}/{{sample}}/{{sample}}.match.fastq.gz'.format(flames_d),
    output:
        tsv = protected('{}/{{sample}}.flames.tsv.gz'.format(eval_d)),
    run:
        outfile = gzip.open(output.tsv, 'wt+')
        for l in gzip.open(input.fastq,'rt'):
            if l[0]!='@':
                continue
            bc,rname = l.rstrip().split()
            bc = bc.lstrip('@')
            rname = rname.split('#')[1]
            outfile.write(f'{rname}\t{bc}\n')
        outfile.close()

rule convert_trie:
    input:
        lr_tsv = '{}/{{sample}}/{{sample}}.lr_bc_matches.TRIE.tsv.gz'.format(extr_d),
    output:
        tsv = protected('{}/{{sample}}.trie.tsv.gz'.format(eval_d)),
    run:
        outfile = gzip.open(output.tsv, 'wt+')
        for l in gzip.open(input.lr_tsv,'rt'):
            l = l.rstrip('\n').split('\t')
            outfile.write(f'{l[0]}\t{l[4]}\n')
        outfile.close()

rule convert_alns:
    input:
        lr_tsv = '{}/{{sample}}/{{sample}}.lr_bc_matches.tsv.gz'.format(extr_d),
    output:
        tsv = protected('{}/{{sample}}.alns.tsv.gz'.format(eval_d)),
    run:
        outfile = gzip.open(output.tsv, 'wt+')
        for l in gzip.open(input.lr_tsv,'rt'):
            l = l.rstrip('\n').split('\t')
            outfile.write(f'{l[0]}\t{l[4]}\n')
        outfile.close()

rule evaluate:
    input:
        truth_tsv = '{}/{{sample}}.truth.tsv.gz'.format(eval_d),
        tool_tsv = '{}/{{sample}}.{{tool}}.tsv.gz'.format(eval_d),
    output:
        stats = protected('{}/{{sample}}.{{tool}}.stats'.format(eval_d)),
    wildcard_constraints:
        sample = 'sim_.+',
        tool = '|'.join([re.escape(s) for s in ['flames','trie','alns']])
    run:
        stats = dict(
            matches = 0,
            ambiguous = 0,
            mismatch = 0,
            skipped = 0,
        )
        rname_to_bc = dict()
        for l in gzip.open(input.truth_tsv, 'rt'):
            rname,bc = l.rstrip('\n').split('\t')
            assert not rname in rname_to_bc, (rname, l) 
            rname_to_bc[rname] = bc
        for l in gzip.open(input.tool_tsv, 'rt'):
            rname,bc = l.rstrip('\n').split('\t')
            assert rname in rname_to_bc, (rname,l)
            if len(bc) == 0:
                continue
            bc = bc.split(',')
            match_count = 0
            for b in bc:
                if b==rname_to_bc[rname] or rev_compl(b)==rname_to_bc[rname]:
                    match_count += 1
            if match_count == 0:
                stats['mismatch'] += 1
            elif match_count == len(bc):
                stats['matches'] += 1
            else:
                stats['ambiguous'] += 1
        stats['skipped'] = len(rname_to_bc) - stats['matches'] - stats['mismatch'] - stats['ambiguous']
        outfile = open(output.stats, 'w+')
        for k,v in stats.items():
            outfile.write(f'{k}:\t{v}\t({v/len(rname_to_bc):.2%})\n')
        outfile.close()

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


rule extract_lr_br:
    input:
        reads = lambda wildcards: config['samples'][wildcards.sample]['reads'],
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
        '{params.script} -r {input.reads} -o {output.tsv} -t {threads} -p {output.plot}'


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
    output:
        tsv = protected('{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d)),
        plot = protected('{}/{{sample}}/{{sample}}.sr_bc.TOP.pdf'.format(extr_d)),
    resources:
        mem  = "8G",
        time = 59,
    shell:
        '{input.script} -i {input.tsv} -o {output.tsv} -p {output.plot}'


rule match_trie:
    input:
        script  = config['exec']['match_trie'],
        lr_tsv = '{}/{{sample}}/{{sample}}.lr_bc.tsv.gz'.format(extr_d),
        sr_tsv = '{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d),
    output:
        lr_tsv = protected('{}/{{sample}}/{{sample}}.lr_bc_matches.TRIE.tsv.gz'.format(extr_d)),
    threads:
        32
    resources:
        mem  = "64G",
        time = 60*5-1,
    shell:
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
        'rm -r {params.outdir} && mkdir -p {params.outdir} && cd {params.outdir} && '
        'cellranger count '
        ' --id={wildcards.sample}'
        ' --chemistry=SC3Pv3'
        ' --transcriptome={input.ref}'
        ' --fastq={input.sr_fastq_dir}'
        ' --sample={params.sample_prefix}'
        ' --localcores {threads}'        
        ' --localmem {params.mem_gb}'


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
clrg_d = f'{outpath}/cellranger-out'
extr_d = f'{outpath}/extract-out'
flames_d = f'{outpath}/flames-out'
simu_d = f'{outpath}/simulation'
eval_d = f'{outpath}/evaluation'

for sample in list(config.get('sim_samples', [])):
    for gcsSL in config['gene_cell_seed_sr_lr']:
        k = f'sim_{sample}_{gcsSL}'
        config['samples'][k] = dict(
            ref = config['samples'][sample]['ref'],
            reads = os.path.abspath(f'{simu_d}/reads/{sample}/{gcsSL}/lr.fq'),
            sr_fastq_dir = os.path.abspath(f'{simu_d}/reads/{sample}/{gcsSL}'),
            sr_fastq_prefix = 'hg_100',
            sr_reads = dict(
                R1 = os.path.abspath(f'{simu_d}/reads/{sample}/{gcsSL}/hg_100_S1_L001_R1_001.fastq.gz'),
                R2 = os.path.abspath(f'{simu_d}/reads/{sample}/{gcsSL}/hg_100_S1_L001_R2_001.fastq.gz'),
            ),
            cell_count = int(gcsSL.split('-')[1]),
        )

rev_compl_l = [chr(i) for i in range(128)]
rev_compl_l[ord('A')] = 'T'
rev_compl_l[ord('C')] = 'G'
rev_compl_l[ord('G')] = 'C'
rev_compl_l[ord('T')] = 'A'

def rev_compl(s):
    return ''.join(rev_compl_l[ord(c)] for c in reversed(s))

rule all:
    input:
        expand('{}/{{sample}}.{{tool}}.stats'.format(eval_d), sample=[s for s in config['samples'] if s.startswith('sim_')], tool=[
            'flames',
            'alns',
            'trie',
        ]),
        expand('{}/{{sample}}.all_vs_all.stats'.format(eval_d), sample=[s for s in config['samples'] if not s.startswith('sim_')]),
        expand('{}/{{sample}}.{{tool}}.tsv.gz'.format(eval_d), sample=[s for s in config['samples']], tool=[
            'flames',
            'alns',
            'trie',
        ]),

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

rule download_real_trim:
    output:
        'data/real_trim/short-read-barcodes/N_trim.sr_bc.TOP.tsv.gz',
        'data/real_trim/short-read-barcodes/NOA1_trim.sr_bc.TOP.tsv.gz',
        'data/real_trim/short-read-barcodes/NOA2_trim.sr_bc.TOP.tsv.gz',
        'data/real_trim/long-read/N_trim.fq',
        'data/real_trim/long-read/NOA1_trim.fq',
        'data/real_trim/long-read/NOA2_trim.fq',
    params:
        url = 'https://figshare.com/ndownloader/files/35062999?private_link=6b618c15911ef480bd92',
    threads:
        64
    shell:
        'wget {params.url} -O data/real_trim.tar && '
        'tar -vxf data/real_trim.tar -C data/real_trim/ && '
        'pigz --decompress -p {threads} data/real_trim/long-read/* && '
        'rm data/real_trim.tar'

rule download_simulation:
    output:
        'data/simulation/short-reads/S/hg_100_S1_L001_R1_001.fastq.gz',
        'data/simulation/short-reads/S/hg_100_S1_L001_R2_001.fastq.gz',
        'data/simulation/long-reads/S_lr.fq',
        'data/simulation/short-reads/M/hg_100_S1_L001_R1_001.fastq.gz',
        'data/simulation/short-reads/M/hg_100_S1_L001_R2_001.fastq.gz',
        'data/simulation/long-reads/M_lr.fq',
        'data/simulation/short-reads/L/hg_100_S1_L001_R1_001.fastq.gz',
        'data/simulation/short-reads/L/hg_100_S1_L001_R2_001.fastq.gz',
        'data/simulation/long-reads/L_lr.fq',
    params:
        url = 'https://figshare.com/ndownloader/files/35076283?private_link=7d7078977a735fef9e57',
    threads:
        64
    shell:
        'wget {params.url} -O data/sim.tar && '
        'tar -vxf data/sim.tar -C data/simulation/ && '
        'pigz --decompress -p {threads} data/simulation/long-reads/*gz && '
        'awk \'BEGIN{{q=""; for (i=0;i<100000;i++)q=q"!"}} NR%2==1{{print "@"substr($0,2)}} NR%2==0 {{print $0; print "+";  print substr(q,1,length($0))}}\' data/simulation/long-reads/S_lr.fa > data/simulation/long-reads/S_lr.fq && '
        'awk \'BEGIN{{q=""; for (i=0;i<100000;i++)q=q"!"}} NR%2==1{{print "@"substr($0,2)}} NR%2==0 {{print $0; print "+";  print substr(q,1,length($0))}}\' data/simulation/long-reads/M_lr.fa > data/simulation/long-reads/M_lr.fq && '
        'awk \'BEGIN{{q=""; for (i=0;i<100000;i++)q=q"!"}} NR%2==1{{print "@"substr($0,2)}} NR%2==0 {{print $0; print "+";  print substr(q,1,length($0))}}\' data/simulation/long-reads/L_lr.fa > data/simulation/long-reads/L_lr.fq && '
        'rm data/simulation/long-reads/*_lr.fa data/sim.tar'

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
        script = config['exec']['extract_lr_bc'],
    shell:
        'GNU_TIME=$(which time) && $GNU_TIME -f "{rule}\\t%e\\t%U\\t%M\\t{threads}" -a -o {input.time} '
        '{params.script} -r {input.reads} -o {output.tsv} -t {threads} -p {output.plot}'

rule extract_sr_br:
    input:
        bam = '{}/{{sample}}/{{sample}}/outs/possorted_genome_bam.bam'.format(clrg_d),
    wildcard_constraints:
        sample = '|'.join([re.escape(s) for s in config['samples'] if 'sr_fastq_dir' in config['samples'][s]]+['^$'])
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

rule select_sr_bc:
    input:
        script = config['exec']['select_sr_bc'],
        tsv = '{}/{{sample}}/{{sample}}.sr_bc.tsv.gz'.format(extr_d),
        time = ancient('{}/{{sample}}/gtime.tsv'.format(bnch_d))
    wildcard_constraints:
        sample = '|'.join([re.escape(s) for s in config['samples'] if 'sr_fastq_dir' in config['samples'][s]]+['^$'])
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
        # sr_tsv = '{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d),
        sr_tsv = lambda wildcards: config['samples'][wildcards.sample].get('top_bc', '{e}/{s}/{s}.sr_bc.TOP.tsv.gz'.format(e=extr_d,s=wildcards.sample)),
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
        # sr_tsv = '{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d),
        sr_tsv = lambda wildcards: config['samples'][wildcards.sample].get('top_bc', '{e}/{s}/{s}.sr_bc.TOP.tsv.gz'.format(e=extr_d,s=wildcards.sample)),
        time = ancient('{}/{{sample}}/gtime.tsv'.format(bnch_d))
    output:
        lr_tsv = protected('{}/{{sample}}/{{sample}}.lr_bc_matches.TRIE.tsv.gz'.format(extr_d)),
    threads:
        32
    resources:
        mem  = "64G",
        time = 60*5-1,
    shell:
        'GNU_TIME=$(which time) && $GNU_TIME -f "{rule}\\t%e\\t%U\\t%M\\t{threads}" -a -o {input.time} '
        '{input.script} -lr {input.lr_tsv} -sr {input.sr_tsv} -o {output.lr_tsv} -t {threads}'

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

rule flames:
    input:
        script = 'extern/FLAMES/src/bin/match_cell_barcode',
        reads  = lambda wildcards: config['samples'][wildcards.sample]['reads'],
        # whitelist = '{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d),
        whitelist = lambda wildcards: config['samples'][wildcards.sample].get('top_bc', '{e}/{s}/{s}.sr_bc.TOP.tsv.gz'.format(e=extr_d,s=wildcards.sample)),
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
        for idx,l in tqdm(enumerate(gzip.open(input.fastq,'rt'))):
            if idx%4!=0:
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

rule evaluate_sim:
    input:
        truth_tsv = lambda wildcards: config['samples'][wildcards.sample].get('truth', '{}/{{sample}}.truth.tsv.gz'.format(eval_d)),
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

rule evaluate_real:
    input:
        flames_tsv = '{}/{{sample}}.flames.tsv.gz'.format(eval_d),
        trie_tsv   = '{}/{{sample}}.trie.tsv.gz'.format(eval_d),
        alns_matches = '{}/{{sample}}/{{sample}}.lr_bc_matches.tsv.gz'.format(extr_d),
    output:
        stats = protected('{}/{{sample}}.all_vs_all.stats'.format(eval_d)),
    wildcard_constraints:
        sample = '|'.join([re.escape(s) for s in config['samples'] if not s.startswith('sim_')]),
    run:
        read_details = dict()
        print(f'Reading {input.alns_matches}')
        for l in tqdm(gzip.open(input.alns_matches, 'tr')):
            rid,dist,count,segment,barcodes = l.rstrip('\n').split('\t')
            dist = float(dist)
            count = int(count)
            barcodes = tuple(sorted(b if b < rev_compl(b) else rev_compl(b) for b in barcodes.split(',')))
            if barcodes == ('',):
                barcodes = tuple()
            assert len(barcodes)==count,(count,barcodes)
            assert not rid in read_details, (rid,l)
            read_details[rid] = dict(
                dist=dist,
                count=count,
                sctagger=tuple(),
                flames=tuple(),
                brute=barcodes,
            )
        for tsv,n in ((input.trie_tsv,'sctagger'),(input.flames_tsv,'flames'),):
            print(f'Reading {tsv}')
            for l in tqdm(gzip.open(tsv, 'tr'), total=len(read_details)):
                rid,barcodes = l.rstrip('\n').split('\t')
                barcodes = tuple(sorted(b if b < rev_compl(b) else rev_compl(b) for b in barcodes.split(',')))
                if barcodes == ('',):
                    barcodes = tuple()
                read_details[rid][n] = barcodes            
        
        c = Counter()
        for details in tqdm(read_details.values(), total=len(read_details)):
            match = (
                'Y' if len(details['brute'])!=0 else 'N',
                'Y' if len(details['brute'])==1 else 'N',        
                'Y' if len(details['sctagger'])!=0 and details['sctagger'] == details['brute'] else 'N',        
                'Y' if len(details['flames'])!=0   and details['flames'] == details['brute'] else 'N',        
                'Y' if len(details['sctagger'])!=0 and details['flames'] == details['brute'] else 'N',        
            )
            c[match]+=1

        outfile = open(output.stats, 'w+')
        record = [
            'Brute_matched?',
            'Brute_match_unique?',
            'scTagger=Brute?',
            'FLAMES=Brute?',
            'scTagger=FLAMES?',
            '%LRs',
            '#LRs',
        ]
        outfile.write('\t'.join(record))
        outfile.write('\n')
        
        for states,count,percent in sorted([(''.join(k),f'{v/sum(c.values()):.1%}', f'{v:,}') for k,v in c.items()], reverse=True):
            record = [X for X in states] + [count,percent]
            outfile.write('\t'.join(record))
            outfile.write('\n')
        outfile.close()
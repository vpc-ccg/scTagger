import pysam
import gzip
from tqdm import tqdm
from collections import Counter

configfile: 'config.yaml'

for k in config.keys():
    if not k.startswith('override_'):
        continue
    keys = k[len('override_'):].split('_')
    top_dict = eval('config{}'.format(''.join(['["{}"]'.format(x) for x in keys[:-1]])))
    assert keys[-1] in top_dict
    top_dict[keys[-1]]=config[k]

outpath = config['outpath'].rstrip('/')
extr_d = f'{outpath}/extract-out'
fred_d = f'{outpath}/freddie-out'

rule all:
    input:
        expand('{}/{{sample}}/{{sample}}.lr_bc_matches.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.sr_bc.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.lr_bc.tsv.gz'.format(extr_d), sample=config['samples']),
        expand('{}/{{sample}}/{{sample}}.sorted.bam'.format(fred_d), sample=config['samples']),
        expand('{}/{{sample}}/freddie.split'.format(fred_d),         sample=config['samples']),
        expand('{}/{{sample}}/freddie.segment'.format(fred_d),       sample=config['samples']),
        expand('{}/{{sample}}/freddie.cluster'.format(fred_d),       sample=config['samples']),
        expand('{}/{{sample}}/freddie.isoforms.gtf'.format(fred_d),  sample=config['samples']),

rule extract_lr_br:
    input:
        script = config['exec']['split'],
        reads = lambda wildcards: config['samples'][wildcards.sample]['reads'],
    output:
        tsv = protected('{}/{{sample}}/{{sample}}.lr_bc.tsv.gz'.format(extr_d)),
    threads:
        16
    resources:
        mem  = "256G",
        time = 59,
    shell:
        '{input.script} -r {input.reads} -o {output.tsv} -t {threads}'

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

rule count:
    input:
        fastq=directory(config["samples"]),
        trans=config["references"]["homo_sapiens"]["transcriptome"],
        respath=directory(config["outpath"])
    output:
        directory("/outs/web_summary.html")
    params:
        samples=config[""]
    shell:
        "cd {input.respath} && cellranger count --id={params.samples} --transcriptome={input.trans} --fastq={input.fastq} --sample={params.samples}"

rule extract_sr_br:
    input:
        bam = lambda wildcards: config['samples'][wildcards.sample]['sr_bam'],
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
        tsv = '{}/{{sample}}/{{sample}}.sr_bc.tsv.gz'.format(extr_d),
    output:
        tsv = protected('{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d)),
    threads:
        1
    resources:
        mem  = "8G",
        time = 59,
    run:
        bc = Counter()
        for l in tqdm(gzip.open(input.tsv)):
            l = l.decode().rstrip().split('\t')
            bc[l[1]]+=1
        total = sum(bc.values())
        del bc['NA']
        bc = sorted(((c,b) for b,c in  bc.items()), reverse=True)
        outfile = gzip.open(output.tsv, 'wt')
        print(f'{len(bc)}, {total}')
        for idx in range(0,len(bc),1000):
            S = sum(c for c,b in bc[idx:idx+1000])
            print(idx, S, S/total)
            if S/total < 0.005:
                break
            for c,b in bc[idx:idx+1000]:
                outfile.write(f'{b}\t{c}\n')
            if bc[idx:idx+1000][-1][0] < 2:
                break
        outfile.close()

rule match_aln:
    input:
        script  = config['exec']['match_aln'],
        lr_tsv = '{}/{{sample}}/{{sample}}.lr_bc.tsv.gz'.format(extr_d),
        sr_tsv = '{}/{{sample}}/{{sample}}.sr_bc.TOP.tsv.gz'.format(extr_d),
    output:
        lr_tsv = protected('{}/{{sample}}/{{sample}}.lr_bc_matches.tsv.gz'.format(extr_d)),
    threads:
        32
    resources:
        mem  = "250G",
        time = 60*6-1,
    shell:
        '{input.script} -lr {input.lr_tsv} -sr {input.sr_tsv} -t {threads} -o {output.lr_tsv}'
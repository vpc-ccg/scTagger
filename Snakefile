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
        expand('{}/{{sample}}/{{sample}}.lr_bc.tsv'.format(extr_d), sample=config['samples']),
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
        tsv = protected('{}/{{sample}}/{{sample}}.lr_bc.tsv'.format(extr_d)),
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
        'envs/minimap2.yml'
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
        'envs/freddie.yml'
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
        'envs/freddie.yml'
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
        'envs/freddie.yml'
    params:
        logs    = directory('{}/{{sample}}/freddie.cluster_logs'.format(fred_d)),
        log     = '{}/{{sample}}/freddie.cluster.log'.format(fred_d),
        timeout = config['gurobi']['timeout'],
    conda:
        'envs/freddie.yml'
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
        'envs/freddie.yml'
    threads:
        8
    resources:
        mem  = "16G",
        time = 359,
    shell:
        '{input.script} -s {input.split} -c {input.cluster} -o {output.isoforms} -t {threads}'
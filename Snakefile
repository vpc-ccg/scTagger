configfile: 'config.yaml'

for k in config.keys():
    if not k.startswith('override_'):
        continue
    keys = k[len('override_'):].split('_')
    top_dict = eval('config{}'.format(''.join(['["{}"]'.format(x) for x in keys[:-1]])))
    assert keys[-1] in top_dict
    top_dict[keys[-1]]=config[k]

outpath = config['outpath'].rstrip('/')

rule all:
    input:
        expand('{}/{{sample}}/{{sample}}.lr_bc.tsv'.format(outpath), sample=config['samples']),

rule extract_lr_br:
    input:
        script = config['exec']['split'],
        reads = lambda wildcards: config['samples'][wildcards.sample]['reads'],
    output:
        tsv = protected('{}/{{sample}}/{{sample}}.lr_bc.tsv'.format(outpath)),
    threads:
        32
    resources:
        mem  = "128G",
        time = 1439,
    shell:
        '{input.script} -r {input.reads} -o {output.tsv} -t {threads}'


import glob
from SCRIPTS import tetyperead
from SCRIPTS import coverage_filter
import pandas as pd
import sys
from Bio import SeqIO


configfile: "samples/samples.yaml"
configfile: "config/config.yaml"
configfile: "config/reference.yaml"


OUTDIR = config['OUTPUT'] + "/"
SAMPLES = config['SAMPLES']
CASETTE = config['CASETTE']['ref1']
STRETCHDEL = [True, False]



def get_length_dna(location_fasta):
    # first do sanity check to confirm that only 1 fasta is present as reference
    if len([x.id for x in SeqIO.parse(location_fasta, 'fasta')]) > 1:
        sys.exit("\n Please provide a reference file with only 1 reference seq")

    #actually check length of fragment
    for i in SeqIO.parse(location_fasta, 'fasta'):
        length = len(i.seq)
    return length

onerror:
    print('an error has occured')
    sys.exit(1)

################################################################################
# rule all
################################################################################

rule all:
  input:
    expand(OUTDIR + "output/{sample}/kma_alignment/KMA.res",  sample = SAMPLES),
    "samples/good_samples.yaml"

################################################################################
# KMA building and mapping
################################################################################

rule build_database:
    input:
        CASETTE
    output:
        comp = "ref/kma_database/reference.comp.b",
        length = "ref/kma_database/reference.length.b",
        name = "ref/kma_database/reference.name",
        seq = "ref/kma_database/reference.seq.b",
    conda:
        "envs/KMA.yaml"
    params:
        name = "ref/kma_database/reference",
    shell:
        "kma index -i {input} -o {params.name}"

rule kma_alignment:
    input:
        forward = lambda wildcards: SAMPLES[wildcards.sample]['forward'],
        rev = lambda wildcards: SAMPLES[wildcards.sample]['reverse'],
        reference = "ref/kma_database/reference.comp.b"
    params:
        reference = "ref/kma_database/reference",
        output = OUTDIR + "output/{sample}/kma_alignment/KMA"
    output:
        aln = OUTDIR + "output/{sample}/kma_alignment/KMA.aln",
        frag = OUTDIR + "output/{sample}/kma_alignment/KMA.frag.gz",
        fsa = OUTDIR + "output/{sample}/kma_alignment/KMA.fsa",
        res =OUTDIR + "output/{sample}/kma_alignment/KMA.res",
    conda:
        "envs/KMA.yaml"
    shell:
        "kma -ipe {input.forward} {input.rev} -o {params.output} -t_db {params.reference}"

rule coverage_filter:
    input:
        res = expand(OUTDIR + "output/{sample}/kma_alignment/KMA.res", sample = SAMPLES),
        yaml = "samples/samples.yaml"
    output:
        yaml = "samples/good_samples.yaml"
    params:
        sample = lambda wildcards: SAMPLES[wildcards.sample]
    run:
        print(f'writing suitable samples to {output.yaml}')
        for i in input.res:
            coverage_filter.selection(i, params.sample)
        coverage_filter.filter(input.yaml, output.yaml)

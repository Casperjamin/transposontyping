import glob
from SCRIPTS import tetyperead
from SCRIPTS import coverage_filter
import pandas as pd
import sys
from Bio import SeqIO
configfile: "samples/samples.yaml"
configfile: "config/config.yaml"
configfile: "config/reference.yaml"
SAMPLES = config['SAMPLES']
CASETTE = config['CASETTE']
STRETCHDEL = [True, False]

cassette = CASETTE["ref1"]

def get_length_dna(location_fasta):
    #first do sanity check to confirm that only 1 fasta is present as reference
    if len([x.id for x in SeqIO.parse(location_fasta, 'fasta')]) > 1:
        sys.exit("\n Please provide a reference file with only 1 reference seq")

    #actually check length of fragment
    for i in SeqIO.parse(location_fasta, 'fasta'):
        length = len(i.seq)
    return length

################################################################################
# rule all
################################################################################

rule all:
  input:
    expand("output/{sample}/kma_alignment/KMA.res",  sample = SAMPLES),
    "samples/good_samples.yaml"

################################################################################
# KMA building and mapping
################################################################################

rule build_database:
    input:
        CASETTE
    output:
        comp = "ref/kma_database/reference.comp.b",
        index = "ref/kma_database/reference.index.b",
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
        reverse = lambda wildcards: SAMPLES[wildcards.sample]['reverse'],
        reference = "ref/kma_database/reference.comp.b"
    params:
        reference = "ref/kma_database/reference",
        output = "output/{sample}/kma_alignment/KMA"
    output:
        aln = "output/{sample}/kma_alignment/KMA.aln",
        frag = "output/{sample}/kma_alignment/KMA.frag.gz",
        fsa = "output/{sample}/kma_alignment/KMA.fsa",
        res ="output/{sample}/kma_alignment/KMA.res",
    conda:
        "envs/KMA.yaml"
    shell:
        "kma -ipe {input.forward} {input.reverse} -o {params.output} -t_db {params.reference}"

rule coverage_filter:
    input:
        res = expand("output/{sample}/kma_alignment/KMA.res", sample = SAMPLES),
        yaml = "samples/samples.yaml"
    output:
        yaml = "samples/good_samples.yaml"
    run:
        for i in input.res:
            coverage_filter.selection(i)
        coverage_filter.filter(input.yaml, output.yaml)

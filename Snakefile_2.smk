import glob
from SCRIPTS import tetyperead, merger, report, coverage_filter
import pandas as pd
import sys
from Bio import SeqIO


configfile: "samples/good_samples.yaml"
configfile: "config/config.yaml"
configfile: "config/reference.yaml"


SAMPLES = config['SAMPLES']
OUTDIR = config['OUTPUT'] + "/"
CASETTE = config['CASETTE']
casette = CASETTE["ref1"]


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
    expand(OUTDIR + "output/{sample}/tetype/{sample}_summary.txt", sample = SAMPLES),
    expand(OUTDIR + "output/{sample}/TEType_convert/{sample}_TEType_df.tsv", sample = SAMPLES),
    expand(OUTDIR + "output/SNPs.csv"),
    expand(OUTDIR + "output/report/pairwise_distance.csv"),
    expand(OUTDIR + "output/report/dendrogram.png"),


################################################################################
# generating figures, pairwise comparisson and bedgraph used for coverage plot
################################################################################

rule generate_heatmap_dendrogram:
    input:
        csv = OUTDIR + "output/SNPs.csv"
    output:
        OUTDIR + "output/report/pairwise_distance.csv",
        OUTDIR + "output/report/dendrogram.png"
    params:
        OUTDIR + "output/report/"
    run:
        os.system(f'mkdir -p {params}')
        report.generate_report(input, params, get_length_dna(casette))



################################################################################
# TETyper, conversion and merging
################################################################################

rule TETyping:
  input:
    forward = lambda wildcards: SAMPLES[wildcards.sample]['forward'],
    rev = lambda wildcards: SAMPLES[wildcards.sample]['reverse']
  output:
    OUTDIR + "output/{sample}/tetype/{sample}_summary.txt",
    OUTDIR + "output/{sample}/tetype/{sample}.bam",
    OUTDIR + "output/{sample}/tetype/{sample}.vcf"
  threads:
    4
  conda:
    "envs/TEType_environment.yml"
  log:
    "logs/tetype/{sample}.txt"
  params:
    OUTDIR + "output/{sample}/tetype/{sample}"
  shell:
    " TETyper.py --outprefix {params}"
    " --ref {casette}"
    " --fq1 {input.forward} "
    " --fq2 {input.rev} "
    " --flank_len 10 "
    " --threads {threads} "
    "2> {log} "

rule convert:
  input:
    loc = OUTDIR + "output/{sample}/tetype/{sample}_summary.txt",
  output:
    tsv = OUTDIR + "output/{sample}/TEType_convert/{sample}_TEType_df.tsv"
  run:
    summary = tetyperead.tetype_result(summary=input.loc, fragment=casette)
    df = summary.result_to_frame()
    df.to_csv(output.tsv, sep = "\t")

rule merge:
    input:
        samples = expand(OUTDIR + "output/{sample}/TEType_convert/{sample}_TEType_df.tsv", sample = SAMPLES)
    output:
        OUTDIR + "output/SNPs.csv"
    run:
        df = merger.frame_merge(list(input))
        df.to_csv(str(output), sep = '\t')

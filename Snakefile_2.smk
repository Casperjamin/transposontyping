import glob
from SCRIPTS import tetyperead
from SCRIPTS import merger
from SCRIPTS import report
from SCRIPTS import coverage_filter
import pandas as pd
import sys
from Bio import SeqIO
configfile: "samples/good_samples.yaml"
configfile: "config/config.yaml"
configfile: "config/reference.yaml"
configfile: "config/score.yaml"
SAMPLES = config['SAMPLES']


CASETTE = config['CASETTE']
SCORE = config['SCORE']
STRETCHDEL = [True, False]

casette = CASETTE["ref1"]
SNP_score = SCORE["minimal_quality"]


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
    expand("output/{sample}/coverage_plots/", sample = SAMPLES),
    expand('output/{sample}/quality_plots/', sample = SAMPLES),
    expand("output/{sample}/tetype/{sample}_summary.txt", sample = SAMPLES),
    expand("output/{sample}/TEType_convert/{sample}_{stretchdel}_TEType_df.tsv", sample = SAMPLES, stretchdel = STRETCHDEL),
    expand("output/SNPs_{stretchdel}.csv", stretchdel = STRETCHDEL),
    expand("output/report/pairwise_distance_{stretchdel}.csv", stretchdel = STRETCHDEL),
    expand("output/report/dendrogram_{stretchdel}.png",stretchdel = STRETCHDEL),

################################################################################
# fastp read trimming
################################################################################
rule trim_reads:
    input:
        forward = lambda wildcards: SAMPLES[wildcards.sample]['forward'],
        reverse = lambda wildcards: SAMPLES[wildcards.sample]['reverse']
    output:
        forward = "input/{sample}/forward_trimmed.fastq.gz",
        reverse = "input/{sample}/reverse_trimmed.fastq.gz",
        report = directory("output/{sample}/fastp_reports/")
    conda:
        "envs/fastp_environment.yml"
    shell:
        "fastp -i {input.forward} -I {input.reverse} -o {output.forward} -O {output.reverse} -j {output.report}fastp.json -h {output.report}fastp.html"


################################################################################
# generating figures, pairwise comparisson and bedgraph used for coverage plot
################################################################################

rule generate_heatmap_dendrogram:
    input:
        csv = "output/SNPs_{stretchdel}.csv"
    output:
        "output/report/pairwise_distance_{stretchdel}.csv",
        "output/report/dendrogram_{stretchdel}.png"
    params:
        "output/report/"
    run:
        os.system(f'mkdir -p {params}')
        report.generate_report(input, params, get_length_dna(casette), wildcards.stretchdel)


rule generate_coverage_matrix:
    input:
        "output/{sample}/tetype/{sample}.bam"
    output:
        "output/{sample}/coverage/{sample}_bedgraph"
    conda:
        "envs/TEType_environment.yml"
    params:
        "output/{sample}/tetype/{sample}"
    shell:
        "samtools sort -o {params}_sorted.bam {params}.bam &&"
        "bedtools genomecov -ibam {params}_sorted.bam -d > {output}"

rule generate_coverage_plot:
    input:
        "output/{sample}/coverage/{sample}_bedgraph"
    output:
        directory("output/{sample}/coverage_plots/")
    run:
        report.generate_plot_coverage(str(input), output, wildcards.sample)
        report.generate_hist_coverage(str(input), output, wildcards.sample)

rule make_sam_files:
    input:
        "output/{sample}/tetype/{sample}.bam"
    output:
        "output/{sample}/tetype/{sample}.sam"
    shell:
        "samtools view -h -o {output} {input}"

rule generate_SNP_quality_plot:
    input:
        "output/{sample}/tetype/{sample}.vcf"
    output:
        directory('output/{sample}/quality_plots/')
    run:
        report.generate_plot_SNP_quality(str(input), output, wildcards.sample)

################################################################################
# TETyper, conversion and merging
################################################################################

rule TETyping:
  input:
    forward = "input/{sample}/forward_trimmed.fastq.gz",
    reverse = "input/{sample}/reverse_trimmed.fastq.gz"
  output:
    "output/{sample}/tetype/{sample}_summary.txt",
    "output/{sample}/tetype/{sample}.bam",
    "output/{sample}/tetype/{sample}.vcf"
  threads:
    4
  conda:
    "envs/TEType_environment.yml"
  log:
    "logs/tetype/{sample}.txt"
  params:
    "output/{sample}/tetype/{sample}"
  shell:
    " TETyper.py --outprefix {params}"
    " --ref {casette}"
    " --fq1 {input.forward} "
    " --fq2 {input.reverse} "
    " --flank_len 10"
    " --threads {threads} "
    "2> {log} "

rule convert:
  input:
    loc = "output/{sample}/tetype/{sample}_summary.txt",
    vcf = "output/{sample}/tetype/{sample}.vcf"
  output:
    tsv = "output/{sample}/TEType_convert/{sample}_{stretchdel}_TEType_df.tsv"
  run:
    summary = tetyperead.tetype_result(f"{input.loc}", f"{casette}")
    df = summary.result_to_frame(input.loc, input.vcf, f"{SNP_score}", wildcards.stretchdel, output.tsv)
    df.to_csv(f"{output.tsv}", sep = "\t")

rule merge:
    input:
        samples = expand("output/{sample}/TEType_convert/{sample}_{stretchdel}_TEType_df.tsv", sample = SAMPLES, stretchdel = STRETCHDEL)
    output:
        "output/SNPs_{stretchdel}.csv"
    run:
        df = merger.frame_merge([f"{input}"], wildcards.stretchdel, f'{output}')

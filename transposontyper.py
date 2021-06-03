#!/usr/bin/env python

from argparse import ArgumentParser
import os
import glob
import sys
import yaml


def obtain_repoloc():
    return os.path.dirname(os.path.abspath(__file__))

def get_absolute_path(path):
    return os.path.abspath(path)
LOCATIONREPO = obtain_repoloc()

def define_ref(ref):
    casette_loc = str(ref)[2:-2]
    casette_dict = {"CASETTE": {}}
    location = glob.glob(casette_loc)[0]

    casette_dict["CASETTE"] = {
    "ref1": get_absolute_path(location)}
    casette_dump = yaml.dump(casette_dict, default_flow_style=False)

    with open(f"{LOCATIONREPO}/config/reference.yaml", 'w+') as f:
        f.write(casette_dump)


def check_presence_fastq(SAMPLES):
    # exit statement if no fastq can be found
    if len(SAMPLES) == 0:
        print("\nNo fastq files found in subdirectories\n")
        print("Please give a input directory subdirectories with fastq files \nexiting now.")
        sys.exit()

def define_input(inputdir):
    inputdir = get_absolute_path(inputdir)
    # check if valid
    SAMPLES = glob.glob(inputdir + "/*/*fastq.gz")
    check_presence_fastq(SAMPLES)
    samplesdict = {"SAMPLES":{}}
    for i in SAMPLES:
        samplename = i.split("/")[-2]
        forward = glob.glob(f"{inputdir}/{samplename}/*1*f*q*")[0]
        reverse = glob.glob(f"{inputdir}/{samplename}/*2*f*q*")[0]

        samplesdict["SAMPLES"][str(samplename)] = {
        "forward": get_absolute_path(forward),
        "reverse": get_absolute_path(reverse)}

    samples_loc = yaml.dump(samplesdict, default_flow_style=False)
    os.system(f"mkdir -p {LOCATIONREPO}/samples")
    with open(f"{LOCATIONREPO}/samples/samples.yaml", 'w+') as f:
        f.write(samples_loc)

def define_output(outputdir):
    outdir = {"OUTPUT":outputdir.strip("/")}
    outdirloc = yaml.dump(outdir, default_flow_style=False)
    os.system(f"mkdir -p {LOCATIONREPO}/config")
    with open(f"{LOCATIONREPO}/config/config.yaml", 'w+') as f:
        f.write(outdirloc)

def define_score(score):
    score_dict = {"SCORE": {}}
    score_dict["SCORE"] = {
    "minimal_quality": int(score)}
    score_dump = yaml.dump(score_dict, default_flow_style=False)
    with open("config/score.yaml", 'w+') as score_file:
        score_file.write(score_dump)

def launch(cores):
    #change to the location of the repo, this will make sure all envs, databases and other stuff sticks in the repocryptic
    os.chdir(f"{LOCATIONREPO}")
    os.system(f"snakemake  --use-conda --cores {cores} --snakefile Snakefile_1.smk")
    os.system(f"snakemake  --use-conda --cores {cores} --snakefile Snakefile_2.smk")

def main(command_line = None):
    #add main parser object
    parser = ArgumentParser(description = "input directory with subdirectories with reads")
    parser.add_argument("-i", required = True, dest = "input_dir")
    parser.add_argument("-o", required = True, dest = "output_dir")
    parser.add_argument("--ref", required = True, dest = "ref", nargs = 1)
    parser.add_argument("--cores", required = False, dest = "cores", default = 8, type = int)
    parser.add_argument("-SNP_score", required = False, dest = "SNP_score", default = 20, type = int)
    args = parser.parse_args(command_line)
    define_input(args.input_dir)
    define_output(args.output_dir)
    #define_ref(args.ref)
    define_ref(args.ref)
    define_score(args.SNP_score)
    launch(args.cores)

if __name__ == "__main__":
    main()

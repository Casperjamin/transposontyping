import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
import sys
import csv
import operator
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
import os
import math

def generate_matrix(SNPfile, lengthfrag):
    """The generate_matrix function takes a merged dataframe of SNP locations for al samples as input.
    This is then converted to a squareform matrix and hamming distance is calculated"""
    df = pd.read_csv(str(SNPfile), sep = "\t", index_col = 0)
    distance = pdist(df, 'hamming')
    namelist = df.index
    hamming_matrix = pd.DataFrame(squareform(distance), columns = namelist, index = namelist)
    SNP_matrix = hamming_matrix * lengthfrag
    return SNP_matrix


def generate_plot_hamming(SNP_matrix, output):
    """The generate_plot_hamming function takes the hamming distance matrix from the generate_matrix function as input.
    this is then used to construct a heatmap and dendrogram indicating phylogeny of compared samples"""
    #make linkage for heatmap and dendrogram
    namelist = SNP_matrix.index
    matrix = squareform(SNP_matrix)
    #print(type(matrix)) #matrix is hier een numpy.ndarray
    Z = linkage(matrix)
    order = leaves_list(Z)
    SNP_matrix = SNP_matrix.iloc[order, order]
    #generate heatmap figure
    plt.figure(figsize = [10,10])
    sns.heatmap(SNP_matrix, cmap = 'plasma')
    plt.savefig(f"{output}/heatmap.png")
    #generate dendrogram figure
    plt.figure(figsize = [10,10])
    boom = dendrogram(Z, orientation = "right", labels = namelist)
    plt.savefig(f'{output}/dendrogram.png')

def determine_fraction(pos, bin):
    """
    takes nucleotide position and binsize
    return a string of the range of the bin
    """
    frac_num = math.ceil(pos/bin)
    frac_top = (frac_num*bin)
    frac_bot = (frac_top - bin)
    return str(frac_bot) + '-' + str(frac_top)

def generate_plot_coverage(bedgraph, output, sample):
    """
    Take the bedgraph as input and generate a coverage plot for each position
    """
    df1 = pd.read_csv(bedgraph, sep='\t', usecols=[1,2], names=["position", "read_count"])
    df1['fraction'] = [determine_fraction(x, 500) for x in df1["position"]]
    plt.figure(figsize = [20,15])
    plot1 = sns.boxplot(data=df1, x=df1["fraction"], y=df1['read_count'])
    plot1.set_xticklabels(df1['fraction'].unique(), rotation=90)
    plt.savefig(f'{output}/{sample}_coverage_plot.png')

def generate_hist_coverage(bedgraph, output, sample):
    df = pd.read_csv(bedgraph, sep='\t')
    list_coverage = list(df.iloc[:,2])
    plt.figure(figsize = [20,15])
    plt.hist(list_coverage, bins=50)
    plt.savefig(f'{output}/{sample}_coverage_hist.png')

# mogelijk een global variable maken die de minimale gewenste snip kwaliteit aangeeft dit is nu handmatig ingevoerd hieronder
def generate_plot_SNP_quality(vcf, output, sample):
    vcf_file = open(vcf, 'r')
    df = pd.DataFrame()
    positions = []
    qualities = []
    for line in vcf_file:
        if str(line)[0] != "#" and int(math.ceil(float(line.split('\t')[5]))) > 1:
            quality = float(line.split('\t')[5])
            position = int(line.split('\t')[1])

            positions.append(int(line.split('\t')[1]))
            qualities.append(float(line.split('\t')[5]))
    df['positions'] = positions
    df['qualities'] = qualities
    df['fraction'] = [determine_fraction(x, 250) for x in df['positions']]
    plt.figure(figsize = [20,15])
    plot1 = sns.boxplot(data=df, x='fraction',  y='qualities')
    plt.savefig(f'{output}/{sample}_SNP_boxplot.png')
    plot2 = sns.swarmplot(data=df,x='fraction', y ='qualities')
    plt.savefig(f'{output}/{sample}_SNP_swarmplot.png')


"""The generate_pairwise_distance function produces a report of pair wise comparison for all analysed samples"""
def generate_pairwise_distance(matrix, outdir):
    df = matrix.unstack().reset_index()
    df.columns = ["sample1", "sample2", "distance"]
    df.to_csv(f"{str(outdir)}/pairwise_distance.csv", sep = '\t', header = True)

"""The generate_report function is used to call upon the functions defined above in this script"""
def generate_report(input, outdir, lengthfrag):
    matrix = generate_matrix(input, lengthfrag)
    generate_pairwise_distance(matrix, outdir)
    generate_plot_hamming(matrix, outdir)

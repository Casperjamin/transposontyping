import pandas as pd
from Bio import SeqIO
import math

class tetype_result:
    def __init__(self, summary, fragment):
        self.summary = summary
        self.fragment = fragment
        self.name = str(self.summary.split("/")[-1].split("_")[0])


    def result_to_frame(self, summary, vcf, score, stretchdel, output, DEL = True, homoSNPs = True):
        result = pd.read_csv(self.summary, sep = '\t')
        for i in SeqIO.parse(self.fragment, 'fasta'):
            lengthfrag = len(i.seq)

        df = pd.DataFrame(columns = [str(x + 1) for x in range(lengthfrag)], index = [self.name]).fillna(0)

        with open(vcf, 'r') as vcf_input:
            approved = []
            snplist = result["SNPs_homozygous"][0].split('|')
            for line in vcf_input:
                if str(line)[0] != "#":
                    if float(line.split('\t')[5]) > int(score):
                        approved.append(int(line.split('\t')[1]))
            for i in snplist:
                if str(i) != "none":
                    if str(int(str(i)[1:-1]) in approved) == 'True':
                        position = int(str(i)[1:-1])
                        SNP = str(i)[-1]
                        if SNP == "A":
                            df.loc[self.name, str(position)] = 1
                        if SNP == "T":
                            df.loc[self.name, str(position)] = 2
                        if SNP == "C":
                            df.loc[self.name, str(position)] = 3
                        if SNP == "G":
                            df.loc[self.name, str(position)] = 4




        if stretchdel == "True":
            listregion = []
            if DEL == True:
                delete = str(result["Deletions"][0])
                delete_list = delete.split('|')
                for i in delete_list:
                    if str(i) != "none":
                        if len(str(i).split("-")) == 2:
                            region = str(i).split("-")
                            for i in range(int(region[0]), int(region[1]) + 1):
                                listregion.append(i)
                            for i in listregion:
                                df.loc[self.name, str(i)] = int(5)
                        if len(str(i).split("-")) == 1:
                            #als we 1 positie deletie hebben dan is krijgt deze de waarde DEL
                            df.loc[self.name, str(i)] = int(5)
                            #als de regio van een deletie groter is dan 1000 word dit bescouwd als een enkele SNP op de eerste positie

        if stretchdel == "False":
            if DEL == True:
                delete = str(result["Deletions"][0])
                delete_list = delete.split('|')
                for i in delete_list:
                    if str(i) != "none":
                        if len(str(i).split("-")) == 2:
                            region = str(i).split("-")
                            position = int(region[0])
                            df.loc[self.name, str(position)] = int(5)
                        if len(str(i).split("-")) == 1:
                            df.loc[self.name, str(i)] = int(5)

        return df

import pandas as pd
from Bio import SeqIO


class tetype_result:
    def __init__(self, summary, fragment):
        self.summary = summary
        self.fragment = fragment
        self.name = str(self.summary.split("/")[-1].split("_")[0])


    def result_to_frame(self, homoSNPs = True, DEL = True):
        result = pd.read_csv(self.summary, sep = '\t')

        #get legnth fragment to make df
        for i in SeqIO.parse(self.fragment, 'fasta'):
            lengthfrag = len(i.seq)
        #generate empty DF
        df = pd.DataFrame(columns = [str(x + 1) for x in range(lengthfrag)], index = [self.name]).fillna(0)
        if homoSNPs:
            for i in result.SNPs_homozygous[0].split("|"):
                if i !='none':
                    position = str(i[1:-1])
                    SNP = i[-1]

                    df.loc[self.name, position] = SNP

        if DEL:
            delete = result.Deletions.iloc[0]
            for i in delete.split("|"):
                if i != "none":
                    # if we have 1 deletion instead of a stretch, this position assigned as deleted
                    if len(i.split("-")) == 1:
                        df.loc[namesample, i] = "DEL"

                    if len(i.split("-")) == 2:
                        region = i.split("-")
                        listregion = [str(x) for x in range(int(region[0]), int(region[1]))]
                        df.loc[self.name, listregion] = "DEL"

        recode = {
            "A":1,
            "T":2,
            "C":3,
            "G":4,
            "DEL":5
        }
        df.replace(recode, inplace = True)
        return df

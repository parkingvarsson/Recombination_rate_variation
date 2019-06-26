
##### Works by typing complement.py <vcf> <agp> #####
##### Switches the negative strands to their complement, as required by #####
##### LDhelmet #####
import sys
import pandas as pd
import numpy as np

### Tries to read the needed columns of your files
def inputs():
    try:
        vcf = pd.read_table(sys.argv[1], header=None, sep = "\t", comment="#", engine="c")
    except:
        raise
        print("{0} cannot be read.".format(sys.argv[1]))
    try:
        agp = pd.read_table(sys.argv[2], header=None, sep = "\t", comment="#", engine="c")
    except:
        raise
        print("{0} cannot be read.".format(sys.argv[2]))
    return vcf, agp

### Complements the negative strand positions within the vcf -file
def complement(vcf, agp):
    agp = agp[agp[8]=="-"]
    for chrom in np.unique(agp[0]):
        agp_sub = agp[agp[0]==chrom]
        for index, row in agp_sub.iterrows():
            vcf.loc[(vcf[0]==chrom) & (vcf[1] >= row[1]) & (vcf[1] <= row[2]) & (vcf[3]=="A"), 3] = 2
            vcf.loc[(vcf[0]==chrom) & (vcf[1] >= row[1]) & (vcf[1] <= row[2]) & (vcf[3]=="T"), 3] = 1
            vcf.loc[(vcf[0]==chrom) & (vcf[1] >= row[1]) & (vcf[1] <= row[2]) & (vcf[3]=="G"), 3] = 4
            vcf.loc[(vcf[0]==chrom) & (vcf[1] >= row[1]) & (vcf[1] <= row[2]) & (vcf[3]=="C"), 3] = 3
            vcf.loc[(vcf[0]==chrom) & (vcf[1] >= row[1]) & (vcf[1] <= row[2]) & (vcf[4]=="A"), 4] = 2
            vcf.loc[(vcf[0]==chrom) & (vcf[1] >= row[1]) & (vcf[1] <= row[2]) & (vcf[4]=="T"), 4] = 1
            vcf.loc[(vcf[0]==chrom) & (vcf[1] >= row[1]) & (vcf[1] <= row[2]) & (vcf[4]=="G"), 4] = 4
            vcf.loc[(vcf[0]==chrom) & (vcf[1] >= row[1]) & (vcf[1] <= row[2]) & (vcf[4]=="C"), 4] = 3
    vcf.loc[vcf[3]== 2, 3] = "T"
    vcf.loc[vcf[3]== 1, 3] = "A"
    vcf.loc[vcf[3]== 4, 3] = "C"
    vcf.loc[vcf[3]== 3, 3] = "G"
    vcf.loc[vcf[4]== 2, 4] = "T"
    vcf.loc[vcf[4]== 1, 4] = "A"
    vcf.loc[vcf[4]== 4, 4] = "C"
    vcf.loc[vcf[4]== 3, 4] = "G"
    return vcf

### Prints out the results as a table
def print_out(vcf):
    name = sys.argv[1].rsplit(".", 1)[0] + ".negative_strand_complemented." + sys.argv[1].rsplit(".", 1)[-1]
    vcf.to_csv(name, header = False, index = None, sep = "\t")

### The main module calling all others
if __name__ == "__main__":
    vcf, agp = inputs()
    vcf = complement(vcf, agp)
    print_out(vcf)

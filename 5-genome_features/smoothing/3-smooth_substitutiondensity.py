##### Script for smoothing the different substitution densities #####

import sys
import pandas as pd

### Config.
### Choose here the window and step size for your smoothing operation.

window_size = 1000000
step_size = 250000

### Actual program
### Tries to read the input
def inputs():
    try:
        smooth = pd.read_table(sys.argv[1], header=0, sep = "\t", comment="#", engine="c")
    except:
        raise
        print("{0} cannot be read.".format(sys.argv[1]))
    smooth.sort_values(["CHROM","BIN_START"], axis=0, ascending=True, inplace=True)
    Chromosome_list = smooth["CHROM"].copy()
    Chromosome_list.drop_duplicates(inplace=True)
    return smooth, Chromosome_list

### Smoothes the file
def smoothe(smooth, Chromosome_list):
    tuples = []
    na = []
    for chromosome in Chromosome_list:
#        print(smooth["Chromosome"])
        subset_smooth = smooth[smooth["CHROM"]==chromosome]
        a = 0
        b = window_size
        while a < max(subset_smooth["BIN_START"]):
            if b < min(subset_smooth["BIN_START"]):
                pass
            elif len(subset_smooth[(a < subset_smooth["BIN_START"]) & (subset_smooth["BIN_START"] <= b)]) == 0:
                c = "NA"
                d = "NA"
                e = "NA"
                na.append([chromosome, a, b, c, d, e])
            elif b > max(subset_smooth["BIN_START"]):
                c = subset_smooth.loc[(a < subset_smooth["BIN_START"]) & (subset_smooth["BIN_START"] <= b), "total_subs_density"].mean(skipna=True)
                d = subset_smooth.loc[(a < subset_smooth["BIN_START"]) & (subset_smooth["BIN_START"] <= b), "old_subs_density"].mean(skipna=True)
                e = subset_smooth.loc[(a < subset_smooth["BIN_START"]) & (subset_smooth["BIN_START"] <= b), "new_subs_density"].mean(skipna=True)
                tuples.append([chromosome, a, b, c, d, e])
                break
            else:
                c = subset_smooth.loc[(a < subset_smooth["BIN_START"]) & (subset_smooth["BIN_START"] <= b), "total_subs_density"].mean(skipna=True)
                d = subset_smooth.loc[(a < subset_smooth["BIN_START"]) & (subset_smooth["BIN_START"] <= b), "old_subs_density"].mean(skipna=True)
                e = subset_smooth.loc[(a < subset_smooth["BIN_START"]) & (subset_smooth["BIN_START"] <= b), "new_subs_density"].mean(skipna=True)
                tuples.append([chromosome, a, b, c, d, e])
            a = a+step_size
            b = b+step_size
    results = pd.DataFrame(tuples)
#    genome_average = len(subs[0])/len(neuts[0])
#    results[3] = results[2]/genome_average
    if len(na) > 0:
        na_values = pd.DataFrame(na)
        na_values[3] = na_values[2]
        results = results.append(na_values)
    else:
        pass
    results.sort_values([0, 1], axis=0, ascending=True, inplace=True)
    return results


### Columnizes and writes the result
def making_the_result(results):
    results.columns = ["chromosome", "win_start", "win_end", "total_subs_density", "old_subs_density", "new_subs_density"]
    name = sys.argv[1].rsplit(".", 1)[0] + ".1Mbp_window." + sys.argv[1].rsplit(".", 1)[-1]
    results.to_csv(name, header = True, index = None, sep = "\t")

### The main module calling all others
if __name__ == "__main__":
    smooth, Chromosome_list = inputs()
    results = smoothe(smooth, Chromosome_list)
    making_the_result(results)

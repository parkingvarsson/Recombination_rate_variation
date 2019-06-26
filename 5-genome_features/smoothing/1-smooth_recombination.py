##### Script for smoothing the two recombination rate estimates #####

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
    smooth.sort_values(["chromosome","position"], axis=0, ascending=True, inplace=True)
    chromosome_list = smooth["chromosome"].copy()
    chromosome_list.drop_duplicates(inplace=True)
#    print(smooth)
    return smooth, chromosome_list

### Smoothes the file
def smoothe(smooth, chromosome_list):
    tuples = []
    na = []
    for chromosome in chromosome_list:
#        print(smooth["chromosome"])
        subset_smooth = smooth[smooth["chromosome"]==chromosome]
        print(chromosome)
        a = 0
        b = window_size
        print(max(subset_smooth["position"]))
        while a < max(subset_smooth["position"]):
            print(b)
            if b < min(subset_smooth["position"]):
                pass
            elif len(subset_smooth[(a < subset_smooth["position"]) & (subset_smooth["position"] <= b)]) == 0:
                c = "NA"
                d = "NA"
                na.append([chromosome, a, b, c, d])
            elif b > max(subset_smooth["position"]):
                c = subset_smooth.loc[(a < subset_smooth["position"]) & (subset_smooth["position"] <= b), "gen_slidingwindow"].mean(skipna=True)
                d = subset_smooth.loc[(a < subset_smooth["position"]) & (subset_smooth["position"] <= b), "seq_slidingwindow"].mean(skipna=True)
                tuples.append([chromosome, a, b, c, d])
                break
            else:
                c = subset_smooth.loc[(a < subset_smooth["position"]) & (subset_smooth["position"] <= b), "gen_slidingwindow"].mean(skipna=True)
                d = subset_smooth.loc[(a < subset_smooth["position"]) & (subset_smooth["position"] <= b), "seq_slidingwindow"].mean(skipna=True)
                tuples.append([chromosome, a, b, c, d])
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
    results.columns = ["chromosome", "win_start", "win_end", "gen_slidingwindow", "seq_slidingwindow"]
    name = sys.argv[1].rsplit(".", 1)[0] + ".1Mbp_window." + sys.argv[1].rsplit(".", 1)[-1]
    results.to_csv(name, header = True, index = None, sep = "\t")

### The main module calling all others
if __name__ == "__main__":
    smooth, chromosome_list = inputs()
    results = smoothe(smooth, chromosome_list)
    making_the_result(results)

##### Script for smoothing the 3 methylation contexts in all 6 individuals. ####

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
    smooth.sort_values(["chromosome","win_start"], axis=0, ascending=True, inplace=True)
    chromosome_list = smooth["chromosome"].copy()
    chromosome_list.drop_duplicates(inplace=True)
    return smooth, chromosome_list

#### Smoothes the file
def smoothe(smooth, chromosome_list):
    tuples = []
    na = []
    for chromosome in chromosome_list:
#        print(smooth["chromosome"])
        subset_smooth = smooth[smooth["chromosome"]==chromosome]
        a = 0
        b = window_size
        while a < max(subset_smooth["win_start"]):
            if b < min(subset_smooth["win_end"]):
                pass
            elif len(subset_smooth[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b)]) == 0:
                c = "NA"
                d = "NA"
                e = "NA"
                f = "NA"
                g = "NA"
                h = "NA"
                i = "NA"
                j = "NA"
                k = "NA"
                l = "NA"
                m = "NA"
                n = "NA"
                o = "NA"
                p = "NA"
                q = "NA"
                r = "NA"
                s = "NA"
                t = "NA"
                na.append([chromosome, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t])
            elif b > max(subset_smooth["win_start"]):
                c = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CpG_swasp001"].mean(skipna=True)
                d = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHG_swasp001"].mean(skipna=True)
                e = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHH_swasp001"].mean(skipna=True)
                f = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CpG_swasp005"].mean(skipna=True)
                g = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHG_swasp005"].mean(skipna=True)
                h = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHH_swasp005"].mean(skipna=True)
                i = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CpG_swasp044"].mean(skipna=True)
                j = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHG_swasp044"].mean(skipna=True)
                k = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHH_swasp044"].mean(skipna=True)
                l = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CpG_swasp046"].mean(skipna=True)
                m = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHG_swasp046"].mean(skipna=True)
                n = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHH_swasp046"].mean(skipna=True)
                o = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CpG_swasp113"].mean(skipna=True)
                p = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHG_swasp113"].mean(skipna=True)
                q = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHH_swasp113"].mean(skipna=True)
                r = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CpG_swasp116"].mean(skipna=True)
                s = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHG_swasp116"].mean(skipna=True)
                t = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHH_swasp116"].mean(skipna=True)
                tuples.append([chromosome, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t])
                break
            else:
#                c = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "SwAsp001_meth%"].mean(skipna=True)
#                d = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "SwAsp005_meth%"].mean(skipna=True)
#                e = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "SwAsp044_meth%"].mean(skipna=True)
#                f = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "SwAsp046_meth%"].mean(skipna=True)
#                g = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "SwAsp113_meth%"].mean(skipna=True)
#                h = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "SwAsp116_meth%"].mean(skipna=True)
                c = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CpG_swasp001"].mean(skipna=True)
                d = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHG_swasp001"].mean(skipna=True)
                e = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHH_swasp001"].mean(skipna=True)
                f = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CpG_swasp005"].mean(skipna=True)
                g = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHG_swasp005"].mean(skipna=True)
                h = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHH_swasp005"].mean(skipna=True)
                i = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CpG_swasp044"].mean(skipna=True)
                j = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHG_swasp044"].mean(skipna=True)
                k = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHH_swasp044"].mean(skipna=True)
                l = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CpG_swasp046"].mean(skipna=True)
                m = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHG_swasp046"].mean(skipna=True)
                n = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHH_swasp046"].mean(skipna=True)
                o = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CpG_swasp113"].mean(skipna=True)
                p = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHG_swasp113"].mean(skipna=True)
                q = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHH_swasp113"].mean(skipna=True)
                r = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CpG_swasp116"].mean(skipna=True)
                s = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHG_swasp116"].mean(skipna=True)
                t = subset_smooth.loc[(a < subset_smooth["win_start"]) & (subset_smooth["win_end"] <= b), "CHH_swasp116"].mean(skipna=True)
                tuples.append([chromosome, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t])
            a = a+step_size
            b = b+step_size
    results = pd.DataFrame(tuples)
#    genome_average = len(subs[0])/len(neuts[0])
#    results[3] = results[2]/genome_average
    if len(na) > 0:
        na_values = pd.DataFrame(na)
#        na_values[3] = na_values[2]
        results = results.append(na_values)
    else:
        pass
    results.sort_values([0, 1], axis=0, ascending=True, inplace=True)
    return results


### Columnizes and writes the result
def making_the_result(results):
    results.columns = ["chromosome", "win_start", "win_end",
 "CpG_swasp001", "CHG_swasp001", "CHH_swasp001",
 "CpG_swasp005", "CHG_swasp005", "CHH_swasp005",
 "CpG_swasp044", "CHG_swasp044", "CHH_swasp044",
 "CpG_swasp046", "CHG_swasp046", "CHH_swasp046",
 "CpG_swasp113", "CHG_swasp113", "CHH_swasp113",
 "CpG_swasp116", "CHG_swasp116", "CHH_swasp116",]
    name = sys.argv[1].rsplit(".", 1)[0] + ".1Mbp_window." + sys.argv[1].rsplit(".", 1)[-1]
    results.to_csv(name, header = True, index = None, sep = "\t")

### The main module calling all others
if __name__ == "__main__":
    smooth, chromosome_list = inputs()
    results = smoothe(smooth, chromosome_list)
    making_the_result(results)

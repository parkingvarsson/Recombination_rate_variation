
##### Script for calculating methylation levels for all three motifs #####
##### in all six individuals in the same file #####

import sys
import pandas as pd

## Config.
# Choose here the window size for your mutation rate estimate.

window_size = 50000

## Actual program
# Tries to read the needed columns of your meths and ends tables
def inputs():
    try:
        meths = pd.read_table(sys.argv[1], header=0, sep = "\t", comment="#", engine="c") #File containing methylation
    except:
        raise
        print("{0} cannot be read.".format(sys.argv[2]))
    try:
        ends = pd.read_table(sys.argv[2], header=None, sep = "\t", comment="#", engine="c") #File containing chromosome names and lengths
    except:
        raise
        print("{0} cannot be read.".format(sys.argv[3]))
    return meths, ends

# Sorts the tables as the way I am going to do it is based upon sorted table.
def sort(meths, ends):
    meths.sort_values(["chromosome", "position"], axis=0, ascending=True, inplace=True)
    ends.sort_values([0, 1], axis=0, ascending=True, inplace=True)
    return meths, ends

# Calculates window averages for all methylation contexts
def methylation(meths, ends):
    tuples = []
    na = []
    for index, row in ends.iterrows():
        subset_meths = meths.loc[(meths["chromosome"] == ends.iloc[index, 0]), :]
        a = 0
        b = window_size
        while a < max(subset_meths["position"]):
            if b < min(subset_meths["position"]):
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
                na.append([ends.iloc[index, 0], a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t])
            elif len(subset_meths[(a < subset_meths["position"]) & (subset_meths["position"] <= b)]) == 0:
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
                na.append([ends.iloc[index, 0], a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t])
            elif b > ends.iloc[index, 1]:
                b = ends.iloc[index, 1]
                c = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CpG_swasp001"].mean(skipna=True)
                d = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHG_swasp001"].mean(skipna=True)
                e = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHH_swasp001"].mean(skipna=True)
                f = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CpG_swasp005"].mean(skipna=True)
                g = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHG_swasp005"].mean(skipna=True)
                h = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHH_swasp005"].mean(skipna=True)
                i = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CpG_swasp044"].mean(skipna=True)
                j = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHG_swasp044"].mean(skipna=True)
                k = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHH_swasp044"].mean(skipna=True)
                l = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CpG_swasp046"].mean(skipna=True)
                m = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHG_swasp046"].mean(skipna=True)
                n = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHH_swasp046"].mean(skipna=True)
                o = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CpG_swasp113"].mean(skipna=True)
                p = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHG_swasp113"].mean(skipna=True)
                q = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHH_swasp113"].mean(skipna=True)
                r = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CpG_swasp116"].mean(skipna=True)
                s = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHG_swasp116"].mean(skipna=True)
                t = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHH_swasp116"].mean(skipna=True)
                tuples.append([ends.iloc[index, 0], a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t])
            else:
                c = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CpG_swasp001"].mean(skipna=True)
                d = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHG_swasp001"].mean(skipna=True)
                e = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHH_swasp001"].mean(skipna=True)
                f = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CpG_swasp005"].mean(skipna=True)
                g = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHG_swasp005"].mean(skipna=True)
                h = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHH_swasp005"].mean(skipna=True)
                i = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CpG_swasp044"].mean(skipna=True)
                j = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHG_swasp044"].mean(skipna=True)
                k = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHH_swasp044"].mean(skipna=True)
                l = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CpG_swasp046"].mean(skipna=True)
                m = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHG_swasp046"].mean(skipna=True)
                n = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHH_swasp046"].mean(skipna=True)
                o = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CpG_swasp113"].mean(skipna=True)
                p = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHG_swasp113"].mean(skipna=True)
                q = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHH_swasp113"].mean(skipna=True)
                r = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CpG_swasp116"].mean(skipna=True)
                s = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHG_swasp116"].mean(skipna=True)
                t = subset_meths.loc[(a < subset_meths["position"]) & (subset_meths["position"] <= b), "CHH_swasp116"].mean(skipna=True)
                tuples.append([ends.iloc[index, 0], a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t])
            a = a+window_size
            b = b+window_size
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


# Columnizes and writes the result
def making_the_result(results):
    results.columns = ["chromosome", "win_start", "win_end",
 "CpG_swasp001", "CHG_swasp001", "CHH_swasp001",
 "CpG_swasp005", "CHG_swasp005", "CHH_swasp005",
 "CpG_swasp044", "CHG_swasp044", "CHH_swasp044",
 "CpG_swasp046", "CHG_swasp046", "CHH_swasp046",
 "CpG_swasp113", "CHG_swasp113", "CHH_swasp113",
 "CpG_swasp116", "CHG_swasp116", "CHH_swasp116",]
    name = sys.argv[1].rsplit(".", 1)[0] + ".windowed." + sys.argv[1].rsplit(".", 1)[-1]
    results.to_csv(name, header = True, index = None, sep = "\t")

# The main module calling all others
if __name__ == "__main__":
    meths, ends = inputs()
    meths, ends = sort(meths, ends)
    results = methylation(meths, ends)
    making_the_result(results)

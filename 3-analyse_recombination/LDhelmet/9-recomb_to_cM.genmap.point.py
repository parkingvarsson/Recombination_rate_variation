##### Works by typing recomb_to_cM.py <recomb> <mareymap> #####

import sys
import pandas as pd

### Actual program ###

### Tries to read the needed columns of your recomb and mareymap tables
def inputs():
    try:
        recomb = pd.read_table(sys.argv[1], header=None, skiprows=3, sep = " ", comment="#", engine="c")
    except:
        raise
        print("{0} cannot be read.".format(sys.argv[1]))
    try:
        mareymap = pd.read_table(sys.argv[2], header=0, skiprows=1, sep = " ", comment="#", engine="c")
    except:
        raise
        print("{0} cannot be read.".format(sys.argv[2]))
    #print(mareymap)
    mareymap = mareymap[mareymap["map"] == sys.argv[3]]
    mareymap_gen_max = max(mareymap["gen"])
    mareymap.sort_values(by=["phys"], axis=0, inplace=True)
    mareymap = mareymap["phys"].tolist()
    return recomb, mareymap, mareymap_gen_max

### Produces a gene map from recombination estimates of 1/Mb
def points(recomb, mareymap, mareymap_gen_max):
    tuples = []
    for pos in mareymap:
        pos_recomb = recomb.loc[((pos < recomb[1]) & (pos >= recomb[0])), 2].values[0]
        print(pos_recomb)
        tuples.append([pos, pos_recomb])
    results = pd.DataFrame(tuples)
    tuples = []
    whole_genome_recomb = sum(results[1])
    scaling_factor = whole_genome_recomb/mareymap_gen_max
    additive_centimorgan = 0
    results[2] = results[1]/scaling_factor
    for index, row in results.iterrows():
        additive_centimorgan = row[2] + additive_centimorgan
        tuples.append(additive_centimorgan)
    results[3] = tuples
    return results


### Columnizes and writes the result
def making_the_result(results):
    name = sys.argv[1].rsplit(".", 1)[0] + ".cM." + sys.argv[1].rsplit(".", 1)[-1]
    results.to_csv(name, header = False, index = None, sep = "\t")

### The main module calling all others
if __name__ == "__main__":
    recomb, mareymap, mareymap_gen_max = inputs()
    results = points(recomb, mareymap, mareymap_gen_max)
    making_the_result(results)

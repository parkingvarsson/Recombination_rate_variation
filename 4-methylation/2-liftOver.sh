#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 03:00:00

###Liftover of bismark cov files.
###Perform step by step.

###Convert to bed

for file in /path/to/unzipped/covs/PE*
do
    output=${file%%.bismark.cov}.bed
    awk '{print $1"\t"$2"\t"$3+1"\t"$4"\t"$5"\t"$6}' $file > $output
done

###bedPlus=3 means that only the first 3 columns are used in the conversion. The rest just tag along.

#Liftover to split scaffolds
liftOver -bedPlus=3 $1 v1.1_Potra01-genome.fa.masked.chain ${1%%.bed}.split.bed ${1%%.bed}.split.unmapped

#Liftover to mixed
liftOver -bedPlus=3 $1 F1.chain ${1%%.split.bed}.chr.bed ${1%%.split.bed}.chr.unmapped

###Back to coverage files and removal on unlifted scaffolds
for file in /path/to/wherever/you/put/the/bed/transformed/covs/PE*.chr.bed
do
    output=${file%%.chr.bed}.chr.bismark.cov
    awk '{print $1"\t"$2"\t"$3-1"\t"$4"\t"$5"\t"$6}' $file | grep "chr" > $output
    echo $output
done

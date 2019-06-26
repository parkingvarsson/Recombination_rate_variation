#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:10:00

##### Depth filtering for bismark cov -files #####
awk -F "\t" '$5+$6 > 4 && $5+$6 < 46 { print $0 }' $1 > ${1%%.bismark.cov}.5-45_filter_cov.bismark.cov

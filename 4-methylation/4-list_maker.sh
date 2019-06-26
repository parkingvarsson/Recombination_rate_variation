#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00

##### Script that produces list of methylated positions #####

filename=${1##*/} # This is one of the context files
zcat $1 | awk '{print $3"\t"$4"\t"$2}' > "/output/folder/"${filename%%.txt.gz}".methylated_list.txt"

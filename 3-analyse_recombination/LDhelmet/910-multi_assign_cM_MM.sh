#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00

##### Looper script with two loops 1) for performing the cM conversion. #####
##### See 8-recomb_to_cM.genmap.point.py for more information. 2) for #####
##### making the mareymaps input files from the cM converted files. See #####
##### 9-MM_input.sh for more info. #####

### Run these separately by commenting one out ###

### For making cM conversion:
file_path="/path/to/readable/ldhelmet/outputs/"
mareymap="/path/to/consensus/mareymap/for/length/" # This file is LMB estimates
for file in $file_path/*.txt; do
    filename=${file##*/}
    chrom=${filename%%.*}
    echo $file $chrom
    sbatch run_recomb_to_cM.genmap.py.sh $file $mareymap $chrom
done

### For making mareymap inputs:
file_path="/path/to/cM/converted/ldhelmet/results/"
for file in $file_path/*.cM.txt; do
    sh MM_input.sh $file
done

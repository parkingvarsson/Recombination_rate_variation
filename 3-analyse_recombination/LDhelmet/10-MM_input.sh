#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00

##### Produces MareyMaps input files for each chromosome from the cM #####
##### converted ldhelmet outputs. #####

cMldhelmet=$1 # The cM converted LDhelmet output file

filename=${cMldhelmet##*/}
chrom=${filename%%.*}

#echo $1 $chrom

### Actual script. produces first the required headers and then inserts the
### data in required format with all the "" and tabs and whatnots.
awk -v chr=$chrom 'BEGIN {print "\"""set""\"""\t""\"""map""\"""\t""\"""mkr""\"""\t""\"""phys""\"""\t""\"""gen""\""} {print "\"""Populus_tremula_consensus_sequence""\"""\t""\""chr"\"""\t""\""chr":"$1"\"""\t"$1"\t"$4}' $cMldhelmet > $chrom.MM_input.txt

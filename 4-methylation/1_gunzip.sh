#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00

##### This script contains a loop for unzipping the bismark cov -files #####

for folder in path/to/bismark/results/SwAsp*
do
    echo $folder
    for file in $folder/new/*bismark.cov.gz
    do
        echo $file
        unzipped=${file##*/}
        named=${unzipped%%.gz}
        gunzip -c $file > $unzipped
        mv $unzipped $named
    done
done

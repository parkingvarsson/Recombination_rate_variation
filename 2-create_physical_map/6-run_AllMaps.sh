##### This file contains the commands for running AllMaps to get #####
##### Potra v1.2 as well as to estimate gap length between the assembled #####
##### scaffolds. #####

#### Run these separately ####
### 1. Run Allmaps as instructed on AllMaps github:
### https://github.com/tanghaibao/jcvi/wiki/ALLMAPS
python -m jcvi.assembly.allmaps merge F1_Female.csv F1_Male.csv -o F1.bed
python -m jcvi.assembly.allmaps path F1.bed F1_SPLIT.fasta

### 2. Estimate gaps with AllMaps
python -m jcvi.assembly.allmaps estimategaps F1.bed
mv F1.estimategaps.agp F1.chr.agp
python -m jcvi.assembly.allmaps build F1.bed F1_SPLIT.fasta

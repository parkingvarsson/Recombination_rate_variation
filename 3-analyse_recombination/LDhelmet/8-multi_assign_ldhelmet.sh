##### Looping multi assigner to run each chromosome separately on cluster. #####
##### Run each assigner loop separately. #####

### For conf files:
file_path="/path/to/fasta_all/" #Contains file each chromosome separately, containing the two haplotypes for each of the 25 individuals. Simple cat -command works to make these files.
for file in $file_path/*.all.fa; do
sbatch run_ldhelmet.sh $file
done

### For	lkuptable files:
file_path="/path/to/conf/"
for file in $file_path/*.conf; do
sbatch run_ldhelmet.sh $file
done

### For pade files:
file_path="/path/to/conf/"
for file in $file_path/*.conf; do
sbatch run_ldhelmet_pade.sh $file
done

### For unsplit rjmcmc run:
file_path="/path/to/ldhelmet/topdirectory/"
for file in $file_path/5_LDhelmet_snps_pos/*.pos; do
filename=${file##*/}
chrom=${filename%%.*}
sbatch run_ldhelmet.sh 3_lkuptable/$chrom.all.fa.conf.lk $file_path/4_pade/$chrom.all.fa.conf.pade $file $file_path/5_LDhelmet_snps_pos/$chrom.ldhelmet.snps
# Edit the sbatch command to meet your directory structure.
done

### For making readables:
file_path="/path/to/rjmcmc/results"
for file in $file_path/*.post; do
echo $file
sbatch run_ldhelmet_to_readable.sh $file
done

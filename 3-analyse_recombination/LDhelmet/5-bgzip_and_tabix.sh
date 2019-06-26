#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

##### Zips and indexes the vcf files #####

FILE_DIR="/path/to/location/with/chr.vcfs/" #path to the location with the chr.vcf files
fastapath="/path/to/all/chr.fasta/" #path to all chr.fasta files


### Saves the header
grep "#" $FILE_DIR/chr1.vcf > header

### Sorts the vcfs, for this the header must be removed and reattached
for filename in	$FILE_DIR/chr*.vcf; do
    cat header > temp
    grep -v	"#" $filename | sort -k2,2n >> temp
    mv temp $filename
done

### Loops through, zips and indexes the vcfs
for filename in $FILE_DIR/chr*.vcf; do
bgzip $filename ;done;
for zipped in $FILE_DIR/chr*.vcf.gz; do
tabix $zipped ;done;

### Loops through and produces the haplotype fasta-files for each individual
declare -a samples=$(zcat $FILE_DIR/chr1.vcf.gz | cut -f 10- | awk 'FNR==10 {print $0}')
for zipped in $FILE_DIR/chr*.vcf.gz; do
    zippedname="${zipped##*/}"
    zippedstart="${zippedname%%.*}"
    fasta="$fastapath/$zippedstart.fa"
    for sample in ${samples[@]} ;do
        echo $sample
        sbatch consensus_vcftools.sh $fasta $zipped $sample
    done
done

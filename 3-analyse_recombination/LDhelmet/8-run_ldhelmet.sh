##### LDhelmet run. Run step by step for each chromosome separately. #####
##### See input file information on "7_multi_assign.sh" #####
##### This script has been performed in ldhelmetv1.10 #####

### Create the configurations from your fasta files
ldhelmet find_confs --num_threads 10 -w 50 -o $1".conf" $1

### Create a lookup table
ldhelmet table_gen --num_threads 19 -c $1 -t 0.0000000375 -r 0.0 0.1 10.0 1.0 100.0 -o $1.lk

### Pade table
ldhelmet pade --num_threads 20 -c $1 -t 0.0000000375 -x 11 -o $1.pade

### Run the actual thing with pos and snps files
ldhelmet rjmcmc --num_threads 18 --lk_file $1 --pade_file $2 -w 50 -b 50.0 --pos_file $3 --snps_file $4 --burn_in 100000 --num_iter 1000000 -o $1.post

### We did not use fasta -files, but left this in as we dabbled in it ###
### Run the actual thing with fasta file
#######ldhelmet rjmcmc --num_threads 5 -l $1 -p $2 -w 50 -b 50.0 -s $3 --burn_in 1000 -n 10000 -o test_fasta.post

### Binary to readable
ldhelmet post_to_text -m -p 0.025 -p 0.50 -p 0.0975 -o $1.txt $1

##### Loop for deduplicate.sh #####
for file in "/path/to/files/to/be/deduplicated/"*
do
  echo $file
  sbatch deduplicate.sh $file
done

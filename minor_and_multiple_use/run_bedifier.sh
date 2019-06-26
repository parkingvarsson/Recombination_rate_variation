##### Loop for bedifier.sh #####

for file in "path/to/name/pos/file"*
  do
  echo $file
  sbatch bedifier.sh $file
done

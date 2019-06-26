##### Loop that goes through the folders containing bismark results and #####
##### loops through the files in those folders submitting them to #####
##### 4-list_maker.sh #####


file_path="/path/to/bismark/results/"
for folder in $file_path/SwAsp*; do
  for file in $folder/new/C*_context_*.deduplicated.txt.gz; do
    echo $file
    sbatch 4-list_maker.sh $file
  done
done

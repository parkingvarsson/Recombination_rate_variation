for file in /path/to/lifted/over/covs/*sorted.bismark.cov
    do
        echo $file
        sh awk_filter.sh $file
    done

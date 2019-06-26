library(BatchMap)
load('../om-base.RData')

args <- commandArgs(trailingOnly = TRUE)

message(args[1])

ovlp <- 25
size <- 50

lg <- testcrosses[[ args[1] ]]

lg <- record.parallel(lg, times = 16, cores = 16)

mp <- map.overlapping.batches(
        lg,
        overlap=ovlp,
        phase.cores=4,
        ripple.cores=32,
        ws=11, min.tries=1,
        size=pick.batch.sizes(lg, size, ovlp),
        fun.order=ripple.ord,
        method='one',
        verbosity=c('batch','order'))

save(mp, file=paste0(args[1],'.rd'))

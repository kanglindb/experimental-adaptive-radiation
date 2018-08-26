### This file contains R command lines of haplotype block reconstruction for selection regime CS
library("haploReconstruct")
CS.data <- sync_to_frequencies(file="CS.sync", base.pops=c(rep(TRUE, 5),rep(FALSE, 10)), header=FALSE)
dat_filtered=initialize_SNP_time_series(chr=CS.data$chr, pos=CS.data$pos, base.freq=CS.data$basePops, lib.freqs=CS.data[,7:ncol(CS.data), with=FALSE], pop.ident=rep(c(rep(1:5)), 3), pop.generation=c(rep(0, 5), rep(21, 5), rep(31, 5)), use.libs=rep(TRUE,15), minfreqchange=0.2, minrepl=2, max.minor.freq=0.1, winsize=1e+06)
dat_reconst.chr2L=reconstruct_hb(dat_filtered, chrom="chr2L")
res <- tryCatch({
  summary(dat_reconst.chr2L)
}, error = function(e) {
  message(e)
  return(NULL)
})
if(!is.null(res)) {
  pdf("CS.chr2L.freq.pdf")
  par(mfrow=c(5,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
  write.table(file="CS.data_reconst.chr2L", summary(dat_reconst.chr2L), quote=F)
  for(i in 1:nrow(summary(dat_reconst.chr2L)[,2])) {
    write.table(file=paste("CS.chr2L", i, "markers", sep="."), markers(dat_reconst.chr2L, i))
  }
  plot_hbr_freq(dat_reconst.chr2L, hbr_id=1, replicate=1, timepoint=c(0,21,31), window=1, xlab = "chr2L [Mb]")
  plot_hbr_freq(dat_reconst.chr2L, hbr_id=1, replicate=2, timepoint=c(0,21,31), window=1, xlab = "chr2L [Mb]")
  plot_hbr_freq(dat_reconst.chr2L, hbr_id=1, replicate=3, timepoint=c(0,21,31), window=1, xlab = "chr2L [Mb]")
  plot_hbr_freq(dat_reconst.chr2L, hbr_id=1, replicate=4, timepoint=c(0,21,31), window=1, xlab = "chr2L [Mb]")
  plot_hbr_freq(dat_reconst.chr2L, hbr_id=1, replicate=5, timepoint=c(0,21,31), window=1, xlab = "chr2L [Mb]")
  dev.off()
}
dat_reconst.chr2R=reconstruct_hb(dat_filtered, chrom="chr2R")
res <- tryCatch({
  summary(dat_reconst.chr2R)
}, error = function(e) {
  message(e)
  return(NULL)
})
if(!is.null(res)) {
  pdf("CS.chr2R.freq.pdf")
  par(mfrow=c(5,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
  write.table(file="CS.data_reconst.chr2R", summary(dat_reconst.chr2R), quote=F)
  for(i in 1:nrow(summary(dat_reconst.chr2R)[,2])) {
    write.table(file=paste("CS.chr2R", i, "markers", sep="."), markers(dat_reconst.chr2R, i))
  }
  plot_hbr_freq(dat_reconst.chr2R, hbr_id=1, replicate=1, timepoint=c(0,21,31), window=1, xlab = "chr2R [Mb]")
  plot_hbr_freq(dat_reconst.chr2R, hbr_id=1, replicate=2, timepoint=c(0,21,31), window=1, xlab = "chr2R [Mb]")
  plot_hbr_freq(dat_reconst.chr2R, hbr_id=1, replicate=3, timepoint=c(0,21,31), window=1, xlab = "chr2R [Mb]")
  plot_hbr_freq(dat_reconst.chr2R, hbr_id=1, replicate=4, timepoint=c(0,21,31), window=1, xlab = "chr2R [Mb]")
  plot_hbr_freq(dat_reconst.chr2R, hbr_id=1, replicate=5, timepoint=c(0,21,31), window=1, xlab = "chr2R [Mb]")
  dev.off()
}
dat_reconst.chr3L=reconstruct_hb(dat_filtered, chrom="chr3L")
res <- tryCatch({
  summary(dat_reconst.chr3L)
}, error = function(e) {
  message(e)
  return(NULL)
})
if(!is.null(res)) {
  pdf("CS.chr3L.freq.pdf")
  par(mfrow=c(5,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
  write.table(file="CS.data_reconst.chr3L", summary(dat_reconst.chr3L), quote=F)
  for(i in 1:nrow(summary(dat_reconst.chr3L)[,2])) {
    write.table(file=paste("CS.chr3L", i, "markers", sep="."), markers(dat_reconst.chr3L, i))
  }
  plot_hbr_freq(dat_reconst.chr3L, hbr_id=1, replicate=1, timepoint=c(0,21,31), window=1, xlab = "chr3L [Mb]")
  plot_hbr_freq(dat_reconst.chr3L, hbr_id=1, replicate=2, timepoint=c(0,21,31), window=1, xlab = "chr3L [Mb]")
  plot_hbr_freq(dat_reconst.chr3L, hbr_id=1, replicate=3, timepoint=c(0,21,31), window=1, xlab = "chr3L [Mb]")
  plot_hbr_freq(dat_reconst.chr3L, hbr_id=1, replicate=4, timepoint=c(0,21,31), window=1, xlab = "chr3L [Mb]")
  plot_hbr_freq(dat_reconst.chr3L, hbr_id=1, replicate=5, timepoint=c(0,21,31), window=1, xlab = "chr3L [Mb]")
  dev.off()
}
dat_reconst.chr3R=reconstruct_hb(dat_filtered, chrom="chr3R")
res <- tryCatch({
  summary(dat_reconst.chr3R)
}, error = function(e) {
  message(e)
  return(NULL)
})
if(!is.null(res)) {
  pdf("CS.chr3R.freq.pdf")
  par(mfrow=c(5,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
  write.table(file="CS.data_reconst.chr3R", summary(dat_reconst.chr3R), quote=F)
  for(i in 1:nrow(summary(dat_reconst.chr3R)[,2])) {
    write.table(file=paste("CS.chr3R", i, "markers", sep="."), markers(dat_reconst.chr3R, i))
  }
  plot_hbr_freq(dat_reconst.chr3R, hbr_id=1, replicate=1, timepoint=c(0,21,31), window=1, xlab = "chr3R [Mb]")
  plot_hbr_freq(dat_reconst.chr3R, hbr_id=1, replicate=2, timepoint=c(0,21,31), window=1, xlab = "chr3R [Mb]")
  plot_hbr_freq(dat_reconst.chr3R, hbr_id=1, replicate=3, timepoint=c(0,21,31), window=1, xlab = "chr3R [Mb]")
  plot_hbr_freq(dat_reconst.chr3R, hbr_id=1, replicate=4, timepoint=c(0,21,31), window=1, xlab = "chr3R [Mb]")
  plot_hbr_freq(dat_reconst.chr3R, hbr_id=1, replicate=5, timepoint=c(0,21,31), window=1, xlab = "chr3R [Mb]")
  dev.off()
}
dat_reconst.chrX=reconstruct_hb(dat_filtered, chrom="chrX")
res <- tryCatch({
  summary(dat_reconst.chrX)
}, error = function(e) {
  message(e)
  return(NULL)
})
if(!is.null(res)) {
  pdf("CS.chrX.freq.pdf")
  par(mfrow=c(5,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
  write.table(file="CS.data_reconst.chrX", summary(dat_reconst.chrX), quote=F)
  for(i in 1:nrow(summary(dat_reconst.chrX)[,2])) {
    write.table(file=paste("CS.chrX", i, "markers", sep="."), markers(dat_reconst.chrX, i))
  }
  plot_hbr_freq(dat_reconst.chrX, hbr_id=1, replicate=1, timepoint=c(0,21,31), window=1, xlab = "chrX [Mb]")
  plot_hbr_freq(dat_reconst.chrX, hbr_id=1, replicate=2, timepoint=c(0,21,31), window=1, xlab = "chrX [Mb]")
  plot_hbr_freq(dat_reconst.chrX, hbr_id=1, replicate=3, timepoint=c(0,21,31), window=1, xlab = "chrX [Mb]")
  plot_hbr_freq(dat_reconst.chrX, hbr_id=1, replicate=4, timepoint=c(0,21,31), window=1, xlab = "chrX [Mb]")
  plot_hbr_freq(dat_reconst.chrX, hbr_id=1, replicate=5, timepoint=c(0,21,31), window=1, xlab = "chrX [Mb]")
  dev.off()
}


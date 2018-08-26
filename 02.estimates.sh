### This file contains the R commands for Ne estimates and s (selection coefficient) estimates for CS regime, other regimes can be processed in the same manner
### script in R
  ### Part 1, Ne estimate
  library("poolSeq")
  gen <- c(0, 0, 0, 0, 0, 21, 21, 21, 21, 21, 31, 31, 31, 31, 31)
  repl <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
  CS.Sync <- read.sync(file="CS.sync", gen=gen, repl=repl)
  CS.Af <- af(CS.Sync, repl=repl, gen=gen)
  CS.Cov <- coverage(CS.Sync, repl=repl, gen=gen)
  Ne.file <- paste("CS", "Ne.est", sep=".")
  CS.NEest.R1 <- estimateNe(p0=CS.Af[,"F0.R1.freq"], pt=CS.Af[,"F31.R1.freq"], cov0=CS.Cov[,"F0.R1.cov"], covt=CS.Cov[,"F31.R1.cov"], t=31, Ncensus=300, poolSize=c(300, 300), method=c("W.planI", "W.planII"))
  write.table(CS.NEest.R1, file=Ne.file, append = TRUE, col.names = FALSE)
  CS.NEest.R2 <- estimateNe(p0=CS.Af[,"F0.R2.freq"], pt=CS.Af[,"F31.R2.freq"], cov0=CS.Cov[,"F0.R2.cov"], covt=CS.Cov[,"F31.R2.cov"], t=31, Ncensus=300, poolSize=c(300, 300), method=c("W.planI", "W.planII"))
  write.table(CS.NEest.R2, file=Ne.file, append = TRUE, col.names = FALSE)
  CS.NEest.R3 <- estimateNe(p0=CS.Af[,"F0.R3.freq"], pt=CS.Af[,"F31.R3.freq"], cov0=CS.Cov[,"F0.R3.cov"], covt=CS.Cov[,"F31.R3.cov"], t=31, Ncensus=300, poolSize=c(300, 300), method=c("W.planI", "W.planII"))
  write.table(CS.NEest.R3, file=Ne.file, append = TRUE, col.names = FALSE)
  CS.NEest.R4 <- estimateNe(p0=CS.Af[,"F0.R4.freq"], pt=CS.Af[,"F31.R4.freq"], cov0=CS.Cov[,"F0.R4.cov"], covt=CS.Cov[,"F31.R4.cov"], t=31, Ncensus=300, poolSize=c(300, 300), method=c("W.planI", "W.planII"))
  write.table(CS.NEest.R4, file=Ne.file, append = TRUE, col.names = FALSE)
  CS.NEest.R5 <- estimateNe(p0=CS.Af[,"F0.R5.freq"], pt=CS.Af[,"F31.R5.freq"], cov0=CS.Cov[,"F0.R5.cov"], covt=CS.Cov[,"F31.R5.cov"], t=31, Ncensus=300, poolSize=c(300, 300), method=c("W.planI", "W.planII"))
  write.table(CS.NEest.R5, file=Ne.file, append = TRUE, col.names = FALSE)

  ### Part 2, s estimate
  Sh.file <- paste("CS", "SH.est", sep=".")
  for (i in seq(1:nrow(CS.Af))) {
    CS.SHest <- estimateSH(CS.Af[i,], Ne=round(107.8714), t=gen, simulate.p.value=TRUE, method="NLS", h=0.5, N.ctraj=1000)
    write.table(file=Sh.file, row.names=F, col.names=F, t(c(row.names(CS.Af)[i], CS.SHest$s, CS.SHest$p.value)), append=T)
  }


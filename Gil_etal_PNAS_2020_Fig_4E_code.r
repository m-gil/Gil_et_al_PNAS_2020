# Code below reproduces main text Fig. 4E
require(doSNOW)
require(foreach)

fmaxXrt.plot.fun.dynf <- function(tf = 5000, f_hi = 0.22, f_lo = 0.1, f_rez = 12, fmax_vec = seq(from = f_lo, to = f_hi, length.out = f_rez), fr_lo = 0, fr_hi = 0.018, fr_rez = 12, frate_vec = seq(from = fr_lo, to = fr_hi, length.out = fr_rez), add0 = 0.11, b = 0.5, a = 1, rA = 2, rC = 0.1, alpha = 0, m = 0.04, mu = 0.02, p = 0, q = 1, s2 = 0){
  
  cl<-makeCluster(6) # number of CPU cores
  registerDoSNOW(cl)
  
  outlist <- foreach(i = 1:length(fmax_vec)) %dopar% {
    require(deSolve)
    source("population_models.r")
    f_hi <- fmax_vec[i]
    SIfwd_mat <- matrix(nrow = length(frate_vec), ncol = 8)
    noSIfwd_mat <- matrix(nrow = length(frate_vec), ncol = 8)
    f <- 0
    
    for (j in 1:length(frate_vec)){
      s <- 1
      add <- add0/(a*b)
      coral_dom0 <- c(1, 0.1, 0.56, 0) # [4] is f
      SIparmsi <- c(add=add, a=a,b=b,q=q,rA=rA,rC=rC,m=m,mu=mu, frate = 0, fmax = 0)
      times <- seq(0, tf, by = 100)
      SIfwd_dynfi <- ode(coral_dom0, times, SI.fun2f, SIparmsi, method = "adams")
      
      prist_H <- as.numeric(SIfwd_dynfi[length(times), 2])
      noSIfwd_pcap <- a*b*(add + prist_H^s/(1 + q*prist_H^s))
      b2 <- noSIfwd_pcap/a
      
      SIcoral_dom <- SIfwd_dynfi[length(times), 2:5]
      SIparms <- c(add=add, a=a,b=b,q=q,rA=rA,rC=rC,m=m,mu=mu, frate = frate_vec[j], fmax = fmax_vec[i])
      SIfwd_dynf <- ode(SIcoral_dom, times, SI.fun2f, SIparms, method = "adams")
      
      Cf1 = SIfwd_dynf[length(SIfwd_dynf[,1]), 4]
      Af1 = SIfwd_dynf[length(SIfwd_dynf[,1]), 3]
      Hf1 = SIfwd_dynf[length(SIfwd_dynf[,1]), 2]
      ff1 = SIfwd_dynf[length(SIfwd_dynf[,1]), 5]
      
      SIfwd_mat[j, ] <- c(1, s, fmax_vec[i], frate_vec[j], round(Cf1, 2), round(Af1, 2), round(Hf1, 2), round(ff1, 2))
      
      s <- 0
      noSIparmsi <- c(add=add, a=a,b=b2,q=q,f=f,rA=rA,rC=rC,m=m,mu=mu, frate = 0, fmax = 0)
      noSIfwd_dynfi <- ode(coral_dom0, times, noSI.fun2f, noSIparmsi, method = "adams")
      noSIcoral_dom <- noSIfwd_dynfi[length(times), 2:5]
      noSIparms <- c(add=add, a=a,b=b2,q=q,rA=rA,rC=rC,m=m,mu=mu, frate = frate_vec[j], fmax = fmax_vec[i])
      noSIfwd_dynf <- ode(noSIcoral_dom, times, noSI.fun2f, noSIparms, method = "adams")
      
      Cf2 = noSIfwd_dynf[length(noSIfwd_dynf[,1]), 4]
      Af2 = noSIfwd_dynf[length(noSIfwd_dynf[,1]), 3]
      Hf2 = noSIfwd_dynf[length(noSIfwd_dynf[,1]), 2]
      ff2 = noSIfwd_dynf[length(noSIfwd_dynf[,1]), 5]
      
      noSIfwd_mat[j, ] <- c(1, s, fmax_vec[i], frate_vec[j], round(Cf2, 2), round(Af2, 2), round(Hf2, 2), round(ff2, 2))
    } 
    fwd_mat <- rbind(SIfwd_mat, noSIfwd_mat)
    colnames(fwd_mat) <- c("fwd", "SI", "fmax", "frate", "Cf", "Af", "Hf", "ff")
    return(fwd_mat)
  }
  stopCluster(cl)
  
  outdf <- data.frame(do.call(rbind, outlist))
  SIoutdf <- subset(outdf, outdf$SI == 1)
  noSIoutdf <- subset(outdf, outdf$SI == 0)
  
  SIoutdf$Cf01 <- as.numeric(SIoutdf$Cf > 0)
  SIoutdf$Hf01 <- as.numeric(SIoutdf$Hf > 0)
  SIoutdf$CfHf_check <- abs(SIoutdf$Cf01 - SIoutdf$Hf01)
  
  noSIoutdf$Cf01 <- as.numeric(noSIoutdf$Cf > 0)
  noSIoutdf$Hf01 <- as.numeric(noSIoutdf$Hf > 0)
  noSIoutdf$CfHf_check <- abs(noSIoutdf$Cf01 - noSIoutdf$Hf01)
  par(mfrow = c(2, 1))
  xrng <- range(SIoutdf$frate)
  yrng <- range(SIoutdf$fmax)
  with(subset(SIoutdf, Cf > 0), plot(fmax ~ frate, xlim = xrng, ylim = yrng, xlab = "fishing increase per year", ylab = "target fishing level", xaxs = "i", col = "chocolate3", pch = 15))
  with(subset(SIoutdf, Cf < 0.01), points(fmax ~ frate, col = "white", pch = 15))
  title(paste("SI, rC=", rC, ", m=", m, ", rA=", rA, ", b=", b, ", a=", round(a, 2), ", mu=", mu, sep = ""))
  with(subset(noSIoutdf, Cf > 0), plot(fmax ~ frate, xlim = xrng, col = "chocolate3", ylim = yrng, xlab = "fishing increase per year", ylab = "target fishing level", pch =15))
  with(subset(noSIoutdf, Cf < 0.01), points(fmax ~ frate, col = "white", pch = 15))
  title("noSI")
  return(SIoutdf)
}

ptm <- proc.time()
rC02xx <- fmaxXrt.plot.fun.dynf(tf = 5000, rC = 0.2, m = 0.08, mu = 0.02, f_lo = 0.17, f_hi = 0.25, fr_hi = 0.035, fr_lo = 0.0, fr_rez = 80, f_rez = 36)
proc.time() - ptm


rC = 0.2; rA = 2; m = 0.08; mu = 0.02; b = 0.5; a = 1

par(mfrow = c(1, 1))
with(subset(rC02xx, Cf > 0), plot(fmax ~ frate, xlab = "fishing increase per year", ylab = "target fishing level", xaxs = "i", yaxs = "i", xlim = c(0, 0.035), ylim = c(0.16, 0.26), col = rgb(219/255, 186/255, 197/255), pch = 15))
with(subset(rC02xx, Cf < 0.01), points(fmax ~ frate, col = "white", pch = 15))
title(paste("SI, rC=", rC, ", m=", m, ", rA=", rA, ", b=", b, ", a=", round(a, 2), ", mu=", mu, sep = ""))


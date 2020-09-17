# Code below reproduces Fig. 2B
require(deSolve)
source("population_models.R")
democurve.fun <- function(X) 0.11 + 0.5*X/(1 + X)

traj.fun.trans.trait = function(alg_dom = c(0.000001, 0.999999, 0.000001, 0.5), frfwd = 0.01, frbwd = -0.01, tf = 5000, tfeqi = 5000, add0 = 0.11, s = 1, b = 0.5, a = 1, rA = 2, rC = 0.02, alpha = 0, m = 0.008, f_lo = 0, f_hi = 0.5, mu = 0.02, p = 0, frez = 160, q = 1, DD = "no", ee = 2, g = 6, plt = "yes") {
  
  SIfwd_out <- matrix(999999, nrow = frez, ncol = 8)
  noSIfwd_out <- SIfwd_out
  SIbwd_out <- SIfwd_out
  noSIbwd_out <- SIfwd_out
  
  SIfwdtrait_out <- matrix(999999, nrow = frez, ncol = 8)
  SIbwdtrait_out <- matrix(999999, nrow = frez, ncol = 8)
  noSIbwdtrait_out <- SIfwd_out
  
  add = add0/(a*b) 
  
  
  coral_dom0 = c(1, 0.1, 0.56, 0)
  fmax_vec = seq(from=f_lo, to=f_hi, length.out=frez)
  q = q; p = p; mu = mu
  
  SIparmsFWD0 = c(fmax = 0, frate = 0, add=add, a=a,b=b,q=q,rA=rA, rC=rC,m=m,mu=mu)
  timesi <- seq(0, tfeqi, by = 100)
 
  SIfwd0 = ode(coral_dom0, timesi, SI.fun2f, SIparmsFWD0, method = "adams")
  prist_H <- SIfwd0[length(SIfwd0[, 1]), 2]
  noSIfwd_pcap <- a*b*(add + prist_H^s/(1 + q*prist_H^s))
  b2fwd <- noSIfwd_pcap/a
  noSIparmsFWD0 <- c(a = a, b = as.numeric(b2fwd), rA=rA, rC=rC, m=m, mu=mu, fmax = 0, frate = 0)
  noSIfwd0 = ode(coral_dom0, timesi, noSI.fun2f, noSIparmsFWD0, method = "adams")
   
  load("bestfitparms_trait.Rda")
  

  d1 <- bestestest_trait_mod[1]
  d2 <- bestestest_trait_mod[2]
  dgrd_H <-alg_dom[1]
  SIparmsFWDtrait0 <- c(a = a, rA=rA,  rC=rC, m=m, mu=mu, fmax = 0, frate = 0, d1 = d1, d2 = d2)
  SIfwd0trait <- ode(coral_dom0, timesi, SI.fun2f.trait, SIparmsFWDtrait0, method = "adams")
  noSIbwd_pcap <- a*b*(add + dgrd_H^s/(1 + q*dgrd_H^s))
  b2bwd <- noSIbwd_pcap/a
  b3bwd <- d1 
  for (i in 1:frez) { 
    fmax <- fmax_vec[i]
    coral_dom <- SIfwd0[length(timesi), 2:5]
    coral_dom_trait <- SIfwd0trait[length(timesi), 2:5]
    SIparmsFWD= c(fmax = fmax, frate = frfwd, add=add, a=a,b=b,q=q,rA=rA,rC=rC,m=m,mu=mu)
    noSIparmsFWD= c(fmax = fmax, frate = frfwd, a=a, b=as.numeric(b2fwd), rA=rA, rC=rC,m=m,mu=mu)
    SIparmsBWD<- c(a = a, b = b, add = add, q=q,rA=rA,  rC=rC, m=m, mu=mu, fmin = fmax, frate = frbwd)
    noSIparmsBWD <- c(a = a, b = as.numeric(b2bwd), rA=rA,  rC=rC, m=m, mu=mu, fmin = fmax, frate = frbwd)
    
    SIparmsFWDtrait <- c(fmax = fmax, frate = frfwd, a=a,rA=rA,rC=rC,m=m,mu=mu, d1 = d1, d2 = d2)
    SIparmsBWDtrait <- c(fmin = fmax, frate = frbwd, a=a,rA=rA,rC=rC,m=m,mu=mu, d1 = d1, d2 = d2)
    noSIparmsBWDtrait <- c(a = a, b = as.numeric(b3bwd), rA=rA,  rC=rC, m=m, mu=mu, fmin = fmax, frate = frbwd)
    
    times = seq(0, tf, by = 100)
    SIfwd = ode(coral_dom, times, SI.fun2f, SIparmsFWD, method = "adams")
    noSIfwd = ode(coral_dom, times, noSI.fun2f, noSIparmsFWD, method = "adams")
    SIbwd = ode(alg_dom, times, SI.fun2fBWD, SIparmsBWD, method = "adams")
    noSIbwd = ode(alg_dom, times, noSI.fun2fBWD, noSIparmsBWD, method = "adams")
    
    SIfwdtrait <- ode(coral_dom_trait, times, SI.fun2f.trait, SIparmsFWDtrait, method = "adams")
    SIbwdtrait <- ode(alg_dom, times, SI.fun2fBWD.trait, SIparmsBWDtrait, method = "adams")
    noSIbwdtrait = ode(alg_dom, times, noSI.fun2fBWD, noSIparmsBWDtrait, method = "adams")
       
    if (sum(abs(SIfwd[length(SIfwd[, 1]), 2:4] - SIfwd[length(SIfwd[, 1]) - 1, 2:4])) > 0.001) print("SI fwd mod failed to converge maybe")
    if (sum(abs(noSIfwd[length(noSIfwd[, 1]), 2:4] - noSIfwd[length(noSIfwd[, 1]) - 1, 2:4])) > 0.001) print("noSI fwd mod failed to converge maybe")
    if (sum(abs(SIbwd[length(SIbwd[, 1]), 2:4] - SIbwd[length(SIbwd[, 1]) - 1, 2:4])) > 0.001) print("SI bwd mod failed to converge maybe")
    if (sum(abs(noSIbwd[length(noSIbwd[, 1]), 2:4] - noSIbwd[length(noSIbwd[, 1]) - 1, 2:4])) > 0.001) print("noSI bwd mod failed to converge maybe")
    
    if (sum(abs(SIfwdtrait[length(SIfwdtrait[, 1]), 2:4] - SIfwdtrait[length(SIfwdtrait[, 1]) - 1, 2:4])) > 0.001) print("SI fwdtrait mod failed to converge maybe")
    if (sum(abs(SIbwdtrait[length(SIbwdtrait[, 1]), 2:4] - SIbwdtrait[length(SIbwdtrait[, 1]) - 1, 2:4])) > 0.001) print("SI bwdtrait mod failed to converge maybe")
    if (sum(abs(noSIbwdtrait[length(noSIbwdtrait[, 1]), 2:4] - noSIbwdtrait[length(noSIbwdtrait[, 1]) - 1, 2:4])) > 0.001) print("noSI bwd mod failed to converge maybe")
    
    SIfwd_out[i, ] <- c(1, 1, fmax, SIfwd[length(SIfwd[, 1]), ])
    noSIfwd_out[i, ] <- c(0, 1, fmax, noSIfwd[length(noSIfwd[, 1]), ])
    SIbwd_out[i, ] <- c(1, 0, fmax, SIbwd[length(SIbwd[, 1]), ])
    noSIbwd_out[i, ] <- c(0, 0, fmax, noSIbwd[length(noSIbwd[, 1]), ])
    
    SIfwdtrait_out[i, ] <- c(1, 1, fmax, SIfwdtrait[length(SIfwdtrait[, 1]), ])
    SIbwdtrait_out[i, ] <- c(1, 0, fmax, SIbwdtrait[length(SIbwdtrait[, 1]), ])
    noSIbwdtrait_out[i, ] <- c(0, 0, fmax, noSIbwdtrait[length(noSIbwdtrait[, 1]), ])
    
    
  }
  colnames(SIfwd_out) <- c("SI", "fwd", "ftarg", "tf", "Hf", "Af", "Cf", "ff")
  colnames(noSIfwd_out) <- c("SI", "fwd", "ftarg", "tf", "Hf", "Af", "Cf", "ff")
  colnames(SIbwd_out) <- c("SI", "fwd", "ftarg", "tf", "Hf", "Af", "Cf", "ff")
  colnames(noSIbwd_out) <- c("SI", "fwd", "ftarg", "tf", "Hf", "Af", "Cf", "ff")
  
  colnames(SIfwdtrait_out) <- c("SI", "fwd", "ftarg", "tf", "Hf", "Af", "Cf", "ff")
  colnames(SIbwdtrait_out) <- c("SI", "bwd", "ftarg", "tf", "Hf", "Af", "Cf", "ff")
  colnames(noSIbwdtrait_out) <- c("SI", "bwd", "ftarg", "tf", "Hf", "Af", "Cf", "ff")
  
  if (plt == "yes") {
    SIfwd_df <- data.frame(SIfwd_out)
    noSIfwd_df <- data.frame(noSIfwd_out)
    SIbwd_df <- data.frame(SIbwd_out)
    noSIbwd_df <- data.frame(noSIbwd_out)
    
    SIfwdtrait_df <- data.frame(SIfwdtrait_out)
    SIbwdtrait_df <- data.frame(SIbwdtrait_out)
    noSIbwdtrait_df <- data.frame(noSIbwdtrait_out)
    
    par(mfrow = c(3, 2))
    with(SIfwd_df, plot(Hf ~ ftarg, type = "l", col = "blue", cex.lab = 1.5, ylim = c(0, 2), xaxs = "i", yaxs = "i", bty = "n", lwd = 2, ylab = "herbivore biomass", xlab = "fishing"))
    axis(side = 1, lwd = 2)
    axis(side = 2, lwd = 2)
    with(noSIfwd_df, points(Hf ~ ftarg, type = "l", lwd = 2, col = "red"))
    with(SIfwdtrait_df, points(Hf ~ ftarg, type = "l", lwd = 2, col = "purple"))
    title(paste("SI (blue) vs. noSI (red)\nFishing pristine reef: ", frfwd, " increase per year\n ", sep = ""))
    
    with(SIbwd_df, plot(Hf ~ ftarg, type = "l", bty = "n", cex.lab = 1.5, col = "blue", ylim = c(0, max(SIbwd_df$Hf, noSIbwd_df$Hf)), xaxs = "i", yaxs = "i", lwd = 2, ylab = "herbivore biomass", xlab = "fishing"))
    axis(side = 1, lwd = 2)
    axis(side = 2, lwd = 2)
    with(noSIbwd_df, points(Hf ~ ftarg, type = "l", col = "red", lwd = 2))
    title(paste("Remediating overfished reef: ", frbwd, " decrease per year", sep = ""))
    
    with(SIfwd_df, plot(Cf ~ ftarg, type = "l", lwd = 2, cex.lab = 1.5, bty = "n", col = "blue", ylim = c(0, 0.6), xaxs = "i", yaxs = "i", ylab = "coral cover", xlab = "fishing"))
    axis(side = 1, lwd = 2)
    axis(side = 2, lwd = 2)
    with(noSIfwd_df, points(Cf ~ ftarg, type = "l", lwd = 2, col = "red"))
    with(SIfwdtrait_df, points(Cf ~ ftarg, type = "l", lwd = 2, col = "purple"))
    
    with(SIbwd_df, plot(Cf ~ ftarg, type = "l", lwd = 2, cex.lab = 1.5, bty = "n", col = "blue", ylim = c(0, 0.6), xaxs = "i", yaxs = "i", ylab = "coral cover", xlab = "fishing"))
    axis(side = 1, lwd = 2)
    axis(side = 2, lwd = 2)
    with(noSIbwd_df, points(Cf ~ ftarg, type = "l", lwd = 2, col = "red"))
    
    with(SIbwdtrait_df, plot(Hf ~ ftarg, type = "l", bty = "n", cex.lab = 1.5, col = "purple", ylim = c(0, max(SIbwd_df$Hf, noSIbwd_df$Hf)), xaxs = "i", yaxs = "i", lwd = 2, ylab = "herbivore biomass", xlab = "fishing"))
    axis(side = 1, lwd = 2)
    axis(side = 2, lwd = 2)
    with(noSIbwdtrait_df, points(Hf ~ ftarg, type = "l", col = "red", lwd = 2))
    title(paste("Remediating overfished reef: ", frbwd, " decrease per year", sep = ""))
    
    with(SIbwdtrait_df, plot(Cf ~ ftarg, type = "l", lwd = 2, cex.lab = 1.5, bty = "n", col = "purple", ylim = c(0, 0.6), xaxs = "i", yaxs = "i", ylab = "coral cover", xlab = "fishing"))
    axis(side = 1, lwd = 2)
    axis(side = 2, lwd = 2)
    with(noSIbwdtrait_df, points(Cf ~ ftarg, type = "l", lwd = 2, col = "red"))
    title(paste("Remediating overfished reef: ", frbwd, " decrease per year", sep = ""))
  }
  adj <- c(prist_H = prist_H, b2fwd = b2fwd, b2bwd = b2bwd)
  return(list(SIfwd_out, noSIfwd_out, SIbwd_out, noSIbwd_out, adj, SIfwdtrait_out, SIbwdtrait_out, noSIbwdtrait_out))
}


ptm <- proc.time()
trajC002A2trait <- traj.fun.trans.trait(rC = 0.02, m = 0.008, rA = 2, alg_dom = c(0.000001, 0.999999, 0.000001, 0.5), tf = 10000, frez = 40)
proc.time() - ptm

# Hfwd SI
HfwdSI <- trajC002A2trait[[1]][, 3][which(trajC002A2trait[[1]][, 5] < 0.01)[1]]
HfwdnoSI <- trajC002A2trait[[2]][, 3][which(trajC002A2trait[[2]][, 5] < 0.01)[1]]
# Cfwd
CfwdSI <- trajC002A2trait[[1]][, 3][which(trajC002A2trait[[1]][, 7] < 0.01)[1]]
CfwdnoSI <- trajC002A2trait[[2]][, 3][which(trajC002A2trait[[2]][, 7] < 0.01)[1]]
# Hbwd
HbwdSI <- trajC002A2trait[[3]][, 3][tail(which(trajC002A2trait[[3]][, 5] >= 1.08), n = 1)] 
HbwdnoSI <- trajC002A2trait[[4]][, 3][tail(which(trajC002A2trait[[4]][, 5] >= 3.635), n = 1)] 
# Cbwd
CbwdSI <- trajC002A2trait[[3]][, 3][tail(which(trajC002A2trait[[3]][, 7] >= 0.277), n = 1)] 
CbwdnoSI <- trajC002A2trait[[4]][, 3][tail(which(trajC002A2trait[[4]][, 7] >= 0.277), n = 1)] 
#trait
HfwdSItrait <- trajC002A2trait[[6]][, 3][which(trajC002A2trait[[6]][, 5] < 0.01)[1]]
CfwdSItrait <- trajC002A2trait[[6]][, 3][which(trajC002A2trait[[6]][, 7] < 0.01)[1]]
HbwdSItrait <- trajC002A2trait[[7]][, 3][tail(which(trajC002A2trait[[7]][, 5] >= 1.08), n = 1)]
HbwdnoSItrait <- trajC002A2trait[[8]][, 3][tail(which(trajC002A2trait[[8]][, 5] >= 3.635), n = 1)] 
CbwdSItrait <- trajC002A2trait[[7]][, 3][tail(which(trajC002A2trait[[7]][, 7] >= 0.277), n = 1)] 
CbwdnoSItrait <- trajC002A2trait[[8]][, 3][tail(which(trajC002A2trait[[8]][, 7] >= 0.277), n = 1)]

dev.new()
trait_df <- data.frame(ttt = c("SI", "noSI", "SI+trait"), Hclps = c(HfwdSI, HfwdnoSI, HfwdSItrait),  Cclps = c(CfwdSI, CfwdnoSI, CfwdSItrait))
par(mfrow = c(1, 2))
with(trait_df, plot(Hclps ~ ttt, log = "y", ylim = c(0.07, 0.42)))
title("rC=0.02, rA=4; H")
with(trait_df, plot(Cclps ~ ttt, log = "y", ylim = c(0.07, 0.42)))
title("C")

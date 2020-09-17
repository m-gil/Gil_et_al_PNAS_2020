# Code below reproduces Fig. 3

require(deSolve)
source("population_models.R")

dynf.plot.fun = function(tf = 5000, xtf = 150, y1max = 1, y2max = 2.5, tstep = 0.1, add0 = 0.11, s = 1, b = 0.5, a = 1, rA = 2, rC = 0.2, alpha = 0, m = 0.08, f = 0, mu = 0.02, p = 0, q = 1, add = add0/(a*b), s2 = 0, alg_dom = c(0.2, 0.95, 0.05, 0), coral_dom = c(1, 0.17, 0.56, 0), fmax = 0.25, frate_vec = c(0.008, 0.01)) {
  SI_dynf <- list()
  noSI_dynf <- list()
  
  times = seq(0, tf, by = tstep)
  parms0 = c(add=add, a=a,b=b,q=q,f=f,rA=rA,rC=rC,m=m,mu=mu, frate = 0, fmax = 0) 
  SI_dynf0 = lsoda(coral_dom, times, SI.fun2f, parms0)
  noSIfwd_pcap = a*b*(add + SI_dynf0[length(times), 2]/(1 + q*SI_dynf0[length(times), 2]))
  b2 = noSIfwd_pcap/a
  parms20 = c(a=a, b=as.numeric(b2), f=f,rA=rA,rC=rC,m=m,mu=mu, frate = 0, fmax = 0)
  
  for (i in 1:length(frate_vec)) {
	frate = frate_vec[i]
    parms = c(add=add, a=a,b=b,q=q,f=f,rA=rA,rC=rC,m=m,mu=mu, frate = frate, fmax = fmax)
    SI_dynf[[i]] <- lsoda(SI_dynf0[length(times), 2:5], times, SI.fun2f, parms)
    noSI_dynf0 <- lsoda(coral_dom, times, noSI.fun2f, parms20)
    parms2 = c(a=a, b=as.numeric(b2), f=f,rA=rA,rC=rC,m=m,mu=mu, frate = frate, fmax = fmax)
    noSI_dynf[[i]] <- lsoda(noSI_dynf0[length(times), 2:5], times, noSI.fun2f, parms2)
  }
  
  par(mar=c(5.1, 6.1, 4.1, 8.1))
  plot(SI_dynf[[1]][ , 1], SI_dynf[[1]][ , 4], axes = F, type = "l", col = "#B5838D", xlim = c(0, xtf), ylim = c(0, 1), xlab ="", ylab = "", xaxs = "i", yaxs = "i", lty = 1, lwd = 2)
  lines(SI_dynf[[1]][ , 1], SI_dynf[[1]][ , 3], col = "#916953", lty = 1, lwd = 2)
  axis(2, ylim = c(0, max(c(SI_dynf[[1]][ , 3], SI_dynf[[1]][ , 4]))), col = "black", lwd = 2)
  mtext(2, text = "coral cover", col = "#B5838D", line = 4)
  
  par(new = T)
  plot(SI_dynf[[1]][ , 1], SI_dynf[[1]][ , 2], type = "l", axes=F, xlab=NA, ylab=NA, ylim = c(0, y2max), xlim = c(0, xtf), xaxs = "i", yaxs = "i", lty = 1, lwd = 2, col = "#628395")
  axis(4, ylim = c(0, 2.5), lwd = 2, line = 2)
  title(paste(fmax, "*biomass/year reached in ", fmax/frate_vec[1], " years", sep = ""))

  axis(1, pretty(range(c(0, xtf)), 10), lwd = 2)
  mtext("time (years)", side = 1, line=2)
  abline(v = fmax/frate_vec[1], lty = 2)
  
  dev.new()
  par(mar=c(5.1, 6.1, 4.1, 8.1))
  plot(SI_dynf[[2]][ , 1], SI_dynf[[2]][ , 4], axes = F, type = "l", lty = 1, lwd = 2, col = "#B5838D", xlim = c(0, xtf), ylim = c(0, 1), xlab ="", ylab = "", xaxs = "i", yaxs = "i")
  lines(SI_dynf[[2]][ , 1], SI_dynf[[2]][ , 3], col = "#916953", lty = 1, lwd = 2)
  axis(2, ylim = c(0, max(SI_dynf[[2]][ , 4])), col = "black", lwd = 2)
  mtext(2, text = "coral cover", col = "#B5838D", line = 4)
  
  par(new = T)
  plot(SI_dynf[[2]][ , 1], SI_dynf[[2]][ , 2], type = "l", axes=F, xlab=NA, ylab=NA, ylim = c(0, y2max), xlim = c(0, xtf), xaxs = "i", yaxs = "i", lty = 1, lwd = 2, col = "#628395")
  
  axis(4, ylim = c(0, 2.5), lwd = 2, line = 2)
  title(paste(fmax, "*biomass/year reached in ", fmax/frate_vec[2], " years", sep = ""))

  axis(1, pretty(range(c(0, xtf)), 10), lwd = 2)
  mtext("time (years)", side = 1, line=2)
  abline(v = fmax/frate_vec[2], lty = 2)
  
  return(list(SI_dynf, noSI_dynf))
}

Fig3plots <- dynf.plot.fun(tf = 250, xtf = 200, y1max = 0.85, y2max = 2.5, tstep = 0.01, rA = 2, b = 0.5, rC = 0.2, m = 0.2*0.4, mu = 0.02, fmax = 0.18, frate_vec = c(0.01, 0.02))#
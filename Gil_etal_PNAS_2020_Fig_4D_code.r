# Code below reproduces Fig. 4D
require(doSNOW)
require(foreach)
require(deSolve)
require(ggplot2)

source("population_models.R")

dynf.plot.fun = function(tf = 5000, xtf = 150, tstep = 0.1, add0 = 0.11, s = 1, b = 0.5, a = 1, rA = 2, rC = 0.2, alpha = 0, m = 0.08, f = 0, mu = 0.02, p = 0, q = 1, add = add0/(a*b), s2 = 0, alg_dom = c(0.2, 0.95, 0.05, 0), coral_dom = c(1, 0.1, 0.56, 0), fmax = 0.25, frate_vec = c(0.008, 0.01)) {
  SI_dynf <- list()
  noSI_dynf <- list()
  times = seq(0, tf, by = tstep)
  parms0 = c(add=add, s2=s2,a=a,b=b,s=s,q=q,f=f,rA=rA,alpha=alpha,rC=rC,m=m,p=p,mu=mu, frate = 0, fmax = 0) 
  SI_dynf0 = lsoda(coral_dom, times, SI.fun2f, parms0)
  noSIfwd_pcap = a*b*(add + SI_dynf0[length(times), 2]/(1 + q*SI_dynf0[length(times), 2]))
  b2 = noSIfwd_pcap/a
  parms20 = c(a=a, b=as.numeric(b2), f=f,rA=rA,alpha=alpha,rC=rC,m=m,p=p,mu=mu, frate = 0, fmax = 0)
  for (i in 1:length(frate_vec)) {
    frate = frate_vec[i]
    parms = c(add=add, s2=s2,a=a,b=b,s=s,q=q,f=f,rA=rA,alpha=alpha,rC=rC,m=m,p=p,mu=mu, frate = frate, fmax = fmax)
    SI_dynf[[i]] <- lsoda(SI_dynf0[length(times), 2:5], times, SI.fun2f, parms)
    noSI_dynf0 <- lsoda(coral_dom, times, noSI.fun2f, parms20)
    parms2 = c(a=a, b=as.numeric(b2), f=f,rA=rA,alpha=alpha,rC=rC,m=m,p=p,mu=mu, frate = frate, fmax = fmax)
    noSI_dynf[[i]] <- lsoda(noSI_dynf0[length(times), 2:5], times, noSI.fun2f, parms2)
  } 
  par(mar=c(5.1, 6.1, 4.1, 8.1))
  plot(SI_dynf[[1]][ , 1], SI_dynf[[1]][ , 4], axes = F, type = "l", col = "chocolate3", xlim = c(0, xtf), ylim = c(0, max(SI_dynf[[1]][ , 4])), xlab ="", ylab = "", xaxs = "i", yaxs = "i", lty = 2, lwd = 3)
  lines(SI_dynf[[2]][ , 1], SI_dynf[[2]][ , 4], lty = 3, lwd = 3, col = "chocolate3")
  axis(2, ylim = c(0, max(SI_dynf[[1]][ , 4])), col = "black", lwd = 2)
  mtext(2, text = "coral cover", col = "chocolate3", line = 4)
  par(new = T)
  plot(SI_dynf[[1]][ , 1], SI_dynf[[1]][ , 2], type = "l", axes=F, xlab=NA, ylab=NA, ylim = c(0, max(SI_dynf[[1]][ , 2])), xlim = c(0, xtf), xaxs = "i", yaxs = "i", lty = 2, lwd = 3)
  lines(SI_dynf[[2]][ , 1], SI_dynf[[2]][ , 2], lty = 3, lwd = 3)
  axis(4, ylim = c(0, max(SI_dynf[[1]][ , 2])), lwd = 2, line = 2)
  title(paste("Target fishing (", fmax, "*biomass/year)\nreached in ", fmax/frate_vec[1], " years", sep = ""))
  axis(1, pretty(range(c(0, xtf)), 10), lwd = 2)
  mtext("time (years)", side = 1, line=2)
  return(list(SI_dynf, noSI_dynf))
}

HC.2D.ASS.plot.fun <- function(fmax = 0.05, lab = "SI.fun2f", plotrez = 0.05, Hi_min = 0, Hi_max = 2, Ci_min = 0, Ci_max = 0.6, Ai = 0.05, tf = 5000, b = 0.5, a = 1, s = 1, s2 = 0, q = 1, rA = 2, alpha = 0, rC = 0.2, m = 0.08, p = 0, mu = 0.02, add0 = 0.11, frate = 0.01, fi = 0){ 

  cl<-makeCluster(6) # number of CPU cores
  registerDoSNOW(cl)
  
  require(ggplot2)
  require(reshape2)
  require(deSolve)

  Hi_vec = seq(from = Hi_min, to = Hi_max, by = plotrez) 
  Ci_vec = seq(from = Ci_min, to = Ci_max, by = plotrez/3) 
  s <- 1
  add <- add0/(a*b)
  SIparms <- c(s = s, s2 = s2, a = a, b = b, q=q, rA=rA, alpha=alpha, rC=rC, m=m, p=p, mu=mu, add = add, frate = frate, fmax = fmax)
  times = seq(0, tf, by = 100)
  
	outlist <- foreach(i = 1:length(Hi_vec)) %dopar% {
		require(deSolve)
		source("population_models.R")
		Hi <- Hi_vec[i]
		SI_Cf_vec <- rep(9999, length(Ci_vec))

		for (j in 1:length(Ci_vec)){
			Ci <- Ci_vec[j]
			SIout = ode(c(Hi, Ai, Ci, fi), times, SI.fun2f, SIparms, method = "adams")
			SI_Cf_vec[j] <- round(SIout[length(times), 4], 3)

		}
		SIoutmat <- cbind(fmax = rep(fmax, length(Ci_vec)), SI = rep(1, length(Ci_vec)), fwd = rep(1, length(Ci_vec)), Hi = rep(Hi, length(Ci_vec)), Ci = Ci_vec, Cf = SI_Cf_vec, rC = rC, rA = rA, frate = rep(fmax, length(Ci_vec)))
		outmat <- SIoutmat 
		return(outmat)
  }
  stopCluster(cl)
  outdf <- data.frame(do.call(rbind, outlist))
  SIoutdf <- outdf
  SIgg <- ggplot(SIoutdf, aes(Ci, Hi, z = Cf)) + geom_contour(size = 0.1, bins = 1) + ggtitle(paste("SI, fmax=", fmax))
  return(list(outdf, SIgg))
}


frate_veca = c(0.005, 0.02)
frate_vecb = c(0.01, 0.01) 

rC <- 0.2; m <- 0.08

test18a <- dynf.plot.fun(xtf = 500, rC = rC, m = m, frate_vec = frate_veca, fmax = 0.18)
test18b <- dynf.plot.fun(xtf = 500, rC = rC, m = m, frate_vec = frate_vecb, fmax = 0.18)

f18_slo_df2 <- data.frame(x = test18a[[1]][[1]][, 4], y = test18a[[1]][[1]][, 2])
f18_mid_df2 <- data.frame(x = test18b[[1]][[2]][, 4], y = test18b[[1]][[2]][, 2])
f18_fst_df2 <- data.frame(x = test18a[[1]][[2]][, 4], y = test18a[[1]][[2]][, 2])
colnames(f18_slo_df2) <- c("Ci", "Hi")
colnames(f18_mid_df2) <- c("Ci", "Hi")
colnames(f18_fst_df2) <- c("Ci", "Hi")

slo018 <- c(x = as.numeric(test18a[[1]][[1]][which(test18a[[1]][[1]][ , 5] > 0.17999 & test18a[[1]][[1]][ , 5] < 0.1801), ][1, 4]), y = as.numeric(test18a[[1]][[1]][which(test18a[[1]][[1]][ , 5] > 0.17999 & test18a[[1]][[1]][ , 5] < 0.1801), ][1, 2]))

mid018 <- c(x = as.numeric(test18b[[1]][[2]][which(test18b[[1]][[2]][ , 5] > 0.1799 & test18b[[1]][[2]][ , 5] < 0.1801), ][1, 4]), y = as.numeric(test18b[[1]][[2]][which(test18b[[1]][[2]][ , 5] > 0.1799 & test18b[[1]][[2]][ , 5] < 0.1801), ][1, 2]))

fst018 <- c(x = as.numeric(test18a[[1]][[2]][which(test18a[[1]][[2]][ , 5] > 0.1799 & test18a[[1]][[2]][ , 5] < 0.1801), ][1, 4]), y = as.numeric(test18a[[1]][[2]][which(test18a[[1]][[2]][ , 5] > 0.1799 & test18a[[1]][[2]][ , 5] < 0.1801), ][1, 2]))

pts018 <- rbind(slo018, mid018, fst018)
colnames(pts018) <- c("Ci", "Hi")

lines18x <- ggplot() +
  geom_path(data = f18_fst_df2[1:(length(f18_fst_df2$Ci) - 4), ], aes(x = Ci, y = Hi), size = 1, linetype = 1) +
  geom_path(data = f18_mid_df2[1:(length(f18_mid_df2$Ci) - 4), ], aes(x = Ci, y = Hi), size = 1, linetype = 1) +
  geom_path(data = f18_slo_df2[1:(length(f18_slo_df2$Ci) - 4), ], aes(x = Ci, y = Hi), size = 1, linetype = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(0.0, 0.6)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.25)) + 
  geom_point(aes(x = as.numeric(pts018[, 1]), y = as.numeric(pts018[, 2])), colour="blue")
dev.new()
lines18x

dev.new()
f018hirezX18 <- HC.2D.ASS.plot.fun(fmax = 0.18, plotrez = 0.01, Hi_min = 0.0, Hi_max = 2.25, Ci_min = 0.0, Ci_max = 0.6, Ai = 0.17, rC = 0.2, m = 0.08, frate = 0.18, fi = 0.18, tf = 5000)

ggplot(f018hirezX18[[1]], aes(Ci, Hi, z = Cf)) + geom_raster(aes(fill = Cf)) + scale_fill_gradientn(colours=c("white", rgb(219/255, 186/255, 197/255))) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.6)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")

# Code below reproduces Fig. 4A-C and S7

rm(list=ls())
require(rootSolve)
require(deSolve)
require(doSNOW)
require(foreach)
source("population_models.R")

nc.eq.fun <- function(f, eq = "Hnc", a = 1, b = 0.5, add = 0.22, rA = 2, mu = 0.02){
  cub <- -(a*b^2*(add + 1)^2)
  sq <- rA*(a*b*add - mu - f + a*b) - 2*a*b^2*add*(add + 1)
  uni <- 2*rA*(a*b*add - mu - f) + rA*a*b - a*b^2*add^2
  scal <- rA*(a*b*add - mu - f)
  Hnc0 <- polyroot(c(scal, uni, sq, cub))
  Hnc <- Re(Hnc0[which(Re(Hnc0) > 0)])
  Anc <- (mu + f)/(a*b*(add + Hnc/(1 + Hnc)))
  if (eq == "Hnc") return(Hnc)
  if (eq == "Anc") return(Anc)
}

in.eq.fun <- function(f, eq = "Cin", a = 1, b = 0.5, add = 0.22, rA = 2, mu = 0.02, rC = 0.2, m = 0.4*rC){ # rC = 0.2
  Hin0 <- polyroot(c(m*rA/rC, (m*rA/rC - b*add), (-b*add - b)))
  Hin <- Re(Hin0[which(Re(Hin0) > 0)])
  Ain <- (mu + f)/(a*b*(add + Hin/(1 + Hin)))
  Cin <- 1 - (mu + f)/(a*b*(add + Hin/(1 + Hin))) - m/rC
  if (eq == "Hin") return(Hin)
  if (eq == "Ain") return(Ain)
  if (eq == "Cin") return(Cin)
}

a <- 1
b <- 0.5
add <- 0.22
rA <- 2
mu <- 0.02
f <- 0
rC <- 0.02
m <- rC*0.4

fvec <- seq(from = 0, to = 0.5, length.out = 2000)

# equilibria:
ad_eqs <- cbind(f = fvec, H = rep(0, length(fvec)), A = rep(0, length(fvec)), C = rep(0, length(fvec))) # for 'all dead'
ao_eqs <- cbind(f = fvec, H = rep(0, length(fvec)), A = rep(1, length(fvec)), C = rep(0, length(fvec))) # for 'alg only'
co_eqs <- cbind(f = fvec, H = rep(0, length(fvec)), A = rep(0, length(fvec)), C = rep((1 - m/rC), length(fvec))) # for 'coral only'
# no coral:
nc_eqs <- cbind(f = fvec, H = as.vector(do.call(rbind, as.list(sapply(fvec, nc.eq.fun, eq = "Hnc", rA = rA)))), A = as.vector(do.call(rbind, as.list(sapply(fvec, nc.eq.fun, eq = "Anc", rA = rA)))), C = rep(0, length(as.vector(do.call(rbind, as.list(sapply(fvec, nc.eq.fun, eq = "Anc")))))))
# internal:
in_eqsH0 <- as.vector(do.call(rbind, as.list(sapply(fvec, in.eq.fun, eq = "Hin", rA = rA, rC = rC, m = m))))
in_eqsA0 <- as.vector(do.call(rbind, as.list(sapply(fvec, in.eq.fun, eq = "Ain", rA = rA, rC = rC, m = m))))
in_eqsC0 <- as.vector(do.call(rbind, as.list(sapply(fvec, in.eq.fun, eq = "Cin", rA = rA, rC = rC, m = m))))
# restricting internal eq case to only cases when coral (and H and A) persist:
in_eqs <- cbind(f = fvec[which(in_eqsC0 > 0)], H = in_eqsH0[which(in_eqsC0 > 0)], A = in_eqsA0[which(in_eqsC0 > 0)], C = in_eqsC0[which(in_eqsC0 > 0)])

# all eqilibria, first 'all dead', then 'algae only', then 'coral only', then, from above, 'no coral' and 'internal':
eqs <- rbind(ad_eqs, ao_eqs, co_eqs, nc_eqs, in_eqs)
stab <- rep(999999, length(eqs[, 1]))
parms0 <- c(s = 1, q = 1, p = 0, s2 = 0, alpha = 0, a = a, b = b, add = add, mu = mu, rA = rA, rC = rC, m = m, f = f)
for (i in 1:length(eqs[, 1])){
  parms0["f"] <- eqs[i, "f"]
  eigs <- eigen(jacobian.full(y = c(H = eqs[i, "H"], A = eqs[i, "A"], C = eqs[i, "C"]), func = SI.fun2, parms = parms0))$values
  lead_eig <- max(Re(eigs[abs(Im(eigs)) < 1e-10]))
  if (lead_eig < 0){
    stab[i] <- 1
  }else{
    stab[i] <- 0
  }
}
eqs <- data.frame(eqs)
eqs$stab <- stab

eqs$type <- c(rep("ad", length(fvec)), rep("ao", length(fvec)), rep("co", length(fvec)), rep("nc", length(nc_eqs[, 1])), rep("in", length(in_eqs[, 1]))) 
sub <- subset(eqs, type == "nc" & f <= 0.28 & H < 0.9)

with(subset(eqs, stab == 0), plot((H/1.8428) ~ f, col = "grey", ylim = c(-0.06, with(subset(eqs, stab == 1), max(H)/1.8428)), pch = 16, ylab = "Normalized fish biomass", xlab = ""))
title(expression(paste(italic("r")[A], " = 2", , sep = "")))
with(subset(eqs, stab == 1), points((H/1.8428) ~ f, col = "black", pch = 16))


bifur.hyst.rt.fun <- function(Ncores = 6, nondynf = "no", scanHi = "no", tf = 5000, tstep = 1, add0 = 0.11, s = 1, b = 0.5, a = 1, rA = 2, rC = 0.2, alpha = 0, m = 0.08, f = 0, mu = 0.02, p = 0, q = 1, add = add0/(a*b), s2 = 0, frate_vec = c(0.01, 0.02), fmax_rez = 4, Ci_rez = 4, fmax_lo = 0.18, fmax_hi = 0.25, fmax_vec = seq(from = fmax_lo, to = fmax_hi, length.out = fmax_rez), Ci_vec = seq(from = 0.001, to = 0.65, length.out = Ci_rez), Hi_vec = seq(from = 0.01, to = 1.85, length.out = Ci_rez), Ci = Ci0, Ai = Ai0, Hi = 1){

  times <- seq(from = 0, to = tf, by = tstep)
  
  cl<-makeCluster(Ncores) # number of CPU cores
  registerDoSNOW(cl)
  
  outlist <- foreach(rr = 1:length(frate_vec)) %dopar% {
    output <- matrix(999999, nrow = Ci_rez*fmax_rez, ncol = 9)
    source("population_models.R")
    require(deSolve)
    for (fff in 1:fmax_rez){
      for (cc in 1:Ci_rez){
        if (nondynf == "yes"){
          if (scanHi == "yes"){
            xi <- c(Hi_vec[cc], Ai, Ci, fmax_vec[fff])
          }else{
            xi <- c(Hi, Ai, Ci_vec[cc], fmax_vec[fff])
          }
        }else{
          if (scanHi == "yes"){
            xi <- c(Hi_vec[cc], Ai, Ci, f)
          }else{
            xi <- c(Hi, Ai, Ci_vec[cc], f)
          }
        }
        parms <- c(add=add, s2=s2,a=a,b=b,s=s,q=q,rA=rA,alpha=alpha,rC=rC,m=m,p=p,mu=mu, frate = frate_vec[rr], fmax = fmax_vec[fff]) 
        SI_dynf <- ode(xi, times, SI.fun2f, parms, method = "ode23")
        Afmax <- SI_dynf[which(SI_dynf[, 5] >= fmax_vec[fff])[1], 3]
        if(abs(SI_dynf[length(times), 2] - SI_dynf[(length(times) - 1), 2]) > 0.0001 | abs(SI_dynf[length(times), 3] - SI_dynf[(length(times) - 1), 3]) > 0.0001 | abs(SI_dynf[length(times), 4] - SI_dynf[(length(times) - 1), 4]) > 0.0001) print(paste("failed to equilibrate @ rr/fff/cc=", rr, fff, cc))
        
        if (scanHi == "yes"){
          output[(cc + Ci_rez*(fff - 1)), ] <- c(frate_vec[rr], fmax_vec[fff], Hi_vec[cc], Ai, Ci, SI_dynf[length(times), 2], SI_dynf[length(times), 3], SI_dynf[length(times), 4], Afmax)
        }else{
          output[(cc + Ci_rez*(fff - 1)), ] <- c(frate_vec[rr], fmax_vec[fff], Hi, Ai, Ci_vec[cc], SI_dynf[length(times), 2], SI_dynf[length(times), 3], SI_dynf[length(times), 4], Afmax)
        }
      }
    }
    return(output)
  }
  stopCluster(cl)
  output2 <- do.call(rbind, outlist)
  if (scanHi == "yes"){
    colnames(output2) <- c("frate", "fmax", "Hi", "Ai", "Ci", "Hf", "Af", "Cf", "Afmax")
  }else{
    colnames(output2) <- c("frate", "fmax", "Hi", "Ai", "Ci", "Hf", "Af", "Cf", "Afmax")
  }
  return(output2)
}

# Fig. 4A-C:
fmax_lo0 <- 0
fmax_hi0 <- 0.3
frate_vec0 <- c(0.001, 0.01, 0.1)
Ci0 <- 0.56
Ai0 <- 0.17
scanHi0 <- "yes"

scanHi = scanHi0; frate_vec = frate_vec0; fmax_lo = fmax_lo0; fmax_hi = fmax_hi0; fmax_rez = 4; Ci_rez = 4; Ci = Ci0; Ai = Ai0
Hscan5 <- bifur.hyst.rt.fun(scanHi = scanHi0, frate_vec = frate_vec0, fmax_lo = fmax_lo0, fmax_hi = fmax_hi0, fmax_rez = 400, Ci_rez = 400, Ci = Ci0, Ai = Ai0, tf = 5000, tstep = 100)

outdf2 <- data.frame(Hscan5)
outdf2$Clive <- as.numeric(outdf2$Cf > 0.006)
outdf2$Hlive <- as.numeric(outdf2$Hf > 0.018)

dev.new()
par(mfrow = c(1, 3))
for (pp in 1:length(frate_vec0)){
	with(subset(outdf2, outdf2$frate == unique(outdf2$frate)[pp] & outdf2$Hlive == 1), plot(Hi ~ fmax, pch = 15, cex = 1, xlim = c(0.08, 0.30), col = "blue", ylim = c(-0.005, 1.85), xaxs = "i", yaxs = "i", xlab = "target fishing level (f @ equilibrium)", ylab = "initial H biomass"))
	with(subset(eqs, stab == 1 & H >= 1.84), lines(H ~ f, lwd = 4))
	with(subset(eqs, stab == 1 & H < 1.84 & H > 0.1), lines(H ~ f, lwd = 4))
	with(subset(eqs, stab == 1 & H < 1), lines(H ~ f, lwd = 4))
	with(subset(eqs, stab == 0 & H > 0 & H < 1.2), lines(H ~ f, lty = 2, lwd = 4))
	title(paste("fishing increase: +", unique(outdf2$frate)[pp], " biomass per yr"))
}

# Fig. S7:
rbPal <- colorRampPalette(c('yellow','brown'))
outdf2$Col <- rbPal(100)[as.numeric(cut(outdf2$Afmax, breaks = 100))]

dev.new()
par(mfrow = c(1, length(frate_vec0)))
for (pp in 1:length(frate_vec0)){
  with(subset(outdf2, outdf2$frate == unique(outdf2$frate)[pp] & outdf2$Hlive == 1), plot(Hi ~ fmax, pch = 15, cex = 1, xlim = c(0.08, 0.30), col = Col, ylim = c(-0.005, 1.85), xaxs = "i", yaxs = "i", xlab = "target fishing level (f @ equilibrium)", ylab = "initial H biomass"))
  with(subset(eqs, stab == 1 & H >= 1.84), lines(H ~ f, lwd = 4))
  with(subset(eqs, stab == 1 & H < 1.84 & H > 0.1), lines(H ~ f, lwd = 4))
  with(subset(eqs, stab == 1 & H < 1), lines(H ~ f, lwd = 4))
  with(subset(eqs, stab == 0 & H > 0 & H < 1.2), lines(H ~ f, lty = 2, lwd = 4))
  title(paste("frate:", unique(outdf2$frate)[pp], ", Ci Ai:,", Ci0, Ai0))
}

dev.new()
my.colors = colorRampPalette(c("yellow", "brown"))
z=matrix(1:100,nrow=1)
x=1
y=seq(range(outdf2$Afmax)[1], range(outdf2$Afmax)[2], len=100)
image(x,y,z,col=my.colors(100),axes=FALSE,xlab="",ylab="")
axis(2)

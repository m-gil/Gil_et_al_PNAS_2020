# Code below reproduces main text Fig. 2A and SI Fig. S1.
# Note: I, O, X and P below correspond to System 1 state variables P, U, E and L in the Main Text, respectively.

require(deSolve)
load("bestfitparms.Rda")
parms <- as.vector(bestfitparms[1:9])
i2 <- bestfitparms[10]
h2 <- bestfitparms[11]
prist_H <- 1.8428
Oentext <- 0.149
democurve.fun <- function(X) 0.11 + 0.5*X/(1 + X)
Hs <- c(seq(from = 0.01, to = 1, length.out = 10), seq(from = 1.01, to = 2, length.out = 10)) 
tt <- seq(from = 0, to = 1000, by = 1)
Ivec <- rep(999999, length(Hs))
arrt <- Ivec
Edwell <- Ivec

behavODE.fun <- function(t, x, pars){
  with(as.list(c(pars, x)),{
    lam1 = pars[1]; alphlam = pars[2]; mu1 = pars[3]; etamu = pars[4]; alphmu = pars[5]; etaphi = pars[6]; tauX = pars[7]; tauY = pars[8]; mindI = pars[9]
    I <- x[1]; O <- x[2]; X <- x[3]; Y <- x[4]
    dI <- O*(lam1 + alphlam*X) + mindI - I*(mu1*I^etamu + alphmu*I^etaphi*Y)    
    dO <- I*(mu1*I^etamu + alphmu*I^etaphi*Y) - (O*(lam1 + alphlam*X) + mindI)
    dX <- -1/tauX*X + O*(lam1 + alphlam*X)
    dY <- -1/tauY*Y + I*(mu1*I^etamu + alphmu*I^etaphi*Y)
    dx <- c(dI, dO, dX, dY)
    return(list(dx))
  })
}

for (h in 1:length(Hs)){
  H <- Hs[h]
  x <- c(H*0.5, H*0.5, H*0.1, H*0.1)
  out <- ode(x, tt, behavODE.fun, parms, method = "lsoda", maxsteps = 50000)
  if (abs(out[length(out[, 1]), 2] - out[length(out[, 1]) - 1, 2]) > 0.0001 | abs(out[length(out[, 1]), 3] - out[length(out[, 1]) - 1, 3]) > 0.0001 | abs(out[length(out[, 1]), 4] - out[length(out[, 1]) - 1, 4]) > 0.0001 | abs(out[length(out[, 1]), 5] - out[length(out[, 1]) - 1, 5]) > 0.0001 | is.na(out[length(out[, 1]), 4]) | is.na(out[length(out[, 1]), 5])){
    print(paste("h=", h, ", failed to equilibrate", sep = ""))
    Ivec[h] <- 9999 
    arrt[h] <- 1 
    if (Ivec[h] >= 4.5 & Ivec[h] <= 6){
      I <- 0; O <- 0; X <- 0; Y <- 0
    }
  }else{
    Ivec[h] <- out[length(out[, 1]), 2]
    arrt[h] <- out[length(out[, 1]), 3]*(parms[1] + parms[2]*out[length(out[, 1]), 4])
    if (Ivec[h] >= 4.5 & Ivec[h] <= 6.5){
      if (h == 1){
        O <- out[length(out[, 1]), 3]; 
        X <- out[length(out[, 1]), 4];
      }
      if (abs(Ivec[h] - 5.6) < abs(Ivec[h - 1] - 5.6)){
        O <- out[length(out[, 1]), 3]; 
        X <- out[length(out[, 1]), 4];
        h2 <- H
        i2 <- Ivec[h]
      }
      Eentext <- O*(lam1 + alphlam*X)  
      RS <- (Oentext - Eentext)^2 
    }
  }
}

if (length(which(Ivec <= 11)) < 1){
  print(paste("no Ivec values <=11", sep = ""))
  RSS <- 99999
}else{
  Hs2 <- Hs[which(Hs <= 2)]
  Epcap <- 0.63*0.856*Ivec[1:length(Hs2)]*0.000171*10.65*(i2/h2)*3600/(Hs2*2) 
  Opcap <- democurve.fun(Hs2)
}

colmod <- rgb(0, 0, 255, max = 255, alpha = 0.5)

par(mfrow = c(1, 1))
plot(Opcap ~ Hs2, ylim = c(0, 0.5), type = "l", lwd = 4, xlab = "H", ylab = "fish per-capita feeding rate (lambda)", yaxs = "i", xaxs = "i", col = rgb(0, 0, 0, max = 255, alpha = 255))
lines(Epcap ~ Hs2, lwd = 4, col = rgb(180, 180, 180, max = 255, alpha = 255))

behavODE.fun_trait <- function(t, x, pars){
  with(as.list(c(pars, x)),{
    
    lam1 = pars[1]; alphlam = pars[2]; mu1 = pars[3]; etamu = pars[4]; alphmu = pars[5]; etaphi = pars[6]; tauX = pars[7]; tauY = pars[8]; mindI = pars[9]; Hf0 = pars[10]
    
    I <- x[1]; O <- x[2]; X <- x[3]; Y <- x[4]
    I = (I >= 0.00001)*I; O = (O >= 0.00001)*O; X = (X >= 0.00001)*X; Y = (Y >= 0.00001)*Y;
    
    dI <- O*(lam1 + alphlam*X) + mindI - I*((mu1 + (mu1/Hf0)*(Hf0 - (I + O)))*I^etamu + alphmu*I^etaphi*Y)    
    dO <- I*((mu1 + (mu1/Hf0)*(Hf0 - (I + O)))*I^etamu + alphmu*I^etaphi*Y) - (O*(lam1 + alphlam*X) + mindI)
    dX <- -1/tauX*X + O*(lam1 + alphlam*X)
    dY <- -1/tauY*Y + I*((mu1 + (mu1/Hf0)*(Hf0 - (I + O)))*I^etamu + alphmu*I^etaphi*Y)
    
    dx <- c(dI, dO, dX, dY)
    return(list(dx))
  })
}

parms[10] <- prist_H 
Hs <- c(seq(from = 0.01, to = 1, length.out = 10), seq(from = 1.01, to = prist_H, length.out = 10))
tt <- seq(from = 0, to = 5000, by = 100)
Ivec <- rep(999999, length(Hs))
arrt <- Ivec
Edwell <- Ivec


lam1 <- parms[1]; alphlam <- parms[2]
RS <- 100

for (h in 1:length(Hs)){
  H <- Hs[h]
  x <- c(H*0.5, H*0.5, H*0.1, H*0.1)
  out <- ode(x, tt, behavODE.fun_trait, parms, method = "lsoda", maxsteps = 50000)
  if (abs(out[length(out[, 1]), 2] - out[length(out[, 1]) - 1, 2]) > 0.0001 | abs(out[length(out[, 1]), 3] - out[length(out[, 1]) - 1, 3]) > 0.0001 | abs(out[length(out[, 1]), 4] - out[length(out[, 1]) - 1, 4]) > 0.0001 | abs(out[length(out[, 1]), 5] - out[length(out[, 1]) - 1, 5]) > 0.0001 | is.na(out[length(out[, 1]), 4]) | is.na(out[length(out[, 1]), 5])){
    print(paste("h=", h, ", failed to equilibrate", sep = ""))
    Ivec[h] <- 9999
    arrt[h] <- 1 

    if (Ivec[h] >= 4.5 & Ivec[h] <= 6){
      I <- 0; O <- 0; X <- 0; Y <- 0
    }
  }else{
    Ivec[h] <- out[length(out[, 1]), 2]
    arrt[h] <- out[length(out[, 1]), 3]*(parms[1] + parms[2]*out[length(out[, 1]), 4])
    if (Ivec[h] >= 4.5 & Ivec[h] <= 6.5){
      if (h == 1){
        O <- out[length(out[, 1]), 3]; 
        X <- out[length(out[, 1]), 4];
      }
      if (abs(Ivec[h] - 5.6) < abs(Ivec[h - 1] - 5.6)){
        O <- out[length(out[, 1]), 3]; 
        X <- out[length(out[, 1]), 4];
      }
      Eentext <- O*(lam1 + alphlam*X)  
      RS <- (Oentext - Eentext)^2 
    }
  }
}
if (length(which(Ivec <= 11)) < 1){
  print(paste("no Ivec values <=11", sep = ""))
  RSS <- 99999
}else{
  Hs2 <- Hs[which(Hs <= 5.6)]
  Epcap <- 0.63*0.856*Ivec[1:length(Hs2)]*0.000171*10.65*(i2/h2)*3600/(Hs2*2)
  Opcap <- democurve.fun(Hs2)
}

lines(Epcap[1:(max(which(Hs2 <= prist_H)))]  ~ Hs2[1:(max(which(Hs2 <= prist_H)))], lwd = 4, col = rgb(0, 0, 0, max = 255, alpha = 255), lty = 1)
mod3b <- lm(Epcap[c(1, (max(which(Hs2 <= prist_H))))] ~ Hs2[c(1, (max(which(Hs2 <= prist_H))))])
abline(a = coef(mod3b)[1], b = coef(mod3b)[2], lwd = 4, col = rgb(180, 180, 180, max = 255, alpha = 255), lty = 1)
abline(h = Epcap[max(which(Hs2 <= prist_H))], lty = 2)
title(paste("prist_H = ", prist_H, sep = ""))
abline(v = 1.8428, lty = 2)

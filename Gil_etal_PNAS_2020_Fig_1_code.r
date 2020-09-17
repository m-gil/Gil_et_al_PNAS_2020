# Code below reproduces Fig. 1B

load("bestfitparms.Rda")
i2 <- bestfitparms[10]
h2 <- bestfitparms[11]

require(deSolve)
Oentext <- 0.149
fig1b.fun <- function(X) 13.0428*(X + 0.5)^0.5

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

Hs <- c(seq(from = 0.01, to = 1, length.out = 10), seq(from = 1, to = 20, length.out = 40))
tt <- seq(from = 0, to = 1000, by = 100)
Ivec <- rep(999999, length(Hs))
arrt <- Ivec
Edwell <- Ivec
parms <- bestfitparms[1:9]
lam1 <- parms[1]; alphlam <- parms[2]

for (h in 1:length(Hs)){
  H <- Hs[h]
  x <- c(H*0.5, H*0.5, H*0.1, H*0.1)
  out <- ode(x, tt, behavODE.fun, parms, method = "lsoda", maxsteps = 50000)
  if (abs(out[length(out[, 1]), 2] - out[length(out[, 1]) - 1, 2]) > 0.0001 | abs(out[length(out[, 1]), 3] - out[length(out[, 1]) - 1, 3]) > 0.0001 | abs(out[length(out[, 1]), 4] - out[length(out[, 1]) - 1, 4]) > 0.0001 | abs(out[length(out[, 1]), 5] - out[length(out[, 1]) - 1, 5]) > 0.0001 | is.na(out[length(out[, 1]), 4]) | is.na(out[length(out[, 1]), 5])){
    print(paste("h=", h, ", failed to equilibrate", sep = ""))
    Ivec[h] <- 9999
    arrt[h] <- 1 
    print(h)
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
  rng <- seq(1:max(which(Ivec <= 12)))
  if (length(rng) > 25){
    drops0 <- which(seq(1:((length(rng) - 25)*2)) %% 2 == 0)
    if (median(drops0) <= 13){
      add <- 13 - median(drops0)
    }else{
      add <- max(rng) - max(drops0)
    }
    drops <- drops0 + add
    rng <- rng[-drops]
  }
  Ivec1 <- Ivec[rng]	
  Edwell <- Ivec1/arrt[rng]
  Odwell <- fig1b.fun(Ivec1)
}

# loading data
dwelltimes_ctrls = read.csv("dwelltimes_data_ctrls.csv")
dwelltimes_ctrls$pcap_bites <- dwelltimes_ctrls$dwelltime*0.63
dwelltimes_ctrls$mgC <- dwelltimes_ctrls$pcap_bites*0.171
dwelltimes_ctrls$avgnum2 <- dwelltimes_ctrls$avgnum/12

# plotting model
plot(dwelltimes_ctrls$mgC ~ dwelltimes_ctrls$avgnum2, cex = 3, ylim = c(0, 8), ylab = "mean per capita algae eaten (mgC/bout)", xlab = "mean density of fish (fish/m^2)")
Ivec2 <- Ivec1/12
lines(Edwell*0.171*0.63 ~ Ivec2, col = "red")

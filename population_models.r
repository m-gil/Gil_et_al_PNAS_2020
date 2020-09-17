# population models

SI.fun2 = function(t,y,parms) {
  with(as.list(c(parms,y)),{
    H=y[1]; A=y[2]; C=y[3];
    dH = a*b*(add + H/(1 + q*H))*A*H - mu*H - f*H
    dA = rA*A*(1 - (A + C)) - b*(add + H/(1 + q*H))*A*H
    dC = rC*C*(1 - (C + A)) - m*C
    dY = c(dH, dA, dC);
    return(list(dY));
  })
}

noSI.fun2 = function(t,y,parms) {
  with(as.list(c(parms,y)), {
    H=y[1]; A=y[2]; C=y[3];
    dH = a*b*A*H - mu*H - f*H
    dA = rA*A*(1 - (A + C)) - b*A*H
    dC = rC*C*(1 - (C + A)) - m*C
    dY = c(dH, dA, dC);
    return(list(dY));
  })
}

SI.fun2f = function(t,y,parms) {
  with(as.list(c(parms,y)),{
    H=y[1]; A=y[2]; C=y[3]; f=y[4]
    if (f < fmax){
      df = frate  
    } else {
      df = 0
    }
    dH = a*b*(add + H/(1 + q*H))*A*H - mu*H - f*H
    dA = rA*A*(1 - (A + C)) - b*(add + H/(1 + q*H))*A*H
    dC = rC*C*(1 - (C + A)) - m*C
    dY = c(dH, dA, dC, df)
    return(list(dY))
  })
}

SI.fun2f.trait <- function(t,y,parms) {
  with(as.list(c(parms,y)),{
    H=y[1]; A=y[2]; C=y[3]; f=y[4]
    if (f < fmax){
      df = frate  
    } else {
      df = 0
    }
    dH = a*(d1 + d2*H)*A*H - mu*H - f*H
    dA = rA*A*(1 - (A + C)) - (d1 + d2*H)*A*H
    dC = rC*C*(1 - (C + A)) - m*C
    dY = c(dH, dA, dC, df)
    return(list(dY))
  })
}

SI.fun2fBWD.trait <- function(t,y,parms) {
  with(as.list(c(parms,y)),{
    H=y[1]; A=y[2]; C=y[3]; f=y[4]
    if (f > fmin){
      df = frate  
    } else {
      df = 0
    }
    dH = a*(d1 + d2*H)*A*H - mu*H - f*H
    dA = rA*A*(1 - (A + C)) - (d1 + d2*H)*A*H
    dC = rC*C*(1 - (C + A)) - m*C
    dY = c(dH, dA, dC, df)
    return(list(dY))
  })
}

noSI.fun2f = function(t,y,parms) {
  with(as.list(c(parms,y)), {
    H=y[1]; A=y[2]; C=y[3]; f=y[4]
    if (f < fmax){
      df = frate  
    } else {
      df = 0
    }
    dH = a*b*A*H - mu*H - f*H
    dA = rA*A*(1 - (A + C)) - b*A*H
    dC = rC*C*(1 - (C + A)) - m*C
    dY = c(dH, dA, dC, df)
    return(list(dY))
  })
}

SI.fun2fBWD = function(t,y,parms) {
  with(as.list(c(parms, y)),{
    H=y[1]; A=y[2]; C=y[3]; f=y[4]
    if (f > fmin){
      df = frate  
    } else {
      df = 0
    }
    dH = a*b*(add + H/(1 + q*H))*A*H - mu*H - f*H
    dA = rA*A*(1 - (A + C)) - b*(add + H/(1 + q*H))*A*H
    dC = rC*C*(1 - (C + A)) - m*C
    dY = c(dH, dA, dC, df)
    return(list(dY))
  })
}

noSI.fun2fBWD = function(t,y,parms) {
  with(as.list(c(parms,y)), {
    H=y[1]; A=y[2]; C=y[3]; f=y[4]
    if (f > fmin){
      df = frate  
    } else {
      df = 0
    }
    dH = a*b*A*H - mu*H - f*H
    dA = rA*A*(1 - (A + C)) - b*A*H
    dC = rC*C*(1 - (C + A)) - m*C
    dY = c(dH, dA, dC, df)
    return(list(dY))
  })
}
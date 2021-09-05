data_generation = function(n, j){
  sx = seq(0, 1, length.out = j)
  sy = seq(0, 1, length.out = j)
  
  np = 5
  mindex = c(1,2,3)
  
  intmat = matrix(0, np+np*(np-1)/2,2);
  k = 1
  for(ii in 1:np){
    for(ij in ii:np){
      intmat[k,1] = ii;
      intmat[k,2] = ij;
      k = k+1;
    }
  }
  
  qindex = c(1,3,10)
  tintmat = intmat[qindex,]
  
  fX = list()
  for(ij in 1:5){
    
    ksi = list()
    for(ik in 1:5){
      if(ij %in% mindex){
        ksi[[ik]] = rnorm(n, 4, sd = (1*ik^(-1)))
      }else{
        ksi[[ik]] = 2*rnorm(n, 0, sd = (4*ik^(-4)))
      }
    }
    
    phi = list()
    for(ik in 1:5){
      if(ij %in% mindex){
        phi[[ik]] = sin(ik * pi * sx) * cos(ik * pi * sx)
      }else{
        phi[[ik]] = 6*cos(ik * pi * sx) * cos(ik * pi * sx)
      }
    }
    
    fX[[ij]] = Reduce("+", lapply(1:5, function(k){ksi[[k]] %*% t(phi[[k]])}))
  }
  
  fBeta = list()
  fBeta[[1]] = function(s,t){
    5*sin(6*pi*s)*cos(pi*t)
  }
  
  fBeta[[2]] = function(s,t){
    4*sin(8*pi*s)*cos(6*pi*t)
  }
  
  fBeta[[3]] = function(s,t){
    2*cos(6*pi*s)*cos(pi*t)
  }
  
  fBeta[[4]] = function(s,t){
    4*(s - 0.5)^2 * (t - 0.5)^2
  }
  
  fBeta[[5]] = function(s,t){
    5 * sqrt(s) * sqrt(t)
  }
  
  vBeta = lapply(1:length(fBeta), function(k){outer(sx, sy, fBeta[[k]])})
  
  
  fGamma = list()
  fGamma[[1]] = function(s1,s2,t){
    8*exp(s1^2-s2^2+t)
  }
  
  fGamma[[2]] = function(s1,s2,t){
    8*(2*s1^2-4*s2^4+t)
  }
  
  
  fGamma[[3]] = function(s1,s2,t)  {
    8*(2*s1-s2+4*t)
  }
  
  vGamma=list()
  for(gj in 1:length(qindex)){
    vGamma[[gj]]=array(0, c(length(sx), length(sx), length(sy)))
    for(gk in 1:length(sy)){
      
      vGamma[[gj]][,,gk]=outer(sx,sx,function(x,y){fGamma[[gj]](x,y,sy[gk])})
    }
  }
  
  fY = Reduce("+", lapply(1:length(mindex), function(k){fX[[mindex[k]]] %*% vBeta[[mindex[k]]] / j}))
  
  fY2 = 0
  
  if(length(qindex)!=0){
    for(ji in 1:length(qindex)){
      fY2 = fY2+sapply(1:length(sy), function(k){diag((fX[[tintmat[ji,1]]]%*%outer(sx,sx,function(x,y){fGamma[[ji]](x,y,sy[k])}))%*%t(fX[[tintmat[ji,2]]]))})/length(sx)^2
    }
  }
  
  
  library(goffda)
  
  err = r_ou(n=n, t = sy, mu=0, alpha = 1, sigma = 1,
             x0=rnorm(n=n, mean = 0, sd = 1/sqrt(2*1)))$data
  
  fYe = fY + fY2 + err
  
  return(list(Y = fYe, X = fX, main_coefs = vBeta, quad_coefs = vGamma))
  
}


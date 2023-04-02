data_generation <- function(n, gpy, gpx, sd.error, out.p = 0){
  
  out.indx = NULL
  Xflist = list()
  
  for(ijk in 1:6){
    Xflist[[ijk]] = function(gp){
      Xp1 = function(gp){
        p1 = matrix(0, 10, length(gp))
        for(i in 1:10)
          p1[i,] = (i^{-2}) * sqrt(2) * sin(i * pi * gp)
        return(p1)
      }
      
      Xp2 = function(gp){
        p2 = matrix(0, 10, length(gp))
        for(i in 1:10)
          p2[i,] = (i^{-2}) * sqrt(2) * cos(i * pi * gp)
        return(p2)
      }
      
      Xi = rnorm(20, 0, 0.5)
      Xii = Xi * rbind(Xp1(gp), Xp2(gp))
      Xp = colSums(Xii)
      return(Xp)
    }
  }
  
  beta = list()
  beta[[1]] = function(t, s) {
    2*sin(2*pi*s) + sin(pi*t)
  }
  beta[[2]] = function(t, s) {
    2*cos(pi*s) + cos(2*pi*t)
  }
  beta[[3]] = function(t, s) {
    cos(pi*s) + sin(2*pi*t)
  }
  beta[[4]] = function(t, s) {
    2*sin(2*pi*s) + sin(2*pi*t)
  }

  
  s = gpx
  t = gpy
  ngpx <- length(gpx)
  ngpy <- length(gpy)
  
  alpha <- function(t) 4*cos(4*pi*t)
  beta.ts = list()
  for(ij in 1:4)
    beta.ts[[ij]] = outer(t, s, beta[[ij]])
  
  Xf = list()
  for(ij in 1:6)
    Xf[[ij]] = t(replicate(n, Xflist[[ij]](s)))
  
  Xbeta.t = list()
  for(ij in 1:4)
    Xbeta.t[[ij]] = Xf[[ij]] %*% t(beta.ts[[ij]]) * (s[2] - s[1])
  
  alpha.t <- matrix(alpha(t), nrow = n, ncol = ngpy, byrow = TRUE)
  
  Xbeta = Reduce("+", Xbeta.t) + alpha.t
  
  np = 6
  mindex = c(1,2,3,4)
  
  intmat = matrix(0, np+np*(np-1)/2,2);
  k = 1
  for(ii in 1:np){
    for(ij in ii:np){
      intmat[k,1] = ii;
      intmat[k,2] = ij;
      k = k+1;
    }
  }
  
  qindex = c(1,4,12,16)
  tintmat = intmat[qindex,]
  
  gamma = list()
  gamma[[1]] = function(s1,s2,t){
    exp((s1^2+s2^2))*t^0.5
  }
  
  gamma[[2]] = function(s1,s2,t){
    2*cos(pi*(s1+s2))*t^0.5
  }
  
  gamma[[3]] = function(s1,s2,t)  {
    (s1*s2+t^2)
  }
  
  gamma[[4]] = function(s1,s2,t)  {
    (s1^2+s2^2)*t
  }
  
  gamma.rst=list()
  for(gj in 1:length(qindex)){
    gamma.rst[[gj]]=array(0, c(length(s), length(s), length(t)))
    for(gk in 1:length(t)){
      
      gamma.rst[[gj]][,,gk]=outer(s,s,function(x,y){gamma[[gj]](x,y,t[gk])})
    }
  }
  
  Yqi = list()
  if(length(qindex)!=0){
    for(ji in 1:length(qindex)){
      Yqi[[ji]] = sapply(1:length(t), 
                         function(k){diag((Xf[[tintmat[ji,1]]] %*% 
                                             outer(s,s,function(x,y){gamma[[ji]](x,y,t[k])}))%*%t(Xf[[tintmat[ji,2]]]))})*(s[2] - s[1])^2
    }
  }
  
  Yq = Reduce("+", Yqi)
  Yt = Xbeta + Yq
  eps <- matrix(rnorm(n * ngpy, 0, sd.error), n, ngpy)
  Yf = Xbeta + Yq + eps
  
  if(out.p > 0){
    
    Xflist_out = list()
    
    for(ijk in 1:2){
      Xflist_out[[ijk]] = function(gp){
        Xp1 = function(gp){
          p1 = matrix(0, 10, length(gp))
          for(i in 1:10)
            p1[i,] = (i^{-2}) * sqrt(6) * cos(i * pi * gp)
          return(p1)
        }
        
        Xp2 = function(gp){
          p2 = matrix(0, 10, length(gp))
          for(i in 1:10)
            p2[i,] = (i^{-2}) * sqrt(6) * cos(i * pi * gp)
          return(p2)
        }
        
        Xi = rnorm(20, 0, 0.5)
        Xii = Xi * rbind(Xp1(gp), Xp2(gp))
        Xp = colSums(Xii)
        return(Xp)
      }
    }
    
    Xf_out1 = list()
    for(ij in 1:2)
      Xf_out1[[ij]] = t(replicate(n, Xflist_out[[ij]](s)))
    
    Xf_out = Xf
    Xf_out[[1]] = Xf_out1[[1]]
    Xf_out[[2]] = Xf_out1[[2]]
    
    Xbeta.t_out = list()
    for(ij in 1:4)
      Xbeta.t_out[[ij]] = Xf_out[[ij]] %*% t(beta.ts[[ij]]) * (s[2] - s[1])
    
    Xbeta_out = Reduce("+", Xbeta.t_out) + alpha.t
    
    Yqi_out = list()
    if(length(qindex)!=0){
      for(ji in 1:length(qindex)){
        Yqi_out[[ji]] = sapply(1:length(t), 
                           function(k){diag((Xf_out[[tintmat[ji,1]]] %*% 
                                               outer(s,s,function(x,y){gamma[[ji]](x,y,t[k])}))%*%t(Xf_out[[tintmat[ji,2]]]))})*(s[2] - s[1])^2
      }
    }
    
    Yq_out = Reduce("+", Yqi_out)
    Yt_out = Xbeta_out + Yq_out
    eps_out <- matrix(rnorm(n * ngpy, 10, sd.error), n, ngpy)
    Yf_out = Xbeta_out + Yq_out + eps_out
    
    nout = round((n-100) * out.p)
    out.indx = sample(1:n, nout)
    
    Yt[out.indx,] = Yt_out[out.indx,]
    Yf[out.indx,] = Yf_out[out.indx,]
    Xf[[1]][out.indx,] = Xf_out[[1]][out.indx,]
    Xf[[2]][out.indx,] = Xf_out[[2]][out.indx,]
  }
  
  all.indx = 1:n
  test.indx = sample(all.indx[!all.indx %in% out.indx], 100)
  train.indx = all.indx[!all.indx %in% test.indx]
  
  Y_train = Yf[-test.indx,]
  Yt_train = Yt[-test.indx,]
  Y_test = Yf[test.indx,]
  Yt_test = Yt[test.indx,]
  X_train = X_test = list()
  for(ix in 1:6){
    X_train[[ix]] = Xf[[ix]][-test.indx,]
    X_test[[ix]] = Xf[[ix]][test.indx,]
  }
  
  
  oindx = which(!match(train.indx, out.indx) %in% NA)
  
  return(list(Y_train = Y_train, Yt_train = Yt_train, Y_test = Y_test, Yt_test = Yt_test,
              X_train = X_train, X_test = X_test, main_coefs = beta.ts, quad_coefs = gamma.rst, out.indx = oindx))
}

tauest = function(X, Y, c1, c2, ka, N){
  
  rhobi = function(u, ccc){
    w = abs(u) <= ccc
    v = (u^2/(2)*(1-(u^2/(ccc^2))+(u^4/(3*ccc^4))))*w +(1-w)*(ccc^2/6)
    v = v*6/ccc^2
    return(v)
  }

  
  psibi = function(u, ccc){
    w = abs(u) <= ccc
    v = w*u*((1-(u/ccc)^2)^2)
    v = v*6/ccc^2
    return(v)
  }

  scale1 = function(X, k, cgl, tol, q){
    bbb = rhobi(k, k)/2
    s0 = median(abs(X))
    if (q==5)  
    {
      s0 = s0/2.1
    }
    if 	(q==2)
    {
      s0 = s0/1.15 
    }
    if (q==1)
    {
      s0 = s0/0.6745
    }
    s = s0
    if (s0>.00001)
    {   
      gamma = 1
    }else
    {
      gamma = 0
    }
    ii = 0
    while ((gamma >tol)|(ii<100))
    {
      ii = ii+1
      s = cgl*s0^2*mean(rhobi(X/s0, k))/bbb
      s = sqrt(s)
      gamma = abs(s-s0)/s0
      s0 = s
    }
    return(s0)
  }

  tau = function(X, c1, c2, ka, cgl, tol, scal){
    if (scal==-1){
      s = scale1(X, c1, cgl, tol, 1)
    }else{
      s = scal
    }
    t = (1/(6*ka))*(s^2)*(c2^2)*mean(rhobi(X/s, c2))
    t = sqrt(t)
    list(t=t, s=s)
  }

  weight = function(d, c1, c2, s){
    n = prod(dim(d))
    d = d/s
    d = pmax(abs(d),.000000001)
    A2 = sum(rhobi(d, c2))
    B1 = sum(psibi(d, c1)*d)
    B2 = sum(psibi(d, c2)*d)
    CC1 = (2*A2-B2)/n
    CC2 = B1/n
    w = CC1*psibi(d, c1)+CC2*psibi(d, c2)
    w = w/d
    list(w=w)
  }
  
  X = as.matrix(X)
  Y = as.matrix(Y)
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)

  n1 = floor(n/2)
  mah = rep(0, n)
  
  for(i in 1:N){   
    ff = 0
    itertest =  0
    while(ff==0 && (itertest<10000)){
      z = runif(n)
      v = sort(z, index.return = TRUE)$ix
      v = v[1:(p+q)]
      X1 = X[v,]
      Y1 = Y[v,]
      z = cbind(X1, Y1)
      co = rcond(z)
      if (co>1e-10){
        ff = 1
      }
      itertest = itertest + 1
      if (itertest==10000) stop("too many degenerate subsamples")
    }
    be = t(Y1)%*%X1%*%ginv(t(X1)%*%X1)
    Ye = X%*%t(be)
    res = as.matrix(Y-Ye)
    cov = t(res[v,])%*%res[v,]
    cov = cov/(det(cov))^(1/q)

    covin=solve(cov)
    for (j in 1:n){
      mah[j]=(res[j,]%*%covin%*%as.matrix(res[j,]))^(1/2);
    }

    so2 = sort(mah, index.return = TRUE)$ix
    X2 = X[so2[1:n1],]
    Y2 = Y[so2[1:n1],]
    
    be = t(Y2)%*%X2%*%ginv(t(X2)%*%X2)
    Ye = X%*%t(be)
    res = as.matrix(Y-Ye)
    cov = t(res[v,])%*%res[v,]
    cov = cov/(det(cov))^(1/q)

    covin = solve(cov)
    for (j in 1:n){
      mah[j]=(res[j,]%*%covin%*%as.matrix(res[j,]))^(1/2);
    }
    if(i>1){
      vv = sum(mah<obj)
    }else{
      vv = n
    }
    if((vv>n/2)){
      obj = median(mah)
      be0 = be
      cov0 = cov
      mah0 = mah
      res0 = res
    }
  }
  beinit = be0

  aa = 0
  ii = 0
  X1 = X
  Y1 = Y
  while(aa==0){
    ii = ii+1
    tscale = tau(mah0, c1, c2, ka, 1, .0001, -1)
    tt = tscale$t
    ss = tscale$s
    u = (tt^2)

    cov0 = u*cov0
    obj = cov0
    mah0 = mah0/sqrt(u)
    ss = ss/sqrt(u)
    WE = weight(mah0, c1, c2, ss)
    ww = WE$w
    for(i in 1:n){
      X1[i,] = sqrt(ww[i])*X[i,]
      Y1[i,] = sqrt(ww[i])*Y[i,]
    }
    be = t(Y1)%*%X1%*%ginv(t(X1)%*%X1)
    g = sum(sum(abs(be-be0)))/sum(sum(abs(be0)))

    Ye = X%*%t(be)
    res = Y-Ye
    res1 = res
    for(j in 1:n){
      res1[j,] = sqrt(ww[j])*res[j,]
    }
    cov0 = t(res1)%*%res1/n
    covin = ginv(cov0)
    for(j in 1:n){
      mah0[j] = ((res[j,])%*%covin%*%(res[j,]))^(1/2)
    }
    if ((g<.0001)|(ii>100)){
      aa = 1
    }
    be0 = be
  }
  
  tscale = tau(mah0, c1, c2, ka, 1, .0001, -1)
  tt = tscale$t
  ss = tscale$s
  u = (tt^2)
  cov0 = u*cov0


  list(B=be0, cov=cov0, Binit=beinit)
}

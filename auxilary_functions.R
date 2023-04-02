#_____________________
#______Packages_______
#_____________________
library('MASS')
library('fda')
library('fda.usc')
library('FRegSigCom')
library('mvtnorm')
library('pcaPP')
library('expm')
library('refund')
library('data.table')
library('mvnTest')
library('ROCR')
library('FDboost')
#_____________________
#_____________________

#_______________________________________________________
#_________________Estimation_function___________________
#_______________________________________________________
est_fun = function(sco_Y, sco_X,
                   emodel = c("classical", "robust")){
  emodel = match.arg(emodel)
  if(emodel == "classical"){
    Bhat = ginv(t(sco_X) %*% sco_X) %*% t(sco_X) %*% sco_Y
  }else if(emodel == "robust"){
    sco_X = as.matrix(sco_X)
    sco_Y = as.matrix(sco_Y)
    tauMod = tauest(X = sco_X, Y = sco_Y, N=dim(sco_X)[1]/2,
                     c1=3, c2=5, ka=1)
    Bhat = as.matrix(t(tauMod$B))
  }
  return(Bhat)
}
#_______________________________________________________
#_______________________________________________________


#________________________________________
#__________Prediction_for_Y______________
#________________________________________
pred_fun = function(comp_Y, sco_X, Bhat){
  ncomp = dim(comp_Y$coefs)[2]
  nest = t(sco_X %*% Bhat)
  if(ncomp == 1){
    nh = nest[1] * comp_Y[1,]
  }else{
    nh = nest[1] * comp_Y[1,]
    for (j in 2:ncomp){
      nh = nh + nest[j] * comp_Y[j,]  
    }
  }
  return(nh)
}
#________________________________________
#________________________________________

#_______________________________________________________________________________________________________
#____________________________________________FPCA_main_effect___________________________________________
#_______________________________________________________________________________________________________
getPCA_main = function(data, nbasis, ncomp = NULL, emodel = c("classical", "robust"), rangeval){
  emodel = match.arg(emodel)
  n = dim(data)[1]
  p = dim(data)[2]
  dimnames(data)=list(as.character(1:n), as.character(1:p))
  grid_points = seq(rangeval[1], rangeval[2], length.out = p)
  bs_basis = create.bspline.basis(rangeval, nbasis = nbasis)
  evalbase = eval.basis(grid_points, bs_basis)
  fdobj = fdPar(bs_basis, int2Lfd(2), lambda=0)
  pcaobj = smooth.basisPar(grid_points, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd

  if(emodel == "classical"){
    mean_coef <- apply(t(pcaobj$coefs), 2, mean)
    sdata <- scale(t(pcaobj$coefs), scale = FALSE)
    dcov <- cov(sdata)
    d.eigen <- eigen(dcov)
    if(is.null(ncomp)){
      var_prop <- cumsum(d.eigen$values)/sum(d.eigen$values)
      ncomp <- which(var_prop>0.9)[1] 
    }
    PCs <- d.eigen$vectors[,1:ncomp]
    PCAcoef <- fd(PCs, bs_basis)
    mean_coef <- fd(as.vector(mean_coef), bs_basis)
    PCAscore <- sdata %*% PCs
  }else if(emodel == "robust"){
    mean_coef <- pcaPP:::l1median(t(pcaobj$coefs), trace = -1)
    sdata <- scale(t(pcaobj$coefs), center = mean_coef, scale = FALSE)
    if(is.null(ncomp)){
      dcov <- covPCAproj(sdata)
      d.eigen <- eigen(dcov$cov)
      var_prop <- cumsum(d.eigen$values)/sum(d.eigen$values)
      ncomp <- which(var_prop>0.9)[1]
    }
    ppur <- PCAproj(sdata, ncomp)
    PCs <- ppur$loadings
    PCAcoef <- fd(PCs, bs_basis)
    mean_coef <- fd(as.vector(mean_coef), bs_basis)
    PCAscore <- (sdata) %*% PCs
  }
  return(list(PCAcoef = PCAcoef, PCAscore = PCAscore, meanScore = mean_coef, evalbase = evalbase,
              bs_basis = bs_basis, gp = grid_points, ncomp = ncomp, emodel = emodel))
}
#_______________________________________________________________________________________________________
#_______________________________________________________________________________________________________



#_______________________________________________________________________________________________________
#_______________________________________FPCA_main_effect_test_data______________________________________
#_______________________________________________________________________________________________________
getTestPCA_main = function(object, data){
  bs_basis <- object$bs_basis
  PCAcoef <- object$PCAcoef
  gp <- object$gp
  mean.tr <- c(object$meanScore$coefs)
  n <- dim(data)[1]
  p <- dim(data)[2]
  dimnames(data) = list(as.character(1:n), as.character(1:p))
  pcaobj <- smooth.basisPar(gp, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
  
  sdata <- scale(t(pcaobj$coefs), center = mean.tr, scale = FALSE)
  PCAscore_test <- sdata %*% object$PCAcoef$coefs
  return(PCAscore.test = PCAscore_test)
}
#_______________________________________________________________________________________________________
#_______________________________________________________________________________________________________


#_____________________________________________________________________________________________________________
#____________________________________________FPCA_quadratic_effect____________________________________________
#_____________________________________________________________________________________________________________
getPCA_quad = function(data, nbasis, ncomp, emodel = c("classical", "robust"), rangeval){
  emodel = match.arg(emodel)
  np = length(data)
  n = dim(data[[1]])[1]
  p = list()
  grid_points = list()
  bs_basis = list()
  w_arg = list()
  w_mat = list()
  for(i in 1:np){
    p[[i]] = dim(data[[i]])[2]
    grid_points[[i]] = seq(rangeval[[i]][1], rangeval[[i]][2], length.out = p[[i]])
    bs_basis[[i]] = create.bspline.basis(rangeval[[i]], nbasis = nbasis[i])
    w_arg[[i]] = matrix(grid_points[[i]], nrow = n, ncol = p[[i]], byrow = T)
    w_mat[[i]] = smooth.basis(argvals = t(w_arg[[i]]), y = t(data[[i]]), fdParobj = bs_basis[[i]])$fd$coefs
  }

  
  int_mat = matrix(, np+np*(np-1)/2,2)
  fk = 1
  for(fi in 1:np){
    for(fj in fi:np){
      int_mat[fk,1] = (1:np)[fi];
      int_mat[fk,2] = (1:np)[fj];
      fk = fk+1;
    }
  }
  
  npq = dim(int_mat)[1]
  qmat = list()
  for(i1 in 1:npq){
    qmat1 = matrix(, nrow = (nbasis[int_mat[i1,1]]*nbasis[int_mat[i1,2]]) , ncol = n)
    for(i2 in 1:n){
      qmat1[,i2] = (w_mat[[int_mat[i1,1]]][,i2] %x% w_mat[[int_mat[i1,2]]][,i2])
    }
    qmat[[i1]] = qmat1
  }
  
  for(i in 1:npq)
    dimnames(qmat[[i]])=list(as.character(1:dim(qmat[[i]])[1]),
                             as.character(1:dim(qmat[[i]])[2]))
  
  mean_coef = list()
  PCAcoef = list()
  PCAscore = list()
  bs_basisiq = list()
  if(emodel == "classical"){
    for(iq in 1:npq){
      bs_basisiq[[iq]] = create.bspline.basis(rangeval[[1]], nbasis = (nbasis[int_mat[iq,1]]*nbasis[int_mat[iq,2]]))
      pcaobj = fd(qmat[[iq]], bs_basisiq[[iq]])
      
      mean_coef1 <- apply(t(pcaobj$coefs), 2, mean)
      sdata <- scale(t(pcaobj$coefs), scale = FALSE)
      dcov <- cov(sdata)
      d.eigen <- eigen(dcov)
      PCs <- d.eigen$vectors[,1:ncomp]
      PCAcoef[[iq]] <- fd(PCs, bs_basisiq[[iq]])
      mean_coef[[iq]] <- fd(as.vector(mean_coef1), bs_basisiq[[iq]])
      PCAscore[[iq]] <- (sdata) %*% PCs
    }
  }else if(emodel == "robust"){
    for(iq in 1:npq){
      bs_basisiq[[iq]] = create.bspline.basis(rangeval[[1]], nbasis = (nbasis[int_mat[iq,1]]*nbasis[int_mat[iq,2]]))
      pcaobj = fd(qmat[[iq]], bs_basisiq[[iq]])
      
      mean_coef1 <- pcaPP:::l1median(t(pcaobj$coefs), trace = -1)
      sdata <- scale(t(pcaobj$coefs), center = mean_coef1, scale = FALSE) 
      ppur <- PCAproj(sdata, ncomp)
      PCs <- ppur$loadings
      PCAcoef[[iq]] <- fd(PCs, bs_basisiq[[iq]])
      mean_coef[[iq]] <- fd(as.vector(mean_coef1), bs_basisiq[[iq]])
      PCAscore[[iq]] <- (sdata) %*% PCs
    }
  }
  return(list(PCAcoef = PCAcoef, PCAscore = PCAscore, bs_basis = bs_basis, bs_basisiq = bs_basisiq,
              gp = grid_points, meanScore = mean_coef, nbasis = nbasis))
}
#_____________________________________________________________________________________________________________
#_____________________________________________________________________________________________________________



#_____________________________________________________________________________________________________________
#_______________________________________FPCA_quadratic_effect_test_data_______________________________________
#_____________________________________________________________________________________________________________
getTestPCA_quad = function(object, data){
  bs_basis <- object$bs_basis
  bs_basisiq <- object$bs_basisiq
  PCAcoef <- object$PCAcoef
  gp <- object$gp
  mean.tr <- object$meanScore
  nbasis = object$nbasis
  np = length(data)
  n = dim(data[[1]])[1]
  p = list()

  w_arg = list()
  w_mat = list()
  for(i in 1:np){
    p[[i]] = dim(data[[i]])[2]
    bs_basis[[i]] = create.bspline.basis(c(gp[[i]][1], gp[[i]][length(gp[[i]])]), nbasis = nbasis[i])
    w_arg[[i]] = matrix(gp[[i]], nrow = n, ncol = p[[i]], byrow = T)
    w_mat[[i]] = smooth.basis(argvals = t(w_arg[[i]]), y = t(data[[i]]), fdParobj = bs_basis[[i]])$fd$coefs
  }
  
  
  int_mat = matrix(, np+np*(np-1)/2,2)
  fk = 1
  for(fi in 1:np){
    for(fj in fi:np){
      int_mat[fk,1] = (1:np)[fi];
      int_mat[fk,2] = (1:np)[fj];
      fk = fk+1;
    }
  }
  
  npq = dim(int_mat)[1]
  qmat = list()
  for(i1 in 1:npq){
    qmat1 = matrix(, nrow = (nbasis[int_mat[i1,1]]*nbasis[int_mat[i1,2]]) , ncol = n)
    for(i2 in 1:n){
      qmat1[,i2] = (w_mat[[int_mat[i1,1]]][,i2] %x% w_mat[[int_mat[i1,2]]][,i2])
    }
    qmat[[i1]] = qmat1
  }
  
  for(i in 1:npq)
    dimnames(qmat[[i]])=list(as.character(1:dim(qmat[[i]])[1]),
                             as.character(1:dim(qmat[[i]])[2]))
  
  PCAscore_test = list()
    for(iq in 1:npq){
      pcaobj = fd(qmat[[iq]], bs_basisiq[[iq]])
      
      sdata <- scale(t(pcaobj$coefs), center = mean.tr[[iq]]$coefs, scale = FALSE)
      PCAscore_test[[iq]] <- (sdata) %*% PCAcoef[[iq]]$coefs
    }
  return(PCAscore_test)
}
#_____________________________________________________________________________________________________________
#_____________________________________________________________________________________________________________



#_____________________________________________________________________________________________________________
#________________________________________________Prediction_function__________________________________________
#_____________________________________________________________________________________________________________
predict_fun = function(object, Xnew){
  fpca_results = object$fpca_results
  model.type = object$model.type
  fmodel = model.type[1]
  fp = object$model.details$fp
  grdp = object$model.details$grdp
  mindex = object$var_used[[1]]
  fq_index = object$model.details$fq_index
  Xnew1 = Xnew[mindex]
  fnp = length(Xnew1)
  fn = dim(Xnew1[[1]])[1]
  fBhat = object$model.details$fBhat
  fBhat_main = object$model.details$fBhat_main
  
  fPCA_Y = fpca_results$Y
  fcomp_Y = fPCA_Y$PCAcoef
  fmean_Y = fPCA_Y$meanScore
  fPCA_X = fpca_results$X
  
  fsco_X = list()
  for(fij in 1:fnp)
    fsco_X[[fij]] = getTestPCA_main(object = fPCA_X[[fij]], data = Xnew1[[fij]])
  
  quad_pca = object$model.details$quad_pca
  if(fmodel == "full" | fmodel == "true"){
    fsco_X_quad1 = getTestPCA_quad(object = quad_pca, data = Xnew1)
  }else{
    fsco_X_quad1 = getTestPCA_quad(object = quad_pca, data = Xnew)
  }
  
  fsco_X_quad = fsco_X_quad1[fq_index]
  
  fYfit = matrix(, nrow = fn, ncol = fp)
  for(fk in 1:fn){
    fXk = cbind(do.call(cbind, fsco_X), do.call(cbind, fsco_X_quad))[fk,]
    fmodel_k = pred_fun(comp_Y = fcomp_Y, sco_X = fXk, Bhat = fBhat) + fmean_Y
    fYfit[fk,] = eval.fd(grdp, fmodel_k)
  }
  
  fYfit_main = matrix(, nrow = fn, ncol = fp)
  for(fk in 1:fn){
    fXk_main = do.call(cbind, fsco_X)[fk,]
    fmodel_k_main = pred_fun(comp_Y = fcomp_Y, sco_X = fXk_main, Bhat = fBhat_main) + fmean_Y
    fYfit_main[fk,] = eval.fd(grdp, fmodel_k_main)
  }
  
  return(list(pred_main = fYfit_main, pred_quad = fYfit))
}
#_____________________________________________________________________________________________________________
#_____________________________________________________________________________________________________________




#__________________________________________________________________
#_________________________BIC_function_____________________________
#__________________________________________________________________
BIC_fun = function(Y, Yfit, ncompX, ncompY, emodel){
  n = dim(Y)[1]
  arg_bic = numeric(n)
  for(m in 1:n)
    arg_bic[m] = t(Y[m,] - Yfit[m,]) %*% (Y[m,] - Yfit[m,])
  
  bic_index = sort.int(arg_bic, decreasing = FALSE,
                       index.return = TRUE)
  ntrim = round(0.8 * n)
  index_trunc = bic_index$ix[1:ntrim]
  
  if(emodel == "classical"){
    BIC_val = n * log(sum(arg_bic) / n) +
      (ncompX * ncompY + 1) * log(n)
  }else if(emodel == "robust"){
    BIC_val = ntrim * log(sum(arg_bic[index_trunc]) / ntrim) +
      (ncompX * ncompY + 1) * log(ntrim)
  }
  return(BIC_val)
}
#__________________________________________________________________
#__________________________________________________________________


#_____________________________________________________________________________________________________________
#________________________________________________BIC_for_ncomp________________________________________________
#_____________________________________________________________________________________________________________
bic_fun_nc = function(Y, X, nbasisY, nbasisX, npca_max_Y, npca_max_X, rangevalY, rangevalX,
                      emodel = c("classical", "robust")){
  emodel = match.arg(emodel)
  np = length(X)
  n = dim(Y)[1]
  p = dim(Y)[2]

  BIC_mat_main = matrix(1e+999, nrow = (npca_max_X+3), ncol = (npca_max_Y+3))
  
  for(i in npca_max_X:(npca_max_X+3)){
    for(j in npca_max_Y:(npca_max_Y+3)){
      
      PCA_Y = getPCA_main(data = Y, nbasis = nbasisY, ncomp = j,
                          rangeval = rangevalY, emodel = emodel)
      
      sco_Y = PCA_Y$PCAscore
      comp_Y = PCA_Y$PCAcoef
      mean_Y = PCA_Y$meanScore
      
      sco_X = list()
      for(ij in 1:np){
        PCA_X = getPCA_main(data = X[[ij]], nbasis = nbasisX[ij], ncomp = i,
                            rangeval = rangevalX[[ij]], emodel = emodel)
        sco_X[[ij]] = PCA_X$PCAscore
      }
      
      Bhat = est_fun(sco_Y = sco_Y, sco_X = do.call(cbind, sco_X),
                     emodel = emodel)
      
      Yhat = matrix(, nrow = n, ncol = p)
      
      for(k in 1:n){
        Xk = do.call(cbind, sco_X)[k,]
        model_k = pred_fun(comp_Y = comp_Y, sco_X = Xk, Bhat = Bhat) + mean_Y
        Yhat[k,] = eval.fd(model_k, seq(rangevalY[1], rangevalY[2], length.out = p))
      }

      BIC_mat_main[(i), (j)] = BIC_fun(Y = Y, Yfit = Yhat, ncompX = i, ncompY = j,
                                       emodel = emodel)
    }
  }
  
  BIC_opt_main = which(BIC_mat_main == min(BIC_mat_main), arr.ind = TRUE)
  ncomp_opt_X = BIC_opt_main[1]
  ncomp_opt_Y = BIC_opt_main[2]

  return(list(ncompX = ncomp_opt_X, ncompY = ncomp_opt_Y, bic = BIC_mat_main[ncomp_opt_X, ncomp_opt_Y]))
}
#_____________________________________________________________________________________________________________
#_____________________________________________________________________________________________________________


#____________________________________________________________________________________________
#_____________________________________Variable_selection_____________________________________
#____________________________________________________________________________________________
var_sel = function(Y, X, nbasisY, nbasisX, ncompX = 2, ncompY = 2, rangevalY, rangevalX,
                   emodel = c("classical", "robust")){
  emodel = match.arg(emodel)
  np = length(X)
  n = dim(Y)[1]
  p = dim(Y)[2]
  
  PCA_Y = getPCA_main(data = Y, nbasis = nbasisY, ncomp = ncompY,
                      rangeval = rangevalY, emodel = emodel)
  
  sco_Y = PCA_Y$PCAscore
  comp_Y = PCA_Y$PCAcoef
  mean_Y = PCA_Y$meanScore
  
  sco_X = list()
  for(ij in 1:np){
    PCA_X = getPCA_main(data = X[[ij]], nbasis = nbasisX[ij], ncomp = ncompX,
                        rangeval = rangevalX[[ij]], emodel = emodel)
    sco_X[[ij]] = PCA_X$PCAscore
  }
  
  sco_Xqq = getPCA_quad(data = X, nbasis = nbasisX, ncomp = ncompX,
                       rangeval = rangevalX, emodel = emodel)
  sco_Xq = sco_Xqq$PCAscore
  
  BIC_individuals = numeric()
  
  for(ind in 1:np){
    Bhat = est_fun(sco_Y = sco_Y, sco_X = sco_X[[ind]],
                   emodel = emodel)
    
    Yhat = matrix(, nrow = n, ncol = p)
    
    for(k in 1:n){
      Xk = sco_X[[ind]][k,]
      model_k = pred_fun(comp_Y = comp_Y, sco_X = Xk, Bhat = Bhat) + mean_Y
      Yhat[k,] = eval.fd(model_k, seq(rangevalY[1], rangevalY[2], length.out = p))
    }
    
    BIC_individuals[ind] = BIC_fun(Y = Y, Yfit = Yhat, ncompX = ncompX,
                                   ncompY = ncompY, emodel = emodel)
  }
  BIC_order = order(BIC_individuals)
  
  main_model_start = sco_X[[BIC_order[1]]]
  BIC_forw = min(BIC_individuals)
  
  X_next = c(BIC_order[1], rep(NA, (length(BIC_order)-1)))
  X_out = c(which.min(BIC_individuals))
  
  for(f1 in 2:np){
    BIC_sel = rbind(subset(BIC_order, !(BIC_order %in% X_out)), NA)
    for(f2 in 1:ncol(BIC_sel)){
      sco_X_forw = cbind(main_model_start, sco_X[[BIC_sel[1,f2]]])
      
      Bhat = est_fun(sco_Y = sco_Y, sco_X = sco_X_forw,
                     emodel = emodel)
      
      Yhat = matrix(, nrow = n, ncol = p)
      
      for(k in 1:n){
        Xk = sco_X_forw[k,]
        model_k = pred_fun(comp_Y = comp_Y, sco_X = Xk, Bhat = Bhat) + mean_Y
        Yhat[k,] = eval.fd(model_k, seq(rangevalY[1], rangevalY[2], length.out = p))
      }
      
      BIC_sel[2,f2] = BIC_fun(Y = Y, Yfit = Yhat, ncompX = ncompX,
                              ncompY = ncompY, emodel = emodel)
    }
    
    BIC_next = BIC_sel[2,][which.min(BIC_sel[2,])]
    BIC_next2 = BIC_sel[1,][which.min(BIC_sel[2,])]
    
    X_out = c(X_out, BIC_next2)
    
    if(BIC_next/BIC_forw < 0.97){
      main_model_start = cbind(main_model_start, sco_X[[BIC_next2]])
      BIC_forw = BIC_next
      X_next[f1] = BIC_next2
    }else if(BIC_next/BIC_forw > 0.97){
      main_model_start = main_model_start
      BIC_forw = BIC_forw
      X_next[f1] = X_next[f1]
    }
  }
  
  selected_main = sort(subset(X_next, !(X_next %in% NA)))
  sco_X_selected = do.call(cbind, sco_X[selected_main])

  quad_mat = matrix(0, np+np*(np-1)/2,2)
  kk = 1
  for(ik in 1:np){
    for(jk in ik:np){
      quad_mat[kk,1] = (1:np)[ik];
      quad_mat[kk,2] = (1:np)[jk];
      kk = kk+1;
    }
  }
  
  np_sel = length(selected_main)
  
  quad_mat_trunc = matrix(0, np_sel+np_sel*(np_sel-1)/2,2)
  kk = 1
  for(ik in 1:np_sel){
    for(jk in ik:np_sel){
      quad_mat_trunc[kk,1] = sort(selected_main)[ik];
      quad_mat_trunc[kk,2] = sort(selected_main)[jk];
      kk = kk+1;
    }
  }
  
  mat_f1 = setkey(data.table(quad_mat_trunc))
  mat_f2 = setkey(data.table(quad_mat))
  trunc_quad_indx = mat_f2[mat_f1, which=TRUE]
  npq = length(trunc_quad_indx)
  
  sco_Xq_truc = sco_Xq[trunc_quad_indx]
  
  BIC_quads = numeric()
  
  for(indq in 1:npq){
    Bhat = est_fun(sco_Y = sco_Y, sco_X = cbind(sco_X_selected, sco_Xq_truc[[indq]]),
                   emodel = emodel)
    
    Yhat = matrix(, nrow = n, ncol = p)
    
    for(k in 1:n){
      Xk = cbind(sco_X_selected, sco_Xq_truc[[indq]])[k,]
      model_k = pred_fun(comp_Y = comp_Y, sco_X = Xk, Bhat = Bhat) + mean_Y
      Yhat[k,] = eval.fd(model_k, seq(rangevalY[1], rangevalY[2], length.out = p))
    }
    
    BIC_quads[indq] = BIC_fun(Y = Y, Yfit = Yhat, ncompX = ncompX,
                                   ncompY = ncompY, emodel = emodel)
  }
  BIC_orderq = order(BIC_quads)
  
  quad_model_start = sco_X_selected
  
  quad_term_selected = rep(NA, npq)
  q_out = numeric()
  
  for(f1 in 1: npq){
    BICq_sel = rbind(subset(BIC_orderq, !(BIC_orderq %in% q_out)), NA)
    for(f2 in 1:ncol(BICq_sel)){
      sco_Xq_forw = cbind(quad_model_start, sco_Xq_truc[[BICq_sel[1,f2]]])
      
      Bhat = est_fun(sco_Y = sco_Y, sco_X = sco_Xq_forw,
                     emodel = emodel)
      
      Yhat = matrix(, nrow = n, ncol = p)
      
      for(k in 1:n){
        Xk = sco_Xq_forw[k,]
        model_k = pred_fun(comp_Y = comp_Y, sco_X = Xk, Bhat = Bhat) + mean_Y
        Yhat[k,] = eval.fd(model_k, seq(rangevalY[1], rangevalY[2], length.out = p))
      }
      
      BICq_sel[2,f2] = BIC_fun(Y = Y, Yfit = Yhat, ncompX = ncompX,
                               ncompY = ncompY, emodel = emodel)
    }
    
    BICq_next = BICq_sel[2,][which.min(BICq_sel[2,])]
    BICq_next2 = BICq_sel[1,][which.min(BICq_sel[2,])]
    q_out = c(q_out, BICq_next2)
    
    if(BICq_next/BIC_forw < 0.995){
      quad_model_start = cbind(quad_model_start, sco_Xq_truc[[BICq_next2]])
      BIC_forw = BICq_next
      quad_term_selected[f1] = BICq_next2
    }else if(BICq_next/BIC_forw > 0.995){
      quad_model_start = quad_model_start
      BIC_forw = BIC_forw
      quad_term_selected[f1] = quad_term_selected[f1]
    }
  }
  
  sco_Xq_selected = subset(quad_term_selected, !(quad_term_selected %in% NA))
  if(is.logical(sco_Xq_selected) == FALSE)
    quad_term_selected = matrix(quad_mat_trunc[sort(sco_Xq_selected),], ncol=2)

  
  return(list(maine = selected_main, quade = quad_term_selected, qterms = trunc_quad_indx))
}
#____________________________________________________________________________________________
#____________________________________________________________________________________________



#___________________________________________________________________________________________________________
#______________________________________Function_for_method_of_Luo_and_Qi____________________________________
#___________________________________________________________________________________________________________

luoqi = function(X, Y, t.y, t.x, main.effect, inter.effect, model = c("full", "true", "selected", "main"),
                 s.n.basis=20, t.n.basis=20, inter.n.basis=20){
  if(model == "full"){
    model = cv.ff.interaction(X, Y, t.x, t.y, main.effect, inter.effect, adaptive=TRUE,
                              s.n.basis=s.n.basis, t.n.basis=t.n.basis, inter.n.basis=inter.n.basis)
    preds = pred.ff.interaction(model, X)
    coefs = getcoef.ff.interaction(model)
    resid = Y - preds
    resid = fdata(resid, argvals = NULL, rangeval = NULL,
                  names = NULL, fdata2d = FALSE)
    fdepth = depth.mode(resid)
  }else if(model == "true"){
    model = cv.ff.interaction(X, Y, t.x, t.y, main.effect, inter.effect, adaptive=TRUE,
                              s.n.basis=s.n.basis, t.n.basis=t.n.basis, inter.n.basis=inter.n.basis)
    preds = pred.ff.interaction(model, X)
    coefs = getcoef.ff.interaction(model)
    resid = Y - preds
    resid = fdata(resid, argvals = NULL, rangeval = NULL,
                  names = NULL, fdata2d = FALSE)
    fdepth = depth.mode(resid)
  }else if(model == "selected"){
    model = step.ff.interaction(X, Y, t.x, t.y,
                                s.n.basis=s.n.basis, t.n.basis=t.n.basis, inter.n.basis=inter.n.basis)
    preds = pred.ff.interaction(model, X)
    coefs = getcoef.ff.interaction(model)
    resid = Y - preds
    resid = fdata(resid, argvals = NULL, rangeval = NULL,
                  names = NULL, fdata2d = FALSE)
    fdepth = depth.mode(resid)
  }else if(model == "main"){
    model = cv.ff.interaction(X, Y, t.x, t.y, main.effect, inter.effect, adaptive=TRUE,
                              s.n.basis=s.n.basis, t.n.basis=t.n.basis, inter.n.basis=inter.n.basis)
    preds = pred.ff.interaction(model, X)
    coefs = getcoef.ff.interaction(model)
    resid = Y - preds
    resid = fdata(resid, argvals = NULL, rangeval = NULL,
                  names = NULL, fdata2d = FALSE)
    fdepth = depth.mode(resid)
  }
  
  return(list(model = model, preds = preds, coefs = coefs, fdepth = fdepth))
  
}
#___________________________________________________________________________________________________________
#___________________________________________________________________________________________________________


#___________________________________________________________________________________________________________
#___________________________________________________REISEE__________________________________________________
#___________________________________________________________________________________________________________

reisee = function(beta, beta_hat, domainY, domainX){
  f = (beta - beta_hat)^2
  nrow.f = dim(f)[1]
  ncol.f = dim(f)[2]
  gap.mat = matrix(diff(domainY), nrow = nrow.f, ncol = (length(domainY)-1), byrow = TRUE)
  r1 = matrix((rowSums(f[,-1]*gap.mat) + rowSums(f[, -ncol.f]*gap.mat))/2, nrow = nrow.f, ncol = 1)
  
  f1 = beta^2
  nrow.f1 = dim(f1)[1]
  ncol.f1 = dim(f1)[2]
  r2 = matrix((rowSums(f1[, -1]*gap.mat) + rowSums(f1[, -ncol.f1]*gap.mat))/2, nrow = nrow.f1, ncol = 1)
  
  arg1 = as.vector(r1)
  len.arg1 = length(arg1)
  s1 = sum((arg1[-1] + arg1[-len.arg1])*diff(domainX))/2
  
  arg2 = as.vector(r2)
  len.arg2 = length(arg2)
  s2 = sum((arg2[-1] + arg2[-len.arg2])*diff(domainX))/2
  
  return(s1/s2)
}
#___________________________________________________________________________________________________________
#___________________________________________________________________________________________________________


#_____________________________________________________
#____________________AUC_function_____________________
#_____________________________________________________
auc_fun = function(depth_values, labels){
  p_class = prediction(depth_values, labels)
  pf_class = performance(p_class, 'auc')
  return(as.numeric(pf_class@y.values))
}
#_____________________________________________________
#_____________________________________________________


#________________________________________________________________________
#________________________Cut_off_for_outlier_detection___________________
#________________________________________________________________________
cutC = function(data, depth, alpha, B){
  Call = numeric()
  for(nb in 1:B){
    prbs = depth/sum(depth)
    indx = sample(1:length(depth), length(depth), replace = TRUE, prob = prbs)
    
    Yboot = data[indx,] + 0.05 * matrix(rmvnorm(dim(data)[1], mean = rep(0, dim(data)[2]), 
                                                sigma = cov(data)), ncol = dim(data)[2])
    boot_depths = depth.mode(fdata(Yboot, argvals = NULL, rangeval = NULL,
                                   names = NULL, fdata2d = FALSE))$dep
    Call[nb] = quantile(boot_depths, probs = alpha)
  }
  
  return(median(Call))
}
#________________________________________________________________________
#________________________________________________________________________





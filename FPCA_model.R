#___________________________________________________FPCA_model__________________________________________________
#_______________________________________________________________________________________________________________

fpca = function(fY, fX, fmodel = c("full", "true", "selected"), emodel = c("classical", "robust"), rangevalY, rangevalX,
                mindex, qindex, fnbasisY, fnbasisX, fnpca_max_Y, fnpca_max_X, fncomp_Y, fncomp_X, BIC = TRUE){
  
  fmodel = match.arg(fmodel)
  emodel = match.arg(emodel)
  if(fmodel == "full"){
    
    fnp = length(fX)
    fn = dim(fY)[1]
    fp = dim(fY)[2]
    
    ntq = fnpca_max_X^2 * (fnp+fnp*(fnp-1)/2)
    ntm = fnp * fnpca_max_X
    ntX = ntq+ntm
    
    if(ntX > fp){
      while(ntX > fp){
        fnpca_max_X = fnpca_max_X-1
        ntq = fnpca_max_X^2 * (fnp+fnp*(fnp-1)/2)
        ntm = fnp * fnpca_max_X
        ntX = ntq+ntm
      } 
    }
    
    if(BIC){
      fbic = bic_fun_nc(Y = fY, X = fX, nbasisY = fnbasisY, nbasisX = fnbasisX, npca_max_X = fnpca_max_X,
                        rangevalY = rangevalY, rangevalX = rangevalX, npca_max_Y = fnpca_max_Y, emodel = emodel)
      ncY = fbic$ncompY
      ncX = fbic$ncompX
      }else{
      ncY = fncomp_Y
      ncX = fncomp_X
    }
  
    fPCA_Y = getPCA_main(data = fY, nbasis = fnbasisY, ncomp = ncY,
                         rangeval = rangevalY, emodel = emodel)
    fsco_Y = fPCA_Y$PCAscore
    fcomp_Y = fPCA_Y$PCAcoef
    fmean_Y = fPCA_Y$meanScore
    evaly = fPCA_Y$evalbase
    
    fsco_X = list()
    evalx = list()
    fcomp_X = list()
    for(fij in 1:fnp){
      fPCA_X = getPCA_main(data = fX[[fij]], nbasis = fnbasisX[[fij]], ncomp = ncX,
                           rangeval = rangevalX[[fij]], emodel = emodel)
      fsco_X[[fij]] = fPCA_X$PCAscore
      evalx[[fij]] = fPCA_X$evalbase
      fcomp_X[[fij]] = fPCA_X$PCAcoef
    }
    
    fnpq = fnp+fnp*(fnp-1)/2
    fq_index = 1:fnpq
    
    qindex = matrix(0, fnp+fnp*(fnp-1)/2,2)
    fkk = 1
    for(fik in 1:fnp){
      for(fjk in fik:fnp){
        qindex[fkk,1] = (1:fnp)[fik];
        qindex[fkk,2] = (1:fnp)[fjk];
        fkk = fkk+1;
      }
    }
    
    fsco_X_quad = getPCA_quad(data = fX, nbasis = fnbasisX, ncomp = ncX^2,
                              rangeval = rangevalX, emodel = emodel)

    fBhat = est_fun(sco_Y = fsco_Y, sco_X = cbind(do.call(cbind, fsco_X),
                                                  do.call(cbind, fsco_X_quad)),
                    emodel = emodel)
    
    fYfit = matrix(, nrow = fn, ncol = fp)
    for(fk in 1:fn){
      fXk = cbind(do.call(cbind, fsco_X), do.call(cbind, fsco_X_quad))[fk,]
      fmodel_k = pred_fun(comp_Y = fcomp_Y, sco_X = fXk, Bhat = fBhat) + fmean_Y
      fYfit[fk,] = eval.fd(fmodel_k, seq(rangevalY[1], rangevalY[2], length.out = fp))
    }
    
    fresids = fY - fYfit
    fdresids = fdata(fresids, argvals = NULL, rangeval = NULL,
                    names = NULL, fdata2d = FALSE)
    fdepth = depth.mode(fdresids)
    
    fBhat_main = est_fun(sco_Y = fsco_Y, sco_X = do.call(cbind, fsco_X),
                         emodel = emodel)
    
    fYfit_main = matrix(, nrow = fn, ncol = fp)
    for(fk in 1:fn){
      fXk_main = do.call(cbind, fsco_X)[fk,]
      fmodel_k_main = pred_fun(comp_Y = fcomp_Y, sco_X = fXk_main, Bhat = fBhat_main) + fmean_Y
      fYfit_main[fk,] = eval.fd(fmodel_k_main, seq(rangevalY[1], rangevalY[2], length.out = fp))
    }
  
    
    fresids_main = fY - fYfit_main
    fdresids_main = fdata(fresids_main, argvals = NULL, rangeval = NULL,
                     names = NULL, fdata2d = FALSE)
    fdepth_main = depth.mode(fdresids_main)
    
    coef_main = list()
    km = 1
    for(im in 1:fnp){
      coef_main[[im]] = evalx[[im]] %*% (fcomp_X[[im]]$coefs %*% fBhat_main[km: (km+ncX-1),] %*% t(fcomp_Y$coefs)) %*% t(evaly)
      km = im*ncX+1
    }
    
    fGamhat = as.matrix(fBhat[-(1:(ncX*fnp)),])
    coef_quad1 = list()
    km = 1
    for(iq in 1:length(fq_index)){
      coef_quad1[[iq]] = (evalx[[qindex[iq,1]]] %x% evalx[[qindex[iq,1]]]) %*% 
        (fcomp_X[[qindex[iq,1]]]$coefs %x% fcomp_X[[qindex[iq,2]]]$coefs %*% fGamhat[km:(km+ncX^2-1),] %*% t(fcomp_Y$coefs)) %*% t(evaly)
      km = iq*ncX^2+1
    }
    
    coef_quad = list()
    for(iq in 1:length(fq_index)){
      cq_array = array(dim = c(fp,fp,fp))
      km = 1
      for(iqj in 1:fp){
        cq_array[,,iqj] = coef_quad1[[iq]][km:(km+fp-1),]
        km = iqj*fp+1
      }
      coef_quad[[iq]] = cq_array
    }

    BIC_main = BIC_fun(Y = fY, Yfit = fYfit_main, ncompX = ncX, ncompY = ncY,
                       emodel = emodel)
    BIC_quad = BIC_fun(Y = fY, Yfit = fYfit, ncompX = ncX, ncompY = ncY,
                       emodel = emodel)
    BIC_values = list(main = BIC_main, int = BIC_quad)
    
    var_index = c("all the main, quadratic, and interaction terms were used in the model")
  }
  
  if(fmodel == "true"){
    fX = fX[mindex]
    
    fnp = length(fX)
    fn = dim(fY)[1]
    fp = dim(fY)[2]
    
    ntq = fnpca_max_X^2 * (fnp+fnp*(fnp-1)/2)
    ntm = fnp * fnpca_max_X
    ntX = ntq+ntm
    
    if(ntX > fp){
      while(ntX > fp){
        fnpca_max_X = fnpca_max_X-1
        ntq = fnpca_max_X^2 * (fnp+fnp*(fnp-1)/2)
        ntm = fnp * fnpca_max_X
        ntX = ntq+ntm
      } 
    }
    
    if(BIC){
      fbic = bic_fun_nc(Y = fY, X = fX, nbasisY = fnbasisY, nbasisX = fnbasisX, npca_max_X = fnpca_max_X,
                        rangevalY = rangevalY, rangevalX = rangevalX, npca_max_Y = fnpca_max_Y, emodel = emodel)
      ncY = fbic$ncompY
      ncX = fbic$ncompX
    }else{
      ncY = fncomp_Y
      ncX = fncomp_X
    }
    
    fPCA_Y = getPCA_main(data = fY, nbasis = fnbasisY, ncomp = ncY,
                         rangeval = rangevalY, emodel = emodel)
    fsco_Y = fPCA_Y$PCAscore
    fcomp_Y = fPCA_Y$PCAcoef
    fmean_Y = fPCA_Y$meanScore
    evaly = fPCA_Y$evalbase
    
    fsco_X = list()
    evalx = list()
    fcomp_X = list()
    for(fij in 1:fnp){
      fPCA_X = getPCA_main(data = fX[[fij]], nbasis = fnbasisX[[fij]], ncomp = ncX,
                           rangeval = rangevalX[[fij]], emodel = emodel)
      fsco_X[[fij]] = fPCA_X$PCAscore
      evalx[[fij]] = fPCA_X$evalbase
      fcomp_X[[fij]] = fPCA_X$PCAcoef
    }
    
    fsco_X_quad = getPCA_quad(data = fX, nbasis = fnbasisX, ncomp = ncX^2,
                              rangeval = rangevalX, emodel = emodel)

    fquad_mat = matrix(0, fnp+fnp*(fnp-1)/2,2)
    fkk = 1
    for(fik in 1:fnp){
      for(fjk in fik:fnp){
        fquad_mat[fkk,1] = (1:fnp)[fik];
        fquad_mat[fkk,2] = (1:fnp)[fjk];
        fkk = fkk+1;
      }
    }
    
    fmat_f1 = setkey(data.table(qindex))
    fmat_f2 = setkey(data.table(fquad_mat))
    fq_index = fmat_f2[fmat_f1, which=TRUE]
    
    fqsco_X_quad = fsco_X_quad[fq_index]

    fBhat = est_fun(sco_Y = fsco_Y, sco_X = cbind(do.call(cbind, fsco_X),
                                                  do.call(cbind, fqsco_X_quad)),
                    emodel = emodel)
    
    fYfit = matrix(, nrow = fn, ncol = fp)
    for(fk in 1:fn){
      fXk = cbind(do.call(cbind, fsco_X), do.call(cbind, fqsco_X_quad))[fk,]
      fmodel_k = pred_fun(comp_Y = fcomp_Y, sco_X = fXk, Bhat = fBhat) + fmean_Y
      fYfit[fk,] = eval.fd(fmodel_k, seq(rangevalY[1], rangevalY[2], length.out = fp))
    }
    
    fresids = fY - fYfit
    fdresids = fdata(fresids, argvals = NULL, rangeval = NULL,
                     names = NULL, fdata2d = FALSE)
    fdepth = depth.mode(fdresids)
    
    fBhat_main = est_fun(sco_Y = fsco_Y, sco_X = do.call(cbind, fsco_X),
                         emodel = emodel)
    
    fYfit_main = matrix(, nrow = fn, ncol = fp)
    for(fk in 1:fn){
      fXk_main = do.call(cbind, fsco_X)[fk,]
      fmodel_k_main = pred_fun(comp_Y = fcomp_Y, sco_X = fXk_main, Bhat = fBhat_main) + fmean_Y
      fYfit_main[fk,] = eval.fd(fmodel_k_main, seq(rangevalY[1], rangevalY[2], length.out = fp))
    }
    
    fresids_main = fY - fYfit_main
    fdresids_main = fdata(fresids_main, argvals = NULL, rangeval = NULL,
                          names = NULL, fdata2d = FALSE)
    fdepth_main = depth.mode(fdresids_main)
    
    coef_main = list()
    km = 1
    for(im in 1:fnp){
      coef_main[[im]] = evalx[[im]] %*% (fcomp_X[[im]]$coefs %*% fBhat_main[km: (km+ncX-1),] %*% t(fcomp_Y$coefs)) %*% t(evaly)
      km = im*ncX+1
    }
    
    fGamhat = as.matrix(fBhat[-(1:(ncX*fnp)),])
    coef_quad1 = list()
    km = 1
    qind = qindex
    for(iq1 in 1:nrow(qind)){
      for(iq2 in 1:ncol(qind)){
        qind[iq1,iq2] = which(qind[iq1,iq2] == mindex)
      }
    }
    for(iq in 1:length(fq_index)){
      coef_quad1[[iq]] = (evalx[[qind[iq,1]]] %x% evalx[[qind[iq,1]]]) %*% 
        (fcomp_X[[qind[iq,1]]]$coefs %x% fcomp_X[[qind[iq,2]]]$coefs %*% fGamhat[km:(km+ncX^2-1),] %*% t(fcomp_Y$coefs)) %*% t(evaly)
      km = iq*ncX^2+1
    }

    coef_quad = list()
    for(iq in 1:length(fq_index)){
      cq_array = array(dim = c(fp,fp,fp))
      km = 1
      for(iqj in 1:fp){
        cq_array[,,iqj] = coef_quad1[[iq]][km:(km+fp-1),]
        km = iqj*fp+1
      }
      coef_quad[[iq]] = cq_array
    }
    
    BIC_main = BIC_fun(Y = fY, Yfit = fYfit_main, ncompX = ncX, ncompY = ncY,
                       emodel = emodel)
    BIC_quad = BIC_fun(Y = fY, Yfit = fYfit, ncompX = ncX, ncompY = ncY,
                       emodel = emodel)
    BIC_values = list(main = BIC_main, int = BIC_quad)
    
    var_index = list(main = mindex, int = qindex)
  }
  
  if(fmodel == "selected"){
    fX1 = fX
    nfp1 = length(fX)
    
    svar = var_sel(Y = fY, X = fX, rangevalY = rangevalY, rangevalX = rangevalX,
                   nbasisY = fnbasisY, nbasisX = fnbasisX, emodel = emodel)
    mindex = svar$maine
    qindex = svar$quade
    
    fX = fX[mindex]
    
    
    fnp = length(fX)
    fn = dim(fY)[1]
    fp = dim(fY)[2]
    
    ntq = fnpca_max_X^2 * (fnp+fnp*(fnp-1)/2)
    ntm = fnp * fnpca_max_X
    ntX = ntq+ntm
    
    if(ntX > fp){
      while(ntX > fp){
        fnpca_max_X = fnpca_max_X-1
        ntq = fnpca_max_X^2 * (fnp+fnp*(fnp-1)/2)
        ntm = fnp * fnpca_max_X
        ntX = ntq+ntm
      } 
    }
    
    if(BIC){
      fbic = bic_fun_nc(Y = fY, X = fX, nbasisY = fnbasisY, nbasisX = fnbasisX[mindex], npca_max_X = fnpca_max_X,
                        rangevalY = rangevalY, rangevalX = rangevalX[mindex], npca_max_Y = fnpca_max_Y, emodel = emodel)
      ncY = fbic$ncompY
      ncX = fbic$ncompX
    }else{
      ncY = fncomp_Y
      ncX = fncomp_X
    }
    
    fPCA_Y = getPCA_main(data = fY, nbasis = fnbasisY, ncomp = ncY,
                         rangeval = rangevalY, emodel = emodel)
    fsco_Y = fPCA_Y$PCAscore
    fcomp_Y = fPCA_Y$PCAcoef
    fmean_Y = fPCA_Y$meanScore
    evaly = fPCA_Y$evalbase
    
    fsco_X = list()
    evalx = list()
    fcomp_X = list()
    for(fij in 1:fnp){
      fPCA_X = getPCA_main(data = fX[[fij]], nbasis = fnbasisX[mindex][[fij]], ncomp = ncX,
                           rangeval = rangevalX[mindex][[fij]], emodel = emodel)
      fsco_X[[fij]] = fPCA_X$PCAscore
      evalx[[fij]] = fPCA_X$evalbase
      fcomp_X[[fij]] = fPCA_X$PCAcoef
    }
    
    if(is.logical(svar$qterms) == FALSE){
      fsco_X_quad = getPCA_quad(data = fX1, nbasis = fnbasisX, ncomp = ncX^2,
                                rangeval = rangevalX, emodel = emodel)

      fquad_mat = matrix(0, nfp1+nfp1*(nfp1-1)/2,2)
      fkk = 1
      for(fik in 1:nfp1){
        for(fjk in fik:nfp1){
          fquad_mat[fkk,1] = (1:nfp1)[fik];
          fquad_mat[fkk,2] = (1:nfp1)[fjk];
          fkk = fkk+1;
        }
      }
      
      fmat_f1 = setkey(data.table(qindex))
      fmat_f2 = setkey(data.table(fquad_mat))
      fq_index = fmat_f2[fmat_f1, which=TRUE]
      
      fqsco_X_quad = fsco_X_quad[fq_index]

      fBhat = est_fun(sco_Y = fsco_Y, sco_X = cbind(do.call(cbind, fsco_X),
                                                    do.call(cbind, fqsco_X_quad)),
                      emodel = emodel)
      
      fYfit = matrix(, nrow = fn, ncol = fp)
      for(fk in 1:fn){
        fXk = cbind(do.call(cbind, fsco_X), do.call(cbind, fqsco_X_quad))[fk,]
        fmodel_k = pred_fun(comp_Y = fcomp_Y, sco_X = fXk, Bhat = fBhat) + fmean_Y
        fYfit[fk,] = eval.fd(fmodel_k, seq(rangevalY[1], rangevalY[2], length.out = fp))
      }
      
      fresids = fY - fYfit
      fdresids = fdata(fresids, argvals = NULL, rangeval = NULL,
                       names = NULL, fdata2d = FALSE)
      fdepth = depth.mode(fdresids)
      
      fBhat_main = est_fun(sco_Y = fsco_Y, sco_X = do.call(cbind, fsco_X),
                           emodel = emodel)
      
      fYfit_main = matrix(, nrow = fn, ncol = fp)
      for(fk in 1:fn){
        fXk_main = do.call(cbind, fsco_X)[fk,]
        fmodel_k_main = pred_fun(comp_Y = fcomp_Y, sco_X = fXk_main, Bhat = fBhat_main) + fmean_Y
        fYfit_main[fk,] = eval.fd(fmodel_k_main, seq(rangevalY[1], rangevalY[2], length.out = fp))
      }
      
      fresids_main = fY - fYfit_main
      fdresids_main = fdata(fresids_main, argvals = NULL, rangeval = NULL,
                            names = NULL, fdata2d = FALSE)
      fdepth_main = depth.mode(fdresids_main)
      
      coef_main = list()
      km = 1
      for(im in 1:fnp){
        coef_main[[im]] = evalx[[im]] %*% (fcomp_X[[im]]$coefs %*% fBhat_main[km: (km+ncX-1),] %*% t(fcomp_Y$coefs)) %*% t(evaly)
        km = im*ncX+1
      }
      
      fGamhat = as.matrix(fBhat[-(1:(ncX*fnp)),])
      coef_quad1 = list()
      km = 1
      qind = qindex
      for(iq1 in 1:nrow(qind)){
        for(iq2 in 1:ncol(qind)){
          qind[iq1,iq2] = which(qind[iq1,iq2] == mindex)
        }
      }
      for(iq in 1:length(fq_index)){
        coef_quad1[[iq]] = (evalx[[qind[iq,1]]] %x% evalx[[qind[iq,1]]]) %*% 
          (fcomp_X[[qind[iq,1]]]$coefs %x% fcomp_X[[qind[iq,2]]]$coefs %*% fGamhat[km:(km+ncX^2-1),] %*% t(fcomp_Y$coefs)) %*% t(evaly)
        km = iq*ncX^2+1
      }
      
      coef_quad = list()
      for(iq in 1:length(fq_index)){
        cq_array = array(dim = c(fp,fp,fp))
        km = 1
        for(iqj in 1:fp){
          cq_array[,,iqj] = coef_quad1[[iq]][km:(km+fp-1),]
          km = iqj*fp+1
        }
        coef_quad[[iq]] = cq_array
      }
    }else{
      fBhat = est_fun(sco_Y = fsco_Y, sco_X = do.call(cbind, fsco_X),
                      emodel = emodel)
      
      fYfit = matrix(, nrow = fn, ncol = fp)
      for(fk in 1:fn){
        fXk = do.call(cbind, fsco_X)[fk,]
        fmodel_k = pred_fun(comp_Y = fcomp_Y, sco_X = fXk, Bhat = fBhat) + fmean_Y
        fYfit[fk,] = eval.fd(fmodel_k, seq(rangevalY[1], rangevalY[2], length.out = fp))
      }
      
      fresids = fY - fYfit
      fdresids = fdata(fresids, argvals = NULL, rangeval = NULL,
                       names = NULL, fdata2d = FALSE)
      
      fYfit_main = fYfit
      fresids_main = fY - fYfit_main
      fdresids_main = fdata(fresids_main, argvals = NULL, rangeval = NULL,
                            names = NULL, fdata2d = FALSE)
      fdepth_main = depth.mode(fdresids_main)
      #fYhat_main = fYhat
      
      coef_main = list()
      km = 1
      for(im in 1:fnp){
        coef_main[[im]] = evalx[[im]] %*% (fcomp_X[[im]]$coefs %*% fBhat[km: (km+ncX-1),] %*% t(fcomp_Y$coefs)) %*% t(evaly)
        km = im*ncX+1
      }
      
      coef_quad = NA
      fdepth = fdepth_main
    }
    
    BIC_main = BIC_fun(Y = fY, Yfit = fYfit_main, ncompX = ncX, ncompY = ncY,
                       emodel = emodel)
    BIC_quad = BIC_fun(Y = fY, Yfit = fYfit, ncompX = ncX, ncompY = ncY,
                       emodel = emodel)
    BIC_values = list(main = BIC_main, int = BIC_quad)
    
    var_index = list(main = mindex, int = qindex)
  }
  
  return(list(main_fit = fYfit_main,
              quad_fit = fYfit,
              depth_main = fdepth_main, depth_quad = fdepth,
              var_used = var_index, BIC = BIC_values,
              main_coeffs = coef_main,
              quad_coeffs = coef_quad))
}

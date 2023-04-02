rm(list=ls())
source("auxilary_functions.R")
source("FPCA_model.R")
source("tauEst.R")
source("dgp.R")

# Number of observations
n = 250
# Outlier percentage
oper = 0.1
# Number of outliers
out_n = n * oper

# Rangeval
rangevalY = c(0, 1)
rangevalX = list(c(0, 1), c(0, 1), c(0, 1),
                 c(0, 1), c(0, 1), c(0, 1))
# Number of basis functions
nbasisY = 20
nbasisX = c(rep(20, 6))

# Data generation
set.seed(12345)
# Generate data with noise level 0.5 and contamination level 10%
simdat = data_generation(n+100, gpy = seq(0, 1, length.out = 100), gpx = seq(0, 1, length.out = 100), sd.error = 0.5, out.p = oper)
# Outlier indices
out_index = simdat$out.indx

Y = simdat$Y_train
Y_test = simdat$Yt_test
X = simdat$X_train
X_test = simdat$X_test

# Required arguments for the method LQ
t.y = seq(0, 1, length.out = 100)
t.x = list()
for(i in 1:6){
  t.x[[i]]=seq(0,1,length.out = 100)
}

main.effect = 1:6
inter.effect = rbind(c(1,1), c(1,2), c(1,3), c(1,4), c(1,5), c(1,6), 
                     c(2,2), c(2,3), c(2,4), c(2,5), c(2,6), 
                     c(3,3), c(3,4), c(3,5), c(3,6), 
                     c(4,4), c(4,5), c(4,6), 
                     c(5,5), c(5,6), 
                     c(6,6))
###

# Full model
# FoFRM-QI
full_pca = fpca(fY = Y, fX = X, fmodel = c("full"), emodel = c("classical"), rangevalY = rangevalY, rangevalX = rangevalX,
                fnbasisY = nbasisY, fnbasisX = nbasisX)

# Proposed method (RFoFRM-QI)
full_rpca = fpca(fY = Y, fX = X, fmodel = c("full"), emodel = c("robust"), rangevalY = rangevalY, rangevalX = rangevalX,
                 fnbasisY = nbasisY, fnbasisX = nbasisX)

# LQ
full_lq = luoqi(X, Y, t.y, t.x, main.effect, inter.effect, model = "full")

# True model
# Required arguments for the LQ
main.effect = c(1,2,3,4)
inter.effect = rbind(c(1,1), c(1,4), c(3,3), c(4,4))
###

# FoFRM-QI
true_pca = fpca(fY = Y, fX = X, fmodel = c("true"), emodel = c("classical"), rangevalY = rangevalY, rangevalX = rangevalX[main.effect],
                mindex = main.effect, qindex = inter.effect, fnbasisY = nbasisY, fnbasisX = nbasisX[main.effect])

# RFoFRM-QI
true_rpca = fpca(fY = Y, fX = X, fmodel = c("true"), emodel = c("robust"), rangevalY = rangevalY, rangevalX = rangevalX[main.effect],
                 mindex = main.effect, qindex = inter.effect, fnbasisY = nbasisY, fnbasisX = nbasisX[main.effect])

# LQ
true_lq = luoqi(X, Y, t.y, t.x, main.effect, inter.effect, model = "true")

# Selected model
# FoFRM-QI
selected_pca = fpca(fY = Y, fX = X, fmodel = c("selected"), emodel = c("classical"), rangevalY = rangevalY, rangevalX = rangevalX,
                    fnbasisY = nbasisY, fnbasisX = nbasisX)

# RFoFRM-QI
selected_rpca = fpca(fY = Y, fX = X, fmodel = c("selected"), emodel = c("robust"), rangevalY = rangevalY, rangevalX = rangevalX,
                     fnbasisY = nbasisY, fnbasisX = nbasisX)

# LQ
selected_lq = luoqi(X, Y, t.y, t.x, model = "selected")

# Main effect model for LQ
main.effect = 1:6
main_lq = luoqi(X, Y, t.y, t.x, main.effect, inter.effect=NULL, model = "main")


# Main effect model for pffr
x1 = X[[1]]
x2 = X[[2]]
x3 = X[[3]]
x4 = X[[4]]
x5 = X[[5]]
x6 = X[[6]]

pffr_model = pffr(Y ~ ff(x1, xind = t.y) + ff(x2, xind = t.y) + ff(x3, xind = t.y) + ff(x4, xind = t.y) +
                    ff(x5, xind = t.y) + ff(x6, xind = t.y), yind = t.y)
pffr_resids = matrix(pffr_model$residuals, nrow = nrow(Y), byrow=T)

pffr_resids = fdata(pffr_resids, argvals = NULL, rangeval = NULL,
                    names = NULL, fdata2d = FALSE)
fdepth_pffr = depth.mode(pffr_resids)

# Main effect model for FDboost
fdboost_list = list()
fdboost_list$Y=Y
fdboost_list$x1=scale(X[[1]], scale = FALSE)
fdboost_list$x2=scale(X[[2]], scale = FALSE)
fdboost_list$x3=scale(X[[3]], scale = FALSE)
fdboost_list$x4=scale(X[[4]], scale = FALSE)
fdboost_list$x5=scale(X[[5]], scale = FALSE)
fdboost_list$x6=scale(X[[6]], scale = FALSE)
fdboost_list$s = seq(0,1, length.out = 100)
fdboost_list$t = seq(0,1, length.out = 100)

mod_fdboost <- FDboost(Y ~ bsignal(x1, s) + bsignal(x2, s) + bsignal(x3, s) + bsignal(x4, s) +
                         bsignal(x5, s) + bsignal(x6, s), timeformula = ~ bbs(t), data = fdboost_list)

fdboost_resids = fdata(resid(mod_fdboost), argvals = NULL, rangeval = NULL,
                       names = NULL, fdata2d = FALSE)
fdepth_fdboost = depth.mode(fdboost_resids)

test_listt = list()
test_listt$x1=scale(X_test[[1]], scale = FALSE)
test_listt$x2=scale(X_test[[2]], scale = FALSE)
test_listt$x3=scale(X_test[[3]], scale = FALSE)
test_listt$x4=scale(X_test[[4]], scale = FALSE)
test_listt$x5=scale(X_test[[5]], scale = FALSE)
test_listt$x6=scale(X_test[[6]], scale = FALSE)
test_listt$s = seq(0,1, length.out = 100)
test_listt$t = seq(0,1, length.out = 100)

# Prediction of test data
predicted_pca_full = predict_fun(object = full_pca, Xnew = X_test)
predicted_rpca_full = predict_fun(object = full_rpca, Xnew = X_test)
predicted_lq_full = pred.ff.interaction(full_lq$model, X_test)


predicted_pca_true = predict_fun(object = true_pca, Xnew = X_test)
predicted_rpca_true = predict_fun(object = true_rpca, Xnew = X_test)
predicted_lq_true = pred.ff.interaction(true_lq$model, X_test)

predicted_pca_selected = predict_fun(object = selected_pca, Xnew = X_test)
predicted_rpca_selected = predict_fun(object = selected_rpca, Xnew = X_test)
predicted_lq_selected = pred.ff.interaction(selected_lq$model, X_test)

predicted_lq_main = pred.ff.interaction(main_lq$model, X_test)
names(X_test) = c("x1","x2","x3","x4","x5","x6")
predicted_pffr_main = predict(pffr_model, newdata = X_test)
predicted_fdboost_main = predict(mod_fdboost, newdata = test_listt)





# MSPE values
# Main effect model
# FDboost
mean((Y_test-predicted_fdboost_main)^2) # 1.86142
# pffr
mean((Y_test-predicted_pffr_main)^2) # 2.754522
# LQ
mean((Y_test-predicted_lq_main)^2) # 2.640711
# FoFRM-QI
mean((Y_test-predicted_pca_full$pred_main)^2) # 2.990115
# RFoFRM-QI
mean((Y_test-predicted_rpca_full$pred_main)^2) # 0.4044299


# Full model
# LQ
mean((Y_test-predicted_lq_full)^2) # 2.219668
# FoFRM-QI
mean((Y_test-predicted_pca_full$pred_quad)^2) # 15.35223
# RFoFRM-QI
mean((Y_test-predicted_rpca_full$pred_quad)^2) # 13.08159


# True model
# LQ
mean((Y_test-predicted_lq_true)^2) # 2.735734
# FoFRM-QI
mean((Y_test-predicted_pca_true$pred_quad)^2) # 9.524432
# RFoFRM-QI
mean((Y_test-predicted_rpca_true$pred_quad)^2) # 0.1518102


# Selected model
# LQ
mean((Y_test-predicted_lq_selected)^2) # 3.280639
# FoFRM-QI
mean((Y_test-predicted_pca_selected$pred_quad)^2) # 3.190253
# RFoFRM-QI
mean((Y_test-predicted_rpca_selected$pred_quad)^2) # 0.2357575


# RISEE values
risee_lq = numeric()
risee_pca = numeric()
risee_rpca = numeric()
for(re in 1:4){
  risee_lq[re] = reisee(simdat$main_coefs[[re]], t(true_lq$coefs$coef_main[[re]]), t.y, t.y)
  risee_pca[re] = reisee(simdat$main_coefs[[re]], true_pca$main_coeffs[[re]], t.y, t.y)
  risee_rpca[re] = reisee(simdat$main_coefs[[re]], true_rpca$main_coeffs[[re]], t.y, t.y)
}

risee_lq # 0.5769584 0.4359497 0.8741444 0.9720636
risee_pca # 3.934227 18.813288  9.349041  4.135128
risee_rpca # 0.23933392 0.04919713 0.14351054 0.28632693

# AUC values
labelsAll = rep(1, n)
labelsAll[out_index] = 0

# Main effect model
# FDboost
auc_fun(fdepth_fdboost$dep, labelsAll) # 0.9957333
# pffr
auc_fun(fdepth_pffr$dep, labelsAll) # 0.9907556
# LQ
auc_fun(full_lq$fdepth$dep, labelsAll) # 0.9792
# FoFRM-QI
auc_fun(full_pca$depth_main$dep, labelsAll) # 0.9962667
# RFoFRM-QI
auc_fun(full_rpca$depth_main$dep, labelsAll) # 0.9955556

# Full model
# LQ
auc_fun(true_lq$fdepth$dep, labelsAll) # 0.9404444
# FoFRM-QI
auc_fun(true_pca$depth_quad$dep, labelsAll) # 0.7701333
# RFoFRM-QI
auc_fun(true_rpca$depth_quad$dep, labelsAll) # 1

# True model
# LQ
auc_fun(selected_lq$fdepth$dep, labelsAll) # 0.9456
# FoFRM-QI
auc_fun(selected_pca$depth_quad$dep, labelsAll) # 0.9696
# RFoFRM-QI
auc_fun(selected_rpca$depth_quad$dep, labelsAll) # 1

# Selected model
# LQ
auc_fun(main_lq$fdepth$dep, labelsAll) # 0.9456
# FoFRM-QI
auc_fun(full_pca$depth_main$dep, labelsAll) # 0.9696
# RFoFRM-QI
auc_fun(full_rpca$depth_main$dep, labelsAll) # 1



# To test the normality of the estimates obtained by the proposed method:
true_rpca = fpca(fY = Y, fX = X, fmodel = c("true"), emodel = c("robust"), rangevalY = rangevalY, rangevalX = rangevalX[main.effect],
                 mindex = main.effect, qindex = inter.effect, fnbasisY = nbasisY, fnbasisX = nbasisX[main.effect],
                 BIC = FALSE, ncompX = 4, ncompY = 2)
tcf_rpca = true_rpca$Bps
AD.test(tcf_rpca[[1]])
AD.test(tcf_rpca[[2]])
AD.test(tcf_rpca[[3]])
AD.test(tcf_rpca[[4]])

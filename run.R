rm(list=ls())
source("auxilary_functions.R")
source("FPCA_model.R")
source("tauEst.R")
source("dgp_1.R")
source("dgp_2.R")

# Note that the data are generated under Scenario-1 in the paper

# Number of observations
n = 250
# Outlier percentage
oper = 0.1
# Number of outliers
out_n = n * oper

# Rangeval
rangevalY = c(0, 1)
rangevalX = list(c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1))
# Number of basis functions
nbasisY = 20
nbasisX = c(rep(20, 5))

# Data generation
set.seed(12345)
simdat1 = data_generation(n = n, j = 101)
simdat2 = data_generation2(n = n, j = 101)
Yt = simdat1$Yt
Y = simdat1$Y
X = simdat1$X

# Outlier indices
out_index = sample(1:n, out_n, replace = F)

# Replace outliers
Y[out_index,] = simdat2$Y[out_index,]
for(i in 1:5)
  X[[i]][out_index,] = simdat2$X[[i]][out_index,]

# Required arguments for the method LQ
t.y = seq(0, 1, length.out = 101)
t.x = list()
for(i in 1:5){
  t.x[[i]]=seq(0,1,length.out = 101)
}

main.effect = 1:5
inter.effect = rbind(c(1,1), c(1,2), c(1,3), c(1,4), c(1,5),
                     c(2,2), c(2,3), c(2,4), c(2,5),
                     c(3,3), c(3,4), c(3,5),
                     c(4,4), c(4,5),
                     c(5,5))
###

# Full model
# FoFRM-QI
full_pca = fpca(fY = Y, fX = X, fmodel = c("full"), emodel = c("classical"), rangevalY = rangevalY, rangevalX = rangevalX,
                fnbasisY = nbasisY, fnbasisX = nbasisX, fnpca_max_Y = 5, fnpca_max_X = 5, BIC = TRUE)

# Proposed method (RFoFRM-QI)
full_rpca = fpca(fY = Y, fX = X, fmodel = c("full"), emodel = c("robust"), rangevalY = rangevalY, rangevalX = rangevalX,
                 fnbasisY = nbasisY, fnbasisX = nbasisX, fnpca_max_Y=5, fnpca_max_X = 5, BIC = TRUE)

# LQ
full_lq = luoqi(X, Y, t.y, t.x, main.effect, inter.effect, model = "full")

# True model
# Required arguments for the LQ
main.effect = 1:3
inter.effect = rbind(c(1,1), c(1,3), c(3,3))
###

# FoFRM-QI
true_pca = fpca(fY = Y, fX = X, fmodel = c("true"), emodel = c("classical"), rangevalY = rangevalY, rangevalX = rangevalX[c(1,2,3)],
                mindex = c(1,2,3), qindex = rbind(c(1,1),c(1,3),c(3,3)), fnbasisY = nbasisY, fnbasisX = nbasisX[c(1,2,3)],
                fnpca_max_Y=5, fnpca_max_X = 5, BIC = TRUE)

# RFoFRM-QI
true_rpca = fpca(fY = Y, fX = X, fmodel = c("true"), emodel = c("robust"), rangevalY = rangevalY, rangevalX = rangevalX[c(1,2,3)],
                 mindex = c(1,2,3), qindex = rbind(c(1,1),c(1,3),c(3,3)), fnbasisY = nbasisY, fnbasisX = nbasisX[c(1,2,3)],
                 fnpca_max_Y=5, fnpca_max_X = 5, BIC = TRUE)

# LQ
true_lq = luoqi(X, Y, t.y, t.x, main.effect, inter.effect, model = "true")

# Selected model
# FoFRM-QI
selected_pca = fpca(fY = Y, fX = X, fmodel = c("selected"), emodel = c("classical"), rangevalY = rangevalY, rangevalX = rangevalX,
                    fnbasisY = nbasisY, fnbasisX = nbasisX, fnpca_max_Y = 5, fnpca_max_X = 5, BIC = TRUE)

# RFoFRM-QI
selected_rpca = fpca(fY = Y, fX = X, fmodel = c("selected"), emodel = c("robust"), rangevalY = rangevalY, rangevalX = rangevalX,
                     fnbasisY = nbasisY, fnbasisX = nbasisX, fnpca_max_Y=5, fnpca_max_X = 5, BIC = TRUE)

# LQ
selected_lq = luoqi(X, Y, t.y, t.x, model = "selected")

# Main effect model for LQ
main.effect = 1:5
main_lq = luoqi(X, Y, t.y, t.x, main.effect, inter.effect=NULL, model = "main")


# Main effect model for pffr
x1 = X[[1]]
x2 = X[[2]]
x3 = X[[3]]
x4 = X[[4]]
x5 = X[[5]]

pffr_model = pffr(Y ~ ff(x1, xind = t.y) + ff(x2, xind = t.y) + ff(x3, xind = t.y) + ff(x4, xind = t.y) + ff(x5, xind = t.y), yind = t.y)
yy <- pffr_model$fitted.values
pffr_fits = matrix(yy, nrow = nrow(Y), byrow=T)
pffr_resids = matrix(pffr_model$residuals, nrow = nrow(Y), byrow=T)

pffr_resids = fdata(pffr_resids, argvals = NULL, rangeval = NULL,
                    names = NULL, fdata2d = FALSE)
fdepth_pffr = depth.mode(pffr_resids)


# FE values
# Main effect model
# pffr
mean((Yt[-out_index,]-pffr_fits[-out_index,])^2) # 0.7665991
# LQ
mean((Yt[-out_index,]-main_lq$preds[-out_index,])^2) # 0.8032643
# FoFRM-QI
mean((Yt[-out_index,]-full_pca$main_fit[-out_index,])^2) # 1.023157
# RFoFRM-QI
mean((Yt[-out_index,]-full_rpca$main_fit[-out_index,])^2) # 0.2117907


# Full model
# LQ
mean((Yt[-out_index,]-full_lq$preds[-out_index,])^2) # 0.8118226
# FoFRM-QI
mean((Yt[-out_index,]-full_pca$quad_fit[-out_index,])^2) # 1.031069
# RFoFRM-QI
mean((Yt[-out_index,]-full_rpca$quad_fit[-out_index,])^2) # 0.285249


# True model
# LQ
mean((Yt[-out_index,]-true_lq$preds[-out_index,])^2) # 0.7809465
# FoFRM-QI
mean((Yt[-out_index,]-true_pca$quad_fit[-out_index,])^2) # 0.8881281
# RFoFRM-QI
mean((Yt[-out_index,]-true_rpca$quad_fit[-out_index,])^2) # 0.1289196


# Selected model
# LQ
mean((Yt[-out_index,]-selected_lq$preds[-out_index,])^2) # 0.8056885
# FoFRM-QI
mean((Yt[-out_index,]-selected_pca$quad_fit[-out_index,])^2) # 1.009029
# RFoFRM-QI
mean((Yt[-out_index,]-selected_rpca$quad_fit[-out_index,])^2) # 0.1894394


# RISEE values
risee_lq = numeric()
risee_pca = numeric()
risee_rpca = numeric()
for(re in 1:3){
  risee_lq[re] = reisee(simdat1$main_coefs[[re]], true_lq$coefs$coef_main[[re]], t.y, t.y)
  risee_pca[re] = reisee(simdat1$main_coefs[[re]], true_pca$main_coeffs[[re]], t.y, t.y)
  risee_rpca[re] = reisee(simdat1$main_coefs[[re]], true_rpca$main_coeffs[[re]], t.y, t.y)
}

risee_lq # 1.6749407 0.9999997 1.0000047
risee_pca # 1.011893 0.999954 1.000110
risee_rpca # 0.9368857 0.9861124 1.0002662

# AUC values
labelsAll = rep(1, n)
labelsAll[out_index] = 0

# Main effect model
# pffr
auc_fun(fdepth_pffr$dep, labelsAll) # 1
# LQ
auc_fun(full_lq$fdepth$dep, labelsAll) # 1
# FoFRM-QI
auc_fun(full_pca$depth_quad$dep, labelsAll) # 0.9973333
# RFoFRM-QI
auc_fun(full_rpca$depth_quad$dep, labelsAll) # 1

# Full model
# LQ
auc_fun(true_lq$fdepth$dep, labelsAll) # 1
# FoFRM-QI
auc_fun(true_pca$depth_quad$dep, labelsAll) # 1
# RFoFRM-QI
auc_fun(true_rpca$depth_quad$dep, labelsAll) # 1

# True model
# LQ
auc_fun(selected_lq$fdepth$dep, labelsAll) # 1
# FoFRM-QI
auc_fun(selected_pca$depth_quad$dep, labelsAll) # 0.9992889
# RFoFRM-QI
auc_fun(selected_rpca$depth_quad$dep, labelsAll) # 1

# Selected model
# LQ
auc_fun(main_lq$fdepth$dep, labelsAll) # 1
# FoFRM-QI
auc_fun(full_pca$depth_main$dep, labelsAll) # 1
# RFoFRM-QI
auc_fun(full_rpca$depth_main$dep, labelsAll) # 1



# To determine outliers in the empirical applications, first determine the cut-off value for the depth values using cutC function
# For example, for the proposed method under the True model
C_true_rpca = cutC(data = Y, depth = true_rpca$depth_quad$dep, alpha = 0.01, B = 200)
# Then, the outliers are;
out_true_rpca = which(true_rpca$depth_quad$dep < C_true_rpca)

out_true_rpca # 8  14  20  40  41  63  64  66  79  80  87  88  99 131 134 161 175 187 204 207 214 215 234 238 244 
# Compare with the true outliers:
sort(out_index) # 8  14  20  40  41  63  64  66  79  80  87  88  99 131 134 161 175 187 204 207 214 215 234 238 244

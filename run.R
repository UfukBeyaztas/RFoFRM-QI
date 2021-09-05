rm(list=ls())
source("auxilary_functions.R")
source("FPCA_model.R")
source("tauEst.R")
source("dgp_1.R")
source("dgp_2.R")

# Note that the data are generated under Scenario-1 in the paper

# Number of observations
n = 200
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
set.seed(123)
simdat1 = data_generation(n = n, j = 101)
simdat2 = data_generation2(n = n, j = 101)
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


# FE values
# Main effect model
# LQ
mean((Y[-out_index,]-main_lq$preds[-out_index,])^2) # 1.281794
# FoFRM-QI
mean((Y[-out_index,]-full_pca$main_fit[-out_index,])^2) # 1.61968
# RFoFRM-QI
mean((Y[-out_index,]-full_rpca$main_fit[-out_index,])^2) # 0.7666512


# Full model
# LQ
mean((Y[-out_index,]-full_lq$preds[-out_index,])^2) # 1.250122
# FoFRM-QI
mean((Y[-out_index,]-full_pca$quad_fit[-out_index,])^2) # 1.905992
# RFoFRM-QI
mean((Y[-out_index,]-full_rpca$quad_fit[-out_index,])^2) # 0.779937


# True model
# LQ
mean((Y[-out_index,]-true_lq$preds[-out_index,])^2) # 1.270857
# FoFRM-QI
mean((Y[-out_index,]-true_pca$quad_fit[-out_index,])^2) # 1.446903
# RFoFRM-QI
mean((Y[-out_index,]-true_rpca$quad_fit[-out_index,])^2) # 0.5269232


# Selected model
# LQ
mean((Y[-out_index,]-selected_lq$preds[-out_index,])^2) # 1.296777
# FoFRM-QI
mean((Y[-out_index,]-selected_pca$quad_fit[-out_index,])^2) # 1.835438
# RFoFRM-QI
mean((Y[-out_index,]-selected_rpca$quad_fit[-out_index,])^2) # 0.5319707


# RISEE values
risee_lq = numeric()
risee_pca = numeric()
risee_rpca = numeric()
for(re in 1:3){
  risee_lq[re] = reisee(simdat1$main_coefs[[re]], true_lq$coefs$coef_main[[re]], t.y, t.y)
  risee_pca[re] = reisee(simdat1$main_coefs[[re]], true_pca$main_coeffs[[re]], t.y, t.y)
  risee_rpca[re] = reisee(simdat1$main_coefs[[re]], true_rpca$main_coeffs[[re]], t.y, t.y)
}

risee_lq # 0.9969676 1.0000019 1.0000000
risee_pca # 4.712007 1.090265 1.332567
risee_rpca # 1.020841 1.000132 0.999809

# AUC values
labelsAll = rep(1, n)
labelsAll[out_index] = 0

# Main effect model
# LQ
auc_fun(full_lq$fdepth$dep, labelsAll) # 1
# FoFRM-QI
auc_fun(full_pca$depth_quad$dep, labelsAll) # 0.9986111
# RFoFRM-QI
auc_fun(full_rpca$depth_quad$dep, labelsAll) # 0.995

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
auc_fun(selected_pca$depth_quad$dep, labelsAll) # 1
# RFoFRM-QI
auc_fun(selected_rpca$depth_quad$dep, labelsAll) # 0.9947222

# Selected model
# LQ
auc_fun(main_lq$fdepth$dep, labelsAll) # 1
# FoFRM-QI
auc_fun(full_pca$depth_main$dep, labelsAll) # 1
# RFoFRM-QI
auc_fun(full_rpca$depth_main$dep, labelsAll) # 0.9988889



# To determine outliers in the empirical applications, first determine the cut-off value for the depth values using cutC function
# For example, for the proposed method under the True model
C_true_rpca = cutC(data = Y, depth = true_rpca$depth_quad$dep, alpha = 0.01, B = 200)
# Then, the outliers are;
out_true_rpca = which(true_rpca$depth_quad$dep < C_true_rpca)

out_true_rpca # 3   8  37  40  44  49  57  62  78  84 112 124 131 135 140 166 171 182 193 195
# Compare with the true outliers:
sort(out_index) # 3   8  37  40  44  49  57  62  78  84 112 124 131 135 140 166 171 182 193 195

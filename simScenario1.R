### SCENARIO 1: TWO TRUE PREDICTORS ###

## This is a simulation for the SWSR paper that compares the performance of SWSR (with optional quantiles and optional data splitting)
# with competing methods (WQS, linear regression) in a setting with 10 correlated predictors, two of which are truly associated with the 
# outcome, one more so than the other.

source("~/Desktop/functions.R")
source("~/Desktop/SWSR/swsr-functions.R") # brings in dplyr, glmnet, parallel

library(mvtnorm,lib.loc="~/RPackages")
library(gWQS,lib.loc="~/RPackages") # Note: Used Version 1.1.0

# read in arguments
## includes job, from 1 to 1000; each job was sent to a separate node for parallel processing
## also includes c (for "condition"), from 1 to 6, which determines R^2 and rho
args = (commandArgs(TRUE))
print(args)
for (i in 1:length(args)) { eval(parse(text=args[[i]])) }

path=paste0("~/Desktop/SWSR/Sim-2True/Cond",c,"/Scripts/")

conditions = data.frame(condition=1:6,true.r2=c(rep(.3,3),rep(.6,3)),true.pairwise.cov=rep(c(.2,.5,.8),2))
true.r2 = conditions$true.r2[c]
true.pairwise.cov = conditions$true.pairwise.cov[c]

# set parameters for data generation
n=300; p=10; my.mu = rep(10,p)

# fix coefficients of true model
beta0 = 10; beta1=5; beta2=3
betas = c(beta1,beta2,rep(0,p-2)) # not including intercept

## determine var.error as a function of var.reg and R^2
get.var.error = function(var.reg,r2) {
  var.error=(var.reg * (1-r2))/r2
  return(var.error)
}

# set number of bootstraps and permutations
boots = 100; N = 1000 # N is the number of permutations

# start the clock
time.start=proc.time()[[3]]

#### SIMULATION CODE ####
set.seed(as.numeric(paste0(c,job)))

# generate data
my.sigma = make.sigma(true.pairwise.cov,p)
x = rmvnorm(n,my.mu,my.sigma)
colnames(x) = paste0("x",1:p)
my.var.reg=var.reg(betas,my.sigma)
sd.error = sqrt(get.var.error(my.var.reg,true.r2))
y = beta0 + beta1*x[,1] + beta2*x[,2] + rnorm(n,0,sd.error)
index = (beta1/sum(betas))*x[,1] + (beta2/sum(betas))*x[,2]
wqs.dataframe = cbind.data.frame(y,x)

# draw training / test indicators
valid = sample(x=1:n,size=floor(.6*n),replace=F)
wqs.dataframe$validation_set = 0 # 0 for training
wqs.dataframe$validation_set[valid] = 1 # 1 for validation/test
x.train = x[-valid,]
y.train = y[-valid]
x.test = x[valid,]
y.test = y[valid]

# get SWSR.ct results - c for "continuous"; t for "training and test"
swsr.ct.fit = tryCatch(cv.swsr(x.train,y.train),error=function(e) {return(NA)})
if (length(swsr.ct.fit)!=1) {
  swsr.ct.lambda = which(swsr.ct.fit$swsr.fit$lambda==swsr.ct.fit$lambda.1se)
  if (swsr.ct.lambda==1) {swsr.ct.lambda=2} # in case all betas=0
  swsr.ct.betas = swsr.ct.fit$swsr.fit$beta.matrix[,swsr.ct.lambda][-1]
  names(swsr.ct.betas) = paste0("swsr.ct.",names(swsr.ct.betas),".beta")
  swsr.ct.select = ifelse(swsr.ct.betas==0,0,1)
  names(swsr.ct.select) = paste0("swsr.ct.",paste0("x",1:p),".select")
  swsr.ct.weights = swsr.ct.fit$swsr.fit$weight.matrix[,swsr.ct.lambda]
  names(swsr.ct.weights) = paste0("swsr.ct.",names(swsr.ct.weights),".weight")
  swsr.ct.size = swsr.ct.fit$swsr.fit$size[swsr.ct.lambda]
  swsr.ct.index = x.test %*% swsr.ct.weights # estimate index in test set for OLS
  swsr.ct.ols = lm(y.test ~ swsr.ct.index)
  swsr.ct.sumbeta = summary(swsr.ct.ols)$coef[2,1]
  swsr.ct.sumbeta.SE = summary(swsr.ct.ols)$coef[2,2]
  swsr.ct.sumbeta.p = summary(swsr.ct.ols)$coef[2,4]
  swsr.ct.index.est = x %*% swsr.ct.weights # now estimate index for all x
  swsr.ct.index.mse = mean((swsr.ct.index.est - index)^2) # get MSE
}

# get SWSR.cp results - c for "continuous"; p for "permutation"
swsr.cp.fit = tryCatch(cv.swsr(x,y),error=function(e) {return(NA)})
if (length(swsr.cp.fit)!=1) {
  swsr.cp.lambda = which(swsr.cp.fit$swsr.fit$lambda==swsr.cp.fit$lambda.1se)
  if (swsr.cp.lambda==1) {swsr.cp.lambda=2} # in case all betas=0
  swsr.cp.betas = swsr.cp.fit$swsr.fit$beta.matrix[,swsr.cp.lambda][-1]
  names(swsr.cp.betas) = paste0("swsr.cp.",names(swsr.cp.betas),".beta")
  swsr.cp.select = ifelse(swsr.cp.betas==0,0,1)
  names(swsr.cp.select) = paste0("swsr.cp.",paste0("x",1:p),".select")
  swsr.cp.weights = swsr.cp.fit$swsr.fit$weight.matrix[,swsr.cp.lambda]
  names(swsr.cp.weights) = paste0("swsr.cp.",names(swsr.cp.weights),".weight")
  swsr.cp.size = swsr.cp.fit$swsr.fit$size[swsr.cp.lambda]
  swsr.cp.index = x %*% swsr.cp.weights # estimate index
  swsr.cp.ols = lm(y ~ swsr.cp.index)
  swsr.cp.sumbeta = summary(swsr.cp.ols)$coef[2,1]
  swsr.perm.results = perm.test(swsr.cp.sumbeta,x,y,nperms = N)
  swsr.cp.sumbeta.p = swsr.perm.results$pvalue
  swsr.cp.index.mse = mean((swsr.cp.index - index)^2) # get MSE
}

# get WQS.t results - t for "training and test"
wqs.t.results = gwqs(y ~ NULL, mix_name = paste0("x",1:p), data = wqs.dataframe, q = 4,
                     valid_var = "validation_set", b = boots, b1_pos = T, b1_constr = F, family = "gaussian", 
                     plots = F, tables = F)
wqs.t.weights = wqs.t.results$final_weights[paste0("x",1:p),]$mean_weight
names(wqs.t.weights) = paste0("wqs.t.",paste0("x",1:p),".weight")
wqs.t.sumbeta = summary(wqs.t.results$fit)$coef[2,1]
wqs.t.sumbeta.SE = summary(wqs.t.results$fit)$coef[2,2]
wqs.t.sumbeta.p = summary(wqs.t.results$fit)$coef[2,4]
wqs.t.index.est = x %*% wqs.t.weights # estimate index for all x
wqs.t.index.mse = mean((wqs.t.index.est - index)^2) # get MSE

# get quantile X
x.train.q = wqs.t.results$data_t %>% dplyr::select(contains("_q")) 
x.test.q = wqs.t.results$data_v %>% dplyr::select(contains("_q")) 
x.q = rbind.data.frame(x.train.q,x.test.q)
x.q = x.q[order(as.numeric(rownames(x.q))),]

# compute true and estimated qindex
qindex = (beta1/sum(betas))*x.q[,1] + (beta2/sum(betas))*x.q[,2]
wqs.t.qindex.est = as.matrix(x.q) %*% wqs.t.weights
wqs.t.qindex.mse = mean((wqs.t.qindex.est - qindex)^2)

# get WQS.n results - n for "no split"
wqs.n.results = gwqs(y ~ NULL, mix_name = paste0("x",1:p), data = wqs.dataframe, q = 4,
                     validation=0, b = boots, b1_pos = T, b1_constr = F, family = "gaussian", 
                     wqs2 = F, plots = F, tables = F)
wqs.n.weights = wqs.n.results$final_weights[paste0("x",1:p),]$mean_weight
names(wqs.n.weights) = paste0("wqs.n.",paste0("x",1:p),".weight")
wqs.n.sumbeta = summary(wqs.n.results$fit)$coef[2,1]
wqs.n.sumbeta.SE = summary(wqs.n.results$fit)$coef[2,2]
wqs.n.sumbeta.p = summary(wqs.n.results$fit)$coef[2,4]
wqs.n.index.est = x %*% wqs.n.weights # estimate index for all x
wqs.n.index.mse = mean((wqs.n.index.est - index)^2) # get MSE
wqs.n.qindex.est = as.matrix(x.q) %*% wqs.n.weights
wqs.n.qindex.mse = mean((wqs.n.qindex.est - qindex)^2)

# get SWSR.qt results - q for "quantile"; t for "training and test"
swsr.qt.fit = tryCatch(cv.swsr(x.train.q,y.train,standardize=F),error=function(e) {return(NA)})
if (length(swsr.qt.fit)!=1) {
  swsr.qt.lambda = which(swsr.qt.fit$swsr.fit$lambda==swsr.qt.fit$lambda.1se)
  if (swsr.qt.lambda==1) {swsr.qt.lambda=2} # in case all betas=0
  swsr.qt.betas = swsr.qt.fit$swsr.fit$beta.matrix[,swsr.qt.lambda][-1]
  names(swsr.qt.betas) = paste0("swsr.qt.",paste0("x",1:p),".beta")
  swsr.qt.select = ifelse(swsr.qt.betas==0,0,1)
  names(swsr.qt.select) = paste0("swsr.qt.",paste0("x",1:p),".select")
  swsr.qt.weights = swsr.qt.fit$swsr.fit$weight.matrix[,swsr.qt.lambda]
  names(swsr.qt.weights) = paste0("swsr.qt.",paste0("x",1:p),".weight")
  swsr.qt.size = swsr.qt.fit$swsr.fit$size[swsr.qt.lambda]
  swsr.qt.index = x.test %*% swsr.qt.weights # estimate index in test set for OLS
  swsr.qt.ols = lm(y.test ~ swsr.qt.index)
  swsr.qt.sumbeta = summary(swsr.qt.ols)$coef[2,1]
  swsr.qt.sumbeta.SE = summary(swsr.qt.ols)$coef[2,2]
  swsr.qt.sumbeta.p = summary(swsr.qt.ols)$coef[2,4]
  swsr.qt.index.est = x %*% swsr.qt.weights # now estimate index for all x
  swsr.qt.index.mse = mean((swsr.qt.index.est - index)^2) # get MSE
  swsr.qt.qindex.est = as.matrix(x.q) %*% swsr.qt.weights
  swsr.qt.qindex.mse = mean((swsr.qt.qindex.est - qindex)^2)
}

# get SWSR.qp results - q for "quantile"; p for "permutation"
swsr.qp.fit = tryCatch(cv.swsr(x.q,y,standardize=F),error=function(e) {return(NA)})
if (length(swsr.qp.fit)!=1) {
  swsr.qp.lambda = which(swsr.qp.fit$swsr.fit$lambda==swsr.qp.fit$lambda.1se)
  if (swsr.qp.lambda==1) {swsr.qp.lambda=2} # in case all betas=0
  swsr.qp.betas = swsr.qp.fit$swsr.fit$beta.matrix[,swsr.qp.lambda][-1]
  names(swsr.qp.betas) = paste0("swsr.qp.",paste0("x",1:p),".beta")
  swsr.qp.select = ifelse(swsr.qp.betas==0,0,1)
  names(swsr.qp.select) = paste0("swsr.qp.",paste0("x",1:p),".select")
  swsr.qp.weights = swsr.qp.fit$swsr.fit$weight.matrix[,swsr.qp.lambda]
  names(swsr.qp.weights) = paste0("swsr.qp.",paste0("x",1:p),".weight")
  swsr.qp.size = swsr.qp.fit$swsr.fit$size[swsr.qp.lambda]
  swsr.qp.index = x %*% swsr.qp.weights # estimate index for all x
  swsr.qp.ols = lm(y ~ swsr.qp.index)
  swsr.qp.sumbeta = summary(swsr.qp.ols)$coef[2,1]
  swsr.perm.results = perm.test(swsr.qp.sumbeta,x,y,standardize = F,nperms=N)
  swsr.qp.sumbeta.p = swsr.perm.results$pvalue
  swsr.qp.index.mse = mean((swsr.qp.index - index)^2) # get MSE
  swsr.qp.qindex.est = as.matrix(x.q) %*% swsr.qp.weights
  swsr.qp.qindex.mse = mean((swsr.qp.qindex.est - qindex)^2)
}

# get one-at-a-time results - separate linear regressions for each component
one.fit = lapply(1:p, function(j) {
  mod = lm(y ~ x[,j])
  return(mod)
})
one.betas = sapply(one.fit,function(mod) {summary(mod)$coef[2,1]})
names(one.betas) = paste0("one.",paste0("x",1:p),".beta")
one.p = sapply(one.fit,function(mod) {summary(mod)$coef[2,4]})
names(one.p) = paste0("one.",paste0("x",1:p),".p")

# get all-included results
all.fit = lm(y~x); all.sum = summary(all.fit)
all.betas = summary(all.fit)$coef[-1,1]
names(all.betas) = paste0("all.",paste0("x",1:p),".beta")
all.p = summary(all.fit)$coef[-1,4]
names(all.p) = paste0("all.",paste0("x",1:p),".p")
all.mix.p = unname(pf(df1=all.sum$fstatistic["numdf"],df2=all.sum$fstatistic["dendf"],q=all.sum$fstatistic["value"],lower.tail=F))

# return results
output = tryCatch(data.frame(iter=job,
                             t(swsr.ct.betas),t(swsr.ct.weights),t(swsr.ct.select),swsr.ct.size,swsr.ct.sumbeta,swsr.ct.sumbeta.SE,swsr.ct.sumbeta.p,swsr.ct.index.mse,
                             t(swsr.cp.betas),t(swsr.cp.weights),t(swsr.cp.select),swsr.cp.size,swsr.cp.sumbeta,swsr.cp.sumbeta.p,swsr.cp.index.mse,
                             t(wqs.t.weights),wqs.t.sumbeta,wqs.t.sumbeta.SE,wqs.t.sumbeta.p,wqs.t.index.mse,wqs.t.qindex.mse,
                             t(wqs.n.weights),wqs.n.sumbeta,wqs.n.sumbeta.SE,wqs.n.sumbeta.p,wqs.n.index.mse,wqs.n.qindex.mse,
                             t(swsr.qt.betas),t(swsr.qt.weights),t(swsr.qt.select),swsr.qt.size,swsr.qt.sumbeta,swsr.qt.sumbeta.SE,swsr.qt.sumbeta.p,swsr.qt.index.mse,swsr.qt.qindex.mse,
                             t(swsr.qp.betas),t(swsr.qp.weights),t(swsr.qp.select),swsr.qp.size,swsr.qp.sumbeta,swsr.qp.sumbeta.p,swsr.qp.index.mse,swsr.qp.qindex.mse,
                             t(one.betas),t(one.p),t(all.betas),t(all.p),all.mix.p,
                             true.r2,true.pairwise.cov), error=function(e) {return(NULL)})
if (!is.null(output)) {write.table(output,file=paste0(path,"results.",job,".txt"),append=F,row.names=FALSE,
                                   col.names=T)}

# stop the clock
time.end=proc.time()[[3]]
time.end-time.start

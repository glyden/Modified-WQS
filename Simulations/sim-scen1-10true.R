### TEN TRUE PREDICTORS ###

## This is a simulation that compares the performance of penalized weighted sum regression 
# with competing methods in a setting with 10 true predictors out of 20, with equal effect sizes.

source("~/Desktop/Github/functions.R") # functions to help create data w/ desired correlation, signal
source("~/Desktop/Github/method-functions.R") # functions to fit the method

library(mvtnorm)
library(gWQS)

# read in args, a two-item vector with condition c (determines correlation, R^2) and job number for simulation (e.g. 1 to 1000)
args = (commandArgs(TRUE))
args = as.numeric(args)
print(args)
c=args[1]
job=args[2]

# set path
path="~/Desktop/Sim-10True/Results/"

# set conditions
conditions = data.frame(condition=1:6,true.r2=c(rep(.1,3),rep(.6,3)),true.pairwise.cov=rep(c(.2,.5,.8),2))
true.r2 = conditions$true.r2[c]
true.pairwise.cov = conditions$true.pairwise.cov[c]

# set parameters for data generation
n=300; p=20; my.mu = rep(10,p)

# fix coefficients of true model
beta=3
betas = c(rep(beta,10),rep(0,10)) # same coefficient for first 10 predictors, not including intercept
beta0 = 10
betas = c(beta0,betas)

## determine var.error as a function of var.reg and R^2
get.var.error = function(var.reg,r2) {
  var.error=(var.reg * (1-r2))/r2
  return(var.error)
}

# set number of iterations
boots = 100; N = 1000 # N = nperms

#### SIMULATION CODE ####
set.seed(as.numeric(paste0(c,job)))

# generate data
my.sigma = make.sigma(true.pairwise.cov,p) # correlation matrix
x = rmvnorm(n,my.mu,my.sigma)
colnames(x) = paste0("x",1:p)
x.design = cbind(rep(1,n),x)
my.var.reg=var.reg(betas[-1],my.sigma)
sd.error = sqrt(get.var.error(my.var.reg,true.r2))
# r2(betas[-1],my.sigma,sd.error^2) # test
y = x.design %*% betas + rnorm(n,0,sd.error)

# standardize x and make data for wqs
x = apply(x,2,function(vect) return(vect/sd(vect)))
wqs.dataframe = cbind.data.frame(y,x)

# draw training / test indicators
valid = sample(x=1:n,size=floor(.6*n),replace=F)
wqs.dataframe$validation_set = 0 # 0 for training
wqs.dataframe$validation_set[valid] = 1 # 1 for validation/test
x.train = x[-valid,]
y.train = y[-valid]
x.test = x[valid,]
y.test = y[valid]

# get pr.t1 results - pr for "penalized regression"; t for "training and test"; 1 for L1
pr.t1.fit = tryCatch(cv.fit(x.train,y.train),error=function(e) {return(NA)})
if (length(pr.t1.fit)!=1) {
  pr.t1.lambda.min = which(pr.t1.fit$all.fit$lambda==pr.t1.fit$lambda.min)
  #pr.t1.lambda.1se = which(pr.t1.fit$all.fit$lambda==pr.t1.fit$lambda.1se)
  if (pr.t1.lambda.min==1) {pr.t1.lambda.min=2} # in case all betas=0
  #if (pr.t1.lambda.1se==1) {pr.t1.lambda.1se=2} # in case all betas=0
  
  # results for lambda.min
  pr.t1.betas = pr.t1.fit$all.fit$beta.matrix[,pr.t1.lambda.min][-1]
  names(pr.t1.betas) = paste0("pr.t1.",names(pr.t1.betas),".beta")
  pr.t1.select = ifelse(pr.t1.betas==0,0,1)
  names(pr.t1.select) = paste0("pr.t1.",paste0("x",1:p),".select")
  pr.t1.weights = pr.t1.fit$all.fit$weight.matrix[,pr.t1.lambda.min]
  names(pr.t1.weights) = paste0("pr.t1.",names(pr.t1.weights),".weight")
  pr.t1.size = pr.t1.fit$all.fit$size[pr.t1.lambda.min]
  pr.t1.index = x.test %*% pr.t1.weights # estimate index in test set for OLS
  pr.t1.ols = lm(y.test ~ pr.t1.index)
  pr.t1.sumbeta = summary(pr.t1.ols)$coef[2,1]
  pr.t1.sumbeta.SE = summary(pr.t1.ols)$coef[2,2]
  pr.t1.sumbeta.p = summary(pr.t1.ols)$coef[2,4]
  
  # # results for lambda.1se
  # pr.t1.betas.1se = pr.t1.fit$all.fit$beta.matrix[,pr.t1.lambda.1se][-1]
  # names(pr.t1.betas.1se) = paste0("pr.t1.",names(pr.t1.betas.1se),".beta.1se")
  # pr.t1.select.1se = ifelse(pr.t1.betas.1se==0,0,1)
  # names(pr.t1.select.1se) = paste0("pr.t1.",paste0("x",1:p),".select.1se")
  # pr.t1.weights.1se = pr.t1.fit$all.fit$weight.matrix[,pr.t1.lambda.1se]
  # names(pr.t1.weights.1se) = paste0("pr.t1.",names(pr.t1.weights.1se),".weight.1se")
  # pr.t1.size.1se = pr.t1.fit$all.fit$size[pr.t1.lambda.1se]
  # pr.t1.index.1se = x.test %*% pr.t1.weights.1se # estimate index in test set for OLS
  # pr.t1.ols.1se = lm(y.test ~ pr.t1.index.1se)
  # pr.t1.sumbeta.1se = summary(pr.t1.ols.1se)$coef[2,1]
  # pr.t1.sumbeta.SE.1se = summary(pr.t1.ols.1se)$coef[2,2]
  # pr.t1.sumbeta.p.1se = summary(pr.t1.ols.1se)$coef[2,4]
  # pr.t1.index.est.1se = x %*% pr.t1.weights.1se # now estimate index for all x
  # pr.t1.index.mse.1se = mean((pr.t1.index.est.1se - index)^2) # get MSE
}

# get pr.p1 results - p for "permutation"; 1 for L1 penalty
pr.p1.fit = tryCatch(cv.fit(x,y),error=function(e) {return(NA)})
if (length(pr.p1.fit)!=1) {
  pr.p1.lambda.min = which(pr.p1.fit$all.fit$lambda==pr.p1.fit$lambda.min)
  #pr.p1.lambda.1se = which(pr.p1.fit$all.fit$lambda==pr.p1.fit$lambda.1se)
  if (pr.p1.lambda.min==1) {pr.p1.lambda.min=2} # in case all betas=0
  #if (pr.p1.lambda.1se==1) {pr.p1.lambda.1se=2} # in case all betas=0
  
  # results for lambda.min
  pr.p1.betas = pr.p1.fit$all.fit$beta.matrix[,pr.p1.lambda.min][-1]
  names(pr.p1.betas) = paste0("pr.p1.",names(pr.p1.betas),".beta")
  pr.p1.select = ifelse(pr.p1.betas==0,0,1)
  names(pr.p1.select) = paste0("pr.p1.",paste0("x",1:p),".select")
  pr.p1.weights = pr.p1.fit$all.fit$weight.matrix[,pr.p1.lambda.min]
  names(pr.p1.weights) = paste0("pr.p1.",names(pr.p1.weights),".weight")
  pr.p1.size = pr.p1.fit$all.fit$size[pr.p1.lambda.min]
  pr.p1.index = x %*% pr.p1.weights # estimate index for all x
  pr.p1.ols = lm(y ~ pr.p1.index)
  pr.p1.sumbeta = summary(pr.p1.ols)$coef[2,1]
  time.start=proc.time()[[3]]
  pr.perm.results = tryCatch(perm.test(pr.p1.sumbeta,x,y,lambda="lambda.min",nperms=N),
                               error=function(e) {return(NA)})
  pr.p1.permTime = proc.time()[[3]] - time.start
  pr.p1.permN = pr.perm.results$N
  pr.p1.sumbeta.p = pr.perm.results$pvalue
  
  # # results for lambda.1se 
  # pr.p1.betas.1se = pr.p1.fit$all.fit$beta.matrix[,pr.p1.lambda.1se][-1]
  # names(pr.p1.betas.1se) = paste0("pr.p1.",names(pr.p1.betas.1se),".beta.1se")
  # pr.p1.select.1se = ifelse(pr.p1.betas.1se==0,0,1)
  # names(pr.p1.select.1se) = paste0("pr.p1.",paste0("x",1:p),".select.1se")
  # pr.p1.weights.1se = pr.p1.fit$all.fit$weight.matrix[,pr.p1.lambda.1se]
  # names(pr.p1.weights.1se) = paste0("pr.p1.",names(pr.p1.weights.1se),".weight.1se")
  # pr.p1.size.1se = pr.p1.fit$all.fit$size[pr.p1.lambda.1se]
  # pr.p1.index.1se = x %*% pr.p1.weights.1se # estimate index for all x
  # pr.p1.ols.1se = lm(y ~ pr.p1.index.1se)
  # pr.p1.sumbeta.1se = summary(pr.p1.ols.1se)$coef[2,1]
  # pr.perm.results.1se = tryCatch(perm.test(pr.p1.sumbeta.1se,x,y,lambda="lambda.1se",nperms=N),
  #                                  error=function(e) {return(NA)})
  # pr.p1.permN.1se = pr.perm.results.1se$N
  # pr.p1.sumbeta.p.1se = pr.perm.results.1se$pvalue
  # pr.p1.index.mse.1se = mean((pr.p1.index.1se - index)^2) # get MSE
}

# get WQS.t results - t for "training and test"
time.start=proc.time()[[3]]
wqs.t.results = gwqs(y ~ wqs, mix_name = paste0("x",1:p), data = wqs.dataframe, q = NULL,
                     valid_var = "validation_set", b = boots, b1_pos = T, b1_constr = F, family = "gaussian")
wqs.t.bootTime = proc.time()[[3]] - time.start
wqs.t.weights = wqs.t.results$final_weights[paste0("x",1:p),]$mean_weight
names(wqs.t.weights) = paste0("wqs.t.",paste0("x",1:p),".weight")
wqs.t.sumbeta = summary(wqs.t.results$fit)$coef[2,1]
wqs.t.sumbeta.SE = summary(wqs.t.results$fit)$coef[2,2]
wqs.t.sumbeta.p = summary(wqs.t.results$fit)$coef[2,4]

# get WQS.n results - n for "no split"
time.start=proc.time()[[3]]
wqs.n.results = gwqs(y ~ wqs, mix_name = paste0("x",1:p), data = wqs.dataframe, q = NULL,
                     validation=0, b = boots, b1_pos = T, b1_constr = F, family = "gaussian")
wqs.n.bootTime = proc.time()[[3]] - time.start
wqs.n.weights = wqs.n.results$final_weights[paste0("x",1:p),]$mean_weight
names(wqs.n.weights) = paste0("wqs.n.",paste0("x",1:p),".weight")
wqs.n.sumbeta = summary(wqs.n.results$fit)$coef[2,1]
wqs.n.sumbeta.SE = summary(wqs.n.results$fit)$coef[2,2]
wqs.n.sumbeta.p = summary(wqs.n.results$fit)$coef[2,4]

# get pr.t2 results - t for "training and test"; 2 for L2 penalty
pr.t2.fit = tryCatch(cv.fit(x.train,y.train,alpha=0),error=function(e) {return(NA)})
if (length(pr.t2.fit)!=1) {
  pr.t2.lambda.min = which(pr.t2.fit$all.fit$lambda==pr.t2.fit$lambda.min)
  # pr.t2.lambda.1se = which(pr.t2.fit$all.fit$lambda==pr.t2.fit$lambda.1se)
  if (pr.t2.lambda.min==1) {pr.t2.lambda.min=2} # in case all betas=0
  # if (pr.t2.lambda.1se==1) {pr.t2.lambda.1se=2} # in case all betas=0
  
  # results for lambda.min
  pr.t2.betas = pr.t2.fit$all.fit$beta.matrix[,pr.t2.lambda.min][-1]
  names(pr.t2.betas) = paste0("pr.t2.",names(pr.t2.betas),".beta")
  pr.t2.select = ifelse(pr.t2.betas==0,0,1)
  names(pr.t2.select) = paste0("pr.t2.",paste0("x",1:p),".select")
  pr.t2.weights = pr.t2.fit$all.fit$weight.matrix[,pr.t2.lambda.min]
  names(pr.t2.weights) = paste0("pr.t2.",names(pr.t2.weights),".weight")
  pr.t2.size = pr.t2.fit$all.fit$size[pr.t2.lambda.min]
  pr.t2.index = x.test %*% pr.t2.weights # estimate index in test set for OLS
  pr.t2.ols = lm(y.test ~ pr.t2.index)
  pr.t2.sumbeta = summary(pr.t2.ols)$coef[2,1]
  pr.t2.sumbeta.SE = summary(pr.t2.ols)$coef[2,2]
  pr.t2.sumbeta.p = summary(pr.t2.ols)$coef[2,4]
  
  # # results for lambda.1se
  # pr.t2.betas.1se = pr.t2.fit$all.fit$beta.matrix[,pr.t2.lambda.1se][-1]
  # names(pr.t2.betas.1se) = paste0("pr.t2.",names(pr.t2.betas.1se),".beta.1se")
  # pr.t2.select.1se = ifelse(pr.t2.betas.1se==0,0,1)
  # names(pr.t2.select.1se) = paste0("pr.t2.",paste0("x",1:p),".select.1se")
  # pr.t2.weights.1se = pr.t2.fit$all.fit$weight.matrix[,pr.t2.lambda.1se]
  # names(pr.t2.weights.1se) = paste0("pr.t2.",names(pr.t2.weights.1se),".weight.1se")
  # pr.t2.size.1se = pr.t2.fit$all.fit$size[pr.t2.lambda.1se]
  # pr.t2.index.1se = x.test %*% pr.t2.weights.1se # estimate index in test set for OLS
  # pr.t2.ols.1se = lm(y.test ~ pr.t2.index.1se)
  # pr.t2.sumbeta.1se = summary(pr.t2.ols.1se)$coef[2,1]
  # pr.t2.sumbeta.SE.1se = summary(pr.t2.ols.1se)$coef[2,2]
  # pr.t2.sumbeta.p.1se = summary(pr.t2.ols.1se)$coef[2,4]
  # pr.t2.index.est.1se = x %*% pr.t2.weights.1se # now estimate index for all x
  # pr.t2.index.mse.1se = mean((pr.t2.index.est.1se - index)^2) # get MSE
}

# get pr.p2 results - p for "permutation"; 2 for L2 penalty
pr.p2.fit = tryCatch(cv.fit(x,y,alpha=0),error=function(e) {return(NA)})
if (length(pr.p2.fit)!=1) {
  pr.p2.lambda.min = which(pr.p2.fit$all.fit$lambda==pr.p2.fit$lambda.min)
  # pr.p2.lambda.1se = which(pr.p2.fit$all.fit$lambda==pr.p2.fit$lambda.1se)
  if (pr.p2.lambda.min==1) {pr.p2.lambda.min=2} # in case all betas=0
  # if (pr.p2.lambda.1se==1) {pr.p2.lambda.1se=2} # in case all betas=0
  
  # results for lambda.min
  pr.p2.betas = pr.p2.fit$all.fit$beta.matrix[,pr.p2.lambda.min][-1]
  names(pr.p2.betas) = paste0("pr.p2.",names(pr.p2.betas),".beta")
  pr.p2.select = ifelse(pr.p2.betas==0,0,1)
  names(pr.p2.select) = paste0("pr.p2.",paste0("x",1:p),".select")
  pr.p2.weights = pr.p2.fit$all.fit$weight.matrix[,pr.p2.lambda.min]
  names(pr.p2.weights) = paste0("pr.p2.",names(pr.p2.weights),".weight")
  pr.p2.size = pr.p2.fit$all.fit$size[pr.p2.lambda.min]
  pr.p2.index = x %*% pr.p2.weights # estimate index for all x
  pr.p2.ols = lm(y ~ pr.p2.index)
  pr.p2.sumbeta = summary(pr.p2.ols)$coef[2,1]
  time.start=proc.time()[[3]]
  pr.perm.results = tryCatch(perm.test(pr.p2.sumbeta,x,y,lambda="lambda.min",nperms=N,alpha=0),
                               error=function(e) {return(NA)})
  pr.p2.permTime = proc.time()[[3]] - time.start
  pr.p2.permN = pr.perm.results$N
  pr.p2.sumbeta.p = pr.perm.results$pvalue
  
  # # results for lambda.1se 
  # pr.p2.betas.1se = pr.p2.fit$all.fit$beta.matrix[,pr.p2.lambda.1se][-1]
  # names(pr.p2.betas.1se) = paste0("pr.p2.",names(pr.p2.betas.1se),".beta.1se")
  # pr.p2.select.1se = ifelse(pr.p2.betas.1se==0,0,1)
  # names(pr.p2.select.1se) = paste0("pr.p2.",paste0("x",1:p),".select.1se")
  # pr.p2.weights.1se = pr.p2.fit$all.fit$weight.matrix[,pr.p2.lambda.1se]
  # names(pr.p2.weights.1se) = paste0("pr.p2.",names(pr.p2.weights.1se),".weight.1se")
  # pr.p2.size.1se = pr.p2.fit$all.fit$size[pr.p2.lambda.1se]
  # pr.p2.index.1se = x %*% pr.p2.weights.1se # estimate index for all x
  # pr.p2.ols.1se = lm(y ~ pr.p2.index.1se)
  # pr.p2.sumbeta.1se = summary(pr.p2.ols.1se)$coef[2,1]
  # pr.perm.results.1se = tryCatch(perm.test(pr.p2.sumbeta.1se,x,y,alpha=0,lambda="lambda.1se",nperms=N),
  #                                  error=function(e) {return(NA)})
  # pr.p2.permN.1se = pr.perm.results.1se$N
  # pr.p2.sumbeta.p.1se = pr.perm.results.1se$pvalue
  # pr.p2.index.mse.1se = mean((pr.p2.index.1se - index)^2) # get MSE
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
output = tryCatch(data.frame(iter=job,true.r2=true.r2,true.pairwise.cov=true.pairwise.cov,
                             t(pr.t1.betas),t(pr.t1.weights),t(pr.t1.select),pr.t1.size,pr.t1.sumbeta,pr.t1.sumbeta.SE,pr.t1.sumbeta.p,#pr.t1.index.mse,
                             pr.p1.permTime,pr.p1.permN,t(pr.p1.betas),t(pr.p1.weights),t(pr.p1.select),pr.p1.size,pr.p1.sumbeta,pr.p1.sumbeta.p,#pr.p1.index.mse,
                             #t(pr.t1.betas.1se),t(pr.t1.weights.1se),t(pr.t1.select.1se),pr.t1.size.1se,pr.t1.sumbeta.1se,pr.t1.sumbeta.SE.1se,pr.t1.sumbeta.p.1se,pr.t1.index.mse.1se,
                             #pr.p1.permN.1se,t(pr.p1.betas.1se),t(pr.p1.weights.1se),t(pr.p1.select.1se),pr.p1.size.1se,pr.p1.sumbeta.1se,pr.p1.sumbeta.p.1se,pr.p1.index.mse.1se,
                             wqs.t.bootTime,t(wqs.t.weights),wqs.t.sumbeta,wqs.t.sumbeta.SE,wqs.t.sumbeta.p,#wqs.t.index.mse,
                             wqs.n.bootTime,t(wqs.n.weights),wqs.n.sumbeta,wqs.n.sumbeta.SE,wqs.n.sumbeta.p,#wqs.n.index.mse,
                             t(pr.t2.betas),t(pr.t2.weights),t(pr.t2.select),pr.t2.size,pr.t2.sumbeta,pr.t2.sumbeta.SE,pr.t2.sumbeta.p,#pr.t2.index.mse,
                             pr.p2.permTime,pr.p2.permN,t(pr.p2.betas),t(pr.p2.weights),t(pr.p2.select),pr.p2.size,pr.p2.sumbeta,pr.p2.sumbeta.p,#pr.p2.index.mse,
                             #t(pr.t2.betas.1se),t(pr.t2.weights.1se),t(pr.t2.select.1se),pr.t2.size.1se,pr.t2.sumbeta.1se,pr.t2.sumbeta.SE.1se,pr.t2.sumbeta.p.1se,pr.t2.index.mse.1se,
                             #pr.p2.permN.1se,t(pr.p2.betas.1se),t(pr.p2.weights.1se),t(pr.p2.select.1se),pr.p2.size.1se,pr.p2.sumbeta.1se,pr.p2.sumbeta.p.1se,pr.p2.index.mse.1se,
                             t(one.betas),t(one.p),t(all.betas),t(all.p),all.mix.p), error=function(e) {return(NULL)})
if (!is.null(output)) {write.table(output,file=paste0(path,"results.",c,".",job,".txt"),append=F,row.names=FALSE,
                                   col.names=T)}

# stop the clock
# time.end=proc.time()[[3]]
# time.end-time.start
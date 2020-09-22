## Functions to fit penalized constrained regression to estimate weights and 
# perform a permutation test to assess significance of the overall mixture effect ##

library(dplyr) #,lib.loc="~/RPackages")
library(glmnet) #,lib.loc="~/RPackages")
library(parallel) #,lib.loc="~/RPackages")

### cv.fit ### 
# This is the function that performs the estimation step. It can be applied to a training set,
# after which weights corresponding to an optimal lambda are extracted and used to compute a weighted sum of mixture 
# components in a held-out test set. Or, the entire dataset can be used to estimate weights, followed by a permutation test.
# For each value of lambda, the tuning parameter, the sign of the overall mixture effect is estimated, along with the 
# mixture component weights and several metrics associated with the K-fold cross-validation, as defined in glmnet.

## inputs:
# x: matrix including (potentially quantiled) chemical mixture variables and non-mixture covariates
# y: continuous outcome vector
# standardize: logical; whether to standardize predictors in the optimization; default is true
# cov.cols: column indices for non-mixture covariates (which are not subject to penalization)
# seed: set a seed to get consistent results from cv.glmnet, which picks K folds randomly
# alpha: default is 1 (L1 penalty); alternately, set to 0 for L2 penalty
# intercept: default is true; set to false if y is centered and intercept is not desired
cv.fit = function(x,y,standardize=T,cov.cols=NULL,seed=123,alpha=1,intercept=T) {
  # set seed
  set.seed(seed)
  
  # make x a matrix
  x = as.matrix(x)
  
  # remove penalty from non-mixture covariates
  penalties = rep(1,dim(x)[2])
  if (!is.null(cov.cols)) {
    penalties[cov.cols] = 0
  }
  
  # remove sign restrictions from non-mixture covariates
  up.lims = c(rep(0,dim(x)[2]))
  lo.lims = c(rep(0,dim(x)[2]))
  if (!is.null(cov.cols)) {
    up.lims[cov.cols] = Inf
    lo.lims[cov.cols] = -Inf
  }
  
  # fit models with both sign constraints
  neg.mod = cv.glmnet(x=x,y=y,alpha=alpha,intercept=intercept,family="gaussian",
                      standardize=standardize,upper.limits=up.lims,penalty.factor=penalties)
  pos.mod = cv.glmnet(x=x,y=y,alpha=alpha,intercept=intercept,family="gaussian",
                      standardize=standardize,lower.limits=lo.lims,penalty.factor=penalties)
  
  # check that the two fits use the same path of 100 lambdas
  if (sum(round(neg.mod$glmnet.fit$lambda,5) != round(pos.mod$glmnet.fit$lambda,5)) != 0) {
    print("Lambda paths do not match.")
  } 

  # for each lambda, choose sign to minimize RSS
  design = cbind(rep(1,dim(x)[1]),x)
  neg.betas = rbind(neg.mod$glmnet.fit$a0,neg.mod$glmnet.fit$beta)
  pos.betas = rbind(pos.mod$glmnet.fit$a0,pos.mod$glmnet.fit$beta)
  neg.rss = sapply(1:100, function(i) {
    rss = sum((y - design %*% neg.betas[,i])^2)
    return(rss)
  })
  pos.rss = sapply(1:100, function(i) {
    rss = sum((y - design %*% pos.betas[,i])^2)
    return(rss)
  })
  signs = ifelse(neg.rss < pos.rss,"neg","pos")

  # get cvm, cvsd, cvup, cvlo, size
  neg.cv = cbind.data.frame(lambda=neg.mod$lambda,neg.cvm=neg.mod$cvm,
                            neg.cvlo=neg.mod$cvlo,neg.cvup=neg.mod$cvup,
                            neg.cvsd=neg.mod$cvsd)
  pos.cv = cbind.data.frame(lambda=pos.mod$lambda,pos.cvm=pos.mod$cvm,
                            pos.cvlo=pos.mod$cvlo,pos.cvup=pos.mod$cvup,
                            pos.cvsd=pos.mod$cvsd)
  lambda.path = data.frame(lambda=neg.mod$glmnet.fit$lambda)
  cv = merge(lambda.path,neg.cv,by="lambda",all.x=T) %>%
    merge(.,pos.cv,by="lambda",all.x=T) %>%
    arrange(desc(lambda))
  cvm = c(); cvsd = c(); cvup = c(); cvlo = c(); size = c()
  cvm[signs=="neg"] = cv$neg.cvm[signs=="neg"]
  cvm[signs=="pos"] = cv$pos.cvm[signs=="pos"]
  cvsd[signs=="neg"] = cv$neg.cvsd[signs=="neg"]
  cvsd[signs=="pos"] = cv$pos.cvsd[signs=="pos"]
  cvup[signs=="neg"] = cv$neg.cvup[signs=="neg"]
  cvup[signs=="pos"] = cv$pos.cvup[signs=="pos"]
  cvlo[signs=="neg"] = cv$neg.cvlo[signs=="neg"]
  cvlo[signs=="pos"] = cv$pos.cvlo[signs=="pos"]

  # get lambda.min
  lambda.min = cv$lambda[which.min(cvm)]
  
  # get lambda.1se
  lambda.1se = max(cv$lambda[cvm < cvup[which.min(cvm)]],na.rm=T)
  
  # all.fit - uses all 100 lambdas from glmnet.fit
  # get beta matrix, weight matrix, sumbeta vector, size vector
  beta.matrix = neg.betas
  beta.matrix[,which(signs=="pos")] = pos.betas[,which(signs=="pos")]
  sumbeta = if (!is.null(cov.cols)) {
    colSums(beta.matrix[-c(1,cov.cols+1),])} else {
      colSums(beta.matrix[-1,])}
  size = if (!is.null(cov.cols)) {
    colSums(beta.matrix[-c(1,cov.cols+1),] != 0)} else {
      colSums(beta.matrix[-1,] != 0)}
  weight.matrix = sapply(1:100,function(j) {
    mix.betas = if (!is.null(cov.cols)) {
      beta.matrix[-c(1,cov.cols+1),j]} else {
        beta.matrix[-1,j]
      }
    weights = if (sumbeta[j]==0) {
      rep(0,dim(x)[2]-length(cov.cols))} else {
        mix.betas / sumbeta[j]}
    return(weights)
  })
  rownames(weight.matrix) = rownames(beta.matrix)[if (!is.null(cov.cols)) {-c(1,cov.cols+1)}
                                                  else {-1}]
  all.fit = list(signs=signs,sumbeta=unname(sumbeta),size=unname(size),lambda=cv$lambda,
                  beta.matrix=beta.matrix,weight.matrix=weight.matrix)

  output = list(all.fit=all.fit,
                cvm=cvm,cvsd=cvsd,cvup=cvup,cvlo=cvlo,
                lambda.min=lambda.min,
                lambda.1se=lambda.1se)
  return(output)
}

### cv.refit ### 
# This function is used repeatedly in the permutation test. It performs the estimation
# and least-squares refitting steps on the same data and returns the estimated overall mixture effect
# from the least-squares regression, for the specified lambda (min or 1se) and penalty (L1 or L2).

## inputs:
# same as before, but with a choice of tuning lambda to either lambda.1se or lambda.min
cv.refit = function(x,y,standardize=T,cov.cols=NULL,seed=123,alpha=1,intercept=T,
                      lambda=c("lambda.1se","lambda.min")) {
  x = as.matrix(x)
  
  # fit model
  mod.fit = cv.fit(x=x,y=y,standardize=standardize,cov.cols=cov.cols,seed=seed,
                   alpha=alpha,intercept=intercept)
  
  # get my lambda
  lambda = match.arg(lambda)
  if (lambda=="lambda.1se") {
    mylambda = which(mod.fit$all.fit$lambda==mod.fit$lambda.1se)
  } else if (lambda=="lambda.min") {
    mylambda = which(mod.fit$all.fit$lambda==mod.fit$lambda.min)
  }
  if (mylambda==1) {mylambda=2} # tune to one value lower if all mixture betas set to 0

  # get my weights at my lambda
  myweights = mod.fit$all.fit$weight.matrix[,mylambda]

  # compute weighted index
  myindex = if (!is.null(cov.cols)) {x[,-cov.cols] %*% myweights} else {
    x %*% myweights
  }
  
  # refit model; OLS of y~index
  if (!is.null(cov.cols)) {
    cov.names = colnames(x)[cov.cols]
    ols.fit = lm(formula(paste0("y~myindex+",paste0(cov.names,collapse="+"))),
                 data = cbind.data.frame(y,x,myindex)) } 
  else {
    ols.fit = lm(y~myindex,data = cbind.data.frame(y,myindex))
  }
  
  # save sumbeta-hat
  sumbeta.hat = summary(ols.fit)$coef[2,1]
  return(sumbeta.hat)
}

### perm.test ###
# This function executes a permutation test for the null hypothesis of no overall mixture effect.
# The entire null distribution is returned, along with a p-value and the number of permutations.

## new inputs:
# teststat: coefficient for the weighted sum from a linear regression in the data used to estimate the weights
# nperms: number of permutations; default is 1000
# ncores: number of cores to be used in parallel processing, if desired; default is no parallel processing
perm.test = function(teststat,x,y,standardize=T,cov.cols=NULL,seed=123,alpha=1,intercept=T,
                     lambda=c("lambda.1se","lambda.min"),
                     nperms=1000,ncores=1) {
  x = as.matrix(x)
  
  # get residuals from linear regressions on covariates
  if (!is.null(cov.cols)) {
    z = x[,cov.cols]
    mod1 = lm(y ~ z); yz = resid(mod1)
    mod2 = lm(x[,-cov.cols] ~ z); xz = resid(mod2)
  }
  
  # run model for each permutation to get a null distribution
  null_dist = mclapply(1:nperms,function(i) {
    set.seed(1992+i)
    perm = sample(nrow(x))
    if (is.null(cov.cols)) {
      perm.x = x[perm,]
      b1 = tryCatch(cv.refit(perm.x,y,standardize=standardize,cov.cols=NULL,seed=seed,
                      alpha=alpha,intercept=intercept,lambda=lambda), # no cov.cols
                    error=function(e) {return(NA)})
    } else {
      perm.xz = xz[perm,]
      b1 = tryCatch(cv.refit(perm.xz,yz,standardize=standardize,cov.cols=NULL,seed=seed,
                      alpha=alpha,intercept=intercept,lambda=lambda), # also no cov.cols
                    error=function(e) {return(NA)})
    }
    return(b1)
  },mc.cores=ncores)
  
  # format results of permutation test
  sumbetas = unlist(null_dist)
  
  # get p-value
  sumbeta.pvalue = ( sum(sumbetas >= abs(teststat), na.rm=T) + 
                       sum(sumbetas <= -abs(teststat), na.rm=T) + 1 ) / 
    (sum(!is.na(sumbetas)) + 1)
  
  return(list(teststat=teststat,null_dist=sumbetas,pvalue=sumbeta.pvalue,N=sum(!is.na(sumbetas))))
}
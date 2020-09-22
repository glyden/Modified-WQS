# Example use of penalized weighted sum regression with optional permutation test

source("~/Desktop/Github/functions.R") # functions to help create data w/ desired correlation, signal
source("~/Desktop/Github/method-functions.R") # functions to fit the method

library(mvtnorm)

# generate data for example
true.pairwise.cov=0.5
n=300; p=20; my.mu = rep(10,p)
betas = c(10,rep(3,10),rep(0,10))
my.sigma = make.sigma(true.pairwise.cov,p)
x = rmvnorm(n,my.mu,my.sigma)
colnames(x) = paste0("x",1:p)
x.design = cbind(rep(1,n),x)
y = x.design %*% betas + rnorm(n,0,20)

# standardize x
x = apply(x,2,function(vect) return(vect/sd(vect)))

# draw training / test indicators
test = sample(x=1:n,size=floor(.6*n),replace=F)
x.train = x[-test,]
y.train = y[-test]
x.test = x[test,]
y.test = y[test]

### first option ###
### fit our method to training data with L1 penalty, then do inference in test data ###
L1.split.fit = cv.fit(x.train,y.train)

# get weights at lambda.min
L1.split.weights = L1.split.fit$all.fit$weight.matrix[,which(L1.split.fit$all.fit$lambda==L1.split.fit$lambda.min)]

# estimate weighted sum in test data
L1.split.index = x.test %*% L1.split.weights

# regress test outcome on weighted sum
L1.split.ols = lm(y.test ~ L1.split.index)
summary(L1.split.ols)

### second option ###
### fit our method to entire dataset with L1 penalty, then do permutation test ###
L1.perm.fit = cv.fit(x,y)

# get weights at lambda.min
L1.perm.weights = L1.perm.fit$all.fit$weight.matrix[,which(L1.perm.fit$all.fit$lambda==L1.perm.fit$lambda.min)]

# estimate weighted sum for entire dataset
L1.perm.index = x %*% L1.perm.weights # estimate index for all x

# regress y on weighted sum
L1.perm.ols = lm(y ~ L1.perm.index)

# get mixture coefficient
L1.perm.sumbeta = summary(L1.perm.ols)$coef[2,1]

# run permutation test
# note: in practice, use larger N, such as 1000
perm.results = perm.test(L1.perm.sumbeta,x,y,lambda="lambda.min",nperms=10)
perm.results$pvalue
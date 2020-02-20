### USEFUL FUNCTIONS ###

# function to make exchangeable correlation matrix
make.sigma = function(pairwise.corr,p) {
  return(matrix(c(rep(c(1, rep(pairwise.corr,p)),p-1),1),nrow=p)) 
}

### For the remainder of this script: ###
# y is a linear combination of the x variables + error
# we wish to compute R^2, variance of Y, and covariance/correlation between y and one predictor

## inputs:
  # betas = vector, such that Y=XB+error
    # can ignore the intercept in betas and xmatrix, as it will not affect covariance/correlation.
  # which.x = integer; which x variable we want to compute cov/corr with
  # xsigma = matrix; covariance matrix for X vector
  # error.sigma2 = double; variance of the error terms

# variance due to regression
var.reg = function(betas,xsigma) {
  var.reg = 0
  for (i in 1:length(betas)) {
    var.reg = var.reg + (betas[i]^2)*xsigma[i,i]
    for (j in 1:length(betas)) {
      if (j>i) {var.reg = var.reg + 2*betas[i]*betas[j]*xsigma[i,j] }
    }
  }
  return(var.reg)
}

# variance of y
var.y = function(betas,xsigma,error.sigma2) {
  y.var = error.sigma2 + var.reg(betas,xsigma)
  return(y.var)
}

# R^2 of the model
r2 = function(betas,xsigma,error.sigma2) {
  y.var = var.y(betas,xsigma,error.sigma2)
  var.reg = y.var - error.sigma2
  my.r2 = var.reg/y.var
  return(my.r2)
}

# covariance b/w y and one of the x's
cov.y.x = function(betas,which.x,xsigma,error.sigma2) {
  y.var = var.y(betas,xsigma,error.sigma2)
  cov = 0
  for (i in 1:length(betas)) {
      cov = cov + betas[i]*xsigma[i,which.x]
  }
  corr = cov / (sqrt(y.var) * sqrt(xsigma[which.x,which.x]))
  return(cbind.data.frame(cov,corr))
}
#' Covariate adjustment with variable selection in randomized clinical trials
#' @param Y Vector of observed outcomes
#' @param A Vector of treatment assignments, which can be binary or multi-level, but values must be categorical or integers; if binary, we suggest to use 0 for the control and 1 for the treated group
#' @param trt.name Value of the active treatment in the treatment vector `A`; must exactly match the value type in `A`, e.g., 1 for number while "1" for character
#' @param ctrl.name Value of the control in the treatment vector `A`; must exactly match the value type in `A`, e.g., 0 for number while "0" for character
#' @param X Covariate vector (only one covariate) or matrix/data frame (multiple covariates)
#' @param conti.out Specify whether the outcome is continous, `TRUE` or `FALSE`, with no default value
#' @param var.sel.method Variable selection method for COADVISE stage I, including the following available options:
#'                       `No` (Use all covariates without selection), `Lasso` (Lasso, which is the default choice), `A.Lasso` (Adaptive Lasso),
#'                       `Corr.k` (Select k covariates with highest correlations with the outcome; k is between 1 and total number of covariates),
#'                       `Corr.xi` (Select all covariates with correlation > xi with the outcome; xi is between 0 and 1),
#'                       `Pre.test` (A preliminary two sample t-test on each single covariate)
#' @param lasso.family GLM model family for Lasso (`Lasso`) variable selection based on outcome type, chosen from `gaussian`, `binomial`, `poisson`, and `multinomial`; default is `gaussian`
#' @param A.lasso.family GLM model family for adaptive Lasso (`A.Lasso`) variable selection based on outcome type, chosen from `gaussian`, `binomial`, `poisson`, and `multinomial`; default is `gaussian`
#' @param out1.model.aipw Model for E(Y(1)|X) in AIPW estimator. We allow the use of the following models:
#'                        (i) For continuous outcome (but can also for other types): `linear` (the linear model);
#'                        (ii) For binary outcome: `logit`, `probit`, `log`, `cloglog`, `identity`;
#'                        (iii) For categorical outcome more than two categories: `poisson`, `multilogi` (multinomial logistic), `ordlogi` (ordinal logistic)
#' @param out0.model.aipw Model for E(Y(0)|X) in AIPW estimator. Same argument choices as `out1.model.aipw`
#' @param k Number of selected covariates, only if the `Corr.k` variable selection method is used; default is `1`
#' @param xi Threshold of marginal correlation, only if the `Corr.xi` variable selection method is used; default is `0.25`
#' @param pre.alpha Confidence (significance) level of preliminary test; default is `0.05` (for a 0.95 confidence interval)
#' @param conf.level Confidence (significance) level of confidence interval of the ATE; default is `0.05` (for a 0.95 confidence interval)
#' @param MI.method Missing data imputation method if there is missing values in any outcome or covariates;
#'                  notice that the treatment is not allowed (or assumed) to be missing in RCT;
#'                  values chosen from: `cc` (complete-case),
#'                  `mice` (multiple imputation by chain equation),
#'                  `missForest` (random forest),
#'                  `missInd` (missingness indicator method) for missing covariates only, and
#'                  `ipw` (inverse probability weighting) for missing outcome only.
#'                  If there is no missing data, no method will be implemented even the value is specified.
#' @param seed Seed for generating random numbers, when using Lasso or adaptive Lasso; default is `4399`
#' @return The function returns a list of three components: (i) `df.fit`: results on ATE estimates, standard errors, CIs and p-values by
#'         simple, ANCOVA, ANHECOVA and AIPW estimators under the chosen variable selection method; (ii) `stage1.covars`: a list of covariates
#'         data frames by the first variable selection step, including the variable selection for all ANCOVA, ANHECOVA and AIPW methods;
#'         (iii) `AIPW.out.means`: the vector of outcome means of the two comparison groups and the covariance matrix.
Coadvise <- function(Y,
                     A,
                     trt.name,
                     ctrl.name,
                     X=NA,
                     conti.out,
                     var.sel.method="Lasso",
                     lasso.family="gaussian",
                     A.lasso.family="gaussian",
                     out1.model.aipw="linear",
                     out0.model.aipw="linear",
                     k=1,
                     xi=0.25,
                     pre.alpha=0.05,
                     conf.level=0.05,
                     MI.method="cc",
                     seed=4399) {

  set.seed(seed=seed)

  if(!require("dplyr")) install.packages("dplyr")
  if(!require("glmnet")) install.packages("glmnet")
  if(!require("MASS")) install.packages("MASS")
  if(!require("nnet")) install.packages("nnet")
  library(dplyr)
  library(glmnet)

  if(!conti.out%in%c(TRUE, FALSE)) {
    stop("Please specify whether the outcome variable is continuous:
         if yes--conti.out=TRUE; if no--conti.out=FALSE")
  }

  if(!conti.out) {
    Y <- as.factor(Y)
  }

  if(sum(is.na(A))!=0) {
    stop("Missing values in the treatment vector is not supported!")
  }

  group.ind <- which(A%in%c(trt.name, ctrl.name))
  A <- as.character(A[group.ind])
  A[A==as.character(trt.name)] <- "1"
  A[A==as.character(ctrl.name)] <- "0"
  A <- as.numeric(A)
  Y <- Y[group.ind]
  X <- X[group.ind,]

  if(sum(is.na(X))!=0 | sum(is.na(Y))!=0) {
    if(!MI.method%in%c("cc", "mice", "missInd", "missForest")) {
      stop("Missing data found in covariates or outcome!
           Please select an imputation method from cc, mice, missInd, missForest,
           or please impute missing values before analysis.")
    }

    if(MI.method=="cc") {
      # complete-case analysis
      dat.miss <- cbind(Y, A, X)
      dat.impu <- na.omit(dat.miss)
      Y <- dat.impu[,1]
      X <- dat.impu[,-c(1,2)]
    }

    if(MI.method=="mice") {
      # multiple imputation by chained equation
      if(!require("mice")) install.packages("mice")
      library(mice)
      dat.miss <- data.frame(Y, A, X)
      dat.impu <- mice(dat.miss, m=5, method="pmm", maxit=50, seed=500)
      dat.impu <- complete(dat.impu, action=1)
      Y.impu <- dat.impu[,1]
      X <- dat.impu[,-c(1,2)]
    }

    if(MI.method=="missForest") {
      # missing data imputation by random forest
      if(!require("missForest")) install.packages("missForest")
      library(missForest)
      dat.miss <- data.frame(Y, A, X)
      dat.impu <- missForest(dat.miss)$ximp
      Y.impu <- dat.impu[,1]
      X <- dat.impu[,-c(1,2)]
    }

    if(MI.method=="ipw") {
      if(sum(is.na(X))!=0) {
        stop("Missing data found in covariates! The IPW method cannot handle missing covariates.
              Please consider other methods or use only complete-covariate data.")
      }
      dat.miss <- data.frame(Y, A, X)
      obs_md <- glm(!is.na(Y)~A+X, family=binomial(), data=dat.miss)
      p_obs <- predict(obs_md, type="response")
      ipw <- ifelse(!is.na(Y), 1/p_obs, NA)
      Y[is.na(Y)] <- weighted.mean(Y, ipw, na.rm=T)
    }

    if(MI.method=="missInd") {
      if(sum(is.na(Y))!=0) {
        stop("Missing data found in outcome! The missInd method cannot handle missing outcomes.
              Please impute the outcome by other methods first or use only complete-outcome data.")
      }
      n.covar <- ncol(X)
      W <- c()
      for(i in 1:n.covar) {
        if(sum(is.na(X[,i]))!=0) {
          W <- cbind(W, as.numeric(is.na(X[,i])))
          X[,i][is.na(X[,i])] <- 0
        }
      }
      X <- cbind(W, X)
    }
  } else {
    Y.impu <- Y
  }

  ### centering covariates to their sample mean
  X <- as.matrix(X)
  X <- apply(X, 2, function(u) u-mean(u))
  AX <- A*X
  X1 <- X[A==1,]
  X0 <- X[A==0,]
  X.names <- colnames(X1) <- colnames(X0) <- colnames(X) <- paste("X.", 1:ncol(X), sep="")
  colnames(AX) <- paste("AX.", 1:ncol(AX), sep="")
  Y <- as.numeric(Y.impu)
  if(sum(levels(Y.impu)==c("0","1"))==2) {
    Y[Y.impu=="0"] <- 0
    Y[Y.impu=="1"] <- 1
  }
  Y1 <- Y[A==1]
  Y0 <- Y[A==0]
  n <- length(A)

  ### Critical value for the confidence interval (quantile of standard normal)
  quant <- qnorm(1-conf.level/2)

  ### Simple estimator (unadjusted)
  tau.simple <- mean(Y[A==1])-mean(Y[A==0])
  se.simple <- sqrt(var(Y[A==1])/sum(A)+var(Y[A==0])/sum(1-A))
  ci.lwr.simple <- tau.simple - quant*se.simple
  ci.upr.simple <- tau.simple + quant*se.simple
  p.simple <- 2*(1-pnorm(abs(tau.simple/se.simple)))

  ### Estimators by different variable selections
  if(var.sel.method=="No") {
    ind <- ind1 <- ind0 <- 1:ncol(X)
  }
  if(var.sel.method=="Lasso") {
    cv.lasso <- cv.glmnet(x=X, y=Y, alpha=1, family=lasso.family)
    lambda <- cv.lasso$lambda.min
    lasso <- glmnet(x=X, y=Y, alpha=1, lambda=lambda, family=lasso.family)
    ind <- which(as.numeric(lasso$beta)!=0)

    cv.lasso1 <- cv.glmnet(x=X1, y=Y1, alpha=1, family=lasso.family)
    lambda <- cv.lasso1$lambda.min
    lasso1 <- glmnet(x=X1, y=Y1, alpha=1, lambda=lambda, family=lasso.family)
    ind1 <- which(as.numeric(lasso1$beta)!=0)

    cv.lasso0 <- cv.glmnet(x=X0, y=Y0, alpha=1, family=lasso.family)
    lambda <- cv.lasso0$lambda.min
    lasso0 <- glmnet(x=X0, y=Y0, alpha=1, lambda=lambda, family=lasso.family)
    ind0 <- which(as.numeric(lasso0$beta)!=0)
  }
  if(var.sel.method=="A.Lasso") {
    cv.ridge <- cv.glmnet(x=X, y=Y, alpha=0, family=A.lasso.family)
    w3 <- 1/abs(matrix(stats::coef(cv.ridge, s=cv.ridge$lambda.min)[,1][2:(ncol(X)+1)]))^1
    w3[w3[,1]==Inf] <- 999999999
    alasso_cv <- cv.glmnet(x=X, y=Y, type.measure="mse", nfold=5, alpha=1, penalty.factor=w3, keep=TRUE, family=A.lasso.family)
    alasso <- stats::coef(alasso_cv, s=alasso_cv$lambda.min, family=A.lasso.family)
    ind <- which(as.numeric(alasso)!=0)[-1]-1

    cv.ridge <- cv.glmnet(x=X1, y=Y1, alpha=0, family=A.lasso.family)
    w3 <- 1/abs(matrix(stats::coef(cv.ridge, s=cv.ridge$lambda.min)[,1][2:(ncol(X1)+1)]))^1
    w3[w3[,1]==Inf] <- 999999999
    alasso_cv1 <- cv.glmnet(x=X1, y=Y1, type.measure="mse", nfold=5, alpha=1, penalty.factor=w3, keep=TRUE, family=A.lasso.family)
    alasso1 <- stats::coef(alasso_cv1, s = alasso_cv1$lambda.min, family=A.lasso.family)
    ind1 <- which(as.numeric(alasso1)!=0)[-1]-1

    cv.ridge <- cv.glmnet(x=X0, y=Y0, alpha=0, family=A.lasso.family)
    w3 <- 1/abs(matrix(stats::coef(cv.ridge, s=cv.ridge$lambda.min)[,1][2:(ncol(X0)+1)]))^1
    w3[w3[,1]==Inf] <- 999999999
    alasso_cv0 <- cv.glmnet(x=X0, y=Y0, type.measure="mse", nfold=5, alpha=1, penalty.factor=w3, keep=TRUE, family=A.lasso.family)
    alasso0 <- stats::coef(alasso_cv0, s = alasso_cv0$lambda.min, family=A.lasso.family)
    ind0 <- which(as.numeric(alasso0)!=0)[-1]-1
  }
  if(var.sel.method=="Corr.k") {
    if(k<1) stop("k should must be a positive integer!")
    if(k>ncol(X)) stop("k is larger than the number of all covariates!")
    corr <- as.numeric(abs(cor(X, Y)))
    ind <- which(rank(corr)>=ncol(X)-(k-1))

    corr1 <- as.numeric(abs(cor(X1, Y1)))
    ind1 <- which(rank(corr1)>=ncol(X1)-(k-1))

    corr0 <- as.numeric(abs(cor(X0, Y0)))
    ind0 <- which(rank(corr0)>=ncol(X0)-(k-1))
  }
  if(var.sel.method=="Corr.xi") {
    if(xi<=0 | xi>=1) stop("xi must be a number within (0,1)!")
    corr <- as.numeric(abs(cor(X, Y)))
    ind <- which(corr>xi)

    corr1 <- as.numeric(abs(cor(X1, Y1)))
    ind1 <- which(corr1>xi)

    corr0 <- as.numeric(abs(cor(X0, Y0)))
    ind0 <- which(corr0>xi)
  }
  if(var.sel.method=="Pre.test") {
    test_x <- c()
    for(j in 1:ncol(X)) {
      x.pre <- X[,j]
      x1 <- x.pre[A==1]
      x0 <- x.pre[A==0]
      x.diff <- mean(x1)-mean(x0)
      sd.x.diff <- sqrt(var(x1)/length(x1)+var(x0)/length(x0))
      test_x <- c(test_x, abs(x.diff/sd.x.diff) > qt(1-pre.alpha/2, df=length(Y)-2))
    }
    ind <- ind1 <- ind0 <- ifelse(sum(test_x)>0, which(test_x), 0)
  }
  fit.ancova <- lm(Y~as.matrix(cbind(A, X[,ind])))
  fit.anhecova <- lm(Y~as.matrix(cbind(A, X[,ind], AX[,ind])))
  tau.ancova <- summary(fit.ancova)$coefficients[2,1]
  tau.anhecova <- summary(fit.anhecova)$coefficients[2,1]
  x.ancova <- cbind(1, A, X[,ind])
  x.anhecova <- cbind(1, A, X[,ind], AX[,ind])
  W <- as.matrix(X[,ind])

  x <- x.ancova
  sigma <- try(solve(t(x)%*%x)%*%t(x)%*%diag((fit.ancova$residuals)^2)%*%x%*%solve(t(x)%*%x))
  if(class(sigma)[1]!="try-error") { sigma <- sigma } else { sigma <- NA }
  if(sum(is.na(sigma))==0) { se.ancova <- sqrt(sigma[2,2]) } else { se.ancova <- NA }
  p.ancova <- 2*(1-pnorm(abs(tau.ancova/(se.ancova))))
  ci.lwr.ancova <- tau.ancova - quant*se.ancova
  ci.upr.ancova <- tau.ancova + quant*se.ancova

  x <- x.anhecova
  sigma <- try(solve(t(x)%*%x)%*%t(x)%*%diag((fit.anhecova$residuals)^2)%*%x%*%solve(t(x)%*%x))
  if(ncol(x)==2) {
    if(class(sigma)[1]!="try-error") { se.anhecova <- sqrt(sigma[2,2]) } else { se.anhecova <- NA }
  }
  if(ncol(x)>2) {
    W1 <- W[A==1,]
    W0 <- W[A==0,]
    fit1 <- summary(lm(Y1~as.matrix(W1)))$coefficients[-1,1]
    fit0 <- summary(lm(Y0~as.matrix(W0)))$coefficients[-1,1]
    sigma.add <- try(t(fit1-fit0)%*%cov(W)%*%(fit1-fit0) / n)
    if(class(sigma)[1]!="try-error" & class(sigma.add)[1]!="try-error") {
      sigma <- sigma
      sigma.add <- sigma.add
    } else {
      sigma <- NA
      sigma.add <- NA }
    if(sum(is.na(c(sigma, sigma.add)))==0) {
      se.anhecova <- sqrt(sigma[2,2]+sigma.add) } else { se.anhecova <- NA }
  }
  p.anhecova <- 2*(1-pnorm(abs(tau.anhecova/(se.anhecova))))
  ci.lwr.anhecova <- tau.anhecova - quant*se.anhecova
  ci.upr.anhecova <- tau.anhecova + quant*se.anhecova

  X1.AIPW <- as.matrix(X1[, ind1])
  colnames(X1.AIPW) <- X.names[ind1]
  X0.AIPW <- as.matrix(X0[, ind0])
  colnames(X0.AIPW) <- X.names[ind0]

  result <- .AIPW(Y1=Y1, Y0=Y0, X1=X1.AIPW, X0=X0.AIPW, A=A, X.pred=X, Y=Y,
                  out0.model=out0.model.aipw, out1.model=out1.model.aipw)
  tau.aipw <- result$ATE$tau
  se.aipw <- result$ATE$se
  p.aipw <- 2*(1-pnorm(abs(tau.aipw/se.aipw)))
  ci.lwr.aipw <- tau.aipw - quant*se.aipw
  ci.upr.aipw <- tau.aipw + quant*se.aipw

  df.fit <- data.frame(method=c("Simple", "ANCOVA", "ANHECOVA", "AIPW"),
                       tau=c(tau.simple, tau.ancova, tau.anhecova, tau.aipw),
                       se=c(se.simple, se.ancova, se.anhecova, se.aipw),
                       ci.lwr=c(ci.lwr.simple, ci.lwr.ancova, ci.lwr.anhecova, ci.lwr.aipw),
                       ci.upr=c(ci.upr.simple, ci.upr.ancova, ci.upr.anhecova, ci.upr.aipw),
                       p=c(p.simple, p.ancova, p.anhecova, p.aipw))

  stage1.covars <- list(ANCOVA=x.ancova[,-c(1,2)],
                        ANHECOVA=x.anhecova[,-c(1,2)],
                        AIPW.ctrl=X0.AIPW, AIPW.trt=X1.AIPW)

  AIPW.out.means <- list(mu=result$mu, sigma=result$sigma)

  return(list(df.fit=df.fit,
              stage1.covars=stage1.covars,
              AIPW.out.means=AIPW.out.means))

}

.AIPW <- function(Y0, Y1, X1, X0, A, X.pred, Y,
                  out0.model="linear", out1.model="linear") {

  n <- length(A)
  pi1 <- mean(A)
  pi0 <- mean(1-A)

  ### if no covariate is selected in treated and/or controls
  if(ncol(X1)==0) {
    m1.h <- rep(0, n)
    m1.h0 <- rep(0, length(Y0))
    m1.h1 <- rep(0, length(Y1))
  }
  if(ncol(X0)==0) {
    m0.h <- rep(0, n)
    m0.h0 <- rep(0, length(Y0))
    m0.h1 <- rep(0, length(Y1))
  }

  ### Y(1) model prediction
  if(ncol(X1)>0) {
    if(out1.model=="linear") {
      fit1 <- lm(Y1~.-Y1, data=data.frame(Y1, X1))
      m1.h <- predict(fit1, data.frame(X.pred))
    }
    if(out1.model%in%c("logit", "probit", "log", "cloglog", "identity")) {
      fit1 <- glm(Y1~.-Y1, data=data.frame(Y1, X1), family=binomial(link=out1.model))
      m1.h <- predict(fit1, data.frame(X.pred), type="response")
    }
    if(out1.model=="poisson") {
      fit1 <- glm(Y1~.-Y1, data=data.frame(Y1, X1), family=poisson(link="log"))
      m1.h <- predict(fit1, data.frame(X.pred), type="response")
    }
    if(out1.model=="multilogi") {
      fit1 <- nnet::multinom(Y1~.-Y1, data=data.frame(Y1, X1))
      m1.h <- predict(fit1, data=data.frame(X.pred))
    }
    if(out1.model=="ordlogi") {
      fit1 <- MASS::polr(Y1~.-Y1, data=data.frame(Y1, X1), method="logistic")
      m1.h <- predict(fit1, data=data.frame(X.pred))
    }
    m1.h1 <- m1.h[A==1]
    m1.h0 <- m1.h[A==0]
  }

  ### Y(0) model prediction
  if(ncol(X0)>0) {
    if(out0.model=="linear") {
      fit0 <- lm(Y0~.-Y0, data=data.frame(Y0, X0))
      m0.h <- predict(fit0, data.frame(X.pred))
    }
    if(out0.model%in%c("logit", "probit", "log", "cloglog", "identity")) {
      fit0 <- glm(Y0~.-Y0, data=data.frame(Y0, X0), family=binomial(link=out0.model))
      m0.h <- predict(fit0, data.frame(X.pred), type="response")
    }
    if(out0.model=="poisson") {
      fit0 <- glm(Y0~.-Y0, data=data.frame(Y0, X0), family=poisson(link="log"))
      m0.h <- predict(fit0, data.frame(X.pred), type="response")
    }
    if(out0.model=="multilogi") {
      fit0 <- nnet::multinom(Y0~.-Y0, data=data.frame(Y0, X0))
      m0.h <- predict(fit0, data=data.frame(X.pred))
    }
    if(out0.model=="ordlogi") {
      fit0 <- MASS::polr(Y0~.-Y0, data=data.frame(Y0, X0), method="logistic")
      m0.h <- predict(fit0, data=data.frame(X.pred))
    }
    m0.h1 <- m0.h[A==1]
    m0.h0 <- m0.h[A==0]
  }

  ### AIPW estimator and standard error calculation
  tau1 <- mean((Y-m1.h)*I(A==1)/pi1 + m1.h)
  tau0 <- mean((Y-m0.h)*I(A==0)/pi0 + m0.h)
  tau <- tau1 - tau0

  v11 <- var(Y1-m1.h1)/pi1 + 2*cov(Y1, m1.h1) - var(m1.h)
  v00 <- var(Y0-m0.h0)/pi0 + 2*cov(Y0, m0.h0) - var(m0.h)
  v10 <- cov(Y1, m0.h1) + cov(Y0, m1.h0) - cov(m1.h, m0.h)
  se <- sqrt(v11-2*v10+v00)/sqrt(n)

  return(list(ATE=data.frame(tau=tau, se=se),
              mu=data.frame(tau0=tau0, tau1=tau1),
              sigma=matrix(c(v00,v10,v10,v11), ncol=2)))
}

# Coadvise

`Coadvise` is an R package/software implementing a flexible framework of using different covariate adjustment methods with a number of variable selection methods to estimate the average treatment effect (ATE) in randomized clinical trials (RCTs). 

## Usage

To install the latest version of the R package from GitHub, please run following commands in R:

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("yiliu1998/Coadvise")
```

As an illustrative example, we generate the following RCT data. We first generate the treatment and 7 covariates:

```{r}
set.seed(1555)
n <- 500
A <- rbinom(n=n, size=1, prob=0.5)
X1 <- rnorm(n, 10)
X2 <- rt(n, df=5)
X3 <- rbeta(n, shape1=0.5, shape2=0.8)
X4 <- X1^3+X2
X5 <- X3*X1
X6 <- X2^2
X7 <- 6*X2*X3
```

We consider the following two outcomes (`Y1` is continuous and `Y2` is binary): 

```r
X <- cbind(X1, X2, X3, X4, X5, X6, X7)
Y1 <- 0.05*(X1+X2+X3+X4+X5+X6+X7-30) + 0.1*A*(X1+X7-X3-10)^2 + rnorm(n=n, sd=3)
Y2 <- rbinom(n, size=1, prob=1/(1+exp(-(5+0.5*(X2+X5+X6+X7-15) + A*(X1+X2-X3-10)))))
mean(Y2) # proportion of 1 in Y2
```

## Package Usage

The main function of the package is `Coadvise()`. The package includes a number of estimation methods for ATE: Simple (unadjusted), ANCOVA, ANHECOVA and AIPW estimators. They will all be output in a data frame with their point estimates, standard errors, confidence intervals, and p-values. 

```r
# load the package
if (!require("devtools")) install.packages("devtools")
library(devtools)
if(!require("Coadvise")) devtools::install_github("yiliu1998/Coadvise")
library(Coadvise)
```

First, we consider using Lasso selection on the continuous outcome `Y1`, with using linear models in any adjusted method. We specify `Lasso` as the value to the argument `var.sel.method`, and since we want to use linear model, the `lasso.family` is specified by `gaussian`. Then, we specify both `out1.model.aipw` and `out0.model.aipw` to be `linear` as well. The following code implements this case: 

```r
lasso.Y1 <- Coadvise(y=Y1, A=A, X=X, 
                     var.sel.method="Lasso", lasso.family="gaussian", 
                     out1.model.aipw="linear", out0.model.aipw="linear",
                     seed=1811) 
print(lasso.Y1)
```

```r
method      tau        se     ci.lwr   ci.upr           p
1   Simple 2.078197 1.4065720 -0.6786333 4.835028 0.139544128
2   ANCOVA 1.004802 0.3223110  0.3730844 1.636520 0.001823961
3 ANHECOVA 1.045475 0.3332517  0.3923133 1.698636 0.001705743
4     AIPW 1.059058 0.2892087  0.4922198 1.625897 0.000250334
```

Second, we consider using selection based on the marginal correlation (using `k=4`) on the binary outcome `Y2`, with using logistic regression models in any adjusted method. Here, `k` means we want to include k covariates having the highest k marginal correlations with the outcome. So, we specify `Corr.k` as the value to the argument `var.sel.method`, and since we want to use binary regression model, the `lasso.family` is specified by `binomial`. Then, we specify both `out1.model.aipw` and `out0.model.aipw` to be `logit` as well. The following code implements this case: 

```r
Corr.k.Y2 <- Coadvise(y=Y2, A=A, X=X, 
                     var.sel.method="Corr.k", k=4, 
                     out1.model.aipw="logit", out0.model.aipw="logit") 
print(Corr.k.Y2)
```

```r
method         tau         se     ci.lwr       ci.upr          p
1   Simple -0.03548776 0.04442372 -0.1225566  0.051581131 0.42437910
2   ANCOVA -0.05789624 0.03509563 -0.1266824  0.010889941 0.09901031
3 ANHECOVA -0.05920252 0.03490272 -0.1276106  0.009205558 0.08984518
4     AIPW -0.07151859 0.03120419 -0.1326777 -0.010359507 0.02190824
```

We also allow the use of adaptive Lasso (`A.Lasso`) and marginal correlation by specifying the threshold of correlation (`Corr.xi`, with specifying `xi` a number within interval (0,1)). We also allow many other choices of outcome models used in the AIPW estimator, besides `linear` and `logit` (logistic for binary outcome) used above, which include the below cases.  

* For binary outcomes: `cloglog`, `log`, `identity`, `probit` used in arguments `out1.model.aipw` and `out0.model.aipw`. 
* For multi-valued categorical (>=2 levels) outcomes: `poission` for Poisson regression (typically for count data), `multilogi` for multinomial logistic regression, and `ordlogi` for ordinal categorical outcomes (e.g., levels of education, grade, etc.). 

The above demonstration can also be found in the Rmd file here ([click here](https://github.com/yiliu1998/Coadvise/tree/main/vignettes)), which can be downloaded and run on your local. 

## Author Information
The R code is maintained by Yi Liu (Please feel free to reach out at yi.liu.biostat@gmail.com, if you have any questions). 

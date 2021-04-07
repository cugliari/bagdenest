
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bagdenest

-----

<!-- badges: start -->

<!-- badges: end -->

The goal of *bagdenest* is to provide density estimation by averaging
classical density estimators such as the histogram, the frequency
polygon and the kernel density estimators obtained over different
bootstrap samples of the original data.

![workflow](https://github.com/cugliari/bagdenest/actions/workflows/build.yml/badge.svg)  
[![cran](https://www.r-pkg.org/badges/version/bagdenest?color=green)](https://cran.r-project.org/package=bagdenest)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

-----

### License

This package is free and open source software, licensed under GPL-3.

-----

### Installation

  - You can install the released version of bagdenest from
    [CRAN](https://CRAN.R-project.org) with:

<!-- end list -->

``` r
install.packages("bagdenest")
```

  - And the development version from [GitHub](https://github.com/) with:

<!-- end list -->

``` r
# install.packages("devtools")
devtools::install_github("cugliari/bagdenest")
```

-----

### Example

This is a basic example which shows you how to display densities
estimator errors using *bagdenest*:

``` r
library(bagdenest)

simulations=function(n = 100, M = 10, B = 150){
  vec = c(1,3,5,7,8,11,10,13,17,19,20,21)
  AA = matrix(0, nrow = length(vec), ncol = 7)
  colnames(AA)=c("H","FP","Kde", "BagHist", "BagFP", "BagKde", "Rash")
  rownames(AA)=c("normal","chi2_10","MIX_1","MIx_2","CLAW","Triangular",
                 "DENS2_rigollet","uniforme","SmooComb","DistBim","DENS1_rigollet","MIX_Unif")
  for(i in 1:M){
    A = matrix(NA, nrow = length(vec), ncol = 7)
    for(ll in 1:length(vec)){
      dd = gendata(vec[ll],n)
      bopt = bropt(dd$train)$opt  
      zz = hist(dd$train,breaks=mybreaks(dd$train,nbr=bopt),plot=F)
      bhist = BagHistfp(xx=dd$train, grid=dd$test, B)
      modelrash = rash(dd$train, grid = dd$test, nbr = bopt, B)
      modelbagkde <- Bagkde(xx = dd$train, grid = dd$test, B)
      A[ll,]=c(error(dd$dobs,predict_hist(zz,dd$test))[1],
               error(dd$dobs,approxfun(x=zz$mids,y=zz$density)(dd$test))[1],
               error(dd$dobs,onekdeucv(dd$train,dd$test))[1],
               error(dd$dobs,bhist$bh)[1],
               error(dd$dobs,bhist$bhfp)[1],
               error(dd$dobs,modelbagkde)[1],
               error(dd$dobs,modelrash)[1])
      }
    AA=AA+A
  }
  AA=AA/M
  print(AA)
}
simulations(n = 100, M = 2, B = 15)
#>                        H         FP        Kde   BagHist     BagFP
#> normal         0.0071200 0.00438000 0.00574500 0.0198000 0.0044650
#> chi2_10        0.0003425 0.00011545 0.00008755 0.0007375 0.0001035
#> MIX_1          0.0187500 0.01236000 0.00360500 0.0227000 0.0216000
#> MIx_2          0.0066900 0.00641000 0.00475900 0.0055950 0.0128000
#> CLAW           0.0399500 0.02400000 0.01970000 0.0275000 0.0243500
#> Triangular     0.0275500 0.02154000 0.01861000 0.1111500 0.0230600
#> DENS2_rigollet 0.0611500 0.06270000 0.06635000 0.0562500 0.0624000
#> uniforme       0.2382000 0.14200000 0.10075000 0.3030000 0.0603000
#> SmooComb       0.0130500 0.00996000 0.00997500 0.0176000 0.0207000
#> DistBim        0.0381500 0.01520000 0.02325000 0.0248500 0.3570000
#> DENS1_rigollet 0.0572000 0.05610000 0.06080000 0.0396000 0.0633000
#> MIX_Unif       0.0200000 0.01141000 0.01860000 0.0237000 0.0240500
#>                  BagKde       Rash
#> normal         0.013760 0.00447000
#> chi2_10        0.000260 0.00011415
#> MIX_1          0.007930 0.00637500
#> MIx_2          0.004155 0.00472000
#> CLAW           0.040650 0.02700000
#> Triangular     0.057700 0.02135000
#> DENS2_rigollet 0.059850 0.06135000
#> uniforme       0.160100 0.10480000
#> SmooComb       0.013055 0.01214500
#> DistBim        0.028350 0.02415000
#> DENS1_rigollet 0.053750 0.05300000
#> MIX_Unif       0.024000 0.01286000
```

-----

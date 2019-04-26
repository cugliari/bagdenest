## File: simulaciones-tubos.r
## Description: Estudio de simulacion para evaluar calidad de predicciones por
##              tubos de los estimadores de densidad H, FP, Kde, BagHistFP, 
##              BagHist, Rash, ASH.

rm(list = ls())
source("functions2.R")

# LAS FUNCIONES ESTAN COMENTADAS PORQUE SE EXPORTAN DEL functions2.R
# tube_hist <- function(data, nbr = 50, B = 10, conf = c(0.05, 0.95)) {
#   # A chaque etape, on prend un nouveau jeu de donn?es, on construit un 
#   # histogramme, on pr?dit et on agr?ge.
# 
#   xx <- data$train
#   grille <- data$test
#   n <- length(xx)
#   mat <- matrix(ncol = length(grille), nrow = B)
# 
#   for (i in 1:B) {
#     xb  = xx[sample(n, replace = TRUE)]
#     mybreaks <- mybreaks(xb, nbr)
#     hs2 = hist(xb, breaks = mybreaks, plot = FALSE)
#     
#     mat[i, ] <- predict.hist(hh = hs2, x = grille)
#   }
# 
#   punctual <- colMeans(mat, na.rm = TRUE)
#   tube     <- baseboot(data$test, punctual, mat, conf = conf)
#   
#   return(list(punctual = punctual, bunch = mat, tube = tube)) 
# }
# 
# tube_fp <- function(data, nbr = 50, B = 10, conf = c(0.05, 0.95)) {
#   # A chaque etape, on prend un nouveau jeu de donn?es, on construit un 
#   # histogramme, on pr?dit et on agr?ge.
# 
#   xx <- data$train
#   grille <- data$test
#   n <- length(xx)
#   mat <- matrix(ncol = length(grille), nrow = B)
# 
#   for (i in 1:B) {
#     xb  = xx[sample(n, replace = TRUE)]
#     mybreaks <- mybreaks(xb, nbr)
#     hs2 = hist(xb, breaks = mybreaks, plot = FALSE)
#     
#     mids <- hs2$mids
#     delta <- diff(mybreaks)[1]
#     # This avoids NA on prediction for test data between the min on grille
#     # and the min obs value on xb
#     newmids <- c(min(mids) - delta, mids, max(mids) + delta) 
# 
#     mat[i, ] <- approxfun(x = newmids, y = c(0, hs2$density, 0))(grille)
#   }
# 
#   # NAs occur when test data falls beyond the FP's support. They may safetly
#   # be replaced by 0s
#   mat[is.na(mat)] <- 0 
# 
#   punctual <- colMeans(mat, na.rm = TRUE)
#   tube     <- baseboot(data$test, punctual, mat, conf = conf)
#   
#   return(list(punctual = punctual, bunch = mat, tube = tube)) 
# }
# 
# tube_KDE <- function(data, nbr = 50, B = 10, conf = c(0.05, 0.95)) {
#   # A chaque etape, on prend un nouveau jeu de donn?es, on construit un 
#   # histogramme, on pr?dit et on agr?ge.
# 
#   xx <- data$train
#   grille <- data$test
#   n <- length(xx)
#   mat <- matrix(ncol = length(grille), nrow = B)
# 
#   for (i in 1:B) {
#     xb  = xx[sample(n, replace = TRUE)]
#     mat[i, ] <- onekdeucv(xx = xb, grille = grille)
#   }
# 
#   punctual <- colMeans(mat, na.rm = TRUE)
#   tube     <- baseboot(data$test, punctual, mat, conf = conf)
# 
#   # usando la escala raiz  => escala invariante al usar quantiles
# #  tube     <- baseboot(data$test, sqrt(punctual), sqrt(mat))
# #  tube[tube < 0] <- 0
# #  tube <- tube^2
# 
#   return(list(punctual = punctual, bunch = mat, tube = tube)) 
# }
# 
# 
# tube_KDEsm <- function(data) { # only at 90% (hard coded)
#   fit <- sm::sm.density(data$train, h = bw.ucv(data$train), display = "none",
#                         eval.grid = TRUE, eval.points = data$test)
#   return(list(punctual = fit$estimate, 
#               tube = rbind(fit$lower, fit$upper)))  
# }
# 
# 
# tube_wass <- function(data, alpha = 0.1, nbr = 50, plot = FALSE, ...) {
#   xb  <- data$train
#   n   <- length(data$train)
#   
#   hs <- hist(xb, breaks = mybreaks(xb, nbr), 
#              warn.unused = FALSE, probability = TRUE, plot = plot, ...)
# 
#   m   <- length(hs$counts)
#   c   <- qnorm(p = 1 - alpha / 2 / m) * sqrt(m / n)
# 
#   ln1    <- c(pmax(sqrt(hs$density) - c, 0)^2, 0)
#   un1    <- c((sqrt(hs$density) + c)^2       , 0)
#     
#   hs2 <- predict.hist(hs, data$test)
# 
#   ln2    <- pmax(sqrt(hs2) - c, 0)^2
#   un2    <- (sqrt(hs2) + c)^2
#   
#   if (plot) {
#     lines(data$test, ln2, type = 's', col = 4, lwd = 2)
#     lines(data$test, un2, type = 's', col = 4, lwd = 2)
#   }
#   
#   return(list(punctual = hs2,
#               tube = rbind(ln2, un2)))
# }
  


## 0. Simulate data ####
#dd <- list(train = rnorm(1000), test = sort(rnorm(1000)))
#oo <- tube_wass(dd, plot = TRUE,ylim = c(0, 2))
#lines(dd$test, oo$punctual, col = 2, lwd = 2)
#lines(dd$test, oo$tube[1, ], type = 's', col = 4, lwd = 2)
#lines(dd$test, oo$tube[2, ], type = 's', col = 4, lwd = 2)
#curve(dchisq(x, 3), add = TRUE, col = 2, lwd = 2)


# error_tube(oo$punctual, oo$tube)
# uno <- gendata(1, 200)
# bis <- tube_KDE(uno)
# error_tube(bis$punctual, bis$tube)
# 
# ter <- BagHistfp(uno, B = 100)
# 
# qter <- tube_wass(uno)
#   
# error_tube(ter$punctual, 
#            baseboot(uno$test, ter$punctual, na.omit(ter$bunch)))
# error_tube(ter$punctual, 
#            npheuristic(uno$test, ter$punctual, na.omit(ter$bunch)))
# error_tube(ter$punctual, 
#            kfwe(uno$test, ter$punctual, na.omit(ter$bunch), k = 5))
# error_tube(qter$punctual, qter$tube)


## By default all the CI are at 90% (hard coded on sm.density)
simulaciones = function(n = 100,    # data size 
                        M = 50,     # number of replicata
                        B = 150){   # number of bootstrapt samples
  estimators <- c("Hist", "FP", "Kde", "Kde_sm", "Kde_M") #, "wass")
  densities  <- c("normal", "chi2", "mix1", "bart", "triangular", "rara1", "rara2", "rara3")

  res <- foreach(i = 1:M, .combine = rbind.data.frame, .packages = "sm") %:% 
    foreach(ll = c(1, 3, 5, 8, 11, 13, 20, 21), .combine = rbind) %dopar% {

      dd   <- gendata(ll, n)
      
      estim_hist  <- tube_hist(dd, B = B) #tubo_Hist(dd, B = B)
      estim_fp    <- tube_fp(dd,  B = B)
      estim_kde   <- tube_KDE(dd, B = B)
      estim_kdesm <- tube_KDEsm(dd)
      estim_kdeM  <- cteboot(dd$test, estim_kde$punctual, estim_kde$bunch)
#      estim_wass  <- tube_wass(dd)

      #tube2 <- baseboot(dd$test, estim_fp$punctual, estim_fp$bunch)
      #tube3 <- npheuristic(dd$test, estim_fp$punctual, estim_fp$bunch)

      uno    <- c( i, ll, 1, error_tube(dd$dobs,  estim_hist$tube))      
      dos    <- c( i, ll, 2, error_tube(dd$dobs,  estim_fp$tube))  
      tres   <- c( i, ll, 3, error_tube(dd$dobs, estim_kde$tube)) 
      cuatro <- c( i, ll, 4, error_tube(dd$dobs, estim_kdesm$tube))
      cinco  <- c( i, ll, 5, error_tube(dd$dobs, estim_kdeM))
#      seis   <- c( i, ll, 6, error_tube(dd$dobs,  estim_wass$tube))
      
      rbind(uno, dos, tres, cuatro, cinco) #, seis)  
    }
  
  colnames(res) <- c("nb_iter", "density", "estimator", "coverage", "min_width", 
                     "mean_width", "max_width")
  rownames(res) <- NULL
  res$density   <- factor(res$density, labels = densities)
  res$estimator <- factor(res$estimator, labels = estimators)
  return(res)
}  


library(foreach)
library(doParallel)
registerDoParallel(6) #detectCores())

#system.time(S <- simulaciones(n = 200, M = 20, B = 10))
system.time(S <- simulaciones(n = 500, M = 100, B = 200))

# save(S, file = "~/simulationcestubos-S.Rdata")

Sres <- aggregate(S[, c(1, 4:7)],
                 by = list(density = S$density, estimator = S$estimator), mean)


xtable::xtable(Sres[order(Sres$density), ], digits = 3)
plot(Sres$coverage, Sres$mean_width, pch = 19, col = as.integer(Sres$estimator))

tab1 <- 100 * xtabs(coverage ~ density + estimator, data = Sres)[, -5]
tab2 <- xtabs(mean_width ~ density + estimator, data = Sres)[, -5]

xtable::xtable(round(cbind(tab1, tab2), 2))

      
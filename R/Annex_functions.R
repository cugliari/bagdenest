
#'@import ks

ind=function(x,a,b)  {
  # returns an indicator value for w being in ]a,b]
  ifelse(x > a & x <= b, 1,0)
}

mel = function(x,prop=0.5,m1=-1,s1=1,m2=1,s2=1){
  # computes a density for a mixture of two gaussians
  z = prop*dnorm(x,m1,s1) + (1-prop)*dnorm(x,m2,s2)
  z
}

melange = function(n,prop=0.5,m1=-1,s1=1,m2=1,s2=1)	{
  x1=rnorm(n*prop,m1,s1)
  x2=rnorm(n-n*prop,m2,s2)
  c(x1,x2)
}



#####################
#Triangle
#####################

#rtriangle
rtriangle=function (n = 1, a = 0, b = 1, c = 0.5)
{
  if (length(n) > 1)
    n <- length(n)
  if (n < 1 | is.na(n))
    stop(paste("invalid argument: n =", n))
  n <- floor(n)
  if (any(is.na(c(a, b, c))))
    return(rep(NaN, times = n))
  if (a > c | b < c)
    return(rep(NaN, times = n))
  if (any(is.infinite(c(a, b, c))))
    return(rep(NaN, times = n))
  p <- runif(n)
  if (a != c) {
    i <- which((a + sqrt(p * (b - a) * (c - a))) <= c)
    j <- which((b - sqrt((1 - p) * (b - a) * (b - c))) >
                 c)
  }
  else {
    i <- which((a + sqrt(p * (b - a) * (c - a))) < c)
    j <- which((b - sqrt((1 - p) * (b - a) * (b - c))) >=
                 c)
  }
  if (length(i) != 0)
    p[i] <- a + sqrt(p[i] * (b - a) * (c - a))
  if (length(j) != 0)
    p[j] <- b - sqrt((1 - p[j]) * (b - a) * (b - c))
  return(p)
}

#dtriangle

dtriangle=function (q, a = 0, b = 1, c = 0.5)
{
  q1 <- length(q)
  a1 <- length(a)
  b1 <- length(b)
  c1 <- length(c)
  dTest <- function(X) {
    if (any(is.na(X))) {
      if (any(is.nan(X)))
        return(NaN)
      else return(NA)
    }
    else if (X[2] > X[4] | X[3] < X[4] | (X[1] == X[2] &
                                          X[2] == X[4])) {
      warning("values required to be  a <= c <= b (at least one strict inequality)")
      return(NaN)
    }
    else if (any(is.infinite(X[2:4]))) {
      return(NaN)
    }
    else if (X[1] <= X[2]) {
      return(0)
    }
    else if (X[2] != X[4] & X[1] < X[4]) {
      return(2 * (X[1] - X[2])/(X[3] - X[2])/(X[4] - X[2]))
    }
    else if (X[4] != X[3] & X[1] >= X[4] & X[1] < X[3]) {
      return(2 * (X[3] - X[1])/(X[3] - X[2])/(X[3] - X[4]))
    }
    else if (X[1] >= X[3]) {
      return(0)
    }
  }
  k <- max(q1, a1, b1, c1)
  if (k == 1)
    return(dTest(c(q, a, b, c)))
  params <- matrix(nrow = k, ncol = 4)
  tryCatch({
    params[, 1] <- q
    params[, 2] <- a
    params[, 3] <- b
    params[, 4] <- c
  }, error = function(X) {
    stop(paste(" -- Argument Lengths: length of q = ", q1,
               ", a = ", a1, ", b = ", b1, ", c = ", c1, " -- ",
               X, sep = ""))
  })
  return(apply(params, 1, dTest))
}


###################
#Erreurs
###################
error = function(obs,prev)	{
  #1) Erreur quadratique (c'est le ISE).
  er1 = mean((obs-prev)^2,na.rm=TRUE)
  #2) Erreur log vraissemblance
  er2 = sum(ifelse(prev>0,log(prev),0))
  #2bis) Distance KL
  #	er22 = mean(ifelse(obs>0,log(obs),0)) - mean(ifelse(prev>0,log(prev),0))
  #3) Erreur de dissimilarit? d'Hellinger
  #	er3 =mean((sqrt(obs)-sqrt(prev))^2)
  signif(c(er1,er2),3)
}


#######################
#Fonctions d'histogrammes
#######################


# ##Grille pour Histogramme greedy
#
# grid=function(x,leaf=10){
#
#   h=delt::eval.greedy(as.matrix(x),leaf)
#   sp=as.vector(h$split)
#   sp=sp[which (sp>0)]
#   sort(sp)
#
#   nx=NULL
#
#   for(i in 1:length(x)-1){
#     nx[i]=(x[i]+x[i+1])/2
#   }
#   g=c(min(x),nx[sp],max(x))
#   g
# }


#mybreaks = function(x = rnorm(100),nbr)
#{
#  mx = min(x)
#  Mx = max(x)
#  seq(mx,Mx, (Mx-mx)/nbr)
#}


predict_hist_x = function(hh,x1)	{
  breaks = hh$breaks
  intens = hh$density
  lims = breaks[length(breaks)]
  if(x1 < breaks[1] || x1 > lims) res= 0
  else {
    pos1= which(x1 < breaks)[1] - 1
    pos2 = which(x1 <= breaks)[1] - 1
    pos = max(pos1,pos2,na.rm=T)
    res=intens[pos]
  }
  #if(x1==lims) res = intens[length(breaks)-1]
  res
}

predict_hist = function(hh,x)	{
  res=NULL
  for(i in 1:length(x))
    res[i] = predict_hist_x(hh,x[i])
  res
}





###############
#Fontions de Kde
###############

onekdenrd = function(xx,grid) {
  # xx	data vector
  # grille	grid for density evaluation
  ks::kde(xx,bw.nrd(xx),eval.points=grid)$estimate
}

onekdenrd0 = function(xx,grid) {
  # xx	data vector
  # grille	grid for density evaluation
  ks::kde(xx,h=bw.nrd0(xx),eval.points=grid)$estimate
}

onekdesj = function(xx,grid) {
  # xx	data vector
  # grille	grid for density evaluation
  ks::kde(xx,h=bw.SJ(xx),eval.points=grid)$estimate
}




ind2=function(x,a,b)  {
  ifelse(x >= a & x <= b, 1,0)
}

ajustkde=function(numodel=1,n=1000,B=10,test=T)	{
  # trace les evaluations de kde construits sur echantillon bootstrap
  # l'evaluaion trac?e peut etre celle de l'ech test ou de l'echantillon d'appe
  dd=gendata(numodel,n)
  xx = dd$train
  if(test) xtest = sort(dd$test) else 	xtest = sort(xx)
  A=matrix(NA,nrow=length(xtest),ncol=B)
  for (i in 1:B)	{
    db= xx[sample(n,replace=T)]
    A[,i]= ks::kde(db,bw.nrd(db),eval.points=xtest)$estimate
  }
  matplot(xtest,A,type="l")
  points(xtest,kde(xx,bw.nrd(xx),eval.points=xtest)$estimate ,type="l",lwd=2)
}


###################
# Dimarzio BoostKde
###################

#Integration de Simpson

sintegral=function (x, fx, n.pts = 256, ret = FALSE)
{
  if (class(fx) == "function")
    fx = fx(x)
  n.x = length(x)
  if (n.x != length(fx))
    stop("Unequal input vector lengths")
  if (n.pts < 64)
    n.pts = 64
  ap = approx(x, fx, n = 2 * n.pts + 1)
  h = diff(ap$x)[1]
  integral = h * (ap$y[2 * (1:n.pts) - 1] + 4 * ap$y[2 * (1:n.pts)] +
                    ap$y[2 * (1:n.pts) + 1])/3
  value=sum(integral)
  #invisible(list(value = sum(integral), cdf = list(x = ap$x[2 *
  #    (1:n.pts)], y = cumsum(integral))))
}





boostkde = function(xlearn,xtest,K=10,h=bw.nrd(xlearn),eps=0.0001)
{
  n = length(xlearn)
  if(is.null(xtest)) xtest = xlearn
  yy = outer(xlearn,xlearn,"-")/h
  noyau = t(exp( -0.5 *(yy^2)) / (h * sqrt(2* pi)) )
  yy = outer(xtest,xlearn,"-")/h
  noyau2 = t(exp( -0.5 *(yy^2)) / (h * sqrt(2* pi)) )
  wk = rep(1/n,n)
  prevapp =  colSums(wk * noyau)
  res  = colSums(wk * noyau2)
  for(k in 2:K)	{
    B =prevapp - diag(wk * noyau)
    maj = ifelse( prevapp > 0 & B >0& wk>0,log(1-wk)+log(prevapp)-log(B),0)
    wk =  wk + maj
    wk = wk /sum(wk)
    prevapp = colSums(wk * noyau)
    prevtest= colSums(wk * noyau2)
    res = res * prevtest
    norm = sintegral(sort(xtest),res[order(xtest)])
    res=  res /(norm)
    res
  }
  norm = sintegral(sort(xtest),res[order(xtest)])
  res /norm
}

############################################################################
###SIMULATIONS
############################################################################





riskhist <- function(obs, m, xlim = c(0, 1)) {
  obs01  <- (obs - xlim[1]) / (xlim[2] - xlim[1])
  h      <- 1 / m
  n      <- length(obs)
  breaks <- seq(0, 1, length.out = m + 1)
  p_hat  <- hist(obs01, plot = FALSE, breaks = breaks, warn.unused = F)$counts / n

  res <- 2 / h / (n - 1) - (n + 1) / (n - 1) / h * sum(p_hat^2)
  return(res) #(m * sum(p_hat^2))
}

riskfp <- function(obs, m, xlim = c(0, 1)) {
  obs01  <- (obs - xlim[1]) / (xlim[2] - xlim[1])
  h      <- 1 / m
  n      <- length(obs)
  breaks <- seq(0, 1, length.out = m + 1)
  p_hat  <- hist(obs01, plot = FALSE, breaks = breaks,warn.unused = F)$counts #/ n

  vs <- cbind(c(0, 0, p_hat), c(0, - 2 * p_hat, 0), c(p_hat, 0, 0))

  res <- 271 / (480 * n * h) + 49 / (2880 * n^2 * h) * sum(rowSums(vs)^2)
  return(res) #(m * sum(p_hat^2))
}

broptfp = function(x){
  Mgrid <-  2:(5 * floor(sqrt(length(x)))) #2:200 modified to reduce computational burden
  J     <- numeric(length(Mgrid))
  for (m in seq_along(Mgrid)) {
    J[m] <- riskfp(obs = x, Mgrid[m], xlim = c(min(x) - 0.5, max(x) + 0.5))
  }
  #  plot.ts(Mgrid, J);  abline(v = Mgrid[which.min(J)], col = "red")
  list(opt = max(5, Mgrid[which.min(J)]))
}


baseboot <- function(x, punctual, bunch, conf = c(0.05, 0.95), plot = FALSE){

  err <- apply(bunch, 1, '-', punctual)
  CI  <- apply(err, 1, quantile, probs = conf)

  Linf <- punctual + CI[1, ]
  Lsup <- punctual + CI[2, ]

  if(plot) plotbunch(x, punctual, t(bunch), Linf, Lsup)

  return(rbind(Linf, Lsup))
}

cteboot <- function(x, punctual, bunch, conf = c(0.05, 0.95), plot = FALSE){

  err <- apply(bunch, 1, '-', punctual)
  Mb  <- apply(abs(err), 1, max)

  M   <- quantile(Mb, probs = conf[2] - conf[1])

  Linf <- pmax(punctual - M, 0)
  Lsup <- punctual + M

  if(plot) plotbunch(x, punctual, t(bunch), Linf, Lsup)

  return(rbind(Linf, Lsup))
}

plotbunch <- function(x, punctual, bunch, inf, sup) {
  matplot(x, bunch, type = 'l', col = "grey", lty = 1, lwd = 2, ylab = '')
  lines(x, punctual, lwd = 2)
  lines(x, inf, lwd = 2, lty = 2)
  lines(x, sup, lwd = 2, lty = 2)
}

error_tube <- function(target, tube) {
  coverage <- na.omit(target <= tube[2, ] & target >= tube[1, ]) # is target in tube?
  width    <- na.omit(tube[2, ] - tube[1, ])

  return(c(mean(coverage), min(width), mean(width), max(width)))
}

tube_hist <- function(data, nbr = 50, B = 10, conf = c(0.05, 0.95)) {
  # A chaque etape, on prend un nouveau jeu de donn?es, on construit un
  # histogramme, on pr?dit et on agr?ge.

  xx <- data$train
  grid <- data$test
  n <- length(xx)
  mat <- matrix(ncol = length(grid), nrow = B)

  for (i in 1:B) {
    xb  = xx[sample(n, replace = TRUE)]
    mybreaks <- mybreaks(xb, nbr)
    hs2 = hist(xb, breaks = mybreaks, plot = FALSE)

    mat[i, ] <- predict_hist(hh = hs2, x = grid)
  }

  punctual <- colMeans(mat, na.rm = TRUE)
  tube     <- baseboot(data$test, punctual, mat, conf = conf)

  return(list(punctual = punctual, bunch = mat, tube = tube))
}

tube_fp <- function(data, nbr = 50, B = 10, conf = c(0.05, 0.95)) {
  # A chaque etape, on prend un nouveau jeu de donn?es, on construit un
  # histogramme, on pr?dit et on agr?ge.

  xx <- data$train
  grid <- data$test
  n <- length(xx)
  mat <- matrix(ncol = length(grid), nrow = B)

  for (i in 1:B) {
    xb  = xx[sample(n, replace = TRUE)]
    mybreaks <- mybreaks(xb, nbr)
    hs2 = hist(xb, breaks = mybreaks, plot = FALSE)

    mids <- hs2$mids
    delta <- diff(mybreaks)[1]
    # This avoids NA on prediction for test data between the min on grille
    # and the min obs value on xb
    newmids <- c(min(mids) - delta, mids, max(mids) + delta)

    mat[i, ] <- approxfun(x = newmids, y = c(0, hs2$density, 0))(grid)
  }

  # NAs occur when test data falls beyond the FP's support. They may safetly
  # be replaced by 0s
  mat[is.na(mat)] <- 0

  punctual <- colMeans(mat, na.rm = TRUE)
  tube     <- baseboot(data$test, punctual, mat, conf = conf)

  return(list(punctual = punctual, bunch = mat, tube = tube))
}

tube_KDE <- function(data, nbr = 50, B = 10, conf = c(0.05, 0.95)) {
  # A chaque etape, on prend un nouveau jeu de donn?es, on construit un
  # histogramme, on pr?dit et on agr?ge.

  xx <- data$train
  grid <- data$test
  n <- length(xx)
  mat <- matrix(ncol = length(grid), nrow = B)

  for (i in 1:B) {
    xb  = xx[sample(n, replace = TRUE)]
    mat[i, ] <- onekdeucv(xx = xb, grid = grid)
  }

  punctual <- colMeans(mat, na.rm = TRUE)
  tube     <- baseboot(data$test, punctual, mat, conf = conf)

  # usando la escala raiz  => escala invariante al usar quantiles
  #  tube     <- baseboot(data$test, sqrt(punctual), sqrt(mat))
  #  tube[tube < 0] <- 0
  #  tube <- tube^2

  return(list(punctual = punctual, bunch = mat, tube = tube))
}





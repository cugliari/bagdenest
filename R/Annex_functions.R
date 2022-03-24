# functions.R
# cr?e en 01/2011

###Actualis? le 13/2/14 par Mathias

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
  er1 = mean((obs-prev)^2,na.rm=T)
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


# rash_var = function(xx, grid,
#                     nbr = 50, B = 10, alpha = 0.05) {
#   # xx  data vector
#   # grille  grid for density evaluation
#   # br	number of break sfor histogram
#   # B		number of histograms to aggregate
#   #Je perturbe chaque histogramme avec une normale et j'agr?ge (sans bootstrap)
#   zz  = hist(xx, breaks = mybreaks(xx, nbr), plot = FALSE,warn.unused = F)$breaks
#   mx  = min(xx)
#   Mx  = max(xx)
#   sig = abs(min(diff(zz)))
#
#   n   <- length(xx)
#   c   <- qnorm(p = 1 - alpha / 2 / nbr) * sqrt(nbr / n)
#
#   lapply(1:B, function(i)  {
#     newb = zz + rnorm(1, 0, sig^3) #runif(1, -.1, .1)#
#     newb=sort(newb)
#     #         if(min(newb) > mx) newb= c(mx,newb)
#     #         if(max(newb) < Mx) newb= c(newb, Mx)
#     if(min(newb) > mx) newb[1] = mx
#     if(max(newb) < Mx) newb[length(newb)] = Mx
#
#     hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
#
#     ln <- c(pmax(sqrt(hs2$density) - c, 0)^2, 0)
#     un <- c((sqrt(hs2$density) + c)^2       , 0)
#     list(hist = hs2, ln = ln, un = un)
#   })
# }



# rashgreedy = function(xx,grid,nbr = 50, B=10) {
#   fin = 0
#   zz=hist(xx,breaks=grid(xx,nbr),plot=F,warn.unused = F)$breaks
#   mx = min(xx)
#   Mx = max(xx)
#   for(i in 1:B)
#   {
#
#     #		blabla=sqrt(abs(min(diff(zz))))
#     #		newb = zz + rnorm(1,0,blabla)
#     sig=abs(min(diff(zz)))
#     newb = zz + rnorm(1,0,sig^3)
#     newb=sort(newb)
#     if(min(newb) > mx) newb= c(mx,newb)
#     if(max(newb) < Mx) newb= c(newb, Mx)
#     hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
#     fin= fin + predict_hist(hs2,grid)
#     #if(i%%20 == 0) cat(i,">>")
#   }
#   #cat("\n")
#   fin/B
# }



# rashcart = function(xx,grid,B=10) {
#   fin = 0
#   d=data.frame(xx)
#   yy=xx
#   d=cbind(d,yy)
#   toto=rpart(yy~xx,data=d)
#   zz=hist(xx,breaks=c(min(xx),toto$splits[,4],max(xx)),plot=F,warn.unused = F)$breaks
#   mx = min(xx)
#   Mx = max(xx)
#   for(i in 1:B)
#   {
#
#     #blabla=sqrt(abs(min(diff(zz))))
#     #newb = zz + rnorm(1,0,blabla)
#     sig=abs(min(diff(zz)))
#     newb = zz + rnorm(1,0,sig^3)
#     newb=sort(newb)
#     if(min(newb) > mx) newb= c(mx,newb)
#     if(max(newb) < Mx) newb= c(newb, Mx)
#     hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
#     fin= fin + predict_hist(hs2,grid)
#     #if(i%%20 == 0) cat(i,">>")
#   }
#   #cat("\n")
#   fin/B
# }



#######################
##ASH
#######################
# avshift=function(xx,aa,nbr,M){
#   mx=min(xx)
#   Mx=max(xx)
#   h=(Mx-mx)/nbr
#   s=seq(mx,Mx,h)
#   pred=0
#   for (i in 1:M){
#     hh=hist(xx,breaks=c(mx-0.5,s+(i-1)*h/M,Mx+0.5),plot=F,warn.unused = F)
#     pred=pred+predict_hist(hh,sort(aa))
#   }
#   pred=pred/M
# }


#########################
#Samarov
###Version Samarov apres suggestionBadih
##########################


# samarov = function(xx,aa,B=10,alpha=1) {
#   jmin=+Inf
#   for(nbr in c(10,20,50)){
#     zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks
#     mx = min(xx)
#     Mx = max(xx)
#     for(i in 1:B)
#     {
#       sig=abs(min(diff(zz)))
#       newb = zz + rnorm(1,0,sig^3)
#       newb=sort(newb)
#       if(min(newb) > mx) newb= c(mx,newb)
#       if(max(newb) < Mx) newb= c(newb, Mx)
#       hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
#       je=-2/(length(aa))*sum(predict_hist(hs2,aa))+sintegral(sort(aa),(hs2$intensities[order(aa)])^2)
#       if (je<jmin) {
#         jmin=je
#         bropt=newb
#       }
#       #if(i%%20 == 0) cat(i,">>")
#     }
#   }
#   #cat("\n")
#   histsam=hist(xx,breaks=bropt,plot=F,warn.unused = F)
#   predsam=predict_hist(histsam,aa)
#   list(bropt=bropt,jmin=jmin,predsam=predsam)
# }


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


################################
# les fonctions pour le stacking
################################

# genereA = function(x,H=c(.1,.2,.3,.1,.2,.3),V=10){
#   # genere la matrice de fhat(E_{-v},E_v) par validation crois?e
#   L = length(H)
#   N = length(x)
#   A = matrix(NA,nrow=N,ncol=L)
#   tailpaquets = N%/%V
#   reste = N%%V
#   indices = rep(1:V,tailpaquets)
#   indices = c(indices,sample(V,reste))
#   indices = sample(indices)
#   for(l in 1:3)
#     for(v in 1:V)		{
#       bloque = which(indices ==v)
#       kk=ks::kde(x[-bloque],h=H[l],eval.points=x[bloque])
#       A[bloque,l] = kk$estimate
#     }
#   for(l in 4:6)
#     for(v in 1:V)		{
#       bloque = which(indices ==v)
#       kk=KDEtr(x[-bloque],x[bloque],h=H[l])
#       A[bloque,l] = kk
#     }
#   A
# }
#
#
# genereB = function(x,xtest=NULL,H=c(.1,.2,.3,.1,.2,.3)){
#   # genere la matrice des densit?s estim?es sur l'?chantillon d'apprentissage
#   # ou sur un echantillon test
#   L = length(H)
#   if(is.null(xtest)) xtest=x
#   N = length(xtest)
#   A = matrix(NA,nrow=N,ncol=L)
#   for(l in 1:3)	{
#     kk=ks::kde(x,h=H[l],eval.points=xtest)
#     A[,l] = kk$estimate
#   }
#   for(l in 4:6)	{
#     kk=KDEtr(x,xtest,h=H[l])
#     A[,l] = kk
#   }
#   A
# }
#
# estimalpha = function(A,eps=0.0001,itermax=3000)
# {
#   # estime les alphas
#   # renvoie les valeurs estim?es de alpha ainsi que la courbe
#   # de la log vraisemblance du m?lange estim?
#   A = A[which(rowSums(A) > 1e-10),]
#   # les lignes de A o? il n'y a que des zeros sont retir?es...
#   L = ncol(A)
#   N = nrow(A)
#   #	print(summary(rowSums(A)))
#   #	alpha = runif(L)
#   #	alpha = alpha / sum(alpha)
#   alpha = rep(1/L,L)
#   iter = 0
#   Vrais = NULL
#   while(iter <itermax)	{
#     Malpha = matrix(alpha,nrow=N,ncol=L,byrow=T)
#     coeffs = A %*% (alpha)
#     Vrais = c(Vrais, mean( ifelse(coeffs>0 ,log(coeffs),0)  ))
#     Mcoeffs = matrix(coeffs,nrow=N,ncol=L)
#     aux = Malpha/Mcoeffs
#     M = aux * A
#     alphanew = colMeans(M)
#     ecart = sum(abs(alpha - alphanew))
#     # l'ecrat quadratique entre deux valeurs successives de alpha
#     # est tr?s faible, voir le papier sur stacking
#     iter = iter +1
#     alpha = alphanew
#     if(ecart < eps) break
#   }
#   list(alpha=alpha,LogL = Vrais)
# }
#
#
# stack_dens = function(don=rnorm(1000),xtest=NULL,pl=F,H=c(.1,.2,.3,.1,.2,.3),V=10,eps=0.0001,itermax=3000)
# {
#
#
#   # Exemple d'utilisation
#   #	stack.dens()
#   # pour un echantillon test constitu? d'une grille r?guli?re
#   # fabriquer une grille r?guli?re pour estimer la densit?
#   #	don= rnorm(1000)
#   #	pas = (max(don) - min(don))/ (length(don)-1)
#   #	xtest =  seq(min(don),max(don),pas)
#   #	stack.dens(don, xtest)
#
#   mat = genereA(don,H,V)
#   toto = estimalpha(mat,eps=eps,itermax=itermax)
#   if(pl) par(mfrow=c(1,2))
#   # pour voir l'?volution de la vraisemblance
#   if(pl) plot(toto$LogL,type="l")
#   if(is.null(xtest)) test = don
#   else test = xtest
#   mat = genereB(don,xtest=test,H)
#   prev = mat%*%toto$alpha
#   prev = prev[order(test)]
#   # pour voir les valeurs de la densit? estim?e
#   if(pl) plot(prev,type= "l")
#   list(prev=prev,alpha = toto$alpha)
# }
#
#
# #######################
# #Stacking d'histogrammes
# ######################
#
# genereAh = function(x,H=c(5,10,20,30,40,50),V=10){
#
#   # genere la matrice de fhat(E_{-v},E_v) par validation crois?e
#   L = length(H)
#   N = length(x)
#   A = matrix(NA,nrow=N,ncol=L)
#   tailpaquets = N%/%V
#   reste = N%%V
#   indices = rep(1:V,tailpaquets)
#   indices = c(indices,sample(V,reste))
#   indices = sample(indices)
#   for(l in 1:L)
#     for(v in 1:V)	{
#       bloque = which(indices ==v)
#       kk=hist(x[-bloque],breaks=mybreaks(x[-bloque],H[l]),plot=F,warn.unused = F)
#       A[bloque,l] = predict_hist(kk,x[bloque])
#     }
#   A
# }
#
#
#
#
# genereBh = function(x,xtest=NULL,H=c(5,10,20,30,40,50))
# {
#   # genere la matrice des densit?s estim?es sur l'?chantillon d'apprentissage
#   # ou sur un echantillon test
#   L = length(H)
#   if(is.null(xtest)) xtest=x
#   N = length(xtest)
#   A = matrix(NA,nrow=N,ncol=L)
#   for(l in 1:L)	{
#     kk=hist(x,breaks=mybreaks(x,H[l]),plot=F,warn.unused = F)
#     A[,l] = predict_hist(kk,xtest)
#   }
#   A
# }
#
#
# estimalpha = function(A,eps=0.0001,itermax=3000)
# {
#   # estime les alphas
#   # renvoie les valeurs estim?es de alpha ainsi que la courbe
#   # de la log vraisemblance du m?lange estim?
#   A = A[which(rowSums(A) > 1e-10),]
#   # les lignes de A o? il n'y a que des zeros sont retir?es...
#   L = ncol(A)
#   N = nrow(A)
#   #	print(summary(rowSums(A)))
#   #	alpha = runif(L)
#   #	alpha = alpha / sum(alpha)
#   alpha = rep(1/L,L)
#   iter = 0
#   Vrais = NULL
#   while(iter <itermax)	{
#     Malpha = matrix(alpha,nrow=N,ncol=L,byrow=T)
#     coeffs = A %*% (alpha)
#     Vrais = c(Vrais, mean( ifelse(coeffs>0 ,log(coeffs),0)  ))
#     Mcoeffs = matrix(coeffs,nrow=N,ncol=L)
#     aux = Malpha/Mcoeffs
#     M = aux * A
#     alphanew = colMeans(M)
#     ecart = sum(abs(alpha - alphanew))
#     # l'ecrat quadratique entre deux valeurs successives de alpha
#     # est tr?s faible, voir le papier sur stacking
#     iter = iter +1
#     alpha = alphanew
#     if(ecart < eps) break
#   }
#   list(alpha=alpha,LogL = Vrais)
# }
#
#
# stack_denshist = function(don=rnorm(1000),xtest=NULL,pl=F,H=c(5,10,20,30,40,50),V=10,eps=0.0001,itermax=3000)
# {
#   mat = genereAh(don,H,V)
#   toto = estimalpha(mat,eps=eps,itermax=itermax)
#   if(pl) par(mfrow=c(1,2))
#   # pour voir l'?volution de la vraisemblance
#   if(pl) plot(toto$LogL,type="l")
#   if(is.null(xtest)) test = don
#   else test = xtest
#   mat = genereBh(don,xtest=test,H)
#   prev = mat%*%toto$alpha
#   prev = prev[order(test)]
#   # pour voir les valeurs de la densit? estim?e
#   if(pl) plot(prev,type= "l")
#   list(prev=prev,alpha = toto$alpha)
# }
#
#
# #################
# # Tsybakov AggPure
# #################
#
# tsybakov = function(dlearn,dtest, S=10,H=c( 0.001,0.005, 0.01,0.05,0.1,0.5))
# {
#   # dd 	:	data
#   # S		:	number of splitting for teh data
#   n = length(dlearn)
#   prev = 0
#   for(s in 1:S) {
#     ind = sample(n,n/2)
#     Bapp = genereB(dlearn[ind],H=H)
#     Btest0 = genereB(dlearn[ind],dlearn[-ind],H=H)
#     Btest = genereB(dlearn[ind],dtest,H=H)
#     alphs = estimalpha(Btest0)$alpha
#     prevs = Btest %*% alphs
#     prev = prev + prevs
#   }
#   prev = prev[order(dtest)]
#   prev/S
# }

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





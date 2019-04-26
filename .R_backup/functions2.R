# functions.R
# cr?e en 01/2011

###Actualis? le 13/2/14 par Mathias


library(caTools) # pour trapz  calcul numerique d'int pour normalisation dans boosting
library(benchden)
library(ks)
library(rpart)
#library(delt)
library(nor1mix)
##################
#Mod?les
##################


gendata= function(nummodel = 1, n = 100)
{
# nummodel : number of the model
# n : sample size
# returns a list with: the training set, the test set (ordered), and the observed densities for the ordered test set
	#NORMALE STANDARD
	if(nummodel == 1)	{
		train <- rnorm(n)
		minm  <- min(train)
		maxm  <- max(train)
		test  <- sort(rnorm(n))
		test  <- test[ (test <= maxm) & (test> minm)]
		dobs  <- dnorm(test)
	}

	#EXPONENTIELLE e^{-x}
	if(nummodel==2) {
		train = rexp(n)
		minm = min(train)
		maxm = max(train)
		test=sort(rexp(n))
		test = test[ (test <= maxm) & (test> minm)]
		dobs = dexp(test)
	}

	#CHI CARRE 10 degr?s de libert?.
	if(nummodel==3) {
		train=rchisq(n,10)
		minm = min(train)
		maxm = max(train)
		test=sort(rchisq(n,10))
		test = test[ (test <= maxm) & (test> minm)]
		dobs=dchisq(test,10)
	}
	
	#t-STUDENT 4 degr?s de libert?.
	if(nummodel==4) {
		train=rt(n,4)
		minm = min(train)
		maxm = max(train)
		test=sort(rt(n,4))
		test = test[ (test <= maxm) & (test> minm)]
		dobs=dt(test,4)
	}
	
	
	#MELANGE NORMALES 0.5N(-1,0.3)+0.5N(1,0.3)
	if(nummodel==5) {
		prop=0.5
		mu1=-1
		sig1=0.3
		mu2=1
		sig2=0.3
		train=melange(n,prop,mu1,sig1,mu2,sig2)
		minm = min(train)
		maxm = max(train)
		test = sort(melange(n,prop,mu1,sig1,mu2,sig2))
		test = test[ (test <= maxm) & (test> minm)]		
		dobs=mel(test,prop,mu1,sig1,mu2,sig2)
	}
	
	if(nummodel==6) {
		train= c(rnorm(n*0.25,-3,0.5),rnorm(n*0.5,0,1),rnorm(n*0.25,3,0.5))
		minm = min(train)
		maxm = max(train)
		test= sort(c(rnorm(n*0.25,-3,0.5),rnorm(n*0.5,0,1),rnorm(n*0.25,3,0.5)))
		test = test[ (test <= maxm) & (test> minm)]
		dobs=0.25*dnorm(test,-3,0.5) + 0.5*dnorm(test,0,1) + 0.25*dnorm(test,3,0.5)
	}
	
	
	#MELANGE NORMALES 0.55N(-3,0.5)+0.35N(0,1)+0.1*N(3,0.5)
	if(nummodel==7) {
		train= c(rnorm(n*0.55,-3,0.5),rnorm(n*0.35,0,1),rnorm(n*0.10,3,0.5))
		minm = min(train)
		maxm = max(train)
		test= sort(c(rnorm(n*0.55,-3,0.5),rnorm(n*0.35,0,1),rnorm(n*0.10,3,0.5)))
		test = test[ (test <= maxm) & (test> minm)]
		dobs=0.55*dnorm(test,-3,0.5) + 0.35*dnorm(test,0,1) + 0.1*dnorm(test,3,0.5)
	}
	
	
	#CLAW DENSITY . IL FAUT TELECHARGER LE PAQUET BENCHDEN, c'est le numero 23.
	if(nummodel==8) {
		train = rberdev(n,23)
		minm = min(train)
		maxm = max(train)
		test = sort(rberdev(n,23))
		test = test[ (test <= maxm) & (test> minm)]
		dobs = dberdev(test,23)
	}

	#SMOOTH COMB DENSITY . IL FAUT TELECHARGER LE PAQUET BENCHDEN, c'est le numero 24.
	if(nummodel==9) {
		train = rberdev(n,24)
		minm = min(train)
		maxm = max(train)
		test = sort(rberdev(n,24))
		test = test[ (test <= maxm) & (test> minm)]
		dobs = dberdev(test,24)
	}

	# DENS1 rigollet. C'est un melange 0.5 N(0,1) + 0.5 \sum \limits_{i=1}^{10} \bold{1}_{\left(\frac{2(i-1)}{T},\frac{2i-1}{T}  \right]}
	#if(nummodel==10) {
	#	mm = rnorm(n*0.5,0,1)
	#	minm = min(mm)
	#	maxm = max(mm)
	#	borneinf = seq(0,1.8,0.2)
	#	bornesup = seq(0.1,1.9,0.2)
	#	for(j in  1:length(borneinf)) mm = c(mm,runif(n*0.05, borneinf[j],bornesup[j]))
	#	test = rnorm(n*0.5,0,1)
	#	for(j in  1:length(borneinf)) test = c(test,runif(n*0.05, borneinf[j],bornesup[j]))
	#	test = sort(test)
	#	test = test[ (test <= maxm) & (test> minm)]
	#	tot = 0
	#	for(j in  1:length(borneinf)) tot = tot +  ind(test,borneinf[j],bornesup[j])
	#	dobs= 0.5* (dnorm(test) + tot)

	# DENS2 rigollet. C'est un melange 0.5 N(0,1) + 0.5 \sum \limits_{i=1}^{14} \bold{1}_{\left(\frac{2(i-1)}{T},\frac{2i-1}{T}  \right]}
	if(nummodel==10) {
		train = rnorm(n*0.5,0,1)
		minm = min(train)
		maxm = max(train)
		borneinf = seq(0,26/14,1/7)
		bornesup = seq(1/14,27/14,1/7)
		for(j in  1:length(borneinf)) train = c(train,runif(n/28, borneinf[j],bornesup[j]))	
		test = rnorm(n*0.5,0,1)
		for(j in  1:length(borneinf)) test = c(test,runif(n/28, borneinf[j],bornesup[j]))	
		test = sort(test)
		test = test[ (test <= maxm) & (test> minm)]
		tot = 0
		for(j in  1:length(borneinf)) tot = tot + ind(test,borneinf[j],bornesup[j]) 		
		dobs = 0.5 * (dnorm(test) + tot)

#		
#		dobs = 0.5*dnorm(sort(test),0,1)+(1/28)*dunif(sort(test),0,1/14)+(1/28)*dunif(sort(test),2/14,3/14)+(1/28)*dunif(sort(test),4/14,5/14)+(1/28)*dunif(sort(test),6/14,7/14)+(1/28)*dunif(sort(test),8/14,9/14)
#		+(1/28)*dunif(sort(test),10/14,11/14)+(1/28)*dunif(sort(test),12/14,13/14)+(1/28)*dunif(sort(test),1,15/14)+(1/28)*dunif(sort(test),16/14,17/14)+(1/28)*dunif(sort(test),18/14,19/14)+(1/28)*dunif(sort(test),20/14,21/14)
#		+(1/28)*dunif(sort(test),22/14,23/14)+(1/28)*dunif(sort(test),24/14,25/14)+(1/28)*dunif(sort(test),26/14,27/14)
		}

	#Triangulaire [0,2] max en 1.
	if(nummodel==11) {
		train = rtriangle(n,0,2,1)
		minm = min(train)
		maxm = max(train)
		test=sort(rtriangle(n,0,2,1))
		test = test[ (test <= maxm) & (test> minm)]
		dobs = dtriangle(test,0,2,1)
	}


#Beta 2,5.
	if(nummodel==12) {
		train = rbeta(n,2,5)
		minm = min(train)
		maxm = max(train)
		test=sort(rbeta(n,2,5))
		test = test[ (test <= maxm) & (test> minm)]
		dobs = dbeta(test,2,5)
	}

#uniforme 0-1  
  if(nummodel==13) {
    train = runif(n,0,1)
    minm = min(train)
    maxm = max(train)
    test=sort(runif(n,0,1))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = dunif(test,0,1)
  }

  #Kurtotic
  if(nummodel==14) {
    train = rnorMix(n, MW.nm4)
    minm = min(train)
    maxm = max(train)
    test = sort(rnorMix(n, MW.nm4))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = dnorMix(test, MW.nm4)
  }  
  #AsimClaw
  if(nummodel==15) {
    train = rnorMix(n, MW.nm12)
    minm = min(train)
    maxm = max(train)
    test = sort(rnorMix(n, MW.nm12))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = dnorMix(test, MW.nm12)
  }  
  #AsimDoubleClaw
  if(nummodel==16) {
    train = rnorMix(n, MW.nm13)
    minm = min(train)
    maxm = max(train)
    test = sort(rnorMix(n, MW.nm13))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = dnorMix(test, MW.nm13)
  }  
  
  #SmooComb
  if(nummodel==17) {
    train = rnorMix(n, MW.nm14)
    minm = min(train)
    maxm = max(train)
    test = sort(rnorMix(n, MW.nm14))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = dnorMix(test, MW.nm14)
  }  
  
  #DistComb
  if(nummodel==18) {
    train = rnorMix(n, MW.nm15)
    minm = min(train)
    maxm = max(train)
    test = sort(rnorMix(n, MW.nm15))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = dnorMix(test, MW.nm15)
  }  
  
  #DistBim
  if(nummodel==19) {
    train = rnorMix(n, MW.nm16)
    minm = min(train)
    maxm = max(train)
    test = sort(rnorMix(n, MW.nm16))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = dnorMix(test, MW.nm16)
  }  
  
  
 #  DENS1 rigollet. C'est un melange 0.5 N(0,1) + 0.5 \sum \limits_{i=1}^{10} \bold{1}_{\left(\frac{2(i-1)}{T},\frac{2i-1}{T}  \right]}
  if(nummodel==20) {
  	train = rnorm(n*0.5,0,1)
  	minm = min(train)
  	maxm = max(train)
  	borneinf = seq(0,1.8,0.2)
  	bornesup = seq(0.1,1.9,0.2)
  	for(j in  1:length(borneinf)) train = c(train,runif(n*0.05, borneinf[j],bornesup[j]))
  	test = rnorm(n*0.5,0,1)
  	for(j in  1:length(borneinf)) test = c(test,runif(n*0.05, borneinf[j],bornesup[j]))
  	test = sort(test)
  	test = test[ (test <= maxm) & (test> minm)]
  	tot = 0
  	for(j in  1:length(borneinf)) tot = tot +  ind(test,borneinf[j],bornesup[j])
  	dobs= 0.5* (dnorm(test) + tot)
  }
  	#Mezcla uniformes
  	if(nummodel==21) {
  	  train= c(runif(n*0.5,-2,-1),runif(n*0.5,1,2))
  	  minm = min(train)
  	  maxm = max(train)
  	  test= sort(c(runif(n*0.5,-2,-1),runif(n*0.5,1,2)))
  	  test = test[ (test <= maxm) & (test> minm)]
  	  dobs=0.5*dunif(test,-2,-1) + 0.5*dunif(test,1,2)
  	}  
  	
  
	list(train=train,test=test,dobs=dobs)
}

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


error2 = function(obs,prev)	{
#1) Erreur quadratique
	er1 = mean(obs * (obs-prev)^2)
#2) Erreur log vraissemblance
	er2 = sum(ifelse(prev>0,obs * log(prev),0))
#3) Erreur de dissimilarit? d'Hellinger
	er3 =mean(obs * (sqrt(obs)-sqrt(prev))^2)
	signif(c(er1,er2,er3),3)
}


#######################
#Fonctions d'histogrammes
#######################


##Grille pour Histogramme greedy

grille=function(x,leaf=10){

h=eval.greedy(as.matrix(x),leaf)
sp=as.vector(h$split)
sp=sp[which (sp>0)]
sort(sp)

nx=NULL

for(i in 1:length(x)-1){
nx[i]=(x[i]+x[i+1])/2
}
g=c(min(x),nx[sp],max(x))
g
}




#########



mybreaks = function(x = rnorm(100),nbr)
{
   mx = min(x)
   Mx = max(x)
   seq(mx,Mx, (Mx-mx)/nbr)
}


predict.hist.x = function(hh,x1)	{
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

predict.hist = function(hh,x)	{
	res=NULL
	for(i in 1:length(x))
		res[i] = predict.hist.x(hh,x[i])
	res
}


AgregHist = function(xx,grille=aa,nbr = 50, B=10, alpha=1) {
# xx	data vector
# grille	grid for density evaluation
# br	number of break sfor histogram
# B		number of histograms to aggregate
#Je perturbe chaque histogramme avec une normale et j'agr?ge (sans bootstrap)
	fin = 0
	zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks	
	mx = min(xx)
	Mx = max(xx)
	for(i in 1:B)
	{
		blabla=sqrt(abs(min(diff(zz))))
		newb = zz + rnorm(length(zz),0,alpha*blabla)
		newb=sort(newb)
		if(min(newb) > mx) newb= c(mx,newb)
		if(max(newb) < Mx) newb= c(newb, Mx)
		hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
		fin= fin + predict.hist(hs2,grille)
		#if(i%%20 == 0) cat(i,">>")
	}
	#cat("\n")
	fin/B
}

AgregHistunif = function(xx,grille=aa,nbr = 50, B=10) {
# xx	data vector
# grille	grid for density evaluation
# br	number of break sfor histogram
# B		number of histograms to aggregate
#Je perturbe chaque histogramme avec une normale et j'agr?ge (sans bootstrap)
	fin = 0
	zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks	
	z=diff(zz)
	h=z[1]
	mx = min(xx)
	Mx = max(xx)
	for(i in 1:B)
	{
		#blabla=sqrt(abs(min(diff(zz))))
		newb = zz + runif(length(zz),0,h)
		#newb=sort(newb)
		#if(min(newb) > mx) newb= c(mx,newb)
		#if(max(newb) < Mx) newb= c(newb, Mx)
		hs2=hist(xx,breaks=c(mx,newb,Mx),plot=F,warn.unused = F)
		fin= fin + predict.hist(hs2,grille)
		#if(i%%20 == 0) cat(i,">>")
	}
	#cat("\n")
	fin/B
}




rash = function(xx,grille=aa,nbr = 50, B=10) {
  # xx	data vector
  # grille	grid for density evaluation
  # br	number of break sfor histogram
  # B		number of histograms to aggregate
  #Je perturbe chaque histogramme avec une normale et j'agr?ge (sans bootstrap)
  fin = 0
  zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks	
  mx = min(xx)
  Mx = max(xx)
  for(i in 1:B)
  {
    blabla=sqrt(abs(min(diff(zz))))
    newb = zz + rnorm(1,0,blabla)
    newb=sort(newb)
    if(min(newb) > mx) newb= c(mx,newb)
    if(max(newb) < Mx) newb= c(newb, Mx)
    hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
    fin= fin + predict.hist(hs2,grille)
    #if(i%%20 == 0) cat(i,">>")
  }
#  cat("\n")
  fin/B
}

rash.var = function(xx, grille, 
                    nbr = 50, B = 10, alpha = 0.05) {
  # xx  data vector
  # grille  grid for density evaluation
  # br	number of break sfor histogram
  # B		number of histograms to aggregate
  #Je perturbe chaque histogramme avec une normale et j'agr?ge (sans bootstrap)
  zz  = hist(xx, breaks = mybreaks(xx, nbr), plot = FALSE,warn.unused = F)$breaks	
  mx  = min(xx)
  Mx  = max(xx)
  sig = abs(min(diff(zz)))

  n   <- length(xx)
  c   <- qnorm(p = 1 - alpha / 2 / nbr) * sqrt(nbr / n)
  
  lapply(1:B, function(i)  {
        newb = zz + rnorm(1, 0, sig^3) #runif(1, -.1, .1)#
        newb=sort(newb)
#         if(min(newb) > mx) newb= c(mx,newb)
#         if(max(newb) < Mx) newb= c(newb, Mx)
        if(min(newb) > mx) newb[1] = mx
        if(max(newb) < Mx) newb[length(newb)] = Mx
        
        hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
        
        ln <- c(pmax(sqrt(hs2$density) - c, 0)^2, 0)
        un <- c((sqrt(hs2$density) + c)^2       , 0)
        list(hist = hs2, ln = ln, un = un)
      }) 
}

# rash.kde = function(xx,grille=aa,nbr = 50, B=10) {
#   # xx  data vector
#   # grille	grid for density evaluation
#   # br	number of break sfor histogram
#   # B		number of histograms to aggregate
#   #Je perturbe chaque histogramme avec une normale et j'agr?ge (sans bootstrap)
#   fin = 0
#   zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F)$breaks	
#   mx = min(xx)
#   Mx = max(xx)
# 
#   lapply(1:B, function(i)  {
#         sig=abs(min(diff(zz)))
#         newb = zz + runif(1, -.1, .1)#rnorm(1,0,sig^3)
#         newb=sort(newb)
#         if(min(newb) > mx) newb= c(mx,newb)
#         if(max(newb) < Mx) newb= c(newb, Mx)
#         hs2=hist(xx,breaks=newb,plot=F)
#         density(unlist(mapply(rep, hs2$mids, hs2$counts)), bw = "ucv")
#       }) 
# }

rashgreedy = function(xx,grille=aa,nbr = 50, B=10) {
	fin = 0
	zz=hist(xx,breaks=grille(xx,nbr),plot=F,warn.unused = F)$breaks
	mx = min(xx)
	Mx = max(xx)
	for(i in 1:B)
	{

#		blabla=sqrt(abs(min(diff(zz))))
#		newb = zz + rnorm(1,0,blabla)
		sig=abs(min(diff(zz)))
		newb = zz + rnorm(1,0,sig^3)	
		newb=sort(newb)
		if(min(newb) > mx) newb= c(mx,newb)
		if(max(newb) < Mx) newb= c(newb, Mx)
		hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
		fin= fin + predict.hist(hs2,grille)
		#if(i%%20 == 0) cat(i,">>")
	}
	#cat("\n")
	fin/B
}



rashcart = function(xx,grille=aa,B=10) {
	fin = 0
	d=data.frame(xx)
	yy=xx
	d=cbind(d,yy)
	toto=rpart(yy~xx,data=d)
	zz=hist(xx,breaks=c(min(xx),toto$splits[,4],max(xx)),plot=F,warn.unused = F)$breaks
	mx = min(xx)
	Mx = max(xx)
	for(i in 1:B)
	{

		#blabla=sqrt(abs(min(diff(zz))))
		#newb = zz + rnorm(1,0,blabla)
		sig=abs(min(diff(zz)))
		newb = zz + rnorm(1,0,sig^3)	
		newb=sort(newb)
		if(min(newb) > mx) newb= c(mx,newb)
		if(max(newb) < Mx) newb= c(newb, Mx)
		hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
		fin= fin + predict.hist(hs2,grille)
		#if(i%%20 == 0) cat(i,">>")
	}
	#cat("\n")
	fin/B
}




AgregHist.err = function(xx,grille=aa,nbr = 50, B= 10,dobs,alpha=1) {
# xx	data vector
# grille	grid for density evaluation
# br	number of break sfor histogram
# B		number of histograms to aggregate
#Je perturbe chaque histogramme avec une normale et j'agr?ge (sans bootstrap)
	fin = 0
	err00=NULL
	zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks	
	mx = min(xx)
	Mx = max(xx)
	for(i in 1:B)
	{
		blabla=sqrt(abs(min(diff(zz))))
		newb = zz + rnorm(length(zz),0,alpha * blabla)
		newb=sort(newb)
		if(min(newb) > mx) newb= c(mx,newb)
		if(max(newb) < Mx) newb= c(newb, Mx)
		hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
		fin= fin + predict.hist(hs2,grille)
		previ=fin/i
		err00=rbind(err00,error(dobs,previ))
		#if(i%%20 == 0) cat(i,">>")
	}
	#cat("\n")
	list(prev=fin/B,erreur=err00)
}

AgregHistunif.err = function(xx,grille=aa,nbr = 50, B= 10,dobs) {
	fin = 0
	err00=NULL
	zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks	
	z=diff(zz)
	h=z[1]
	mx = min(xx)
	Mx = max(xx)
	for(i in 1:B)
	{
		newb = zz + runif(length(zz),0,h)
		hs2=hist(xx,breaks=c(mx,newb,Mx),plot=F,warn.unused = F)
		fin= fin + predict.hist(hs2,grille)
		previ=fin/i
		err00=rbind(err00,error(dobs,previ))
		#if(i%%20 == 0) cat(i,">>")
	}
	#cat("\n")
	list(prev=fin/B,erreur=err00)
}

rash.err = function(xx,grille=aa,nbr = 50, B= 10,dobs) {
	fin = 0
	err00=NULL
	zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks	
	mx = min(xx)
	Mx = max(xx)
	for(i in 1:B)
	{
		blabla=sqrt(abs(min(diff(zz))))
		newb = zz + rnorm(1,0,blabla)
		newb=sort(newb)
		if(min(newb) > mx) newb= c(mx,newb)
		if(max(newb) < Mx) newb= c(newb, Mx)
		hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
		fin= fin + predict.hist(hs2,grille)
		previ=fin/i
		err00=rbind(err00,error(dobs,previ))
#		if(i%%20 == 0) cat(i,">>")
	}
#	cat("\n")
	list(prev=fin/B,erreur=err00)
}








BagHist = function(xx,grille=aa, B= 10) {
# A chaque etape, on prend un nouveau jeu de donn?es, on construit un histogramme, on pr?dit et on agr?ge.
#
   n = length(xx)
   fin = 0
   for(i in 1:B)    {
       xb = xx[sample(n,replace=TRUE)]
       nbr=bropt(xb)$opt
       hs2=hist(xb,breaks=mybreaks(xb,nbr),plot=F,warn.unused = F)
       fin= fin + predict.hist(hs2,grille)
       #if(i%%20 == 0) cat(i,">>")
   }
   #cat("\n")
   fin/B
}

BagHist.err = function(xx,grille=aa,B= 10,dobs) {
# A chaque etape, on prend un nouveau jeu de donn?es, on construit un histogramme, on pr?dit et on agr?ge.
#
	n = length(xx)
	err00=NULL
	fin = 0
	mx = min(xx)
	Mx = max(xx)
	for(i in 1:B)    {
		xb = xx[sample(n,replace=T)]
		nbr=bropt(xb)$opt
		hs2=hist(xb,breaks=mybreaks(xb,nbr),plot=F,warn.unused = F)
		fin= fin + predict.hist(hs2,grille)
		previ=fin/i
		err00=rbind(err00,error(dobs,previ))
		#if(i%%20 == 0) cat(i,">>")
	}
#	cat("\n")
	list(prev=fin/B,erreur=err00)
}



# BagHistfp = function(xx,grille=aa,nbr = 50, B= 10) {
#   # A chaque etape, on prend un nouveau jeu de donn?es, on construit un histogramme, on pr?dit et on agr?ge.
#   #
#   n = length(xx)
#   fin = 0
#   fin2=0
#   for(i in 1:B)    {
#     xb = xx[sample(n,replace=TRUE)]
#     hs2=hist(xb,breaks=mybreaks(xb,nbr),plot=F)
#     fin= fin + approxfun(x=hs2$mids,y=hs2$density)(grille)
#     fin2=fin2+ predict.hist(hs2,grille)
#   }
#   cat("\n")
#   list(bhfp=fin/B,bh=fin2/B)
# }

BagHistfp = function(xx,grille=aa, B= 10) {
  # A chaque etape, on prend un nouveau jeu de donn?es, on construit un histogramme, on pr?dit et on agr?ge.
  #
  n = length(xx)
  fin = 0
  fin2=0
  for(i in 1:B)    {
    xb = xx[sample(n,replace=TRUE)]
    nbr=bropt(xb)$opt
    nbrfp=broptfp(xb)$opt
    hs=hist(xb,breaks=mybreaks(xb,nbr),plot=F,warn.unused = F)
    hs2=hist(xb,breaks=mybreaks(xb,nbrfp),plot=F,warn.unused = F)
    m <- hs2$mids
    h <- m[2] - m[1]
    m <- c(m[1] - h, m, m[length(m)] + h)
    d <- c(0, hs2$density, 0)
    fin  = fin  + approxfun(x = m, y = d, yright = 0, yleft = 0)(grille)
    fin2 = fin2 + predict.hist(hs,grille)
  }
 # cat("\n")
  list(bhfp=fin/B,bh=fin2/B)
}


BagHistfp.err = function(xx,grille=aa, B= 10,dobs) {
  # A chaque etape, on prend un nouveau jeu de donn?es, on construit un histogramme, on pr?dit et on agr?ge.
  #
  n = length(xx)
  err00=NULL
  err002=NULL
  fin = 0
  fin2=0
  mx = min(xx)
  Mx = max(xx)
  for(i in 1:B)    {
    xb = xx[sample(n,replace=T)]
    nbr=bropt(xb)$opt
    nbrfp=broptfp(xb)$opt
    hs=hist(xb,breaks=mybreaks(xb,nbr),plot=F,warn.unused = F)
    hs2=hist(xb,breaks=mybreaks(xb,nbrfp),plot=F,warn.unused = F)
    m <- hs2$mids
    h <- m[2] - m[1]
    m <- c(m[1] - h, m, m[length(m)] + h)
    d <- c(0, hs2$density, 0)
    fin  = fin  + approxfun(x = m, y = d, yright = 0, yleft = 0)(grille)
    
    fin2=fin2+ predict.hist(hs,grille)
    previ=fin/i
    previ2=fin2/i
    err00=rbind(err00,error(dobs,previ))
    err002=rbind(err002,error(dobs,previ2))
    #if(i%%20 == 0) cat(i,">>")
  }
#  cat("\n")
  list(prev=fin/B,erreurbhfp=err00,erreurbh=err002)
}


#######################
##ASH
#######################
avshift=function(xx,aa,nbr,M){
	mx=min(xx)
	Mx=max(xx)
	h=(Mx-mx)/nbr
	s=seq(mx,Mx,h)
	pred=0
	for (i in 1:M){
		hh=hist(xx,breaks=c(mx-0.5,s+(i-1)*h/M,Mx+0.5),plot=F,warn.unused = F)
		pred=pred+predict.hist(hh,sort(aa))
	}
	pred=pred/M
}


ash.err = function(xx,aa,nbr, B= 10,dobs) {
	mx=min(xx)
	Mx=max(xx)
	h=(Mx-mx)/nbr
	s=seq(mx,Mx,h)
	pred=0
	err0=NULL
	for (i in 1:B){
		hh=hist(xx,breaks=c(mx-0.5,s+(i-1)*h/B,Mx+0.5),plot=F,warn.unused = F)
		pred=pred+predict.hist(hh,sort(aa))
		predi=pred/i
		err0=rbind(err0,error(dobs,predi))
	}
	list(prev=predi,erreur=err0)
}


#########################
#Samarov
###Version Samarov apres suggestionBadih
##########################


samarov = function(xx,aa,B=10,alpha=1) {
    jmin=+Inf
	for(nbr in c(10,20,50)){
		zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks	
		mx = min(xx)
		Mx = max(xx)
		for(i in 1:B)
		{
			sig=abs(min(diff(zz)))
			newb = zz + rnorm(1,0,sig^3)	
			newb=sort(newb)
			if(min(newb) > mx) newb= c(mx,newb)
			if(max(newb) < Mx) newb= c(newb, Mx)
			hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
			je=-2/(length(aa))*sum(predict.hist(hs2,aa))+sintegral(sort(aa),(hs2$intensities[order(aa)])^2)
			if (je<jmin) {
				jmin=je 
				bropt=newb
			}
			#if(i%%20 == 0) cat(i,">>")
		}
	}
	#cat("\n")
	histsam=hist(xx,breaks=bropt,plot=F,warn.unused = F)
	predsam=predict.hist(histsam,aa)
    list(bropt=bropt,jmin=jmin,predsam=predsam)
}
	

###############
#Fontions de Kde
###############

onekdenrd = function(xx,grille=aa) {
# xx	data vector
# grille	grid for density evaluation
	kde(xx,bw.nrd(xx),eval.points=grille)$estimate
}

onekdenrd0 = function(xx,grille=aa) {
# xx	data vector
# grille	grid for density evaluation
	kde(xx,h=bw.nrd0(xx),eval.points=grille)$estimate
}

onekdeucv = function(xx,grille=aa) {
# xx	data vector
# grille	grid for density evaluation
	kde(xx,h=bw.ucv(xx),eval.points=grille)$estimate
}

onekdesj = function(xx,grille=aa) {
# xx	data vector
# grille	grid for density evaluation
	kde(xx,h=bw.SJ(xx),eval.points=grille)$estimate
}


KDE = function(dlearn,dtest,w=NULL,h=0.1)
{
	n = length(dlearn)
	if(is.null(w)) w=rep(1/n,n)
	yy = outer(dtest,dlearn,"-")/h
	aux = t(exp( -0.5 *(yy^2)) / sqrt(2* pi)) 
	colSums(w * aux) / h
}

ind2=function(x,a,b)  {
		ifelse(x >= a & x <= b, 1,0)
}

KDEtr = function(dlearn,dtest,w=NULL,h=0.1)
{
	n = length(dlearn)
	if(is.null(w)) w=rep(1/n,n)
	yy = outer(dtest,dlearn,"-")/h
	noyau = t((1- abs(yy))*ind2(yy,-1,1) )
	colSums(w * noyau) / h
}


Bagkde <- function(xx, grille = aa, dobs, B = 10) {
# xx	data vector
# grille	grid for density evaluation
# B		number of histograms to aggregate
	n   = length(xx)
	kk=0
  for (i in 1:B){
    		xb = xx[sample(n, replace = TRUE)]
		    kk0 = kde(xb, h=bw.ucv(xb), eval.points = grille)$estimate
		    kk=kk+kk0
        }
  kk/B
}

# Bagkde <- function(xx, grille = aa, B = 10, plerr = F) {
#   # xx  data vector
#   # grille	grid for density evaluation
#   # B		number of histograms to aggregate
#   # plerr	to plot error evolution
#   fin = 0
#   n   = length(xx)
#   logerr = numeric(B)
#   for(i in 1:B) {
#     xb = xx[sample(n,replace=T)]+rnorm(1,0,1)
#     kk = kde(xb,bw.SJ(xb),eval.points=grille)   
#     fin = fin + kk$estimate    
#     logerr[i] = mean(log(fin/i))
#     if(i%%20 == 0) cat(i,">>")
#   }
#   cat("\n")
#   if(plerr) plot(logerr,type="l")
#   fin/B
# }

BagKDE.err = function(xx,grille=aa,B= 10,dobs) {
# xx	data vector
# grille	grid for density evaluation
# B		number of histograms to aggregate
# plerr	to plot error evolution
	fin = 0
	err00=NULL
	n = length(xx)
	for(i in 1:B) {
		xb = xx[sample(n,replace=T)]
		kk = kde(xb,h=bw.ucv(xb),eval.points=grille)   
		fin = fin + kk$estimate    
		previ=fin/i
		err00=rbind(err00,error(dobs,previ))
#		if(i%%20 == 0) cat(i,">>")
	}
#	cat("\n")
	list(prev=fin/B,erreur=err00)
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
		A[,i]= kde(db,bw.nrd(db),eval.points=xtest)$estimate   
	}
	matplot(xtest,A,type="l")
	points(xtest,kde(xx,bw.nrd(xx),eval.points=xtest)$estimate ,type="l",lwd=2)
}


################################
# les fonctions pour le stacking
################################

genereA = function(x,H=c(.1,.2,.3,.1,.2,.3),V=10)
{
# genere la matrice de fhat(E_{-v},E_v) par validation crois?e
	L = length(H)
	N = length(x)
	A = matrix(NA,nrow=N,ncol=L)
	tailpaquets = N%/%V
	reste = N%%V
	indices = rep(1:V,tailpaquets)
	indices = c(indices,sample(V,reste))
	indices = sample(indices)
	for(l in 1:3)
		for(v in 1:V)		{
			bloque = which(indices ==v)
			kk=kde(x[-bloque],h=H[l],eval.points=x[bloque])   
			A[bloque,l] = kk$estimate
		}
	for(l in 4:6)
		for(v in 1:V)		{
			bloque = which(indices ==v)
			kk=KDEtr(x[-bloque],x[bloque],h=H[l])   
			A[bloque,l] = kk
		}
	A
}

genereB = function(x,xtest=NULL,H=c(.1,.2,.3,.1,.2,.3))
{
# genere la matrice des densit?s estim?es sur l'?chantillon d'apprentissage
# ou sur un echantillon test
	L = length(H)
	if(is.null(xtest)) xtest=x
	N = length(xtest)
	A = matrix(NA,nrow=N,ncol=L)
	for(l in 1:3)	{
		kk=kde(x,h=H[l],eval.points=xtest)   
		A[,l] = kk$estimate
	}
	for(l in 4:6)	{
		kk=KDEtr(x,xtest,h=H[l])   
		A[,l] = kk
	}
	A
}

estimalpha = function(A,eps=0.0001,itermax=3000)
{
# estime les alphas 
# renvoie les valeurs estim?es de alpha ainsi que la courbe
# de la log vraisemblance du m?lange estim?
	A = A[which(rowSums(A) > 1e-10),]
# les lignes de A o? il n'y a que des zeros sont retir?es...
	L = ncol(A)
	N = nrow(A)
#	print(summary(rowSums(A)))
#	alpha = runif(L)
#	alpha = alpha / sum(alpha)
	alpha = rep(1/L,L)
	iter = 0
	Vrais = NULL
	while(iter <itermax)	{
		Malpha = matrix(alpha,nrow=N,ncol=L,byrow=T)
		coeffs = A %*% (alpha)
		Vrais = c(Vrais, mean( ifelse(coeffs>0 ,log(coeffs),0)  ))
		Mcoeffs = matrix(coeffs,nrow=N,ncol=L)
		aux = Malpha/Mcoeffs
		M = aux * A
		alphanew = colMeans(M)
		ecart = sum(abs(alpha - alphanew))
# l'ecrat quadratique entre deux valeurs successives de alpha
# est tr?s faible, voir le papier sur stacking
		iter = iter +1
		alpha = alphanew
		if(ecart < eps) break
	}
	list(alpha=alpha,LogL = Vrais)
}


stack.dens = function(don=rnorm(1000),xtest=NULL,pl=F,H=c(.1,.2,.3,.1,.2,.3),V=10,eps=0.0001,itermax=3000)
{


# Exemple d'utilisation
#	stack.dens()
# pour un echantillon test constitu? d'une grille r?guli?re
# fabriquer une grille r?guli?re pour estimer la densit?
#	don= rnorm(1000)
#	pas = (max(don) - min(don))/ (length(don)-1)
#	xtest =  seq(min(don),max(don),pas)
#	stack.dens(don, xtest)

	mat = genereA(don,H,V)
	toto = estimalpha(mat,eps=eps,itermax=itermax)
	if(pl) par(mfrow=c(1,2))
# pour voir l'?volution de la vraisemblance
	if(pl) plot(toto$LogL,type="l")
	if(is.null(xtest)) test = don
	else test = xtest
	mat = genereB(don,xtest=test,H)
	prev = mat%*%toto$alpha 
	prev = prev[order(test)]
# pour voir les valeurs de la densit? estim?e
	if(pl) plot(prev,type= "l")
	list(prev=prev,alpha = toto$alpha)	
}


#######################
#Stacking d'histogrammes
######################

genereAh = function(x,H=c(5,10,20,30,40,50),V=10){

# genere la matrice de fhat(E_{-v},E_v) par validation crois?e
	L = length(H)
	N = length(x)
	A = matrix(NA,nrow=N,ncol=L)
	tailpaquets = N%/%V
	reste = N%%V
	indices = rep(1:V,tailpaquets)
	indices = c(indices,sample(V,reste))
	indices = sample(indices)
	for(l in 1:L)
		for(v in 1:V)	{
			bloque = which(indices ==v)
			kk=hist(x[-bloque],breaks=mybreaks(x[-bloque],H[l]),plot=F,warn.unused = F)   
			A[bloque,l] = predict.hist(kk,x[bloque])
		}
	A
}




genereBh = function(x,xtest=NULL,H=c(5,10,20,30,40,50))
{
# genere la matrice des densit?s estim?es sur l'?chantillon d'apprentissage
# ou sur un echantillon test
	L = length(H)
	if(is.null(xtest)) xtest=x
	N = length(xtest)
	A = matrix(NA,nrow=N,ncol=L)
	for(l in 1:L)	{
		kk=hist(x,breaks=mybreaks(x,H[l]),plot=F,warn.unused = F)   
		A[,l] = predict.hist(kk,xtest)
	}
	A
}


estimalpha = function(A,eps=0.0001,itermax=3000)
{
# estime les alphas 
# renvoie les valeurs estim?es de alpha ainsi que la courbe
# de la log vraisemblance du m?lange estim?
	A = A[which(rowSums(A) > 1e-10),]
# les lignes de A o? il n'y a que des zeros sont retir?es...
	L = ncol(A)
	N = nrow(A)
#	print(summary(rowSums(A)))
#	alpha = runif(L)
#	alpha = alpha / sum(alpha)
	alpha = rep(1/L,L)
	iter = 0
	Vrais = NULL
	while(iter <itermax)	{
		Malpha = matrix(alpha,nrow=N,ncol=L,byrow=T)
		coeffs = A %*% (alpha)
		Vrais = c(Vrais, mean( ifelse(coeffs>0 ,log(coeffs),0)  ))
		Mcoeffs = matrix(coeffs,nrow=N,ncol=L)
		aux = Malpha/Mcoeffs
		M = aux * A
		alphanew = colMeans(M)
		ecart = sum(abs(alpha - alphanew))
# l'ecrat quadratique entre deux valeurs successives de alpha
# est tr?s faible, voir le papier sur stacking
		iter = iter +1
		alpha = alphanew
		if(ecart < eps) break
	}
	list(alpha=alpha,LogL = Vrais)
}


stack.denshist = function(don=rnorm(1000),xtest=NULL,pl=F,H=c(5,10,20,30,40,50),V=10,eps=0.0001,itermax=3000)
{
	mat = genereAh(don,H,V)
	toto = estimalpha(mat,eps=eps,itermax=itermax)
	if(pl) par(mfrow=c(1,2))
# pour voir l'?volution de la vraisemblance
	if(pl) plot(toto$LogL,type="l")
	if(is.null(xtest)) test = don
	else test = xtest
	mat = genereBh(don,xtest=test,H)
	prev = mat%*%toto$alpha 
	prev = prev[order(test)]
# pour voir les valeurs de la densit? estim?e
	if(pl) plot(prev,type= "l")
	list(prev=prev,alpha = toto$alpha)	
}


#################
# Tsybakov AggPure
#################

tsybakov = function(dlearn,dtest, S=10,H=c( 0.001,0.005, 0.01,0.05,0.1,0.5))
{
# dd 	:	data
# S		:	number of splitting for teh data	
	n = length(dlearn)
	prev = 0
	for(s in 1:S) {
		ind = sample(n,n/2)
		Bapp = genereB(dlearn[ind],H=H)
		Btest0 = genereB(dlearn[ind],dlearn[-ind],H=H)
		Btest = genereB(dlearn[ind],dtest,H=H)
		alphs = estimalpha(Btest0)$alpha
		prevs = Btest %*% alphs 
		prev = prev + prevs
	}
	prev = prev[order(dtest)]
	prev/S
}

###################
# Dimarzio BoostKde
###################

#maj.poids = function(aux,w=NULL,h=0.1)
#{
# condamn?e au 16 03 2012
# aux : valeur du noyau K((x - xj) /h)
#
#	if(is.null(w)) w=rep(1/n,n)
#	Ai = (rowSums(w*aux) - diag(w*aux)) /h  
#	rapp = (1-w) * (1 +  (w/ (h * Ai * sqrt(2*pi))) )
#	log(rapp)
#}


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

simulations = function(nummodel=1,m=5,n=100,K=5,nbr,dmp=T)
{
# dmp :        If T save results
    nberr = 2
 resHist=resASH=resAggHistu=resRASH=resBagHist=resSama=resKdenrd=resKdenrd0=resKdeucv=resKdesj=resStack=resStackH=resBoost=resTsyb=matrix(NA,nrow=m,ncol=nberr)
#    res=matrix(NA,nrow=m,ncol=nberr)
    ficdon =paste("data/data-",n,"-",nummodel,".txt",sep="")
    if(!file.exists(ficdon)) stop("Data file does not exist")
    source(ficdon)
    for(i in 1:m)
    {
        dd = ldata[[i]]
#HISTOGRAMME
        zz=hist(train,breaks=mybreaks(train,nbr[1]),plot=F,warn.unused = F)
        predictionssurHISTO=predict.hist(zz,dd$test)
        resHist[i,] = error(dd$dobs,predictionssurHISTO)
#ASH
        ash = avshift(train,dd$test,nbr[2],nbr[3])
        resASH[i,] = error(dd$dobs,ash)
#AGGHISTUNI
        agghistu = AgregHistunif(train,dd$test,nbr[4],K)
        resAggHistu[i,] = error(dd$dobs,agghistu)
#RASH
        finrash = rash(train,dd$test,nbr[5],K)
        resRASH[i,] = error(dd$dobs,finrash)
#BAGHIST
        BH = BagHist(train,dd$test,nbr[6],K)
        resBagHist[i,] = error(dd$dobs,BH)
#kdenrd
        finkdenrd = onekdenrd(train,dd$test)
        resKdenrd[i,] = error(dd$dobs,finkdenrd)
#kdenrd0
        finkdenrd0 = onekdenrd0(train,dd$test)
        resKdenrd0[i,] = error(dd$dobs,finkdenrd0)
#kdeucv
        finkdeucv = onekdeucv(train,dd$test)
        resKdeucv[i,] = error(dd$dobs,finkdeucv)
#kdesj
        finkdesj = onekdesj(train,dd$test)
        resKdesj[i,] = error(dd$dobs,finkdesj)
#STACKING
        stacking=stack.dens(train,xtest=dd$test)
        resStack[i,]=error(dd$dobs,stacking$prev)
#Stacking Histogrammes
        stackh= stack.denshist(train,dd$test,pl=F,H=c(5,10,20,30,40,50),V=10,eps=0.0001,itermax=3000)
        resStackH[i,]=error(dd$dobs,stackh$prev)
#boosting
        boost =boostkde(train,xtest=dd$test,K=5,h=bw.nrd(train))
        resBoost[i,]=error(dd$dobs,boost)
#Tsybakov
        tsyb =tsybakov(train,dtest=dd$test)
        resTsyb[i,]=error(dd$dobs,tsyb)
#Samarov
        sama =samarov(train,dd$test,K,1)
        resSama[i,]=error(dd$dobs,sama$predsam)
    }
    fic = paste("resultats/",paste(nummodel,m,n,K,sep="-"),".txt",sep="")
    lres = list(resHist =resHist, resASH =resASH ,resAggHistu = resAggHistu, resRASH = resRASH,resBagHist = resBagHist,resKdenrd= resKdenrd,resKdenrd0 = resKdenrd0,
 resKdeucv=resKdeucv,resKdesj=resKdesj,resStack=resStack,resStackH=resStackH,resBoost=resBoost,resTsyb=resTsyb,resSama=resSama)
    if(dmp) dump("lres",file= fic)
    erreurs  = rbind(colMeans(resHist),colMeans(resASH),colMeans(resAggHistu),colMeans(resRASH),colMeans(resBagHist),colMeans(resKdenrd),colMeans(resKdenrd0),colMeans(resKdeucv),colMeans(resKdesj),colMeans(resStack),colMeans(resStackH),colMeans(resBoost),colMeans(resTsyb),colMeans(resSama))
#    var.err  = rbind(colVars(resHist),colVars(resAggHist),colVars(resBagHist),colVars(resKdenrd0),colVars(resKdeucv),colVars(resStack),colVars(resStackH),colVars(resBoost),colVars(resTsyb))
#    if(dmp) dump("res",file= fic)
 rownames(erreurs)=c("Histogramme","ASH","AgregHist","RASH","BagHist","Nrd","Nrd0","UCV","SJ","Stack","StackH","BoostKde", "Aggpure","Samarov")
    colnames(erreurs)=c("MISE", "LogL")
#    erreurs  = colMeans(res)
    erreurs
}

simulationsAHBHRASH= function(nummodel=1,m=5,n=100,K=5,nbr,dmp=T,alpha) 
{
# dmp :		If T save results
	nberr = 2
	#resHist=resAggHist=resAggHistu=resAggHisttr=resBagHist=resKdenrd0=matrix(NA,nrow=m,ncol=nberr)
	resHist=resAggHist=resAggHistu=resAggHisttr=resBagHist=matrix(NA,nrow=m,ncol=nberr)

	for(i in 1:m)
	{
		dd = gendata(nummodel,n)
		
#HISTOGRAMME
		zz=hist(train,breaks=mybreaks(train,nbr[1]),plot=F,warn.unused = F)
		predictionssurHISTO=predict.hist(zz,dd$test)
		resHist[i,] = error(dd$dobs,predictionssurHISTO)
		
#AGGREG HISTOGRAMMES NORMAL
		agghist = AgregHist(train,dd$test,nbr[2],K,alpha) 
		resAggHist[i,] = error(dd$dobs,agghist)

#AGGREG HISTOGRAMMES UNIF
		agghistunif = AgregHistunif(train,dd$test,nbr[3],K) 
		resAggHistu[i,] = error(dd$dobs,agghistunif)
		
#RASH
		finrash = rash(train,dd$test,nbr[4],K) 
		resAggHisttr[i,] = error(dd$dobs,finrash)		

#BAGHISTOGRAMMES
		baghist = BagHist(train,dd$test,nbr[5],K) 
		resBagHist[i,] = error(dd$dobs,baghist)		
		
#kdenrd0
#		finkdenrd0 = onekdenrd0(train,dd$test)
#		resKdenrd0[i,] = error(dd$dobs,finkdenrd0)	
		
	}
	fic = paste("resultats/",paste(nummodel,m,n,K,sep="-"),".txt",sep="")
#	if(dmp) dump(c("resHist","resAggHist","resAggHistu","resAggHisttr","resBagHist","resKdenrd0"),file= fic)
	if(dmp) dump(c("resHist","resAggHist","resAggHistu","resAggHisttr","resBagHist"),file= fic)

#	erreurs  = rbind(colMeans(resHist),colMeans(resAggHist),colMeans(resAggHistu),colMeans(resAggHisttr),colMeans(resBagHist),colMeans(resKdenrd0))
	erreurs  = rbind(colMeans(resHist),colMeans(resAggHist),colMeans(resAggHistu),colMeans(resAggHisttr),colMeans(resBagHist))

	#	var.err  = rbind(colVars(resHist),colVars(resAggHist),colVars(resBagHist),colVars(resKdenrd0),colVars(resKdeucv),colVars(resStack),colVars(resStackH),colVars(resBoost),colVars(resTsyb))
#	if(dmp) dump("res",file= fic)
	rownames(erreurs)=c("Hist", "AggregHist","AggregHistu","AggHistTr","BagHist")
	colnames(erreurs)=c("MISE", "LogL")
#	erreurs  = colMeans(res)	
	erreurs
}





colVars = function(mat)
{
	apply(mat,2,var,na.rm=T)

}


moyerreur=function(nummodel=1,n=100,K=50,m=5,breaks =50,alpha)
{
	dd = gendata(nummodel,n)
	ErrmoyAH = ErrmoyBH= ErrmoyBH2 = ErrmoyBH3 = 0
	
	for(kk in 1:m)
	{
#AGREGATION HISTOGRAMMES
		errcAH=AgregHist.err(dd$train,dd$test,breaks,K,dd$dobs,alpha)$erreur
		ErrmoyAH=ErrmoyAH+errcAH
#BAGGING HISTOGRAMMES
		errcBH=BagHist.err(dd$train,dd$test,breaks,K,dd$dobs)$erreur
		ErrmoyBH=ErrmoyBH+errcBH
#BAGGING KDE
		errcBH2=BagHist2.err(dd$train,dd$test,breaks,K,dd$dobs,alpha)$erreur
		ErrmoyBH2=ErrmoyBH2+errcBH2
#BAGGING KDE
		errcBH3=BagHist3.err(dd$train,dd$test,breaks, K,dd$dobs,alpha)$erreur
		ErrmoyBH3=ErrmoyBH3+errcBH3
	}
	MSE= cbind(ErrmoyAH[,1],ErrmoyBH[,1],ErrmoyBH2[,1],ErrmoyBH3[,1]) / m
	LogL= cbind(ErrmoyAH[,2],ErrmoyBH[,2],ErrmoyBH2[,2],ErrmoyBH3[,2]) / m
	colnames(MSE) = colnames(LogL) =  c("AgregHist","BagHist","BagHist2","BagHist3")
	list(MSE=MSE,LogL=LogL)
}



moyerreurRASH=function(nummodel=1,n=100,K=50,m=5,brRASH)
{
	dd = gendata(nummodel,n)
	ErrmoyRASH=0
	
	for(kk in 1:m)
	{
	#RASH
		errcRASH=rash.err(dd$train,dd$test,brRASH,K,dd$dobs)$erreur
		ErrmoyRASH=ErrmoyRASH+errcRASH
		}
	ErrmoyRASH/ m
#	LogL=ErrmoyRASH / m
#	colnames(MISE) = colnames(LogL) =  c("AggregHist","BagHist")
#	list(Err=Err)
}



moyerreurAHBHRASH=function(nummodel=1,n=100,K=50,m=5,brAH,brBH,brRASH)
{
	dd = gendata(nummodel,n)
	ErrmoyAH = ErrmoyBH=ErrmoyRASH=0
	
	for(kk in 1:m)
	{
#AGREGATION HISTOGRAMMES
		errcAH=AgregHistunif.err(dd$train,dd$test,brAH,K,dd$dobs)$erreur
		ErrmoyAH=ErrmoyAH+errcAH
#BAGGING HISTOGRAMMES
		errcBH=BagHist.err(dd$train,dd$test,brBH,K,dd$dobs)$erreur
		ErrmoyBH=ErrmoyBH+errcBH
#RASH
		errcRASH=AgregHisttr.err(dd$train,dd$test,brRASH,K,dd$dobs)$erreur
		ErrmoyRASH=ErrmoyRASH+errcRASH		
		}
	MISE= cbind(ErrmoyAH[,1],ErrmoyBH[,1],ErrmoyRASH[,1]) / m
	LogL= cbind(ErrmoyAH[,2],ErrmoyBH[,2],ErrmoyRASH[,2]) / m
	colnames(MISE) = colnames(LogL) =  c("AggregHist","BagHist","RASH")
	list(MISE=MISE,LogL=LogL)
}



moyash.err = function(nummodel=1,n=100,B=50,m=5,brASH=10) {

	dd = gendata(nummodel,n)
	ErrmoyASH=0
	
	for(kk in 1:m)
	{
		#ASH
		errash=ash.err(dd$train,dd$test,brASH,B,dd$dobs)$erreur
		ErrmoyASH=ErrmoyASH+errash
		}
	ErrmoyASH/m
#	colnames(MISE) = colnames(LogL) 
#	list(MISE=MISE,LogL=LogL)
}


simulationsash = function(nummodel=1,m=5,n=100,K=5,nbr,dmp=T) 
{
  ficdon =paste("data/data-",n,"-",nummodel,".txt",sep="")
    if(!file.exists(ficdon)) stop("Data file does not exist")
    source(ficdon)
	resASH=matrix(NA,nrow=m,ncol=2)
	for(i in 1:m)
	{
	     dd = ldata[[i]]	
#ASH
		ashift = avshift(dd$train,dd$test,nbr,K) 
		resASH[i,] = error(dd$dobs,ashift)
	}
	fic = paste("resultats/ASH-",paste(nummodel,m,n,sep="-"),".txt",sep="")
	if(dmp) dump("resASH",file= fic)
	lres = resASH
    if(dmp) dump("lres",file= fic)
	erreurs  = colMeans(resASH)
	erreurs
}


simulations2 = function(nummodel=1,m=5,n=100,K=5,nbr,dmp=T)
{
#nbr sera le vecteur de breaks et/ou de feuilles pour les histogrammes
# dmp :        If T save results
    nberr = 2
 resHist=resHistG=resHistC=resASH=resRASH=resRASHG=resRASHC=resSama=resKdenrd=resKdenrd0=resKdeucv=resKdesj=resKdenrdt=resKdenrd0t=resKdeucvt=resKdesjt=matrix(NA,nrow=m,ncol=nberr)
#    res=matrix(NA,nrow=m,ncol=nberr)
    ficdon =paste("data200/data200-",n,"-",nummodel,".txt",sep="")
    if(!file.exists(ficdon)) stop("Data file does not exist")
    source(ficdon)
    for(i in 1:m)	{
        dd = ldata[[i]]
		train= dd$train
		test = dd$test
		obs = dd$obs
#HISTOGRAMME
        zz=hist(train,breaks=mybreaks(train,nbr[1]),plot=F,warn.unused = F)
        predictionssurHISTO=predict.hist(zz,dd$test)
        resHist[i,] = error(dd$dobs,predictionssurHISTO)
#HISTOGRAMME Greedy		
		zz2=hist(train,breaks=grille(train,nbr[2]),plot=F,,warn.unused = F)
        predictionssurHISTO2=predict.hist(zz2,dd$test)
        resHistG[i,] = error(dd$dobs,predictionssurHISTO2)
#Histogramme avec CART
		d=data.frame(train)
		yy=train
		d=cbind(d,yy)
		toto=rpart(yy~train,data=d)
		zz3=hist(train,breaks=c(min(train),toto$splits[,4],max(train)),plot=F,warn.unused = F)
		predictionssurHISTO3=predict.hist(zz3,dd$test)
		resHistC[i,] = error(dd$dobs,predictionssurHISTO3)				
#ASH
		ash = avshift(train,dd$test,nbr[3],K)
		resASH[i,] = error(dd$dobs,ash)
#RASH
		finrash = rash(train,dd$test,nbr[4],K)
		resRASH[i,] = error(dd$dobs,finrash)
#RASH greedy (finlandais)
		finrashg = rashgreedy(train,dd$test,nbr[5],K)
		resRASHG[i,] = error(dd$dobs,finrashg)
#RASH avec Cart
		finrashc = rashcart(train,dd$test,K)
		resRASHC[i,] = error(dd$dobs,finrashc)
#kdenrd
		finkdenrd = onekdenrd(train,dd$test)
		resKdenrd[i,] = error(dd$dobs,finkdenrd)
#kdenrd0
		finkdenrd0 = onekdenrd0(train,dd$test)
		resKdenrd0[i,] = error(dd$dobs,finkdenrd0)
#kdeucv
		finkdeucv = onekdeucv(train,dd$test)
		resKdeucv[i,] = error(dd$dobs,finkdeucv)
#kdesj
		finkdesj = onekdesj(train,dd$test)
		resKdesj[i,] = error(dd$dobs,finkdesj)
#KdeNrdTriangulaire
		a=KDEtr(train,sort(dd$test),h=bw.nrd(train))
		resKdenrdt[i,] = error(dd$dobs,a)
#KdeNrd0Triangulaire
		b=KDEtr(train,sort(dd$test),h=bw.nrd0(train))
		resKdenrd0t[i,] = error(dd$dobs,b)	
#KdeUcvTriangulaire
		c=KDEtr(train,sort(dd$test),h=bw.ucv(train))
		resKdeucvt[i,] = error(dd$dobs,c)	
#KdeSJTriangulaire
		d=KDEtr(train,sort(dd$test),h=bw.SJ(train))
		resKdesjt[i,] = error(dd$dobs,d)		
#Samarov
		sama =samarov(train,dd$test,K,1)
		resSama[i,]=error(dd$dobs,sama$predsam)
    }
	fic = paste("resultats/",paste(nummodel,m,n,K,sep="-"),".txt",sep="")
	lres = list(resHist =resHist,resHistG=resHistG,resHistC=resHistC, resASH =resASH , resRASH = resRASH, resRASHG=resRASHG,resRASHC = resRASHC,resKdenrd= resKdenrd,resKdenrd0 = resKdenrd0,
	resKdeucv=resKdeucv,resKdesj=resKdesj,resKdenrdt= resKdenrdt,resKdenrd0t = resKdenrd0t,
	resKdeucvt=resKdeucvt,resKdesjt=resKdesjt,resSama=resSama)
	if(dmp) dump("lres",file= fic)
	erreurs  = rbind(colMeans(resHist),colMeans(resHistG),colMeans(resHistC),colMeans(resASH),colMeans(resRASH),colMeans(resRASHG),colMeans(resRASHC),colMeans(resKdenrd),colMeans(resKdenrd0),colMeans(resKdeucv),colMeans(resKdesj),colMeans(resKdenrdt),colMeans(resKdenrd0t),colMeans(resKdeucvt),colMeans(resKdesjt),colMeans(resSama))
	rownames(erreurs)=c("Hist","Histgreedy","Histcart","ASH","RASH","RASHGreedy","RASHCart","KdeNrd","KdeNrd0","KdeUCV","KdeSJ","KdeNrdT","KdeNrd0T","KdeUCVT","KdeSJT","Samarov")
	colnames(erreurs)=c("MISE", "LogL")
    erreurs
}


optbreaks= function(nummodel=1,m=5,n=100,K=5,nbr,dmp=T) 
{
# dmp :		If T save results
	nberr = 2
	resHist=resHistG=resASH=resRASH=resRASHG=matrix(NA,nrow=m,ncol=nberr)

	for(i in 1:m)
	{
		dd = gendata(nummodel,n)
#HISTOGRAMME
        zz=hist(train,breaks=mybreaks(train,nbr),plot=F,warn.unused = F)
        predictionssurHISTO=predict.hist(zz,dd$test)
        resHist[i,] = error(dd$dobs,predictionssurHISTO)

#HISTOGRAMME Greedy		
		zz2=hist(train,breaks=grille(train,nbr),plot=F,warn.unused = F)
        predictionssurHISTO2=predict.hist(zz2,dd$test)
        resHistG[i,] = error(dd$dobs,predictionssurHISTO2)				
		
#ASH
        ash = avshift(train,dd$test,nbr,K)
        resASH[i,] = error(dd$dobs,ash)
#RASH
        finrash = rash(train,dd$test,nbr,K)
        resRASH[i,] = error(dd$dobs,finrash)

#RASH greedy (finlandais)
		finrashg = rashgreedy(train,dd$test,nbr,K)
        resRASHG[i,] = error(dd$dobs,finrashg)
	}
	fic = paste("resultats/",paste(nummodel,m,n,K,sep="-"),".txt",sep="")

	if(dmp) dump(c("resHist","resHistG","resASH","resRASH","resRASHG"),file= fic)
	erreurs  = rbind(colMeans(resHist),colMeans(resHistG),colMeans(resASH),colMeans(resRASH),colMeans(resRASHG))
	rownames(erreurs)=c("Hist","Histgreedy","ASH","RASH","RASHGreedy")
	colnames(erreurs)=c("MISE", "LogL")

	erreurs
}



rashfp.err = function(xx,grille=aa,nbr = 50, B= 10,dobs,alpha=1) {
  fin = 0
  err00=NULL
  zz = hist(xx,breaks=mybreaks(xx,nbr),plot=F,warn.unused = F)$breaks	
  mx = min(xx)
  Mx = max(xx)
  for(i in 1:B)
  {
    blabla=sqrt(abs(min(diff(zz))))
    newb = zz + rnorm(length(zz),0,alpha * blabla)
    newb=sort(newb)
    if(min(newb) > mx) newb= c(mx,newb)
    if(max(newb) < Mx) newb= c(newb, Mx)
    hs2=hist(xx,breaks=newb,plot=F,warn.unused = F)
    fin= fin + approxfun(x=hs2$mids,y=hs2$density)(grille)
    previ=fin/i
    err00=rbind(err00,error(dobs,previ))
    #if(i%%20 == 0) cat(i,">>")
  }
  #cat("\n")
  list(prev=fin/B,erreur=err00[,1])
}


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


bropt=function(x){
  Mgrid <- 2:(5 * floor(sqrt(length(x))))
  J     <- numeric(length(Mgrid))
  for(m in seq_along(Mgrid)) {
    J[m] <- riskhist(obs=x, Mgrid[m], xlim = c(min(x)-0.5, max(x)+0.5))
  #plot(1 / Mgrid, J, type = 'l')
  }
  list(opt=Mgrid[which.min(J)])
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
  grille <- data$test
  n <- length(xx)
  mat <- matrix(ncol = length(grille), nrow = B)
  
  for (i in 1:B) {
    xb  = xx[sample(n, replace = TRUE)]
    mybreaks <- mybreaks(xb, nbr)
    hs2 = hist(xb, breaks = mybreaks, plot = FALSE)
    
    mat[i, ] <- predict.hist(hh = hs2, x = grille)
  }
  
  punctual <- colMeans(mat, na.rm = TRUE)
  tube     <- baseboot(data$test, punctual, mat, conf = conf)
  
  return(list(punctual = punctual, bunch = mat, tube = tube)) 
}

tube_fp <- function(data, nbr = 50, B = 10, conf = c(0.05, 0.95)) {
  # A chaque etape, on prend un nouveau jeu de donn?es, on construit un 
  # histogramme, on pr?dit et on agr?ge.
  
  xx <- data$train
  grille <- data$test
  n <- length(xx)
  mat <- matrix(ncol = length(grille), nrow = B)
  
  for (i in 1:B) {
    xb  = xx[sample(n, replace = TRUE)]
    mybreaks <- mybreaks(xb, nbr)
    hs2 = hist(xb, breaks = mybreaks, plot = FALSE)
    
    mids <- hs2$mids
    delta <- diff(mybreaks)[1]
    # This avoids NA on prediction for test data between the min on grille
    # and the min obs value on xb
    newmids <- c(min(mids) - delta, mids, max(mids) + delta) 
    
    mat[i, ] <- approxfun(x = newmids, y = c(0, hs2$density, 0))(grille)
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
  grille <- data$test
  n <- length(xx)
  mat <- matrix(ncol = length(grille), nrow = B)
  
  for (i in 1:B) {
    xb  = xx[sample(n, replace = TRUE)]
    mat[i, ] <- onekdeucv(xx = xb, grille = grille)
  }
  
  punctual <- colMeans(mat, na.rm = TRUE)
  tube     <- baseboot(data$test, punctual, mat, conf = conf)
  
  # usando la escala raiz  => escala invariante al usar quantiles
  #  tube     <- baseboot(data$test, sqrt(punctual), sqrt(mat))
  #  tube[tube < 0] <- 0
  #  tube <- tube^2
  
  return(list(punctual = punctual, bunch = mat, tube = tube)) 
}


tube_KDEsm <- function(data) { # only at 90% (hard coded)
  fit <- sm::sm.density(data$train, h = bw.ucv(data$train), display = "none",
                        eval.grid = TRUE, eval.points = data$test)
  return(list(punctual = fit$estimate, 
              tube = rbind(fit$lower, fit$upper)))  
}


tube_wass <- function(data, alpha = 0.1, nbr = 50, plot = FALSE, ...) {
  xb  <- data$train
  n   <- length(data$train)
  
  hs <- hist(xb, breaks = mybreaks(xb, nbr), 
             warn.unused = FALSE, probability = TRUE, plot = plot, ...)
  
  m   <- length(hs$counts)
  c   <- qnorm(p = 1 - alpha / 2 / m) * sqrt(m / n)
  
  ln1    <- c(pmax(sqrt(hs$density) - c, 0)^2, 0)
  un1    <- c((sqrt(hs$density) + c)^2       , 0)
  
  hs2 <- predict.hist(hs, data$test)
  
  ln2    <- pmax(sqrt(hs2) - c, 0)^2
  un2    <- (sqrt(hs2) + c)^2
  
  if (plot) {
    lines(data$test, ln2, type = 's', col = 4, lwd = 2)
    lines(data$test, un2, type = 's', col = 4, lwd = 2)
  }
  
  return(list(punctual = hs2,
              tube = rbind(ln2, un2)))
}

#mymat <- matrix(rnorm(2000), nrow = 80)
#mypunctual <- colMeans(mymat)
#myx   <- seq(length(mypunctual))

#uno  <- baseboot(myx, mypunctual, mymat, plot = TRUE)
#dos  <- npheuristic(myx, mypunctual, mymat, plot = TRUE)
#tres <- kfwe(myx, mypunctual, mymat, k = 5, plot = TRUE)

#error_tube(mypunctual, uno)
#error_tube(mypunctual, dos)
#error_tube(mypunctual, tres)

#' @title Data generation whith gendata  function .
#'
#' @description This function provides the Training set, the Test set and density values
#' from a distribution law.
#' @param nummodel number to choose a dristribution low.
#' @param n sample size .
#'
#' @return the training set, the test set and density values .
#'
#' @export
#' @import stats  nor1mix
#' @examples
#' gendata(10,100)
#'
gendata= function(nummodel = 1, n = 100)
{
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
    train = benchden::rberdev(n,23)
    minm = min(train)
    maxm = max(train)
    test = sort(benchden::rberdev(n,23))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = benchden::dberdev(test,23)
  }

  #SMOOTH COMB DENSITY . IL FAUT TELECHARGER LE PAQUET BENCHDEN, c'est le numero 24.
  if(nummodel==9) {
    train = benchden::rberdev(n,24)
    minm = min(train)
    maxm = max(train)
    test = sort(benchden::rberdev(n,24))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = benchden::dberdev(test,24)
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
    train = nor1mix::rnorMix(n, MW.nm4)
    minm = min(train)
    maxm = max(train)
    test = sort(nor1mix::rnorMix(n, MW.nm4))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = nor1mix::dnorMix(test, MW.nm4)
  }
  #AsimClaw
  if(nummodel==15) {
    train = nor1mix::rnorMix(n, MW.nm12)
    minm = min(train)
    maxm = max(train)
    test = sort(nor1mix::rnorMix(n, MW.nm12))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = nor1mix::dnorMix(test, MW.nm12)
  }
  #AsimDoubleClaw
  if(nummodel==16) {
    train = nor1mix::rnorMix(n, MW.nm13)
    minm = min(train)
    maxm = max(train)
    test = sort(nor1mix::rnorMix(n, MW.nm13))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = nor1mix::dnorMix(test, MW.nm13)
  }

  #SmooComb
  if(nummodel==17) {
    train = nor1mix::rnorMix(n, MW.nm14)
    minm = min(train)
    maxm = max(train)
    test = sort(nor1mix::rnorMix(n, MW.nm14))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = nor1mix::dnorMix(test, MW.nm14)
  }

  #DistComb
  if(nummodel==18) {
    train = nor1mix::rnorMix(n, MW.nm15)
    minm = min(train)
    maxm = max(train)
    test = sort(nor1mix::rnorMix(n, MW.nm15))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = nor1mix::dnorMix(test, MW.nm15)
  }

  #DistBim
  if(nummodel==19) {
    train = nor1mix::rnorMix(n, MW.nm16)
    minm = min(train)
    maxm = max(train)
    test = sort(nor1mix::rnorMix(n, MW.nm16))
    test = test[ (test <= maxm) & (test> minm)]
    dobs = nor1mix::dnorMix(test, MW.nm16)
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

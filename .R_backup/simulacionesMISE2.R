rm(list=ls())
source("functions2.R")


simulaciones=function(n = 100, M = 10, B = 150){

  vec = c(1,3,5,7,8,11,10,13,17,19,20,21)
  AA = matrix(0, nrow = length(vec), ncol = 7)
  #AA=matrix(0,nrow=2,ncol=7)
  colnames(AA)=c("H","FP","Kde", "BagHist", "BagFP", "BagKde", "Rash")
  #rownames(AA)=c("normal","chi2","mezcla1","mezcla2","bart","triangular")
  
  
  for(i in 1:M){
    A = matrix(NA, nrow = length(vec), ncol = 7)
    #A=matrix(NA,nrow=2,ncol=7)
  for(ll in 1:length(vec)){
    print(ll)
    dd = gendata(vec[ll],n)
    bopt = bropt(dd$train)$opt  # needed for Hist & RASH
    zz = hist(dd$train,breaks=mybreaks(dd$train,nbr=bopt),plot=F)
    
    bhist = BagHistfp(xx=dd$train, grille=dd$test, B)
    modelrash = rash(dd$train, grille = dd$test, nbr = bopt, B)
    #modelavshift(dd$train,dd$test,nbr=bopt,M=B)
    modelbagkde <- Bagkde(xx = dd$train, grille = dd$test, B)
    
  A[ll,]=c(error(dd$dobs,predict.hist(zz,dd$test))[1],
           error(dd$dobs,approxfun(x=zz$mids,y=zz$density)(dd$test))[1],    
           error(dd$dobs,onekdeucv(dd$train,dd$test))[1],
           error(dd$dobs,bhist$bh)[1],
           error(dd$dobs,bhist$bhfp)[1], 
           #error(dd$dobs,modelavshift)[1]
           error(dd$dobs,modelbagkde)[1],
           error(dd$dobs,modelrash)[1]) 
  }
    AA=AA+A
  }
  AA=AA/M
}


#system.file(
#res <- lapply(c(50, 100, 200, 500, 1000),
#              function(n) simulaciones(n = n, M = 100, B = 200))
#)

library(parallel)
vars2export <- c("BagHistfp", "Bagkde", "bropt", "broptfp", "dtriangle", 
                 "dberdev", "error", "ind", "rnorMix", "MW.nm14", "dnorMix",
                 "MW.nm16",
                 "gendata", "kde", "mel", "melange", "mybreaks", "onekdeucv",
                 "predict.hist", "predict.hist.x",   "rash", "rberdev", 
                 "riskhist", "riskfp","rtriangle", "simulaciones")   
cls <- makeCluster(detectCores() - 1)
clusterExport(cls, vars2export)
system.file(
  res <- parLapplyLB(cls, rev(c(20, 50, 100, 200, 500, 1000, 2000)),
                     function(n) simulaciones(n = n, M = 100, B = 200))
)
stopCluster(cls)


save(res, file = "resultados_mise-testjc.Rdata")

load("resultados_mise-testjc.Rdata")

vec=c(1,2,3,5,6,8,11,12)
library(xtable)
i=2
mise=(res[[i]]*100)[vec,]
rownames(mise)=c("Normal","Chi2","Mixture","Claw","Triangular","Uniform","Tsybakov","Uniform Mixt.")
xtable(mise,digits=4)





#save(res, file = "resultados_mise_new.Rdata")
#save(res,file="resultados_mise_modelos10y13.Rdata")
#save(res,file="resultados_mise_modelosMandW.Rdata")

# test para ver si la funcion riskfp esta bien
#rm(list = ls())
#source("functions2.R")
#vec <- 1 #13 
#M <- 3
#B <- 200
#n <- 2000
#dd <- gendata(vec, n)
#bopt <- bropt(dd$train)$opt  # needed for Hist & RASH
#zz <- hist(dd$train,breaks=mybreaks(dd$train,nbr=bopt),plot=F)
#bhist <- BagHistfp(xx=dd$train, grille=dd$test, B)
#error(dd$dobs,bhist$bhfp)[1]

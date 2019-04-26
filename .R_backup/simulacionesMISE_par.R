library(parallel)

source("functions.r")

simulaciones=function(n=100,M=10,B=150){
 
  require(ks)
  require(benchden)
  AA=matrix(0,nrow=6,ncol=7)
  colnames(AA)=c("H","FP","Kde", "BagHist", "BagHistFP", "Bagkde", "Rash")
  rownames(AA)=c("normal","chi2","mezcla1","mezcla2","bart","triangular")
  for(i in 1:M){
    
   
   A=matrix(NA,nrow=6,ncol=7)
   j=1
  for(ll in c(1,3,5,7,8,11)){
    
    dd=gendata(ll,n)
    bopt=bropt(dd$train)$opt
    zz=hist(dd$train,breaks=mybreaks(dd$train,nbr=bopt),plot=F)
    
    bhist=BagHistfp(xx=dd$train,grille=dd$test,nbr = bopt, B)
    modelrash=rash(dd$train,grille=dd$test,nbr = bopt, B)
    #modelavshift(dd$train,dd$test,nbr=bopt,M=B)
    modelbagkde <- Bagkde(xx = dd$train, grille = dd$test, B)
    
  A[j,]=c(error(dd$dobs,predict.hist(zz,dd$test))[1],
           error(dd$dobs,approxfun(x=zz$mids,y=zz$density)(dd$test))[1],    
           error(dd$dobs,onekdeucv(dd$train,dd$test))[1],
           error(dd$dobs,bhist$bh)[1],
           error(dd$dobs,bhist$bhfp)[1], 
           #error(dd$dobs,modelavshift)[1]
           error(dd$dobs,modelbagkde)[1],
           error(dd$dobs,modelrash)[1] )
  
  #MISE[j,]=(dd$dobs-bhist$bhfp)^2 +
  j=j+1
    }
  AA=AA+A
  }
  AA=AA/M
}

cls <- makeCluster(detectCores() - 1)
vars2export <- c("BagHistfp", "Bagkde", "bropt", "dtriangle", "error", "gendata",         
                "mel", "melange", "mybreaks", "onekdeucv",
                "predict.hist", "predict.hist.x",   "rash", 
                "riskhist", "rtriangle", "simulaciones")   
clusterExport(cls, vars2export)

system.file(
res <- parLapplyLB(cls, rev(c(20, 50, 100, 200, 500, 1000, 2000)),
              function(n) simulaciones(n = n, M = 100, B = 200))
)

stopCluster(cls)

save(res, file = "resultados_mise-par2.Rdata")

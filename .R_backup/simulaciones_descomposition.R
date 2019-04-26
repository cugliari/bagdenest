
rm(list = ls())
source("functions2.R")

sesgo_par = function(n = 100, M = 5, K = 10, B = 10){

    res0 <- foreach::foreach(ll = c(1, 3, 5, 8, 11, 13, 20, 21), .combine = rbind, 
                            .export = c("gendata", "bropt", "riskhist", "mybreaks", 
                                        "broptfp", "riskfp", "ind",
                                        "predict.hist", "predict.hist.x", "onekdeucv",
                                        "kde", "BagHistfp", "melange", "mel", "rberdev",
                                        "dberdev", "rtriangle", "dtriangle", "Bagkde", "rash"),
                           .packages = "ks") %dopar% {
                             
  res = matrix(0, nrow = 1, ncol = 21)
   for (m in 1:M) {
     dd=gendata(ll,n)
     bopt=bropt(dd$train)$opt
     zz=hist(dd$train,breaks=mybreaks(dd$train,nbr=bopt),plot=F)
     h=predict.hist(zz,dd$test)
     fp=approxfun(x=zz$mids,y=zz$density)(dd$test)
     kde=onekdeucv(dd$train,dd$test)
     estim1=BagHistfp(xx=dd$train,grille=dd$test, B)$bhfp
     estim2=BagHistfp(xx=dd$train,grille=dd$test, B)$bh
     estim3=Bagkde(xx=dd$train,grille=dd$test,B)
     estim4=rash(xx=dd$train,grille=dd$test, nbr=bopt,B)
     
     finalh=0
     finalh_b=0
     
     finalfp=0
     finalfp_b=0
     
     finalbh=0
     finalbh_b=0
     
     finalbhfp=0
     finalbhfp_b=0
     
     finalkde=0
     finalkde_b=0
     
     finalbagkde=0
     finalbagkde_b=0
     
     finalrash=0
     finalrash_b=0
     
     for (k in 1:K){
       #if(k%%20 == 0) cat(k,">>")
       ddk=gendata(ll,n)
       bopt=bropt(ddk$train)$opt
       # -- Histogram
       zz=hist(ddk$train,breaks=mybreaks(ddk$train,nbr=bopt),plot=F)
       h=predict.hist(zz,dd$test)
       # -- FP
       fp=approxfun(x=zz$mids,y=zz$density)(dd$test)
       # KDE
       kde=onekdeucv(ddk$train,dd$test)
       
       # BagHist
       estim2l=BagHistfp(xx=ddk$train,grille=dd$test,B)$bh
       
       # BagFP
       estim1l=BagHistfp(xx=ddk$train,grille=dd$test,B)$bhfp
     
       #BagKde
       estim3l=Bagkde(xx=ddk$train,grille=dd$test,B)
       
       #Rash
       estim4l=rash(xx=ddk$train,grille=dd$test, nbr=bopt,B)
       
       finalh      = finalh      + (h - dd$dobs)^2
       finalfp     = finalfp     + (fp - dd$dobs)^2 
       finalkde    = finalkde    + (kde - dd$dobs)^2
       
       finalbh   = finalbh   + (estim2l - dd$dobs)^2 
       finalbhfp   = finalbhfp   + (estim1l - dd$dobs)^2 
       finalbagkde = finalbagkde + (estim3l-dd$dobs)^2
       finalrash = finalrash + (estim4l-dd$dobs)^2
       
       finalbh_b = finalbh_b + estim2l
       finalbagkde_b=finalbagkde_b + estim3l
       finalbhfp_b = finalbhfp_b + estim1l
       finalrash_b = finalrash_b + estim4l
       finalkde_b  = finalkde_b  + kde
       finalh_b    = finalh_b    + h
       finalfp_b   = finalfp_b   + fp 
       
       #print(finalbhfp)
    }
    
     finalh      = finalh/K
     finalbh    = finalbh/K
     
     finalfp     = finalfp/K
     finalbhfp   = finalbhfp/K
     finalkde    = finalkde/K
     finalbagkde= finalbagkde/K
     finalrash = finalrash/K
     
     finalh_b    = finalh_b/K
     finalbh_b = finalbh_b/K
     finalfp_b   = finalfp_b/K
     finalbhfp_b = finalbhfp_b/K
     finalkde_b  = finalkde_b/K
     finalbagkde_b  = finalbagkde_b/K
     finalrash_b = finalrash_b/K
     
    res = res + c(mean(finalh                   , na.rm = T),
                  mean((finalh_b - dd$dobs)^2   , na.rm = T),
                  mean((h-finalh_b)^2           , na.rm = T),
                  
                   
                  mean(finalfp                  , na.rm = T),
                  mean((finalfp_b - dd$dobs)^2  , na.rm = T),
                  mean((fp - finalfp_b)^2       , na.rm = T),
                 
                  mean(finalkde                 , na.rm = T),
                  mean((finalkde_b - dd$dobs)^2 , na.rm = T),
                  mean((kde - finalkde_b)^2     , na.rm = T),
                  
                  mean(finalbh                   , na.rm = T),
                  mean((finalbh_b - dd$dobs)^2   , na.rm = T),
                  mean((estim2-finalbh_b)^2      , na.rm = T),
                  
                  mean(finalbhfp                , na.rm = T),
                  mean((finalbhfp_b - dd$dobs)^2, na.rm = T),
                  mean((estim1 - finalbhfp_b)^2  , na.rm = T),
                  
                  mean(finalbagkde                 , na.rm = T),
                  mean((finalbagkde_b - dd$dobs)^2 , na.rm = T),
                  mean((estim3 - finalbagkde_b)^2  , na.rm = T),
                  
                  mean(finalrash                 , na.rm = T),
                  mean((finalrash_b - dd$dobs)^2 , na.rm = T),
                  mean((estim4 - finalrash_b)^2  , na.rm = T) )
   }
  
  #res0[j,]=res/M
   # write.table(res0,"resbb_par.txt")
   #j=j+1
  res / M
  } 

  colnames(res0) = c("errorH",   "sesgoH",   "varH", 
                     "errorFP",  "sesgoFP",  "varFP",
                     "errorKde",  "sesgoKde",  "varKde",
                     "errorBH",   "sesgoBH",   "varBH",
                     "errorBFP", "sesgoBFP", "varBFP",
                     "errorBKde", "sesgoBKde", "varBKde",
                     "errorRash", "sesgoRash", "varRash")
  rownames(res0) = c("normal", "chi2", "mix1", "bart", "triangular", 
                     "rara1", "rara2", "rar3")
  
  return(res0)
}

library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)
system.time(aa <- sesgo_par(n = 500, B = 200, K = 100, M = 100))
#system.time(aa <- sesgo_par(n = 500, B = 10, K = 5, M = 100))

save(aa, file = "resultados_descomposicion.Rdata")


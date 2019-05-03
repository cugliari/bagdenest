
library(bagdenest)

simulaciones=function(n = 100, M = 10, B = 150){

  vec = c(1,3,5,7,8,11,10,13,17,19,20,21)
  AA = matrix(0, nrow = length(vec), ncol = 7)
  #AA=matrix(0,nrow=2,ncol=7)
  colnames(AA)=c("H","FP","Kde", "BagHist", "BagFP", "BagKde", "Rash")
  rownames(AA)=c("normal","chi2 10","MIX_1","MIx_2","CLAW","Triangulaire",
                 "DENS2 rigollet","uniforme","SmooComb","DistBim","DENS1 rigollet","Mezcla uniformes")


  for(i in 1:M){
    A = matrix(NA, nrow = length(vec), ncol = 7)
    #A=matrix(NA,nrow=2,ncol=7)
    for(ll in 1:length(vec)){
      #print(ll)
      dd = gendata(vec[ll],n)
      bopt = bropt(dd$train)$opt  # needed for Hist & RASH
      zz = hist(dd$train,breaks=mybreaks(dd$train,nbr=bopt),plot=F)

      bhist = BagHistfp(xx=dd$train, grille=dd$test, B)
      modelrash = rash(dd$train, grille = dd$test, nbr = bopt, B)
      #modelavshift(dd$train,dd$test,nbr=bopt,M=B)
      modelbagkde <- Bagkde(xx = dd$train, grille = dd$test, B)

      A[ll,]=c(error(dd$dobs,predicthist(zz,dd$test))[1],
               error(dd$dobs,approxfun(x=zz$mids,y=zz$density)(dd$test))[1],
               error(dd$dobs,onekdeucv(dd$train,dd$test))[1],
               error(dd$dobs,bhist$bh)[1],
               error(dd$dobs,bhist$bhfp)[1],
               #error(dd$dobs,modelavshift)[1]
               error(dd$dobs,modelbagkde)[1],
               error(dd$dobs,modelrash)[1])
    }
    AA=AA+A6
  }
  AA=AA/M
  print(AA)
}


simulaciones(n = 100, M = 2, B = 15)

######################################################################################


system.file(
  res <- lapply(c(50, 100, 200, 500, 1000),
                function(n) simulaciones(n = n, M = 2, B = 15))
)
res

######################################################################################

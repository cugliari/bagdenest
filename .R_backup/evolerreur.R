
rm(list=ls())
setwd("~/Dropbox/Trabajo/BD_pour_hist/codes")
source("functions2.r")


#####################
#Evol. error

evol.error=function(modele=1,n=100,B=200,M=5){
    
    A=matrix(NA,nrow=M,ncol=B)
    C=matrix(NA,nrow=M,ncol=B)
    D=matrix(NA,nrow=M,ncol=B)
    E=matrix(NA,nrow=M,ncol=B)
    for(kk in 1:M){
      print(kk)
      dd=gendata(modele,n)
      bopt=bropt(dd$train)$opt
      
      bhfp=BagHistfp.err(xx=dd$train,grille=dd$test,B=B,dobs=dd$dobs)
      A[kk,]=bhfp$erreurbh[,1]
      C[kk,]=bhfp$erreurbhfp[,1]
      D[kk,]=BagKDE.err(dd$train, dd$test,B,dd$dobs)$erreur[,1]
      E[kk,]=rash.err(dd$train,dd$test,nbr=bopt, B,dobs=dd$dobs)$erreur[,1]
    }
list(A=A,C=C,D=D,E=E)
    }

system.file(
  res <- lapply(c(1, 3, 5, 8, 11,13,20,21),
                function(modele) evol.error(modele = modele, n=500, M = 100, B = 200))
)

save(res, file = "res_evolerror_new.Rdata")



rm(list = ls())
source("functions2.R")
source("simulaciones-tubos.R")
df <- read.csv("datos2plot_tubos.csv", sep = ";", dec = ",", row.names = 1)

#xtabs(coverage ~ density + estimator, data = df)

## Plotear tubos
ll <- gendata(7, 500)
tt <- tube_KDE(ll)

b1 <- baseboot(ll$test, tt$punctual, tt$bunch)  #booststrapkde
#b2 <- cteboot(ll$test, tt$punctual, tt$bunch)
b3 <- tube_KDEsm(ll)$tube   #smdensity

matplot(ll$test, t(tt$bunch), type = 'l', col = "grey", lty = 1, 
        ylim = c(0, max(ll$dobs) *1.5),xlab="x", ylab="y")
lines(ll$test, tt$punctual, col = 2, lwd = 2)
lines(ll$test, ll$dobs, lwd = 2)

matlines(ll$test, t(b1), col = 4, lty = 1, lwd = 2)  #boostrapkde
#matlines(ll$test, t(b2), col = "violet", lty = 1, lwd = 2)
matlines(ll$test, t(b3), col = "darkgreen", lty = 1, lwd = 2) #smdensity
legend(x="top",inset=c(0, 0),cex=0.7,horiz=T, c("Target","Estimation","Boostrap","Sm.Density"),col=c("black","red","blue","darkgreen"),lty=c(1,1,1,1))

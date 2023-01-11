
###JUNE 2022
###re-doing figures

#...................................................................................................

#not sure what packages I need need, so here's all of the ones I use ??...

library(ggplot2)
library(lme4)
library(lmerTest)
library(dplyr)
library(performance)
library(corrplot)
library("PerformanceAnalytics")
library(ggrepel)
library(reshape2)
library(gridExtra)

##BUT I THINK THESE ARE THE IMPORTANT PCKGS: 

library(geoR)

#or

library(gstat)

#bring in all the data

vario.data <- read.csv("C:/Users/robin/Documents/RESEARCH/PAPERS/Spatial Autocorrelation of Reproduction/Amphiprion percula/Updated_datasheet_AperculaPNG_20220210.csv", na.strings="")

#and then I subset this data to make varigorams for each metric (see below)

#...................................................................................................


#get the residuals from linear regression bwtween metrics 

lm.FSL.AEA<-lm(vario.data$Female_SL~vario.data$Anemone_Ellipse)
vario.data$res.FSL.AEA<-rstandard(lm.FSL.AEA)

lm.CCF.FSL<-lm(vario.data$Clutch1_CountF~vario.data$Female_SL)
vario.data$res.CCF.FSL<-rstandard(lm.CCF.FSL)

#...................................................................................................



#VARIOGRAM SET UP

par(mar=c(5,5,4,1)+.1)

## Binding necessary columns
V.DF <- vario.data[,c("Lat", "Lon", "Anemone_Ellipse", "Female_SL", "Clutch1_CountF")]

## Binding necessary columns
V.res.DF <- vario.data[,c("Lat", "Lon", "res.FSL.AEA", "res.CCF.FSL")]


#...................................................................................................


#ANEMONE SIZE

## Assigning Anemone Ellipse Area AEA data as geodata (column 3)
V.AEA <- as.geodata(V.DF, coords.col=1:2, data.col= 3)

summary(V.AEA)

## Building and plotting variogram
Var.V.AEA <- variog(V.AEA)

plot(Var.V.AEA, 
     pch=19, cex=1.5,
     xlab="Distance", ylab= "Semivariance", 
     cex.lab=2.5, cex.axis=1.5)

# 1 dec deg = 111 km (aprrox at the equator) 

#re-write distance labels... 

plot(Var.V.AEA, 
     pch=19, cex=1.5,
     xlab="Distance (km)", ylab= "Semivariance", 
     cex.lab=2.5, cex.axis=1.5, xaxt = "n")
axis(1, at = c(0.00, 0.01, 0.02, 0.03, 0.04),
     labels = c("0", "1", "2", "3", "4"), 
     cex.axis=2)

#Fitting variogram 
AEA.ini.vals <- expand.grid(seq(0, 10000000, by = 2000000), seq(0, 0.05, by = 0.01))
AEA.ols <- variofit(Var.V.AEA, ini = AEA.ini.vals, fix.nug = TRUE, wei = "equal")
AEA.ols.nofix <- variofit(Var.V.AEA, ini = AEA.ini.vals, fix.nug = FALSE, wei = "equal")
lines(AEA.ols, lwd=2)
lines(AEA.ols.nofix, lty=2, lwd=2)

summary(AEA.ols) 
# practical Range = 0.002657386
# sigma sq = 6.485328e+06
# Sum of Squares = 4.105845e+13 ###Better model 

summary(AEA.ols.nofix) 
# practical Range = 0.0009429926
# nugget = 1.450225e-11 
# sigma sq = 6.000000e+06
# Sum of Squares = 4.475234e+13 

#proportion of the variance that is spatially structured?
#(sill - nugget)/max semivariance
summary(Var.V.AEA$v)
#max semivariance = 9647413
(6.485328e+06 - 0)/9647413
# = 0.6722349

#.................................................................

#trying to change bins 
#(still using unprojected dec deg coordinate data though...) 

#to compare: 

Var.V.AEA <- variog(V.AEA)
plot(Var.V.AEA)
AEA.ols <- variofit(Var.V.AEA, ini = AEA.ini.vals, fix.nug = TRUE, wei = "equal")
AEA.ols.nofix <- variofit(Var.V.AEA, ini = AEA.ini.vals, fix.nug = FALSE, wei = "equal")
lines(AEA.ols, lwd=2)
lines(AEA.ols.nofix, lty=2, lwd=2)
summary(AEA.ols)

distance.bin <- seq(0, 0.045, by = 0.0005) 
  
Var.V.AEA.b <- variog(V.AEA, uvec=distance.bin, breaks=distance.bin )
plot(Var.V.AEA.b)

distance.bin <- seq(0, 0.045, by = 0.001) 

Var.V.AEA.b <- variog(V.AEA, uvec=distance.bin, breaks=distance.bin )
plot(Var.V.AEA.b)

distance.bin <- seq(0, 0.045, by = 0.005) 

Var.V.AEA.b <- variog(V.AEA, uvec=distance.bin, breaks=distance.bin )
plot(Var.V.AEA.b)

distance.bin <- seq(0, 0.045, by = 0.0025) 

Var.V.AEA.b <- variog(V.AEA, uvec=distance.bin, breaks=distance.bin )
plot(Var.V.AEA.b)

#ughhh i d k.... by 0.005?

distance.bin <- seq(0, 0.045, by = 0.005) 

Var.V.AEA.b <- variog(V.AEA, uvec=distance.bin, breaks=distance.bin )
plot(Var.V.AEA.b)

#Fitting variogram

AEA.ini.vals <- expand.grid(seq(0, 10000000, by = 2000000), seq(0, 0.05, by = 0.01))

#og lines from default breaks 
AEA.ols <- variofit(Var.V.AEA, ini = AEA.ini.vals, fix.nug = TRUE, wei = "equal")
AEA.ols.nofix <- variofit(Var.V.AEA, ini = AEA.ini.vals, fix.nug = FALSE, wei = "equal")
lines(AEA.ols, lwd=2)
lines(AEA.ols.nofix, lty=2, lwd=2)

#new lines based on breaks by 0.005 (9 bins?)
AEA.b.ols <- variofit(Var.V.AEA.b, ini = AEA.ini.vals, fix.nug = TRUE, wei = "equal")
AEA.b.ols.nofix <- variofit(Var.V.AEA.b, ini = AEA.ini.vals, fix.nug = FALSE, wei = "equal")
lines(AEA.b.ols, lwd=2, col="red")
lines(AEA.b.ols.nofix, lty=2, lwd=2, col="red")

#different!!! ughhhh

#og line
summary(AEA.ols) 
# practical Range = 0.002657386
# sigma sq = 6.485328e+06
# Sum of Squares = 4.105845e+13 

#new line
summary(AEA.b.ols) 
# practical Range = 0.004578714
# sigma sq = 6.000000e+06
# Sum of Squares = 2.266013e+13 ###Better model

####################################################################################

#okay more importantly- the residuals 

#RESIDUALS OF FEMALE SIZE from FSL ~ AEA 

## Assigning Residual FSL data as geodata (column 3)
V.resFSL <- as.geodata(V.res.DF, coords.col=1:2, data.col= 3)

## Building and plotting variogram
Var.V.resFSL <- variog(V.resFSL)

plot(Var.V.resFSL, 
     pch=19, cex=1.5,
     xlab="Distance (km)", ylab= "Semivariance", 
     cex.lab=2.5, cex.axis=1.5)

#changing bins...

#to compare 

Var.resFSL <- variog(V.resFSL)
plot(Var.resFSL)
resFSL.ini.vals <- expand.grid(seq(0, 1.4, by = 0.2), seq(0, 0.05, by = 0.01))
resFSL.ols <- variofit(Var.resFSL, ini = resFSL.ini.vals, fix.nug = TRUE, wei = "equal")
resFSL.ols.nofix <- variofit(Var.resFSL, ini = resFSL.ini.vals, fix.nug = FALSE, wei = "equal")
lines(resFSL.ols, lwd=2)
lines(resFSL.ols.nofix, lty=2, lwd=2)
summary(resFSL.ols)

distance.bin <- seq(0, 0.045, by = 0.0005) 

Var.resFSL.b <- variog(V.resFSL, uvec=distance.bin, breaks=distance.bin )
plot(Var.resFSL.b)

distance.bin <- seq(0, 0.045, by = 0.001) 

Var.resFSL.b <- variog(V.resFSL, uvec=distance.bin, breaks=distance.bin )
plot(Var.resFSL.b)

distance.bin <- seq(0, 0.045, by = 0.005) 

Var.resFSL.b <- variog(V.resFSL, uvec=distance.bin, breaks=distance.bin )
plot(Var.resFSL.b)

#Fitting variogram #og lines
resFSL.ini.vals <- expand.grid(seq(0, 1.4, by = 0.2), seq(0, 0.05, by = 0.01))
resFSL.ols <- variofit(Var.V.resFSL, ini = resFSL.ini.vals, fix.nug = TRUE, wei = "equal")
resFSL.ols.nofix <- variofit(Var.V.resFSL, ini = resFSL.ini.vals, fix.nug = FALSE, wei = "equal")
lines(resFSL.ols, lwd=2)
lines(resFSL.ols.nofix, lty=2, lwd=2)

#new lines based on breaks at 0.005
resFSL.ini.vals <- expand.grid(seq(0, 2, by = 0.5), seq(0, 0.05, by = 0.01))
resFSL.b.ols <- variofit(Var.resFSL.b, ini = resFSL.ini.vals, fix.nug = TRUE, wei = "equal")
resFSL.b.ols.nofix <- variofit(Var.resFSL.b, ini = resFSL.ini.vals, fix.nug = FALSE, wei = "equal")
lines(resFSL.b.ols, lwd=2, col="red")
lines(resFSL.b.ols.nofix, lty=2, lwd=2, col="red")

#comapring models 

summary(resFSL.ols) 
# practical Range = 0.003311736
# sigma sq = 0.999946586
# Sum of Squares = 1.478769 

summary(resFSL.b.ols) 
# practical Range = 0.004474421
# sigma sq = 0.986413161
# Sum of Squares = 0.5776653 #better model !!!

0.004474421 * 111
# = 0.49666 km 
# or 497 meters   


#############


#RESIDUAL CLUTCH SIZE 

#RESIDUALS OF FEMALE SIZE from FSL ~ AEA 

## Assigning Residual FSL data as geodata (column 3)
V.resCS <- as.geodata(V.res.DF, coords.col=1:2, data.col= 4)

## Building and plotting variogram
Var.V.resCS <- variog(V.resCS)

plot(Var.V.resCS, 
     pch=19, cex=1.5,
     xlab="Distance (km)", ylab= "Semivariance", 
     cex.lab=2.5, cex.axis=1.5)

#changing bins...

#to comapre
Var.resCS <- variog(V.resCS)
plot(Var.resCS)

distance.bin <- seq(0, 0.045, by = 0.0005) 

Var.resCS.b <- variog(V.resCS, uvec=distance.bin, breaks=distance.bin )
plot(Var.resCS.b)

distance.bin <- seq(0, 0.045, by = 0.001) 

Var.resCS.b <- variog(V.resCS, uvec=distance.bin, breaks=distance.bin )
plot(Var.resCS.b)

distance.bin <- seq(0, 0.045, by = 0.005) 

Var.resCS.b <- variog(V.resCS, uvec=distance.bin, breaks=distance.bin )
plot(Var.resCS.b)

#Fitting variogram #og lines
resCS.ini.vals <- expand.grid(seq(0, 1.4, by = 0.2), seq(0, 0.05, by = 0.01))
resCS.ols <- variofit(Var.V.resCS, ini = resCS.ini.vals, fix.nug = TRUE, wei = "equal")
resCS.ols.nofix <- variofit(Var.V.resCS, ini = resCS.ini.vals, fix.nug = FALSE, wei = "equal")
lines(resCS.ols, lwd=2)
lines(resCS.ols.nofix, lty=2, lwd=2)

#new lines based on breaks at 0.005
resCS.ini.vals <- expand.grid(seq(0, 1.2, by = 0.2), seq(0, 0.05, by = 0.01))
resCS.b.ols <- variofit(Var.resCS.b, ini = resCS.ini.vals, fix.nug = TRUE, wei = "equal")
resCS.b.ols.nofix <- variofit(Var.resCS.b, ini = resCS.ini.vals, fix.nug = FALSE, wei = "equal")
lines(resCS.b.ols, lwd=2, col="red")
lines(resCS.b.ols.nofix, lty=2, lwd=2, col="red")

 #comapring models 

summary(resCS.ols) 
# practical Range = 0.0001159668
# sigma sq = 0.868026
# Sum of Squares = 1.38116 

summary(resCS.b.ols) 
# practical Range = 0.0001159668
# sigma sq = 0.9007021
# Sum of Squares = 0.5883112 #better model !!!

0.0001159668 * 111
# = 0.0128 km 
# or 1.28 meters   



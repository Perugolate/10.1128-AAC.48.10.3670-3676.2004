###############################################################################
#                                                                             #
# THIS DEMO SHOWS HOW TO FIT THE PHARMACODYNAMIC FUNCTION TO                  #
# TIME-KILL DATA IN THE STATISTICAL COMPUTING LANGUAGE R                      #
#                                                                             #
# Copyright (C) 2003 Roland R. Regoes (rregoes@emory.edu)                     #
#                                                                             #
# This demo is free software; you can redistribute it and/or                  #
# modify it under the terms of the GNU General Public License                 #
# as published by the Free Software Foundation; either version 2              #
# of the License, or (at your option) any later version.                      #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program; if not, write to the Free Software                 #
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA. #
#                                                                             #
###############################################################################

###############################################################################
# FIRST READ THE SAMPLE DATA FROM OUR WEB-SITE:
# (CHANGE THE read.table COMMAND IF YOU WANT TO READ YOUR OWN DATASET)

#data<-read.table(
#      "http://userwww.service.emory.edu/~rregoes/pharmacodynamics/ampic.data",
#      header=TRUE)
data <- read.table("ampic.data", header = TRUE)
attach(data)
antibiotic.conc<-factor(antibiotic.conc)


###############################################################################
# LET'S LOOK AT THE DATA:

par(mfrow=c(1,2),mai=par("mai")+c(2,0,1,0))->op # define useful plot parameters
plot(c(0,5),c(1,1e10),type="n",log="y",
     xlab="Hours",ylab="CFU/ml",main="Time-kill curves") 
for(idose in 1:nlevels(antibiotic.conc)){
  dose<-as.numeric(levels(antibiotic.conc)[idose])
  points(hours[antibiotic.conc==dose],CFU.per.ml[antibiotic.conc==dose],
         col=idose,pch=idose)
  lines(hours[antibiotic.conc==dose],CFU.per.ml[antibiotic.conc==dose],
        col=idose)
}


###############################################################################
# NOW LET'S CALCULATE THE BACTERIAL NET GROWTH RATES FOR
# EACH LEVEL OF ANTIBIOTIC CONCENTRATION

tmin<-0 # lower bound of the time interval over which to calculate the decline
tmax<-1 # upper bound of the time interval over which to calculate the decline

# define a vector for the decline rates:
slopes<-vector("numeric",length=nlevels(antibiotic.conc))


# NET GROWTH RATES ARE CALCULATED BY LINEARLY REGRESSING 
# BACTERIAL CONCENTRATIONS OVER THE CHOSEN TIME INTERVAL
# THE DECLINE RATE IS THE COEFFICIENT OF REGRESSION

for(idose in 1:nlevels(antibiotic.conc)){
  dose<-as.numeric(levels(antibiotic.conc)[idose])
  data2regress<-data.frame(time=hours[hours>=tmin & hours<=tmax &
                                      antibiotic.conc==dose],
                           log.bacteria.rel=log(CFU.per.ml[hours>=tmin &
                             hours<=tmax & antibiotic.conc==dose]/
                             CFU.per.ml[hours==tmin&antibiotic.conc==dose],10))
  slopes[idose]<-coef(lm(log.bacteria.rel~time-1,data2regress))
}


# PLOT THE DECLINES RATES
# (WE PLOT A SHIFT THE CONCENTRATIONS SLIGHTLY,
# BY A TENTH OF THE LOWEST ANTIBIOTIC CONCENTRATION,
# TO ALLOW CONCENTRATION ZERO TO BE PLOTTED ON THE LOG-SCALE)

plot(as.numeric(levels(antibiotic.conc))+
     as.numeric(levels(antibiotic.conc))[2]/10,
     slopes,
     log="x",
     xlab="Antibiotic concentration",ylab="Net bacterial growth rate",
     main="Pharmacodynamic relationship",
     xlim=c(as.numeric(levels(antibiotic.conc))[2]/10,
       2*max(as.numeric(levels(antibiotic.conc)))),
     ylim=c(-4,2))


###############################################################################
# THE PHARMACODYNAMIC FUNCTION:

psi<-function(psimax,psimin,kappa,MIC,dose) {
  psimax-(psimax-psimin)*(dose/MIC)^kappa/((dose/MIC)^kappa-psimin/psimax)
}

###############################################################################
# FITTING THE PHARMACODYNAMIC FUNCTION:

library(nls) # load the nls-package

# GUESS AND DEFINE START PARAMETERS FOR THE FITTING ROUTINE
# (WHEN FITTING YOUR OWN DATA, ADAPT THESE GUESSES TO YOUR DATA SET)
MIC.guess<-3
psimin.guess<--10

# start parameters of the fitting routine
startpar<-c(psimax=1,psimin=psimin.guess,kappa=1,MIC=MIC.guess) 

data2fit<-data.frame(dose=as.numeric(levels(antibiotic.conc)),growth=slopes)


# FITTING PROCEDURE

fittedmodel<-nls(growth~psi(psimax,psimin,kappa,MIC,dose),
                 data=data2fit,
                 start=startpar)

print(summary(fittedmodel)) # print a summary of the fit


# PLOT THE FITTED FUNCTION ONTO THE PREVIOUS PLOT
# (AS BEFORE, WE PLOT A CURVE SHIFTED SLIGHTLY,
# BY A TENTH OF THE LOWEST ANTIBIOTIC CONCENTRATION)

curve(psi(coef(fittedmodel)[[1]],
          coef(fittedmodel)[[2]],
          coef(fittedmodel)[[3]],
          coef(fittedmodel)[[4]],
          x-as.numeric(levels(antibiotic.conc))[2]/10),
      add=TRUE,col=2)


par(op) # reset plotting parameters
rm(op,idose,dose,pu) # remove dummy variables
detach(data)


# TO CLEAN UP YOUR R WORKSPACE OUTCOMMENT THE FOLLOWING LINES
#rm(data,antibiotic.conc,tmin,tmax,data2regress,slopes,psi,MIC.guess,
#psimin.guess,startpar,data2fit,fittedmodel)

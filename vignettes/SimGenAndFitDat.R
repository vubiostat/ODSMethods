setwd("/Users/schildjs/rsch/Git/ODS4LDA/tests/")
library(lme4)
library(mvtnorm)
library(mitools)
library(ODS4LDA)
source("SimGenDatFns.R")
source("ImputationFns.R")

#source(paste(path, "Functions4.R", sep=""))

# set.seed(ITER)
#
#
# start.tm <- date()
#
# RunMods <- function(params,
#                     NsPerStratumUniv,
#                     N,
#                     prev.grp,
#                     n.imp,
#                     p.central,
#                     conf.param,
#                     ni){

    params=c(75, -1, -.5, -2, -.5, log(9),log(1.25),0,log(3.5))
    NsPerStratumUniv=c(150,100,150)
    NsPerStratumBiv=c(100,300)

    N=2000
    prev.grp=0.3
    n.imp=50
    p.central=0.8
    conf.param = c(-0.25, 0.5)
    ni = c(4, 6)


    inits            <- params
    NsRand           <- sum(NsPerStratumUniv)
    ##################################################
    ## Generate cutoffs using a population of 25000
    dat.tmp   <- GenerateX(N=25000, n=ni[2], prev.grp=prev.grp, c.parm=conf.param)
    dat.tmp$y <- GenerateY(X=cbind(1, dat.tmp$time, dat.tmp$grp, dat.tmp$conf, dat.tmp$time*dat.tmp$grp), Z=cbind(1, dat.tmp$time), id=dat.tmp$id,
                           beta=inits[1:5], sig.b0 = exp(inits[6]), sig.b1 =exp(inits[7]), rho = inits[8], sig.e = exp(inits[9]), RanefDist="Gaussian", ErrorDist="Gaussian")

    ## Calculate cutoffs that define central regions for intercept, slope, and bivariate sampling
    cutoffs   <- est.cutoffs(Y=dat.tmp$y, time=dat.tmp$time, id=dat.tmp$id, PropInCentralRegion=p.central)

    ## Generate a full cohort dataframe
    dat      <- GenerateX(N=N, n=ni[2], prev.grp=prev.grp, c.parm=conf.param)
    dat$y    <- GenerateY(X=cbind(1, dat$time, dat$grp, dat$conf, dat$time*dat$grp), Z=cbind(1, dat$time), id=dat$id,
                          beta=inits[1:5], sig.b0 = exp(inits[6]), sig.b1 =exp(inits[7]), rho = inits[8], sig.e = exp(inits[9]), RanefDist="Gaussian", ErrorDist="Gaussian")
    dat      <- dat[order(dat$id, dat$time),]
    droptime <- rep(sample(seq(ni[1], ni[2]), N, replace=TRUE), each=ni[2]) - 1
    dat      <- dat[dat$time<=droptime,]

    ## Calculate subject specific intercepts and slopes
    IntSlps <- CalcSSIntSlp( Y=dat$y, time=dat$time, id=dat$id)
    dat$Int <- IntSlps[[1]]
    dat$Slp <- IntSlps[[2]]

    ## identify stratum membership
    dat$StratInt <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="intercept", cutpoints=cutoffs$IntCutUniv, Int=dat$Int, Slp=dat$Slp)
    dat$StratSlp <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="slope",     cutpoints=cutoffs$SlpCutUniv, Int=dat$Int, Slp=dat$Slp)
    dat$StratBiv <- identify.stratum(Y=dat$y, time=dat$time, id=dat$id, w.function="bivar",     cutpoints=c(cutoffs$IntCutBiv, cutoffs$SlpCutBiv), Int=dat$Int, Slp=dat$Slp)

    ## identify the sample along with individual sampling probs and stratum sampling probs
    SampledRan <- random.sampling(id.long=dat$id, n=NsRand)
    SampledInt <- ods.sampling(id.long=dat$id, stratum.long=dat$StratInt, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumUniv)
    SampledSlp <- ods.sampling(id.long=dat$id, stratum.long=dat$StratSlp, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumUniv)
    SampledBiv <- ods.sampling(id.long=dat$id, stratum.long=dat$StratBiv, SamplingStrategy="IndepODS", NsPerStratum=NsPerStratumBiv)

    ## add who will be sampled for each design to dat
    dat$SampledRan <- SampledRan
    dat$SampledInt <- SampledInt[["Sampled"]]
    dat$SampledSlp <- SampledSlp[["Sampled"]]
    dat$SampledMix1 <- dat$SampledInt*(dat$id <= N/3) + dat$SampledSlp*(dat$id > N/3)
    dat$SampledMix2 <- dat$SampledInt*(dat$id <= 2*N/3) + dat$SampledSlp*(dat$id > 2*N/3)
    dat$SampledBiv <-  SampledBiv[["Sampled"]]

    ## Added subject specific sampling weights for each design if doing a weighted likelihood analysis
    dat$WeightsInt <- 1/SampledInt[["SampProbi"]]
    dat$WeightsSlp <- 1/SampledSlp[["SampProbi"]]
    dat$WeightsMix1 <- 1/(SampledInt[["SampProbi"]]*(dat$id <= N/3) + SampledSlp[["SampProbi"]]*(dat$id > N/3))
    dat$WeightsMix2 <- 1/(SampledInt[["SampProbi"]]*(dat$id <= 2*N/3) + SampledSlp[["SampProbi"]]*(dat$id > 2*N/3))
    dat$WeightsBiv <- 1/SampledBiv[["SampProbi"]]
    ## for acml.lmem (not acml.lmem2)
    dat$SampProbiInt <- SampledInt[["SampProbi"]]
    dat$SampProbiSlp <- SampledSlp[["SampProbi"]]
    dat$SampProbiBiv <- SampledBiv[["SampProbi"]]

    ## This is used if we are not doing sampling weights
    dat$NoWeighting <-1

    dat$IntProbLow  <- SampledInt[["SampProbs"]][1]
    dat$IntProbMid  <- SampledInt[["SampProbs"]][2]
    dat$IntProbHigh <- SampledInt[["SampProbs"]][3]
    dat$SlpProbLow  <- SampledSlp[["SampProbs"]][1]
    dat$SlpProbMid  <- SampledSlp[["SampProbs"]][2]
    dat$SlpProbHigh <- SampledSlp[["SampProbs"]][3]
    dat$Mix1ProbLow  <- SampledInt[["SampProbs"]][1]*(dat$id <= N/3) + SampledSlp[["SampProbs"]][1]*(dat$id > N/3)
    dat$Mix1ProbMid  <- SampledInt[["SampProbs"]][2]*(dat$id <= N/3) + SampledSlp[["SampProbs"]][2]*(dat$id > N/3)
    dat$Mix1ProbHigh <- SampledInt[["SampProbs"]][3]*(dat$id <= N/3) + SampledSlp[["SampProbs"]][3]*(dat$id > N/3)
    dat$Mix2ProbLow  <- SampledInt[["SampProbs"]][1]*(dat$id <= 2*N/3) + SampledSlp[["SampProbs"]][1]*(dat$id > 2*N/3)
    dat$Mix2ProbMid  <- SampledInt[["SampProbs"]][2]*(dat$id <= 2*N/3) + SampledSlp[["SampProbs"]][2]*(dat$id > 2*N/3)
    dat$Mix2ProbHigh <- SampledInt[["SampProbs"]][3]*(dat$id <= 2*N/3) + SampledSlp[["SampProbs"]][3]*(dat$id > 2*N/3)
    dat$BivProbInner <- SampledBiv[["SampProbs"]][1]
    dat$BivProbOuter <- SampledBiv[["SampProbs"]][2]

    dat$IntCutoff1 <- cutoffs$IntCutUniv[1]
    dat$IntCutoff2 <- cutoffs$IntCutUniv[2]
    dat$SlpCutoff1 <- cutoffs$SlpCutUniv[1]
    dat$SlpCutoff2 <- cutoffs$SlpCutUniv[2]
    dat$Mix1Cutoff1 <- dat$IntCutoff1*(dat$id <= N/3) + dat$SlpCutoff1*(dat$id > N/3)
    dat$Mix1Cutoff2 <- dat$IntCutoff2*(dat$id <= N/3) + dat$SlpCutoff2*(dat$id > N/3)
    dat$Mix2Cutoff1 <- dat$IntCutoff1*(dat$id <= 2*N/3) + dat$SlpCutoff1*(dat$id > 2*N/3)
    dat$Mix2Cutoff2 <- dat$IntCutoff2*(dat$id <= 2*N/3) + dat$SlpCutoff2*(dat$id > 2*N/3)
    dat$BivCutoffInt1 <- cutoffs$IntCutBiv[1]
    dat$BivCutoffInt2 <- cutoffs$IntCutBiv[2]
    dat$BivCutoffSlp1 <- cutoffs$SlpCutBiv[1]
    dat$BivCutoffSlp2 <- cutoffs$SlpCutBiv[2]

    ## Define w.function within dat for each design
    dat$Int.w <- "intercept"
    dat$Slp.w <- "slope"
    dat$Mix1.w <- ifelse(dat$id <= N/3, dat$Int.w, dat$Slp.w)
    dat$Mix2.w <- ifelse(dat$id <= 2*N/3, dat$Int.w, dat$Slp.w)
    dat$Biv.w <- "bivar"

    ## Stratum Sampling probabilities
    SampProbRan <- c(1,1,1)
    SampProbInt <- SampledInt[["SampProbs"]]
    SampProbSlp <- SampledSlp[["SampProbs"]]
    SampProbBiv <- SampledBiv[["SampProbs"]]

    ## Datasets for sampled subjects
    datRan  <- dat[dat$SampledRan==1,]
    datInt  <- dat[dat$SampledInt==1,]
    datSlp  <- dat[dat$SampledSlp==1,]
    datMix1 <- dat[dat$SampledMix1==1,]
    datMix2 <- dat[dat$SampledMix2==1,]
    datBiv  <- dat[dat$SampledBiv==1,]

    cutpointsRan=cbind(datRan$IntCutoff1, datRan$IntCutoff2)  ## just need these numbers for the function, not used
    SampProbRan=matrix(1, ncol=3, nrow=length(datRan[,1]))
    w.functionRan=datRan$Int.w                                ## just need these numbers for the function, not used

    cutpointsInt=cbind(datInt$IntCutoff1, datInt$IntCutoff2)
    SampProbInt=cbind(datInt$IntProbLow, datInt$IntProbMid, datInt$IntProbHigh)
    w.functionInt=datInt$Int.w

    cutpointsSlp=cbind(datSlp$SlpCutoff1, datSlp$SlpCutoff2)
    SampProbSlp=cbind(datSlp$SlpProbLow, datSlp$SlpProbMid, datSlp$SlpProbHigh)
    w.functionSlp=datSlp$Slp.w

    cutpointsMix1=cbind(datMix1$Mix1Cutoff1, datMix1$Mix1Cutoff2)
    SampProbMix1=cbind(datMix1$Mix1ProbLow, datMix1$Mix1ProbMid, datMix1$Mix1ProbHigh)
    w.functionMix1=datMix1$Mix1.w

    cutpointsMix2=cbind(datMix2$Mix2Cutoff1, datMix2$Mix2Cutoff2)
    SampProbMix2=cbind(datMix2$Mix2ProbLow, datMix2$Mix2ProbMid, datMix2$Mix2ProbHigh)
    w.functionMix2=datMix2$Mix2.w

    cutpointsBiv=cbind(datBiv$BivCutoffInt1, datBiv$BivCutoffInt2, datBiv$BivCutoffSlp1, datBiv$BivCutoffSlp2)
    SampProbBiv=cbind(datBiv$BivProbInner, datBiv$BivProbOuter)
    w.functionBiv=datBiv$Biv.w

    ## Datasets for unsampled subjects
    datNotRan  <- dat[dat$SampledRan==0,]
    datNotInt  <- dat[dat$SampledInt==0,]
    datNotSlp  <- dat[dat$SampledSlp==0,]
    datNotMix1 <- dat[dat$SampledMix1==0,]
    datNotMix2 <- dat[dat$SampledMix2==0,]
    datNotBiv  <- dat[dat$SampledBiv==0,]

    cutpointsNotRan=cbind(datNotRan$IntCutoff1, datNotRan$IntCutoff2)  ## just need these numbers for the function, not used
    SampProbNotRan=matrix(1, ncol=3, nrow=length(datNotRan[,1]))
    w.functionNotRan=datNotRan$Int.w                                ## just need these numbers for the function, not used

    cutpointsNotInt=cbind(datNotInt$IntCutoff1, datNotInt$IntCutoff2)
    SampProbNotInt=cbind(datNotInt$IntProbLow, datNotInt$IntProbMid, datNotInt$IntProbHigh)
    w.functionNotInt=datNotInt$Int.w

    cutpointsNotSlp=cbind(datNotSlp$SlpCutoff1, datNotSlp$SlpCutoff2)
    SampProbNotSlp=cbind(datNotSlp$SlpProbLow, datNotSlp$SlpProbMid, datNotSlp$SlpProbHigh)
    w.functionNotSlp=datNotSlp$Slp.w

    cutpointsNotMix1=cbind(datNotMix1$Mix1Cutoff1, datNotMix1$Mix1Cutoff2)
    SampProbNotMix1=cbind(datNotMix1$Mix1ProbLow, datNotMix1$Mix1ProbMid, datNotMix1$Mix1ProbHigh)
    w.functionNotMix1=datNotMix1$Mix1.w

    cutpointsNotMix2=cbind(datNotMix2$Mix2Cutoff1, datNotMix2$Mix2Cutoff2)
    SampProbNotMix2=cbind(datNotMix2$Mix2ProbLow, datNotMix2$Mix2ProbMid, datNotMix2$Mix2ProbHigh)
    w.functionNotMix2=datNotMix2$Mix2.w

    cutpointsNotBiv=cbind(datNotBiv$BivCutoffInt1, datNotBiv$BivCutoffInt2, datNotBiv$BivCutoffSlp1, datNotBiv$BivCutoffSlp2)
    SampProbNotBiv=cbind(datNotBiv$BivProbInner, datNotBiv$BivProbOuter)
    w.functionNotBiv=datNotBiv$Biv.w

    CalcDiffs <- function(x,y) {
        out <- list(x[["coefficients"]]-y[["coefficients"]],
                    x[["covariance"]]-y[["covariance"]],
                    x[["robcov"]]-y[["robcov"]],
                    x[["loglik"]]-y[["loglik"]])
        out}


    simplify2array(Fit.ran, Fit.ran2)
    ####### Complete case analyses
    print("Univariate ACML, WL, and ML Analyses")
    print(date())
    Fit.ran     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datRan, InitVals=inits, ProfileCol=NA,
                              cutpoints=cutpointsRan, SampProb=SampProbRan, Weights=NoWeighting, w.function=w.functionRan)
    print(date())

    datRan$SampProbiNoWt <- 1
    Fit.ran2     <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datRan, InitVals=inits, ProfileCol=NA,
                              cutpoints=cutpointsRan[1,], SampProb=SampProbRan[1,], SampProbiWL=SampProbiNoWt, w.function=w.functionRan[1])
    print(date())

    CalcDiffs(Fit.ran,Fit.ran2)

    Fit.int     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datInt, InitVals=inits, ProfileCol=NA,
                              cutpoints=cutpointsInt, SampProb=SampProbInt, Weights=NoWeighting, w.function=w.functionInt)
    print(date())
    datInt$SampProbiNoWt <- 1
    Fit.int2     <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datInt, InitVals=inits, ProfileCol=NA,
                              cutpoints=cutpointsInt[1,], SampProb=SampProbInt[1,], SampProbiWL=SampProbiNoWt, w.function=w.functionInt[1])
    print(date())

    CalcDiffs(Fit.int,Fit.int2)

    Fit.int.wl     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datInt, InitVals=inits, ProfileCol=NA,
                                 cutpoints=cutpointsInt, SampProb=matrix(1, nrow=nrow(datInt), ncol=3), Weights=WeightsInt, w.function=w.functionInt)
    print(date())

    Fit.int.wl2     <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datInt, InitVals=inits, ProfileCol=NA,
                                 cutpoints=cutpointsInt[1,], SampProb=matrix(1, nrow=nrow(datInt), ncol=3)[1,], SampProbiWL =SampProbiInt, w.function=w.functionInt[1])

    print(date())

    CalcDiffs(Fit.int.wl,Fit.int.wl2)


    Fit.slp     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datSlp, InitVals=inits, ProfileCol=NA,
                              cutpoints=cutpointsSlp, SampProb=SampProbSlp, Weights=NoWeighting, w.function=w.functionSlp)
    print(date())
    datSlp$SampProbiNoWt <- 1
    Fit.slp2     <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datSlp, InitVals=inits, ProfileCol=NA,
                              cutpoints=cutpointsSlp[1,], SampProb=SampProbSlp[1,], SampProbiWL=SampProbiNoWt, w.function=w.functionSlp[1])
    print(date())
    CalcDiffs(Fit.slp,Fit.slp2)

    Fit.slp.wl     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datSlp, InitVals=inits, ProfileCol=NA,
                                 cutpoints=cutpointsSlp, SampProb=matrix(1, nrow=nrow(datSlp), ncol=3), Weights=WeightsSlp, w.function=w.functionSlp)
    print(date())

    Fit.slp.wl2     <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datSlp, InitVals=inits, ProfileCol=NA,
                                 cutpoints=cutpointsSlp[1,], SampProb=matrix(1, nrow=nrow(datSlp), ncol=3), SampProbiWL=SampProbiSlp, w.function=w.functionSlp[1])
    print(date())
    CalcDiffs(Fit.slp.wl,Fit.slp.wl2)

    Fit.mix1     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datMix1, InitVals=inits, ProfileCol=NA,
                                  cutpoints=cutpointsMix1, SampProb=SampProbMix1, Weights=NoWeighting, w.function=w.functionMix1)
    print(date())

    Fit.mix2     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datMix2, InitVals=inits, ProfileCol=NA,
                               cutpoints=cutpointsMix2, SampProb=SampProbMix2, Weights=NoWeighting, w.function=w.functionMix2)
    print(date())

    Fit.mix1.wl     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datMix1, InitVals=inits, ProfileCol=NA,
                              cutpoints=cutpointsMix1, SampProb=matrix(1, nrow=nrow(datMix1), ncol=3), Weights=WeightsMix1, w.function=w.functionMix1)
    print(date())

    Fit.mix2.wl     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datMix2, InitVals=inits, ProfileCol=NA,
                                 cutpoints=cutpointsMix2, SampProb=matrix(1, nrow=nrow(datMix2), ncol=3), Weights=WeightsMix2, w.function=w.functionMix2)
    print(date())

    Fit.biv     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datBiv, InitVals=inits, ProfileCol=NA,
                              cutpoints=cutpointsBiv, SampProb=SampProbBiv, Weights=NoWeighting, w.function=w.functionBiv)
    print(date())

    Fit.biv2     <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datBiv, InitVals=inits, ProfileCol=NA,
                              cutpoints=c(cutoffs$IntCutBiv, cutoffs$SlpCutBiv), SampProb=SampProbBiv[1,], SampProbiWL=NoWeighting, w.function="bivar")
    print(date())

    CalcDiffs(Fit.biv,Fit.biv2)

    Fit.biv.wl     <- acml.lmem2(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datBiv, InitVals=inits, ProfileCol=NA,
                              cutpoints=cutpointsBiv, SampProb=matrix(1, nrow=nrow(datBiv), ncol=2), Weights=WeightsBiv, w.function=w.functionBiv)
    print(date())

    Fit.biv.wl2     <- acml.lmem(formula.fixed=y~time*grp+conf, formula.random= ~time, id=id, data=datBiv, InitVals=inits, ProfileCol=NA,
                              cutpoints=c(cutoffs$IntCutBiv, cutoffs$SlpCutBiv), SampProb=c(1,1), SampProbiWL=SampProbiBiv, w.function="bivar")
    print(date())

    CalcDiffs(Fit.biv.wl,Fit.biv.wl2)


    ########## Indirect Imputation Analysis
    print("Indirect Imputation Analyses"); print(date())
    Fit.ran.mi <- IndirectImputation(acml.fit=Fit.ran, datSampled=datRan, datNotSampled=datNotRan,
                                     cutpointsNotSampled=cutpointsNotRan, w.functionNotSampled=w.functionNotRan, SampProbNotSampled=SampProbNotRan, n.imp=n.imp)
    print(date())
    Fit.int.mi <- IndirectImputation(acml.fit=Fit.int, datSampled=datInt, datNotSampled=datNotInt,
                                     cutpointsNotSampled=cutpointsNotInt, w.functionNotSampled=w.functionNotInt, SampProbNotSampled=SampProbNotInt, n.imp=n.imp)
    print(date())
    Fit.slp.mi <- IndirectImputation(acml.fit=Fit.slp, datSampled=datSlp, datNotSampled=datNotSlp,
                                     cutpointsNotSampled=cutpointsNotSlp, w.functionNotSampled=w.functionNotSlp, SampProbNotSampled=SampProbNotSlp, n.imp=n.imp)
    print(date())
    Fit.mix1.mi <- IndirectImputation(acml.fit=Fit.mix1, datSampled=datMix1, datNotSampled=datNotMix1,
                                     cutpointsNotSampled=cutpointsNotMix1, w.functionNotSampled=w.functionNotMix1, SampProbNotSampled=SampProbNotMix1, n.imp=n.imp)
    print(date())

    Fit.mix2.mi <- IndirectImputation(acml.fit=Fit.mix2, datSampled=datMix2, datNotSampled=datNotMix2,
                                     cutpointsNotSampled=cutpointsNotMix2, w.functionNotSampled=w.functionNotMix2, SampProbNotSampled=SampProbNotMix2, n.imp=n.imp)
    print(date())

    out <- list(
    rs.est    = Fit.ran$coefficients,
    rs.mi.est = Fit.ran.mi[["coefficients"]],

    int.est    = Fit.int$coefficients,
    int.wl.est    = Fit.int.wl$coefficients,
    int.mi.est = Fit.int.mi[["coefficients"]],

    slp.est    = Fit.slp$coefficients,
    slp.wl.est    = Fit.slp.wl$coefficients,
    slp.mi.est = Fit.slp.mi[["coefficients"]],

    mix1.est    = Fit.mix1$coefficients,
    mix1.wl.est    = Fit.mix1.wl$coefficients,
    mix1.mi.est = Fit.mix1.mi[["coefficients"]],

    mix2.est    = Fit.mix2$coefficients,
    mix2.wl.est    = Fit.mix2.wl$coefficients,
    mix2.mi.est = Fit.mix2.mi[["coefficients"]],

    rs.cov    = Fit.ran$covariance,
    rs.mi.cov = Fit.ran.mi[["covariance"]],

    int.cov    = Fit.int$covariance,
    int.wl.cov = Fit.int.wl$robcov,
    int.mi.cov = Fit.int.mi[["covariance"]],

    slp.cov    = Fit.slp$covariance,
    slp.wl.cov = Fit.slp.wl$robcov,
    slp.mi.cov = Fit.slp.mi[["covariance"]],

    mix1.cov    = Fit.mix1$covariance,
    mix1.wl.cov = Fit.mix1.wl$robcov,
    mix1.mi.cov = Fit.mix1.mi[["covariance"]],

    mix2.cov    = Fit.mix2$covariance,
    mix2.wl.cov = Fit.mix2.wl$robcov,
    mix2.mi.cov = Fit.mix2.mi[["covariance"]])

}


fit1 <- RunMods(params=c(75, -1, -.5, -2, -.5, log(9),log(1.25),0,log(3.5)),
                NsPerStratumUniv=c(150,100,150),
                N=2000,
                prev.grp=0.3,
                n.imp=50,
                p.central=0.8,
                conf.param = c(-0.25, 0.5),
                ni = c(4, 6))

## Bigger confounder effect
fit2 <- RunMods(params=c(75, -1, -.5, -6, -.5, log(9),log(1.25),0,log(3.5)),
                NsPerStratumUniv=c(150,100,150),
                N=2000,
                prev.grp=0.3,
                n.imp=50,
                p.central=0.8,
                conf.param = c(-0.25, 0.5),
                ni = c(4, 6))

## Smaller variance components (less response dependence)
fit3 <- RunMods(params=c(75, -1, -.5, -2, -.5, log(3),log(.5),0,log(3.5)),
                NsPerStratumUniv=c(150,100,150),
                N=2000,
                prev.grp=0.3,
                n.imp=50,
                p.central=0.8,
                conf.param = c(-0.25, 0.5),
                ni = c(4, 6))

## Bigger snp and snpXtime interaction effects
fit4 <- RunMods(params=c(75, -2.5, -.5, -2, -1.25, log(9),log(1.25),0,log(3.5)),
                NsPerStratumUniv=c(150,100,150),
                N=2000,
                prev.grp=0.3,
                n.imp=50,
                p.central=0.8,
                conf.param = c(-0.25, 0.5),
                ni = c(4, 6))

## negative random intercept and slope correlation
fit5 <- RunMods(params=c(75, -1, -.5, -2, -.5, log(9),log(1.25),-0.25,log(3.5)),
                NsPerStratumUniv=c(150,100,150),
                N=2000,
                prev.grp=0.3,
                n.imp=50,
                p.central=0.8,
                conf.param = c(-0.25, 0.5),
                ni = c(4, 6))

## stronger snp-confounder relationship
fit6 <- RunMods(params=c(75, -1, -.5, -2, -.5, log(9),log(1.25),0,log(3.5)),
                NsPerStratumUniv=c(150,100,150),
                N=2000,
                prev.grp=0.3,
                n.imp=50,
                p.central=0.8,
                conf.param = c(-.75, 1.5),
                ni = c(4, 6))

## Sample more at the extremes
fit7 <- RunMods(params=c(75, -1, -.5, -2, -.5, log(9),log(1.25),0,log(3.5)),
                NsPerStratumUniv=c(175,50,175),
                N=2000,
                prev.grp=0.3,
                n.imp=50,
                p.central=0.8,
                conf.param = c(-0.25, 0.5),
                ni = c(4, 6))

## Use a 60 percent central region
fit8 <- RunMods(params=c(75, -1, -.5, -2, -.5, log(9),log(1.25),0,log(3.5)),
                NsPerStratumUniv=c(150,100,150),
                N=2000,
                prev.grp=0.3,
                n.imp=50,
                p.central=0.6,
                conf.param = c(-0.25, 0.5),
                ni = c(4, 6))

sim_out <- list(fit1=fit1, fit2=fit2, fit3=fit3, fit4=fit4, fit5=fit5, fit6=fit6, fit7=fit7, fit8=fit8)
save(sim_out, file=paste(path,'/output/i',ITER,'.Rdata', sep=""))

# Saving
# save(out, file=paste(path, "/output/i", ITER, ".RData", sep=""))

# Summarizing

#Cases            <- rbind(c(0,log(9)), c(-1, log(9)),c(-5,log(9)), c(-1, log(4)) )

#params           <- t(inits.sim)
#params.tmp       <- split(inits.sim, seq(nrow(inits.sim)))
#NewParms         <- mapply(c, params.tmp, NsCentral)
#Uniqueparams.tmp <- unique(t(NewParms))
#L <- length(NewParms[1,])
#row.keep <- NULL
#for (i in 1:L){ if (all(NewParms[,i]==Uniqueparams.tmp[1,])==TRUE) row.keep <- c(row.keep, i)}

#
# est.parms <- slp.est
# est.cov   <- slp.cov
# params <- NewParms

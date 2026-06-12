pkgname <- "ODSMethods"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ODSMethods')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("WL")
### * WL

flush(stderr()); flush(stdout())

### Name: WL
### Title: Fit model using ascertainment corrected likelihood model (ACML)
### Aliases: WL

### ** Examples

data(gbti)
design <- ods(Response ~ Month|Patient, 'intercept', p_sample=c(1, 0.25, 1),
              data=gbti, quantiles=c(0.1, 0.9), weights="pSample")
est <- WL(Response ~ Month, design)
est
summary(est)



cleanEx()
nameEx("acml")
### * acml

flush(stderr()); flush(stdout())

### Name: acml
### Title: Fit model using ascertainment corrected likelihood model (ACML)
### Aliases: acml

### ** Examples

data(gbti)
design <- ods(Response ~ Month|Patient, 'intercept', p_sample=c(1, 0.25, 1),
              data=gbti, quantiles=c(0.1, 0.9))
est <- acml(Response ~ Month*Genotype, design)
est
summary(est)



cleanEx()
nameEx("bds")
### * bds

flush(stderr()); flush(stdout())

### Name: bds
### Title: Specify a given design for BLUP Dependent Sampling (BDS)
### Aliases: bds

### ** Examples

data(gbti)

odsd <- bds(Response ~ Month|Patient, 'mean', p_sample=c(1, 0.25, 1),
            data=gbti, quantiles=c(0.1, 0.9))
summary(odsd)
plot(odsd)

odsd <- bds(Response ~ Month|Patient, 'intercept', p_sample=c(1, 0.25, 1),
            data=gbti, quantiles=c(0.1, 0.9))
summary(odsd)
plot(odsd)

odsd <- bds(Response ~ Month|Patient, 'slope', p_sample=c(1, 0.25, 1),
            data=gbti, quantiles=c(0.1, 0.9))
summary(odsd)
plot(odsd)

odsd <- bds(Response ~ Month|Patient, 'bivariate', p_sample=c(0.25, 1),
            data=gbti, quantiles=0.8)
summary(odsd)
plot(odsd)




cleanEx()
nameEx("coef")
### * coef

flush(stderr()); flush(stdout())

### Name: fixef
### Title: Extract parameters from fitted models
### Aliases: fixef ranef coef.WL calc_D coef.acml

### ** Examples

data(gbti)
design <- ods(Response ~ Month|Patient, 'intercept', p_sample=c(1, 0.25, 1),
              data=gbti, quantiles=c(0.1, 0.9), weights="pSample")
est <- WL(Response ~ Month, design)
coef(est)
ranef(est)
fixef(est)
data(gbti)
design <- ods(Response ~ Month|Patient, 'intercept', p_sample=c(1, 0.25, 1),
              data=gbti, quantiles=c(0.1, 0.9))
est <- acml(Response ~ Month*Genotype, design)
coef(est)
ranef(est)
fixef(est)



cleanEx()
nameEx("gbti")
### * gbti

flush(stderr()); flush(stdout())

### Name: gbti
### Title: Simulated Group by Time Interaction (gbti) dataset
### Aliases: gbti

### ** Examples

  data(gbti)



cleanEx()
nameEx("ods")
### * ods

flush(stderr()); flush(stdout())

### Name: ods
### Title: Specify a given design for Outcome Dependent Sampling (ODS)
### Aliases: ods

### ** Examples

data(gbti)

odsd <- ods(Response ~ Month|Patient, 'mean', p_sample=c(1, 0.25, 1),
            data=gbti, quantiles=c(0.1, 0.9))
summary(odsd)
plot(odsd)

odsd <- ods(Response ~ Month|Patient, 'intercept', p_sample=c(1, 0.25, 1),
            data=gbti, quantiles=c(0.1, 0.9))
summary(odsd)
plot(odsd)

odsd <- ods(Response ~ Month|Patient, 'slope', p_sample=c(1, 0.25, 1),
            data=gbti, quantiles=c(0.1, 0.9))
summary(odsd)
plot(odsd)

odsd <- ods(Response ~ Month|Patient, 'bivariate', p_sample=c(0.25, 1),
            data=gbti, quantiles=0.8)
summary(odsd)
plot(odsd)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

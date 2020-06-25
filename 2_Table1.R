library(tidyverse)
library(nlme)
library(lme4)
source("neighborhood_functions.R")
#model specification
#response is number of species
#ln.NAtm.NAdd = log(Added N +1)
#Year linear covariate of time
#means.Sr.4_mean times series mean of the mean of the number of species\n
#in the four adjacent plots
#a separate variance term is used for each treatment
dat <- read.csv("2020-02-03-e001.csv")
Table1 <- lme(sr~
                ln.NAtm.NAdd+Year+
                means.Sr.4_mean,
              random=~1|Plot,
              weights=varIdent(form=~1|NTrt),
              correlation=corARMA(form=~Year|Plot,p=2,q=1),
              control=list(maxIter = 10000,msMaxIter=10000,niterEM=10000,msMaxEval = 10000),
              method="REML",
              data=dat)
summary(Table1)
intervals(Table1,which="fixed")
#
#print out summary table
out <- summary(Table1)$tTable
out[,1:4] <- apply(out[,1:4],2,function(x)round(x,2))
out <- as.data.frame(out)
out$Parameter <- row.names(out)
out <- out %>% select(Parameter,everything())
write.csv(x=out,file = "Tables/Table1_neighe001.csv",row.names = FALSE)

#
SuppTable1 <- lme(sr~
                    ln.NAtm.NAdd+Year+
                    means.Sr.4,
                  random=~1|Plot,
                  weights=varIdent(form=~1|NTrt),
                  correlation=corARMA(form=~Year|Plot,p=2,q=1),
                  control=list(maxIter = 10000,msMaxIter=10000,niterEM=10000,msMaxEval = 10000),
                  method="REML",
                  data=dat)
summary(SuppTable1)
#print out summary table
SuppTab <- summary(SuppTable1)$tTable
SuppTab[,1:4] <- apply(SuppTab[,1:4],2,function(x)round(x,2))
SuppTab <- as.data.frame(SuppTab)
SuppTab$Parameter <- row.names(SuppTab)
SuppTab <- SuppTab %>% select(Parameter,everything())
write.csv(x=SuppTab,file = "Tables/SuppTable1_neighe001.csv",row.names = FALSE)

#

#line 172
base <- dat %>% filter(Year==1982) %>% summarise(sr=mean(sr))
0.3/base$sr
#test with type III
car::Anova(Table1,type="III")
#
#test with edge effects
#edge effect is a parameter that indicates if a plot is on the edge
Table1 <- update(Table1, method="ML")
Table1 <- MASS::stepAIC(Table1)#backwards selection
#
Table1.edge <- update(Table1,.~. +edgeeffect)
#
# summary(Table1.edge)
anova(Table1,Table1.edge)
aicfun(Table1,Table1.edge)#custom function to print delta AIC and LRR
#
Table1.edge <- update(Table1.edge, method="REML")
summary(Table1.edge)
#############
#more detailed tests of the model specifications
#
Table1.0 <- lme(sr~
                  ln.NAtm.NAdd+Year,
                random=~1|Plot,
                control=list(maxIter = 10000,msMaxIter=10000,niterEM=10000,msMaxEval = 10000),
                method="ML",
                data=dat)
#add in variance weights and autocorrelation functions
Table1.0.ntrt <- update(Table1.0,weights = varIdent(form=~1|NTrt))
Table1.0.ntrt.cor <- update(Table1.0.ntrt,correlation=corARMA(form=~Year|Plot,p=2,q=1))
#
Table1.1 <- lme(sr~
                  ln.NAtm.NAdd+Year+
                  means.Sr.4_mean,
                random=~1|Plot,
                control=list(maxIter = 10000,msMaxIter=10000,niterEM=10000,msMaxEval = 10000),
                method="ML",
                data=dat)
# summary(Table1.1)
#
Table1.1.ntrt <- update(Table1.1,weights = varIdent(form=~1|NTrt))
Table1.1.ntrt.cor <- update(Table1.1.ntrt,
                            correlation=corARMA(form=~Year|Plot,p=2,q=1))
#
#SupplementalTable1
# anova(Table1.0,Table1.1)#compare neighborhood effect 
# anova(Table1.0.ntrt,Table1.1.ntrt)#compare with variance weights
# anova(Table1.0.ntrt.cor,Table1.1.ntrt.cor)#compare with ACF
#
Table1.1 <- update(Table1.1,method="REML")
Table1.1.ntrt <- update(Table1.1.ntrt,method="REML")
Table1.1.ntrt.cor <- update(Table1.1.ntrt.cor,method="REML")

##############
#other models 
#try with time as cubic polynomial 
Table1.cub <- lme(sr~
                    ln.NAtm.NAdd+poly(Year,3),
                  random=~1|Plot,
                  weights = varIdent(form=~1|NTrt),
                  correlation=corARMA(form=~Year|Plot,p=2,q=1),
                  control=list(maxIter = 10000,msMaxIter=10000,niterEM=10000,msMaxEval = 10000),
                  method="ML",
                  data=dat)
#
Table1.cub.1 <- lme(sr~
                      ln.NAtm.NAdd+means.Sr.4_mean+poly(Year,3),
                    random=~1|Plot,
                    weights = varIdent(form=~1|NTrt),
                    correlation=corARMA(form=~Year|Plot,p=2,q=1),
                    control=list(maxIter = 10000,msMaxIter=10000,niterEM=10000,msMaxEval = 10000),
                    method="ML",
                    data=dat)
# anova(Table1.cub,Table1.cub.1)
Table1.cub.1 <- update(Table1.cub.1,method="REML")
# car::Anova(Table1.cub.1,type="III")
# summary(Table1.cub.1)
########
#try with a poisson as data is count
#perhaps one should use poisson for integers, but it is more complicated and it does not change the results
Table1.glmer.0 <- glmer(sr~
                          ln.NAtm.NAdd_c+Year_c+
                          (1|Plot),
                        family="poisson",
                        data=dat)
#
dat$means.Sr.4_mean_c <- scale_this(dat$means.Sr.4_mean)
Table1.glmer.1 <- glmer(sr~
                          ln.NAtm.NAdd_c+means.Sr.4_mean_c+
                          Year_c+(1|Plot),
                        family="poisson",
                        data=dat)
# summary(Table1.glmer.1)
# anova(Table1.glmer.0,Table1.glmer.1)
#
# plot(Table1.1.ntrt.cor)
# plot(Table1.glmer.1)
#####
Table1.0.ml <- update(Table1.0,method="ML")
Table1.1.ml <- update(Table1.1,method="ML")
Table1.0.ntrt.ml  <- update(Table1.0.ntrt,method="ML")
Table1.1.ntrt.ml  <- update(Table1.1.ntrt,method="ML")
Table1.0.ntrt.cor.ml  <- update(Table1.0.ntrt.cor,method="ML")
Table1.1.ntrt.cor.ml  <- update(Table1.1.ntrt.cor,method="ML")
Table1.cub.ml  <- update(Table1.cub,method="ML")
Table1.cub.1.ml  <- update(Table1.cub.1,method="ML")
#Data for supplementable Table 1
SuppTable1 <- data.frame(
  neighbordiversity = c(
    summary(Table1.1)$tTable["means.Sr.4_mean",1],
    summary(Table1.1.ntrt)$tTable["means.Sr.4_mean",1],
    summary(Table1.1.ntrt.cor)$tTable["means.Sr.4_mean",1],
    summary(Table1.cub.1)$tTable["means.Sr.4_mean",1],
    summary(Table1.glmer.1)$coefficients["means.Sr.4_mean_c",1]),
  SE =c(
    summary(Table1.1)$tTable["means.Sr.4_mean",2],
    summary(Table1.1.ntrt)$tTable["means.Sr.4_mean",2],
    summary(Table1.1.ntrt.cor)$tTable["means.Sr.4_mean",2],
    summary(Table1.cub.1)$tTable["means.Sr.4_mean",2],
    summary(Table1.glmer.1)$coefficients["means.Sr.4_mean_c",2]),
  LRR_Chi=c(
    anova(Table1.0.ml,Table1.1.ml)$L.Ratio[2],
    anova(Table1.0.ntrt.ml,Table1.1.ntrt.ml)$L.Ratio[2],
    anova(Table1.0.ntrt.cor.ml,Table1.1.ntrt.cor.ml)$L.Ratio[2],
    anova(Table1.cub.ml,Table1.cub.1.ml)$L.Ratio[2],
    anova(Table1.glmer.0,Table1.glmer.1)$Chisq[2]),
  p.value=c(
    anova(Table1.0.ml,Table1.1.ml)$`p-value`[2],
    anova(Table1.0.ntrt.ml,Table1.1.ntrt.ml)$`p-value`[2],
    anova(Table1.0.ntrt.cor.ml,Table1.1.ntrt.cor.ml)$`p-value`[2],
    anova(Table1.cub.ml,Table1.cub.1.ml)$`p-value`[2],
    anova(Table1.glmer.0,Table1.glmer.1)$`Pr(>Chisq)`[2]))
write.csv(x = SuppTable1,file = "Tables/SuppTable2_neighe001.csv",row.names = FALSE)
#################
#try the same model but using the last 10 years only
dat2004 <- dat %>% filter(Year>1994) %>% select(-means.Sr.4_mean)
print(length(unique(dat2004$Year)))
dat2004_m <- dat2004 %>% group_by(Plot) %>% 
  summarise(means.Sr.4_mean=mean(means.Sr.4))
dat2004 <- left_join(dat2004,dat2004_m)
final.2004 <- lme(sr~
                    ln.NAtm.NAdd+Year+
                    means.Sr.4_mean,
                  random=~1|Plot,
                  weights = varIdent(form=~1|NTrt),
                  correlation=corARMA(form=~Year|Plot,p=2,q=1),
                  control=list(maxIter = 10000,msMaxIter=10000,niterEM=10000,msMaxEval = 10000),
                  method="ML",
                  data=dat2004)
summary(final.2004)
car::Anova(final.2004,type="III")
#
final2004.edge <- update(final.2004,.~. +edgeeffect)
summary(final2004.edge)
anova(final.2004,final2004.edge)
aicfun(final.2004,final2004.edge)
final.2004 <- update(final.2004,method="REML")
summary(final.2004)
final2004.edge <- update(final2004.edge,method="REML")
summary(final2004.edge)
#######
library(tidyverse)
library(nlme)
library(lme4)
source("helperfunctions_e001Cneigh.R")
#model specification
#response is number of species
#ln.NAtm.NAdd = log(Added N +1)
#Year linear covariate of time
#means.Sr.4_mean times series mean of the mean of the number of species\n
#in the four adjacent plots
#a separate variance term is used for each treatment
dat <- read.csv("2021-06-25-e001.csv")
dat$Year <- dat$Year-1982
Table1 <- lme(sr~
                ln.NAtm.NAdd+Year+
                means.Sr.4_mean,
              random=~1|Plot,
              weights=varIdent(form=~1|NTrt),
              correlation=corARMA(form=~Year|Plot,p=2,q=1),
              control=list(maxIter = 10000,msMaxIter=10000,niterEM=10000,msMaxEval = 10000),
              method="REML",
              data=dat)
plot(Table1)
summary(Table1)
length(unique(dat$Year))
intervals(Table1,which="fixed")
###print out summary table
out <- summary(Table1)$tTable
out[,1:4] <- apply(out[,1:4],2,function(x)round(x,2))
out <- as.data.frame(out)
out$Parameter <- row.names(out)
out <- out %>% select(Parameter,everything())
#####
library(flextable)
library(officer)
outtab <- as.data.frame(out)
# outtab$coef <- rownames(outtab)
#
outtab$`p-value` <- ifelse(outtab$`p-value` ==0,
                      "<0.001",
                      ifelse(outtab$`p-value` <0.01,
                             "<0.01",
                             round(outtab$`p-value` ,3)))
#
outtab1 <- flextable(outtab) %>% fontsize(size=12)
set_table_properties(outtab1, width = 1, layout = "autofit")
outtab1 <- font(outtab1,fontname = "Times")
doc1 <- read_docx() %>% 
  body_add_par(value = "Supplemental Table 1", style = "Normal") %>%
  body_add_flextable(value = outtab1)
#
print(doc1, target = "Tables/SupplementalTable1.docx")

#####
# write.csv(x=out,file = "Tables/SupplementalTable1_neighe001.csv",row.names = FALSE)
#line in results comparing effectsize to base
#test with type III
car::Anova(Table1,type="III")
########
#end of paper results
########
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
# anova(Table1.0.ntrt.cor,Table1.1.ntrt.cor)#compare 
#
Table1.1 <- update(Table1.1,method="REML")
Table1.1.ntrt <- update(Table1.1.ntrt,method="REML")
Table1.1.ntrt.cor <- update(Table1.1.ntrt.cor,method="REML")
#######
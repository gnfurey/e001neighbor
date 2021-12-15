####
source("helperfunctions_e001Cneigh.R")
library(gridExtra)
library(tidyverse)
library(nlme)
library(spdep)
library(viridis)
library(gstat)#see package details for class that is required
e001 <- read_csv("2021-06-25-e001.csv")
library(broom)
dat_f2 <- e001 %>% filter(NTrt=="I-0.00")#true controls
#
dat_f2$Year <- dat_f2$Year-1982#to make intercept understandable
############
#run a model across time with an interaction between time \n
# and neighborhood biodiversity 
mod0 <- lm(sr~Year*means.Sr.4,
                  data=dat_f2)
#
summary(mod0)
car::Anova(mod0,type="III")
#
SuppTab2.0 <- lme(sr~
                    means.Sr.4+Year,
                  random=~1|Plot,
                  correlation=corARMA(form=~Year|Plot,p=2,q=1),
                  control=list(maxIter = 10000,
                               msMaxIter=10000,
                               niterEM=10000,
                               msMaxEval = 10000),
                  method="ML",
                  data=dat_f2)
summary(SuppTab2.0)
plot(SuppTab2.0)
#test interaction with likelihood ratio test
SuppTab2.1 <- update(SuppTab2.0,~Year*means.Sr.4)
summary(SuppTab2.1)
#
anova(SuppTab2.0,SuppTab2.1)
anova(SuppTab2.1,type="marginal")
car::Anova(SuppTab2.1,"III")
#
SuppTab2.1 <- update(SuppTab2.1,method="REML")
summary(SuppTab2.1)
#
#
plot(SuppTab2.1)
#
anova(SuppTab2.1,type="marginal")
#
library(flextable)
library(officer)
outtab <- as.data.frame(anova(SuppTab2.1,type="marginal"))
outtab$coef <- rownames(outtab)
#
outtab$`F-value` <- round(outtab$`F-value`,2)
outtab$pval <- p.adjust(outtab$`p-value`,"fdr")
outtab$`p-value` <- NULL
outtab$pval <- ifelse(outtab$pval==0,
                        "<0.001",
                        ifelse(outtab$pval<0.01,
                               "<0.01",
                               round(outtab$pval,3)))
#
outtab <- outtab %>% select(coef,numDF:pval)
outtab1 <- flextable(outtab) %>% fontsize(size=12)
set_table_properties(outtab1, width = 1, layout = "autofit")
outtab1 <- font(outtab1,fontname = "Times")
doc1 <- read_docx() %>% 
  body_add_par(value = "Supplemental Table 2", style = "Normal") %>%
  body_add_flextable(value = outtab1)
#
print(doc1, target = "Tables/SupplementalTable2.docx")
#
summary(SuppTab2.1)$tTable
outcoef <- as.data.frame(summary(SuppTab2.1)$tTable)
outcoef
outcoef$coef <- rownames(outcoef)
#
outcoef$Value <- round(outcoef$Value,3)
outcoef$Std.Error <- round(outcoef$Std.Error,2)
outcoef$`t-value` <- round(outcoef$`t-value`,2)
outcoef$`p-value` <- ifelse(outcoef$`p-value`==0,
                      "<0.001",
                      ifelse(outcoef$`p-value`<0.01,
                             "<0.01",
                             round(outcoef$`p-value`,3)))
#
outcoef <- outcoef %>% select(coef,Value:`p-value`)
#
outtab2 <- flextable(outcoef) %>% fontsize(size=12)
set_table_properties(outtab2, width = 1, layout = "autofit")
outtab2 <- font(outtab2,fontname = "Times")
doc2 <- read_docx() %>% 
  body_add_par(value = "Supplemental Table 3", style = "Normal") %>%
  body_add_flextable(value = outtab2)
#
print(doc2, target = "Tables/SupplementalTable3.docx")
#####
# write.csv(x = anova(SuppTab2.1,type="marginal"),
#           file="Tables/SuppTable2_neighe001.csv")
####
#build a matrix for predictions
a <- mean(dat_f2$means.Sr.4)+2*sd(dat_f2$means.Sr.4)
mean <- mean(dat_f2$means.Sr.4)
b <- mean(dat_f2$means.Sr.4)-2*sd(dat_f2$means.Sr.4)
summary(SuppTab2.1)
list1 <- list(means.Sr.4=c(a,mean,b),
              Year=seq(from=0,to=22))
##
library(emmeans)
out <- emmeans(object = SuppTab2.1, 
               specs = ~Year*means.Sr.4,
               var="Year",
               df=5,
               type="response",
               sigmaAdjust = TRUE,
               at=list1,
               adjust="none")
out
out1 <- as.data.frame(out)
out1$emmean
p1 <- ggplot(out1,aes(x=Year,y=emmean,
                      col=as.factor(means.Sr.4),
                      fill=as.factor(means.Sr.4)))+
  geom_ribbon(aes(ymax=upper.CL,ymin=lower.CL),
              col= NA,
              alpha=0.4)+
  geom_line(size=1)+
  theme_bw()+
  ylab("Species Richness")+
  theme(legend.position = "bottom",
        legend.margin = margin(0,0,0,0),
        legend.key.height = unit(0.2,"cm"),
        legend.key.width = unit(0.3,"cm"),
        legend.text = element_text(size=8),
        legend.box.spacing = unit(0.1,"cm"),
        legend.box.margin = margin(0,0,0,0,"cm"))+
  scale_color_viridis(begin = 0,end = 0.8,
                      name="NeighborSr",
                      labels=c("-2SD","Mean","+2SD"),
                      discrete = TRUE,
                      option = "A")+
  scale_fill_viridis(begin = 0,end = 0.8,
                     name="NeighborSr",
                     labels=c("-2SD","Mean","+2SD"),
                     discrete = TRUE,
                      option = "A")+
  guides(col=guide_legend(nrow=2,byrow=TRUE))
p1
#
p1 <- p1+
  annotate(x = 26,y = 18,geom = "label",
           label="mu+2*SD",
           parse=TRUE)+
  # annotate(x = 2008,y = 15,geom = "text",
  #          label="mu+1*SD",
  #          parse=TRUE)+
  annotate(x = 26,y = 11.5,geom = "label",
           label="mu",
           parse=TRUE)+
  # annotate(x = 2008,y = 8,geom = "text",
  #          label="mu-1*SD",
  #          parse=TRUE)+
  annotate(x = 26,y = 5.3,geom = "label",
           label="mu-2*SD",
           parse=TRUE)+
  scale_x_continuous(breaks = c(0,10,20),
                     limits=c(0,26.2),
                     expand = expansion(mult=c(0.03,0.15)))+
  scale_y_continuous(limits=c(0,22),
                     breaks=c(0,5,10,15,20,25))
p1
#
library(colorBlindness)
cvdPlot(p1)
#
ggsave(filename = "Figures/Figure3_neighe001.eps",
       plot = p1,
       device = cairo_ps,
       dpi=600,
       height=2.5,width=3,unit="in")


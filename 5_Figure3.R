####
source("neighborhood_functions.R")
library(gridExtra)
library(tidyverse)
library(nlme)
library(spdep)
library(viridis)
library(gstat)#see package details for class that is required
#ln.NAtm.NAdd = log(Added N +1)
e001 <- read_csv("2021-06-25-e001.csv")
library(broom)
dat_f2 <- e001 %>% filter(NTrt=="I-0.00")#true controls
#
############
#run a model across time with an interaction between time \n
# and neighborhood biodiversity 
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
SuppTab2.1 <- update(SuppTab2.0,~means.Sr.4*Year)
summary(SuppTab2.1)
#
anova(SuppTab2.0,SuppTab2.1)
anova(SuppTab2.1,type="marginal")
car::Anova(SuppTab2.1,"III")
aicfun(SuppTab2.0,SuppTab2.1)
#
#
SuppTab2.1 <- update(SuppTab2.1,method="REML")
summary(SuppTab2.1)
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
# write.csv(x = anova(SuppTab2.1,type="marginal"),
#           file="Tables/SuppTable2_neighe001.csv")
####
#build a matrix for predictions
a <- mean(dat_f2$means.Sr.4)+2*sd(dat_f2$means.Sr.4)
b <- mean(dat_f2$means.Sr.4)+1*sd(dat_f2$means.Sr.4)
mean <- mean(dat_f2$means.Sr.4)
c <- mean(dat_f2$means.Sr.4)-1*sd(dat_f2$means.Sr.4)
d <- mean(dat_f2$means.Sr.4)-2*sd(dat_f2$means.Sr.4)
##critical value
crit <- mean(dat_f2$means.Sr.4)+1.37#9.4
#compare to baseline of 12 sr
9.4/12
##
means.Sr.4 <- c(d,c,mean,b,a)
means.Sr.4 <- round(means.Sr.4,2)
mean <- round(mean,2)
Year <- unique(dat_f2$Year)
#build matrix
newdat <- expand_grid(means.Sr.4,Year)
#predictions
newdat$fits <- predict(SuppTab2.1,newdata = newdat,level = 0)
#
parm <- c("-2SD","-1SD","Mean","+1SD","+2SD")
match <- data.frame(means.Sr.4=means.Sr.4,
                    NeighborSr = parm)
#
match$NeighborSr <- factor(match$NeighborSr,levels = parm)
#
newdat <- left_join(newdat,match)
#make an interaction plot
#
p1 <- ggplot(newdat,aes(x=Year,y=fits,
                  col=NeighborSr))+
  geom_line(size=2)+
  theme_classic()+
  geom_hline(yintercept = 13,linetype=2)+
  ylab("Species Richness")+
  theme(legend.position = "bottom",
        legend.margin = margin(0,0,0,0),
        legend.key.height = unit(0.2,"cm"),
        legend.key.width = unit(0.3,"cm"),
        legend.text = element_text(size=8),
        legend.box.spacing = unit(0.1,"cm"),
        legend.box.margin = margin(0,0,0,0,"cm"))+
  scale_color_viridis(begin = 0,end = 0.8,discrete = TRUE,
                      option = "A")+
  annotate(x=2000,y=13,geom="label",label=9.4,size=3)+
  guides(col=guide_legend(nrow=2,byrow=TRUE))
p1
#
p1 <- p1+
  annotate(x = 2008,y = 18,geom = "text",
           label="mu+2*SD",
           parse=TRUE)+
  annotate(x = 2008,y = 15,geom = "text",
           label="mu+1*SD",
           parse=TRUE)+
  annotate(x = 2005.3,y = 11.5,geom = "text",
           label="mu",
           parse=TRUE)+
  annotate(x = 2008,y = 8,geom = "text",
           label="mu-1*SD",
           parse=TRUE)+
  annotate(x = 2008,y = 5,geom = "text",
           label="mu-2*SD",
           parse=TRUE)+
  scale_x_continuous(breaks = c(1982,1992,2002),
                     expand = expansion(mult=c(0.03,0.15)))
#
ggsave(filename = "Figures/Figure3_neighe001.eps",
       plot = p1,
       device = cairo_ps,
       dpi=600,
       height=2.5,width=3,unit="in")


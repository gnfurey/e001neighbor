####
source("neighborhood_functions.R")
library(gridExtra)
library(tidyverse)
library(nlme)
library(spdep)
library(viridis)
library(gstat)#see package details for class that is required
#ln.NAtm.NAdd = log(Added N +1)
e001 <- read_csv("2020-02-03-e001.csv")
library(broom)
dat_f2 <- e001 %>% filter(NTrt=="I-0.00")#true controls
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
SuppTab2.1 <- update(SuppTab2.0,~means.Sr.4*Year)
anova(SuppTab2.0,SuppTab2.1)
aicfun(SuppTab2.0,SuppTab2.1)
#
SuppTab2.1.edge <- update(SuppTab2.1,.~. +edgeeffect)
anova(SuppTab2.1,SuppTab2.1.edge)#minor edge effect
aicfun(SuppTab2.1,SuppTab2.1.edge)
SuppTab2.1.edge <- update(SuppTab2.1.edge,method="REML")
summary(SuppTab2.1.edge)
SuppTab2.1 <- update(SuppTab2.1,method="REML")
summary(SuppTab2.1)
anova(SuppTab2.1,type="marginal")
write.csv(x = anova(SuppTab2.1,type="marginal"),
          file="Tables/SuppTable3_neighe001.csv")
# write.csv(anova(SuppTab2.1,type="marginal"),row.names = FALSE)
dat_f2 <- as.data.frame(dat_f2)
#run a regression of neighborhood effects in each year
mods <- dat_f2 %>% group_by(Year) %>% do(mod1=lm(sr~means.Sr.4,data=.))
mods.t <- tidy(mods,mod1) %>% filter(term=="means.Sr.4")
#
p1 <- ggplot(mods.t,aes(x=Year,y=estimate))+
  geom_smooth(span=1,se=FALSE)+
  geom_point()+
  theme_bw()+
  ylab("Focal Plot Number of Species \n ~ f(Neighborhood Number of Species)")+
  geom_hline(yintercept=0,linetype=2)+
  geom_errorbar(aes(ymax=estimate+std.error,ymin=estimate-std.error),width=0.1)+
  theme(panel.grid = element_blank(),
        plot.margin=margin(t = 0,r = 1,b = 0,l = 1),
        title = element_text(size=8))+
  scale_x_continuous(limits=c(1982,2004),breaks=c(1982,1990,1998,2004))+
  ggtitle("A")
p1
#get the fitted values
aug <- augment(mods,mod1)
#
p2 <- ggplot(aug,aes(x=means.Sr.4,y=sr,col=Year,group=Year))+
  geom_smooth(method="lm",se=FALSE)+
  theme_bw()+ylab("Focal Plot Number of Species")+
  xlab("Neighborhood Number of Species")+
  scale_color_viridis(option="A",begin=0,end = 0.8)+
  scale_x_continuous(breaks=seq(from=3,to=16,by=1))+
  scale_y_continuous(breaks=seq(from=0,to=25,by=3))+
  theme(legend.position = "bottom",
        plot.title = element_text(margin = margin(0,0,0,0)),
        plot.margin=margin(t = 0,r = 1,b = 0,l = 1),
        # panel.grid = element_blank(),
        legend.key.height =  unit(0.20, "cm"),
        legend.key.width = unit(1.6,"cm"),
        title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6),
        legend.box.spacing = unit(0,"cm"),
        legend.box.margin = margin(0,0,0,0,"cm"))+
  ggtitle("B")
grid.arrange(p1,p2)
#
ggsave(filename = "Figures/Figure2_neighe001.eps",
       device = cairo_ps,
       dpi=600,
       plot = grid.arrange(p1,p2),
       height=105,width=82,unit="mm")

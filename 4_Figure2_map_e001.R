#
library(tidyverse)
library(mgcv)
library(sp)
library(viridis)
library(spdep)
library(gstat)
#
source("neighborhood_functions.R")
e001 <- read_csv("2021-06-25-e001.csv")
e001$NTrt <- as.factor(e001$NTrt)
e001$NTrt <- factor(e001$NTrt,levels(e001$NTrt)[c(9,1,2:8)])
#mean of last 10 years
meansr <- e001 %>% 
  filter(Year<2005) %>%
  # filter(Year>1994)%>%
  group_by(Plot,NAdd,ln.NAtm.NAdd,row,col,NTrt,edgeeffect) %>%
  summarise_at(vars(means.Sr.4,sr,means.NAdd.4,Evenness,ShanWinr,EtoH),
               funs(m=mean,se=my.stand))
#in caption
round(range(meansr$sr_m),1)
##reverse columns to have the plots face north
col2 <- data.frame(col2=c(1,2,3,4,5,6,7,8,9),col=c(9,8,7,6,5,4,3,2,1))
meansr <- left_join(meansr,col2)
row2 <- data.frame(row2=c(1,2,3,4,5,6),row=c(6,5,4,3,2,1))
meansr <- left_join(meansr,row2)
#plot n treatment in 2d space
nadd.map <- ggplot(meansr,aes(x=row2,y=col2,fill=NTrt,label=Plot))+
  geom_tile(col="black")+
  scale_fill_viridis(discrete=TRUE,option="A",direction=-1,begin=0.3,end = 1)+
  ggtitle("b")+
  ylab("")+
  xlab("")+
  geom_text()+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6))+
  theme_classic()+
  theme(legend.justification = c(0,-1),
        legend.position = "bottom",
        legend.key.height =  unit(0.05, "cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.text = element_text(size=7),
        legend.box.spacing = unit(-0.4,"cm"),
        legend.box.margin = margin(0,0,0,0,"cm"),
        plot.margin = margin(0,0,0,0),
        axis.ticks = element_blank(),
        axis.text = element_blank())
nadd.map
#
#plot sr in 2d
sr.map <- ggplot(meansr,aes(x=row2,y=col2,fill=sr_m))+
  geom_tile(col="black")+
  scale_fill_viridis(option="D",direction=1,begin=0.0,end = 1,
                     name="SR")+
  ggtitle("a")+
  ylab("")+
  xlab("")+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6))+
  theme_classic()+
  theme(legend.justification = c(0,-1),
        legend.position = "bottom",
        legend.key.height =  unit(0.2, "cm"),
        legend.key.width = unit(1.2,"cm"),
        legend.text = element_text(size=6),
        legend.box.spacing = unit(-0.4,"cm"),
        legend.box.margin = margin(0,0,0,0,"cm"),
        plot.margin = margin(0,0,0,0),
        axis.ticks = element_blank(),
        axis.text = element_blank())
sr.map
#
o1 <- gridExtra::grid.arrange(sr.map,nadd.map,ncol=1)
ggsave(filename = "Figures/Figure2_map.eps",
       plot = o1,
       device = cairo_ps,
       dpi=600,
       width=3,height=4)
################
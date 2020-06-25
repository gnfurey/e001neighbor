source("neighborhood_functions.R")
library(gridExtra)
library(tidyverse)
library(nlme)
library(spdep)
library(viridis)
library(gstat)#see package details for class that is required
#ln.NAtm.NAdd = log(Added N +1)
e001 <- read_csv("2020-02-03-e001.csv")
#get mean of last 10 years 
meansr_I <- e001 %>% 
  filter(NTrt=="I-0.00") %>% 
  filter(Year<2005) %>%
  filter(Year>1994)%>%
  group_by(Plot,NAdd,ln.NAtm.NAdd,row,col,NTrt,edgeeffect) %>%
  summarise_at(vars(means.Sr.4,sr,means.NAdd.4,Evenness,ShanWinr,EtoH,Agrorepe.neighbor,Agrorepe),funs(m=mean,se=my.stand))
#######
#run a simple OLS regression using only control plots
mod1.lm <- lm(sr_m~means.Sr.4_m,data=meansr_I)
summary(mod1.lm)
mod1.edge <- lm(sr_m~means.Sr.4_m+edgeeffect,data=meansr_I)
summary(mod1.edge)
#
#Figure3
justI <- e001 %>% filter(NTrt=="I-0.00")#only control
ggplot(justI,aes(x=Year,y=sr))+
  geom_point()+
  facet_wrap(~Plot)
#run a regression for each plot through time
mods <- e001 %>% filter(NTrt=="I-0.00") %>%
  group_by(Plot) %>% do(mod1=lm(sr~Year,data=.))
mods.t <- tidy(mods,mod1) %>% filter(term=="Year")
mods.t$Plot <- as.factor(mods.t$Plot)
#get fitted values
p.adjust(mods.t$p.value,method="fdr")
fits <- augment(mods,mod1)
#merge
justI <- left_join(justI,fits)
#
# ggplot(justI,aes(x=Year,y=sr))+geom_point()+
#   facet_wrap(~Plot)
#
dat_f3 <- e001 %>% 
  filter(NTrt=="I-0.00") %>% 
  filter(Year<2005) %>%
  # filter(Year>1994)%>%
  # filter(Year<1987) %>%
  group_by(Plot,NAdd,ln.NAtm.NAdd,row,col,NTrt,edgeeffect) %>%
  summarise_at(vars(means.Sr.4,sr,means.NAdd.4,Evenness,ShanWinr,EtoH,Agrorepe.neighbor,Agrorepe),funs(m=mean,se=my.stand))
#
mods.t$Plot <- as.character(mods.t$Plot)
dat_f3$Plot <- as.character(dat_f3$Plot)
#
avg <- left_join(mods.t,dat_f3)
#
mods.t$Plot <- as.factor(mods.t$Plot)
justI$Plot <- as.factor(justI$Plot)
justI <- left_join(justI,mods.t)
justI$Plot <- paste("Plot = ",justI$Plot,sep="")
##############
#show control plot slope through time
# fits$Plot <- paste("Plot = ",fits$Plot,sep="")
fits$Plot <- as.factor(fits$Plot)

fits$yr <- avg$estimate[match(fits$Plot,avg$Plot)]
#
un <- fits %>% ungroup() %>%  select(Plot,yr) %>% unique()
un
#
p1 <- ggplot(fits,aes(x=Year,y=.fitted))+
  geom_line()+
  ylab("Focal Plot Number of Species")+
  theme_bw()+
  geom_label(aes(x=1984,y=23.5,label=Plot),data=un,
             label.padding = unit(0.1, "lines"))+
  facet_wrap(~Plot)+
  geom_line(aes(y=.fitted+.se.fit),linetype=2,alpha=0.4)+
  geom_line(aes(y=.fitted-.se.fit),linetype=2,alpha=0.4)+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(margin =margin(0,0,0,0)),
        plot.margin=margin(t = 0,r = 1,b = 0,l = 1),
        title = element_text(size=8),
        axis.text = element_text(size=8))+
  scale_y_continuous(limits=c(3,25),
                     breaks=seq(from=3,to=24,by=8))+
  scale_x_continuous(limits=c(1982,2004),breaks=c(1982,1992,2002))+
  geom_point(aes(x=Year,y=sr),shape=21)+
  ggtitle("A")
p1
#############
#show slope through time compared to neighborhood biodiversity
p2 <- ggplot(avg,aes(x=means.Sr.4_m,
                     y=estimate,label=Plot,
                     col=p.value<0.05))+
  geom_text(size=5)+
  ylab("Focal Plot Number of Species \n~ f(Year)")+
  # geom_abline(slope =0.038,intercept = -0.48)+
  theme_bw()+
  scale_color_brewer(palette = "Dark2")+
  geom_hline(yintercept = 0,linetype=2)+
  geom_errorbar(aes(ymax=estimate+std.error,
                    ymin=estimate-std.error),
                width=0.1,alpha=0.3)+
  xlab(expression(paste("Neighborhood Number of Species ","(",mu["[1982-2004]"],")")))+
  theme(plot.title = element_text(margin =margin(0,0,0,0)),
        axis.title.x = element_text(margin=margin(0,0,0,0)),
        legend.key.height =  unit(0.05, "cm"),
        legend.box.spacing = unit(0,"cm"),
        plot.margin=margin(t = 0,r = 1,b = 0,l = 1),
        legend.margin=margin(c(0,0,0,0)),
        legend.box.margin=margin(0,0,0,0),
        legend.position = "bottom",
        title = element_text(size=8))+
  ggtitle("B")
grid.arrange(p1,p2)
#
# theme(legend.justification = c(0,-1),
#       legend.position = "bottom",
#       legend.key.height =  unit(0.05, "cm"),
#       legend.key.width = unit(0.4,"cm"),
#       legend.text = element_text(size=6),
#       legend.box.spacing = unit(-0.8,"cm"),
#       legend.box.margin = margin(0,0,0,0,"cm"))
mod1 <- lm(estimate~means.Sr.4_m,data=avg)
summary(mod1)
#
mod2 <- lm(estimate~means.Sr.4_m+edgeeffect,data=avg)
summary(mod2)
#
ggsave(filename = "Figures/Figure3_neighe001.eps",
       device = cairo_ps,
       dpi=600,
       plot = grid.arrange(p1,p2,
                           heights=c(unit(70, c("mm")),
                                     unit(50,c("mm")))),
       height=120,width=82,unit="mm")#############

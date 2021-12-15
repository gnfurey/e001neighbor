#
library(tidyverse)
library(mgcv)
library(sp)
library(viridis)
library(spdep)
library(gstat)
#
source("helperfunctions_e001Cneigh.R")
e001 <- read_csv("2021-06-25-e001.csv")
e001$NTrt <- as.factor(e001$NTrt)
e001$NTrt <- factor(e001$NTrt,levels(e001$NTrt)[c(9,1,2:8)])
#
meansr <- e001 %>% 
  group_by(Plot,NAdd,ln.NAtm.NAdd,row,col,NTrt) %>%
  summarise_at(vars(means.Sr.4,sr),
               funs(m=mean,se=my.stand))
#in caption
round(range(meansr$sr_m),1)
##reverse columns to have the plots face north
col2 <- data.frame(col2=c(1,2,3,4,5,6,7,8,9),col=c(9,8,7,6,5,4,3,2,1))
meansr <- left_join(meansr,col2)
row2 <- data.frame(row2=c(1,2,3,4,5,6),row=c(6,5,4,3,2,1))
meansr <- left_join(meansr,row2)
#plot n treatment in 2d space
#
meansr$control <- ifelse(meansr$NTrt=="I-0.00","Control","Treat")
#
meansr$control <- ifelse(meansr$NTrt=="I-0.00","Control","Treat")
meansr$NTrt_let <- gsub(pattern = "\\d|\\.|\\-",replacement = "",x = meansr$NTrt)
meansr$NTrt_let
meansr$plot_let <- paste(meansr$Plot,meansr$NTrt_let,sep="-")
nadd.map <- ggplot(meansr,aes(x=row2,y=col2,fill=NTrt))+
  geom_tile(col="black")+
  geom_text(aes(label=plot_let,col=NAdd<9),size=2.2,
            nudge_x=-0.28,nudge_y=0.4,
            label.padding = unit(0.07, "lines"))+
  scale_color_manual(values=c("Grey","Black"))+
  guides(col="none")+
  scale_fill_viridis(
    discrete = TRUE,
    option="A",begin=0,end=1,direction=-1,
    labels=
      c("I-Control",
        "A-0 N",
        "B-1.02 N",
        "C-2.04 N",
        "D-3.40 N",
        "E-5.44 N",
        "F-9.52 N",
        "G-17.0 N",
        "H-27.2 N"))+
  ggtitle("b")+
  ylab("")+
  xlab("")+
  scale_y_continuous(breaks = 1, expand = c(0,0)) + 
  scale_x_continuous(breaks = 1, expand = c(0,0)) + 
  theme_classic()+
  theme(legend.justification = c(0,0),
        legend.position = "bottom",
        legend.key.height =  unit(0.1, "cm"),
        legend.key.width = unit(0.55,"cm"),
        legend.text = element_text(size=7),
        legend.box.spacing = unit(-.23,"cm"),
        legend.background = element_rect(color="Black"),
        legend.box.margin = margin(t = 0,r = 0,b = 0.2,l = 0,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(r = 2,t = 0.3,b = 0.3,l = 0.3),
        plot.title = element_text(margin = margin(0,0,0,0)),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank())+
  guides(fill=guide_legend(nrow=3))
nadd.map
#get build
p <- ggplot_build(nadd.map)
#get hex
cols <- p$data[[1]]
cols <-  cols %>% mutate(row2=x,col2=y)
labcol <- meansr %>% ungroup() %>% select(NTrt,row2,col2,Plot)
labcol <- left_join(labcol,cols) %>% 
  select(Plot,NTrt,fill)# used these in sketchup for Figure S1
#
#plot sr in 2d
sr.map <- ggplot(meansr,aes(x=row2,y=col2,fill=sr_m))+
  geom_tile(col="black")+
  geom_text(aes(label=round(sr_m,1),col=sr_m>8),size=2.2,
             nudge_x=-0.32,nudge_y=0.4,
             label.padding = unit(0.07, "lines"))+
  scale_color_manual(values=c("Grey","Black"))+
  guides(col="none")+
  scale_fill_viridis(option="A",begin=0,end=1,direction=1,
                     name="Sr",
                     breaks = c(2,4,6,8,12,14),
                     limits=c(3.43,14.4))+
  ggtitle("a")+
  ylab("")+
  xlab("")+
  theme_classic()+
  scale_y_continuous(breaks = 1, expand = c(0,0)) + 
  scale_x_continuous(breaks = 1, expand = c(0,0)) + 
  theme(
        legend.justification = c(0,0),
        legend.position = "bottom",
        legend.key.height =  unit(0.68, "cm"),
        legend.key.width = unit(0.632,"cm"),
        legend.text = element_text(size=8),
        legend.box.spacing = unit(-0.22,"cm"),
        legend.spacing.x = unit(0.2,"cm"), 
        legend.background = element_rect(color="Black",),
        legend.box.margin = margin(t = 0,l = 0,b = 0.2,r = 0,"cm"),
        plot.margin = margin(r = 0.3,t = 0.3,b = 0.3,l = 0.3),
        plot.title = element_text(margin = margin(0,0,0,0)),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank())+
  guides(fill=guide_legend(ncol=5),)
sr.map
#
#
#
o1 <- gridExtra::grid.arrange(sr.map,nadd.map,ncol=2,
                              widths=list(unit(3,"in"),unit(3,"in")))
library(colorBlindness)
cvdPlot(o1)
ggsave(filename = "Figures/Figure2_map.pdf",
       plot = o1,
       device = cairo_ps,
       dpi=600,
       width=6,height=4)
################
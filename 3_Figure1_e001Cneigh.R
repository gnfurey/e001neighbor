source("helperfunctions_e001Cneigh.R")
library(gridExtra)
library(tidyverse)
library(nlme)
library(spdep)
library(viridis)
library(gstat)#see package details for class that is required
#ln.NAtm.NAdd = log(Added N +1)
e001 <- read_csv("2021-06-25-e001.csv")
e001$NTrt <- as.factor(e001$NTrt)
e001$NTrt <- factor(e001$NTrt,levels(e001$NTrt)[c(9,1,2:8)])
#get mean data set
meansr <- e001 %>% 
  filter(Year<2005) %>%
  # filter(Year>1994)%>%
  group_by(Plot,NAdd,ln.NAtm.NAdd,row,col,NTrt) %>%
  summarise_at(vars(means.Sr.4,sr),
               funs(m=mean,se=my.stand))
#get unique number of years
yrs <- e001 %>% 
  # filter(Year>1994)
  filter(Year<2005)

length(unique(yrs$Year))
##
#run simple OLS model
mod1.lm <- lm(sr_m~ln.NAtm.NAdd,data=meansr)
summary(mod1.lm)
#add in neighborhood effect
mod2.lm <- lm(sr_m~ln.NAtm.NAdd+means.Sr.4_m,data=meansr)
summary(mod2.lm)
#get r2 difference
summary(mod1.lm)$r.squared-summary(mod2.lm)$r.squared
#anova using type III
car::Anova(mod2.lm,type = "III")
#gls model with variance model
mod1.ntrt <-  gls(sr_m~ln.NAtm.NAdd,
                  weights=varIdent(form=~1|NTrt),
                  data=meansr)
#
meansr$r.1 <- rstandard(mod1.lm)
meansr$r.2 <- residuals(mod1.ntrt,type="normalized")
#gls with neighborhood
mod2 <- gls(sr_m~means.Sr.4_m+ln.NAtm.NAdd,
            data=meansr)
#gls with neighborhood and variance function
mod2.ntrt <- gls(sr_m~
                   means.Sr.4_m+
                   ln.NAtm.NAdd,
                 method="REML",
                 weights=varIdent(form=~1|NTrt),
                 data=meansr)  
summary(mod2.ntrt)#line 268
car::Anova(mod2.ntrt,type="III")
#
meansr$r.3 <- residuals(mod2,type="normalized")
meansr$r.4 <- residuals(mod2.ntrt,type="normalized")
##################
#base model fit without neighborhood effects
dat.eff.2 <- data.frame(ln.NAtm.NAdd=
                          c(seq(from=min(e001$ln.NAtm.NAdd),
                                to=3.2,by=0.1),3.339))
#
mod2.ntrt.fit <- AICcmodavg::predictSE(mod1.ntrt,newdata=dat.eff.2)
#
dat.eff.2$fit <- mod2.ntrt.fit$fit
dat.eff.2$se <- mod2.ntrt.fit$se.fit
#
dat.eff.2$NAdd <- exp(dat.eff.2$ln.NAtm.NAdd)-1
#
nadd.p <- ggplot(dat.eff.2,aes(x=NAdd,y=fit))+
  # ylab(expression(paste("Number of Species","(",mu["[1995-2004]"],")")))+
  ylab("Species Richness")+
  
  # xlab(expression(paste("log"["e"],"(", "Added Nitrogen + Deposition ","(g" %.% "m"^-2,")",")")))+
  xlab(expression(paste("Added Nitrogen ","(g" %.% "m"^-2 %.% "yr"^-1,")")))+
  scale_y_continuous(breaks=c(seq(from=2,to=15,by=3)),limits=c(2,15))+
  scale_x_continuous(breaks=round(unique(e001$NAdd),1),
                     labels=round(unique(e001$NAdd),1))+
  theme_classic()+
  theme(axis.title = element_text(size=8),
        axis.text  = element_text(size=6),
        panel.grid.minor = element_blank(),
        plot.title = element_text(margin = margin(0,0,0,0)))+
  ggtitle("a")+
  geom_point(aes(x=NAdd,y=sr_m),data=meansr,size=1.5,shape=21,
             col="Black",fill="Grey")+
  geom_line()+
  geom_line(aes(y=fit+se),linetype=2)+
  geom_line(aes(y=fit-se),linetype=2)
nadd.p
#######3
#run a regression on residuals with neighborhood effects
mod.r <- lm(r.2~means.Sr.4_m,data=meansr)
mod.r.conf <- predict(mod.r,interval = "confidence")
colnames(mod.r.conf) <- c("fit.mod.r","lwr.mod.r","upr.mod.r")
meansr <- as.data.frame(meansr)
meansr <- cbind(meansr,mod.r.conf)
#
mean.resid.p<- ggplot(meansr,aes(x=means.Sr.4_m,y=r.2,label=Plot))+
  geom_hline(yintercept = 0)+
  geom_point(size=2,shape=21,
             fill="Grey",
             col="Black")+
  ylab("Normalized Residuals \n SR ~ Ln[N+D]")+
  scale_y_continuous(breaks=c(-2.5,-1.5,-0.5,0,0.5,1.5,2.5,3.5,4.5))+
  # xlab(expression(paste("Neighborhood Number of Species ","(",mu["[1995,2004]"],")")))+
  xlab("Neighborhood Species Richness")+
  ggtitle("b")+
  theme_classic()+
  theme(title = element_text(size=8),
        plot.title = element_text(margin = margin(0,0,0,0)))+
  geom_line(aes(y=fitted(mod.r)))+
  geom_line(aes(y=lwr.mod.r),linetype=2)+
  geom_line(aes(y=upr.mod.r),linetype=2)
mean.resid.p
########3
#test for autocorrelation
#get row/column associations
coor1 <- meansr %>% select(row,col)
#create neighborhood matrix
neighbor <- cell2nb(nrow=9,ncol = 6,type = "rook",torus = FALSE)
#plot neighborhood
plot(neighbor,coords=coor1)
#
spa <- function(x,rand){
  x <- as.vector(x)
  print(shapiro.test(x))
  #test for autocorrelation using Moran's I
  spcor <- sp.correlogram(neighbor, x, order = 3,
                          method = "I",style=spatweight,
                          zero.policy = TRUE,randomisation = rand)#assumption of normality
  out <- as.data.frame(print(spcor))
  out$dist <- 1:nrow(out)
  out$se <- 2*sqrt(spcor$res[,3])#copied from documentation of plot correlogram
  return(out)
}
spatweight="S"#select edge correcting weights
out1 <- spa(meansr$r.2,FALSE)
out1$mod <- "Ln[N+D]"
out1
#
out2 <- spa(meansr$r.4,FALSE)
out2$mod <- "Ln[N+D]+NeighborSR"
out2
#
out <- rbind(out1,out2)
out$`Pr(I) two sided` <-p.adjust(out$`Pr(I) two sided`,"fdr")
out
#test with added nitrogen
nadd <- spa(meansr$NAdd,rand=TRUE)
nadd
##
c1 <- ggplot(out,aes(x=dist,y=estimate,fill=mod,col=mod))+
  geom_line()+
  geom_hline(yintercept = 0)+
  scale_x_continuous(breaks=c(1,2,3))+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  guides(color=FALSE)+
  geom_errorbar(aes(ymax=estimate+se,
                    ymin=estimate-se,
                    width=0.1,col=mod),
                position = position_dodge(width=0.14))+
  ylab("Moran's I")+
  xlab("Distance (Number of Plots Away)")+
  theme_classic()+
  annotate("text",x=0.97,y=0.55,label="**",col="#1b9e77",size=7)+
  theme(title = element_text(size=8),
        legend.position = c(0.70,0.85),
        plot.title = element_text(margin = margin(0,0,0,0)),
        legend.title = element_blank(),
        legend.key.height =  unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.text=element_text(size=6),
        legend.margin=margin(c(0,0,0,0)),
        legend.box.margin=margin(0,0,0,0))+
  geom_point(shape=21,col="Black",size=3,position = position_dodge(width=0.14))+
  scale_y_continuous(limits=c(-0.3,0.6))+
  ggtitle("c")
c1
############
#run an additional tests 
#not included in the paper
col.W<- nb2listw(neighbor,style="S")#create neighborhood weights
#run moran's I test
moran.test(as.vector(meansr$r.1),listw = col.W,randomisation = FALSE)
t1 <- lm.morantest(mod1.lm,listw=col.W,alternative = "two.sided",resfun=weighted.residuals)
t1
t2 <- lm.morantest(mod2.lm,listw=col.W,alternative = "two.sided",resfun=weighted.residuals)
t2#neighborhood effect addresses autocorrelation
#############
#
# #impose residuals onto a 2d map
# r.map <- ggplot(meansr,aes(x=row,y=col,fill=r.2,label=Plot))+geom_tile(col="black")+
#   geom_text()+
#   scale_fill_viridis(begin=0,name="Resid")+
#   ylab("")+xlab("")+
#   geom_text()+
#   scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9))+
#   theme_bw()+
#   theme(legend.justification = c(0,-1),
#         legend.position = "bottom",
#         legend.key.height =  unit(0.1, "cm"),
#         legend.key.width = unit(1.3,"cm"),
#         legend.text = element_text(size=8),
#         legend.box.spacing = unit(-0.8,"cm"),
#         legend.box.margin = margin(0,0,0,0,"cm"))+
#   scale_x_continuous(breaks=c(1,2,3,4,5,6))+
#   ggtitle("B")
# # coord_fixed(ratio=1)
# r.map
#
nadd.p <- nadd.p+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=8),
        plot.title = element_text(size=10,family = "Helvetica"))
mean.resid.p <- mean.resid.p+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=8),
        plot.title = element_text(size=10,family = "Helvetica"))
c1 <- c1+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=8),
        plot.title = element_text(size=10,family = "Helvetica"))

ggsave(filename = "Figures/Figure1_neighe001.eps",
       device = cairo_ps,
       dpi=600,
       plot = arrangeGrob(nadd.p,mean.resid.p,c1),
       height=140,width=82,unit="mm")


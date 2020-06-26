source("neighborhood_functions.R")
library(gridExtra)
library(tidyverse)
library(nlme)
library(spdep)
library(viridis)
library(gstat)#see package details for class that is required
#ln.NAtm.NAdd = log(Added N +1)
e001 <- read_csv("2020-02-03-e001.csv")
e001$NTrt <- as.factor(e001$NTrt)
e001$NTrt <- factor(e001$NTrt,levels(e001$NTrt)[c(9,1,2:8)])
meansr <- e001 %>% 
  filter(Year<2005) %>%
  filter(Year>1994)%>%
  group_by(Plot,NAdd,ln.NAtm.NAdd,row,col,NTrt,edgeeffect) %>%
  summarise_at(vars(means.Sr.4,sr,means.NAdd.4,Evenness,ShanWinr,EtoH),
               funs(m=mean,se=my.stand))
yrs <- e001 %>% 
  filter(Year<2005) %>%
  filter(Year>1994)
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
summary(mod2.ntrt)#line 192
car::Anova(mod2.ntrt,type="III")
#
meansr$r.3 <- residuals(mod2,type="normalized")
meansr$r.4 <- residuals(mod2.ntrt,type="normalized")
#######
#compare edge effect
#add in edge effect
mod2.ntrt.edge <- update(mod2.ntrt,.~.+edgeeffect)
summary(mod2.ntrt.edge)
car::Anova(mod2.ntrt.edge,type="III")
#
mod2.ntrt.ml <- update(mod2.ntrt,method="ML")
mod2.ntrt.edge.ml <- update(mod2.ntrt.edge,method="ML")
anova(mod2.ntrt.ml,mod2.ntrt.edge.ml)
aicfun(mod2.ntrt.ml,mod2.ntrt.edge.ml)
###########
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
dat.eff.2$ln.NAtm.NAdd
dat.eff.2$NAdd <- exp(dat.eff.2$ln.NAtm.NAdd)-1
#
nadd.p <- ggplot(dat.eff.2,aes(x=NAdd,y=fit))+
  geom_line(aes(y=fit+se),linetype=2,alpha=0.2)+
  geom_line(aes(y=fit-se),linetype=2,alpha=0.2)+
  ylab(expression(paste("Number of Species","(",mu["[1995-2004]"],")")))+
  # xlab(expression(paste("log"["e"],"(", "Added Nitrogen + Deposition ","(g" %.% "m"^-2,")",")")))+
  xlab(expression(paste("Added Nitrogen ","(g" %.% "m"^-2 %.% "yr"^-1,")")))+
  scale_y_continuous(breaks=c(seq(from=2,to=15,by=3)),limits=c(2,15))+
  scale_x_continuous(breaks=round(unique(e001$NAdd),1),
                     labels=round(unique(e001$NAdd),1))+
  theme_bw()+
  theme(title = element_text(size=8),
        panel.grid.minor = element_blank(),
        plot.title = element_text(margin = margin(0,0,0,0)))+
  ggtitle("A")+
  geom_line()+
  geom_point(aes(x=NAdd,y=sr_m),data=meansr,size=1,shape=21)
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
  geom_point(size=1,shape=21)+geom_hline(yintercept = 0,linetype=2)+
  ylab("Normalized Residuals \n NumSp ~ f(Ln[N+D])")+
  scale_y_continuous(breaks=c(-2.5,-1.5,-0.5,0,0.5,1.5,2.5,3.5,4.5))+
  xlab(expression(paste("Neighborhood Number of Species ","(",mu["[1995,2004]"],")")))+
  ggtitle("B")+
  theme_bw()+
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
spa <- function(x){
  x <- as.vector(x)
  #test for autocorrelation using Moran's I
  spcor <- sp.correlogram(neighbor, x, order = 3,
                          method = "I",style=spatweight,
                          zero.policy = TRUE,randomisation = FALSE)
  out <- as.data.frame(print(spcor))
  out$dist <- 1:nrow(out)
  out$se <- 2*sqrt(spcor$res[,3])#copied from documentation of plot correlogram
  return(out)
}
spatweight="S"#select edge correcting weights
out1 <- spa(meansr$r.2)
out1$mod <- "Ln[N+D]"
#
out2 <- spa(meansr$r.4)
out2$mod <- "Ln[N+D]+NeighborNumSp"
out <- rbind(out1,out2)
#
spa(meansr$NAdd)
c1 <- ggplot(out,aes(x=dist,y=estimate,fill=mod,col=mod))+
  geom_line()+
  geom_hline(yintercept = 0)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  guides(color=FALSE)+
  geom_errorbar(aes(ymax=estimate+se,
                    ymin=estimate-se,
                    width=0.05,col=mod),alpha=0.4)+
  ylab("Moran's I")+
  xlab("Distance (Number of Plots Away)")+
  theme_bw()+
  annotate("text",x=1,y=0.5,label="**",col="#1b9e77",size=7)+
  theme(title = element_text(size=8),
        legend.position = c(0.70,0.85),
        plot.title = element_text(margin = margin(0,0,0,0)),
        legend.title = element_blank(),
        legend.key.height =  unit(0.1, "cm"),
        legend.box.spacing = unit(0.1,"cm"),
        legend.text=element_text(size=6),
        legend.margin=margin(c(0,0,0,0)),
        legend.box.margin=margin(0,0,0,0))+
  geom_point(shape=21,col="Black",size=3)+
  scale_y_continuous(limits=c(-0.3,0.6))+
  ggtitle("C")
c1
############
#run an additional tests 
col.W<- nb2listw(neighbor,style="S")#create neighborhood weights
#run moran's I test
moran.test(as.vector(meansr$r.1),listw = col.W,randomisation = FALSE)
t1 <- lm.morantest(mod1.lm,listw=col.W,alternative = "two.sided",resfun=weighted.residuals)
t1
t2 <- lm.morantest(mod2.lm,listw=col.W,alternative = "two.sided",resfun=weighted.residuals)
t2#neighborhood effect addresses autocorrelation
#############
#test for different correlation structures
mod3.exp <- update(mod1.ntrt,correlation=corExp(form=~row+col,nugget = TRUE))
meansr$r.5 <- residuals(mod3.exp,type="normalized")
mod3.Sph <- update(mod1.ntrt,correlation=corSpher(form=~row+col,nugget = TRUE)) 
mod3.Rat <- update(mod1.ntrt,correlation=corRatio(form=~row+col,nugget = TRUE))
mod3.Gau <- update(mod1.ntrt,correlation=corGaus(form=~row+col,nugget = TRUE))
####
col.W<- nb2listw(neighbor,style="S")#create neighborhodo weights
#run moran's I test
moran.test(as.vector(meansr$r.5),listw = col.W,randomisation = FALSE)
####
anova(mod1.ntrt,mod3.exp)
aicfun(mod1 = mod1.ntrt,mod2 = mod3.exp)
summary(mod3.exp)
AIC(mod3.exp,mod3.Sph,mod3.Rat,mod3.Gau)
#
#set coordinates
coordinates(meansr) <- c("row","col")
#create a variogram to data
v <- variogram(r.2~1,meansr,cutoff=5,width=1)
#
vario <- as.data.frame(v)
#fit a variogram
v1 <- fit.variogram(v,model = vgm("Exp"))
v1
#
fitted.vario <- variogramLine(v1,maxdist = 5)
#make a plot
vario1plot <- ggplot(vario,aes(x=dist,y=gamma))+geom_point()+
  geom_line(aes(x=dist,y=gamma),data=fitted.vario)+
  xlab("Distance (Number of Plots)")+
  theme_bw()+
  ylab("Semivariance \n NumSp~f(Ln[N+D])")+ggtitle("A")
vario1plot
#
meansr <- as.data.frame(meansr)
#impose residuals onto a 2d map
r.map <- ggplot(meansr,aes(x=row,y=col,fill=r.2,label=Plot))+geom_tile(col="black")+
  geom_text()+
  scale_fill_viridis(begin=0,name="Resid")+
  ylab("")+xlab("")+
  geom_text()+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9))+
  theme_bw()+
  theme(legend.justification = c(0,-1),
        legend.position = "bottom",
        legend.key.height =  unit(0.1, "cm"),
        legend.key.width = unit(1.3,"cm"),
        legend.text = element_text(size=8),
        legend.box.spacing = unit(-0.8,"cm"),
        legend.box.margin = margin(0,0,0,0,"cm"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6))+
  ggtitle("B")
# coord_fixed(ratio=1)
r.map
#create a map of nitrogen treatments
un <- e001 %>%
  filter(Year==1982) %>%
  select(Plot,NTrt,NAdd,Year,row,col)
nadd.map <- ggplot(un,aes(x=row,y=col,fill=NTrt,label=Plot))+
  geom_tile(col="black")+
  scale_fill_viridis(discrete=TRUE,option="A",direction=-1,begin=0.3,end = 1)+
  ggtitle("C")+
  ylab("")+xlab("")+
  geom_text()+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6))+
  theme_bw()+
  theme(legend.justification = c(0,-1),
        legend.position = "bottom",
        legend.key.height =  unit(0.05, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.text = element_text(size=6),
        legend.box.spacing = unit(-0.8,"cm"),
        legend.box.margin = margin(0,0,0,0,"cm"))
nadd.map
#
nadd.p <- nadd.p+theme(axis.title = element_text(size=7),
                       axis.text  = element_text(size=7))
mean.resid.p <- mean.resid.p+theme(axis.title = element_text(size=7),
                                   axis.text  = element_text(size=7))            
c1 <- c1+theme(axis.title = element_text(size=7),
               axis.text = element_text(size=7))
#
#
ggsave(filename = "Figures/Figure1_neighe001.eps",
       device = cairo_ps,
       dpi=600,
       plot = arrangeGrob(nadd.p,mean.resid.p,c1),
       height=140,width=82,unit="mm")
ggsave(filename = "Figures/SupplementalFigure1.eps",
       plot = grid.arrange(vario1plot,r.map,nadd.map),
       dpi=600,
       height=6,width=4,unit="in")

source("neighborhood_functions.R")
library(gridExtra)
library(tidyverse)
library(nlme)
library(spdep)
library(viridis)
library(gstat)#see package details for class that is required
library(ggrepel)
library(grid)
#ln.NAtm.NAdd = log(Added N +1)
e001 <- read_csv("2020-02-03-e001.csv")
e001.l <- read.csv("2020-02-03-e001.l.csv")
"%ni%" <- Negate("%in%")
e001.l$func <- get_funcgroup(e001.l$Specid)# get functional groups and remove woody plants
ranks <- e001.l %>%
  filter(func!="W") %>%
  group_by(Specid)  %>% 
  filter(Year<2005) %>% 
  summarise(Biomass=mean(Biomass)) %>% filter(Biomass>1) %>% 
  filter(Specid %ni% c("Misclitt","Mosslich","Calacana","Heligiga","Rumeacet")) %>% 
  arrange(desc(Biomass))
ranks$rank <- seq(from=1,to=NROW(ranks))
ranksshort <- ranks 
abundance <- ggplot(ranks,aes(x=rank,y=Biomass,label=Specid))+geom_text_repel()+geom_point()+
  theme_bw()+
  scale_x_continuous(breaks=seq(from=0,to=20))
abundance
splist <- ranks$Specid[1:11]
splist <- splist[-9]#remove Apoccann as it is mostly zeros with a few large values and has a skewed mean
###################
#second function to print out the params
regfun <- function(species){
  # species="Euphcoro"
  species <- species
  # print(species)
  dat=e001.l
  response <- "Biomass"
  dat <- dat %>% filter(Year<2005)
  #get species neighborhood 
  neigh <- paste(species,"neighbor",sep=".")
  dat$NAdd.fac <- as.factor(dat$NAdd)
  dat <- dat %>% filter(Specid==species)
  #calculate mean values across time series 
  dat1 <- dat %>% filter(Specid==species) %>% 
    group_by(Plot,NAdd,NTrt,ln.NAtm.NAdd,NAdd.fac,edgeeffect) %>% 
    summarise(Biomass=mean(Biomass),
              neighbor=mean(!!sym(neigh)))
  ggplot(dat1,aes(x=NAdd,y=Biomass))+
    geom_point()
  #fit three models with nitrogen,neighborhood,and their interaction
  tmp <- paste(response,"~ ",sep=" ")
  parms <- c("ln.NAtm.NAdd")
  form =as.formula(paste(tmp, paste(parms, collapse= "+")))
  #
  parms1 <- c("ln.NAtm.NAdd","neighbor")
  form1 =as.formula(paste(tmp, paste(parms1, collapse= "+")))
  #
  int <- paste("neighbor","ln.NAtm.NAdd",sep="*")
  parms2 <- c(int)
  form2 =as.formula(paste(tmp, paste(parms2, collapse= "+")))
  #
  #https://stackoverflow.com/questions/7666807/anova-test-fails-on-lme-fits-created-with-pasted-formula
  mod1 <- do.call("gls", args = list(form,
                                     method="ML",
                                     weights=varIdent(form=~1|NTrt),data=dat1))
  mod2 <- do.call("gls", args = list(form1,
                                     method="ML",
                                     weights=varIdent(form=~1|NTrt),data=dat1))
  tmp <- gls(form,weights = varIdent(form=~1|NTrt),
             method="ML",data=dat1)
  coef(tmp)
  mod3 <- do.call("gls", args = list(form2,
                                     method="ML",
                                     weights=varIdent(form=~1|NTrt),data=dat1))
  #all lines below pull various statistics into a table
  p.aic <- round(diff(anova(mod2,mod3)$AIC),1)
  p.lrr <- round(anova(mod2,mod3)$L.Ratio[2],1)
  p.int <- anova(mod2,mod3)$`p-value`[2]
  #
  p.neigh <- anova(mod1,mod2)$`p-value`[2]#neighborhood effect
  mod1 <- update(mod1,method="REML")
  mod2 <- update(mod2,method="REML")
  mod3 <- update(mod3,method="REML")
  #
  out1 <- as.data.frame(summary(mod1)$tTable)#nitrogen
  out1 <- out1[rownames(out1)=="ln.NAtm.NAdd",c(1,2,4)]
  colnames(out1) <- c("Coef_N","Std.Err_N","p.value_N")
  out1$parm1 <- rownames(out1)
  out2 <- as.data.frame(summary(mod2)$tTable)
  out2.N <- out2[rownames(out2)=="ln.NAtm.NAdd",c(1,2,4)]
  out2.neigh <- out2[rownames(out2)=="neighbor",c(1,2)]#neighbor
  out2.neigh$`p-value` <- p.neigh#use likelihood ratio test p-value
  colnames(out2.N) <- c("Coef_N_neigh",
                        "Std.Err_N_neigh",
                        "p.value_N_Neigh")
  out2.N$parm2 <- "ln.NAtm_NAdd_neigh"
  colnames(out2.neigh) <- c("Coef_neigh","Std.Err_neigh","p.value_Neigh")
  out2.neigh$parm3 <- "Neigh"
  out3 <- data.frame(parm4="N*neighbor",p.aic,p.lrr,p.int)
  colnames(out1)
  colnames(out2.N)
  colnames(out2.neigh)
  colnames(out3)
  outs <- cbind(out1,out2.N,out2.neigh,out3)
  #
  outs <- outs %>% select(parm1,Coef_N:p.value_N,
                          parm2,Coef_N_neigh:p.value_N_Neigh,
                          parm3,Coef_neigh:p.value_Neigh,
                          parm4,everything())
  
}
#
names(splist) <- set_names(splist)
coef1 <- map_df(.x = splist,.f = regfun,.id = "Species")#run on species
#
#
# modtestfun <- function(spec){
#   # spec="Poaprate"
#   neigh <- paste(spec,"neighbor",sep=".")
#   #calculate mean values across time series 
#   test1 <- e001.l %>% filter(Specid==spec) %>%
#     group_by(Plot,NAdd,NTrt,ln.NAtm.NAdd,edgeeffect) %>% 
#     summarise(Biomass=mean(Biomass),
#               neighbor=mean(!!sym(neigh)))
#   test_mod1 <- gls(Biomass~ln.NAtm.NAdd,weights = varIdent(form=~1|NTrt),
#                    method="REML",data=test1)
#   fun_out <- coef1 %>% filter(Species==spec)
#   #
#   summary(test_mod1)
#   Coef_N <- fun_out$Coef_N==summary(test_mod1)$tTable[2,1]#test 1 passed Nitrogen
#   Coef_N
#   Std.Err_N <- fun_out$Std.Err_N==summary(test_mod1)$tTable[2,2]
#   p.value_N <- fun_out$p.value_N==summary(test_mod1)$tTable[2,4]
#   #
#   test_mod2 <- gls(Biomass~ln.NAtm.NAdd+neighbor,weights = varIdent(form=~1|NTrt),
#                    method="REML",data=test1)
#   #
#   Coef_neigh <- fun_out$Coef_neigh==summary(test_mod2)$tTable[3,1]#passed
#   Std.Err_neigh <- fun_out$Std.Err_neigh==summary(test_mod2)$tTable[3,2]#passed
#   #
#   test_mod1_ml <- update(test_mod1,method="ML")#passed
#   test_mod2_ml <- update(test_mod2,method="ML")#passed
#   p.value_Neigh <- fun_out$p.value_Neigh==anova(test_mod1_ml,test_mod2_ml)$`p-value`[2]#passed
#   #
#   test_mod3 <- gls(Biomass~ln.NAtm.NAdd*neighbor,weights = varIdent(form=~1|NTrt),
#                    method="REML",data=test1)
#   test_mod3_ml <- update(test_mod3,method="ML")
#   p.int <- round(fun_out$p.int,3)==round(anova(test_mod3_ml,test_mod2_ml)$`p-value`[2],3)
#   p.aic <- fun_out$p.aic==round(diff(anova(test_mod2_ml,test_mod3_ml)$AIC),1)
#   p.lrr <- fun_out$p.lrr==round(anova(test_mod2_ml,test_mod3_ml)$L.Ratio[2],1)
#   out <- c(
#     Coef_N,Std.Err_N,p.value_N,
#     Coef_neigh,Std.Err_neigh,p.value_Neigh,
#     p.int,p.aic,p.lrr)
#   out
#   return(all(out)==TRUE)
# }
# #
# tests_spec <- map_df(.x = splist,.f = modtestfun,.id = "Species")#run on species
# all(apply(tests_spec[2:11], 2, isTRUE))==TRUE#passed
# #
#adjust p-values 
coef2 <- coef1
coef2$p.value_N <- p.adjust(coef2$p.value_N,method="fdr")
coef2$p.value_N_Neigh <- p.adjust(coef2$p.value_N_Neigh,method="fdr")
coef2$p.value_Neigh <- p.adjust(coef2$p.value_Neigh,method="fdr")
coef2$p.int <- p.adjust(coef2$p.int,method="fdr")
#
colnames(coef2)
#round coef and SE
coef2[,c(3,4,7,8,11,12)] <- apply(coef2[,c(3,4,7,8,11,12)],
                                  2,function(x)round(as.numeric(as.character(x)),2))
#
coef2[,c(5,9,13,17)] <- apply(coef2[,c(5,9,13,17)],
                              2,function(x)round(as.numeric(as.character(x)),2))
#
coef2[,c(5,9,13,17)] <- apply(coef2[,c(5,9,13,17)],
                              2,
                              function(x)ifelse(x==0,
                                                "<0.001",
                                                ifelse(x<0.01,
                                                       "<0.01",
                                                       x)))
#
coef2 <- coef2 %>% arrange(desc(Coef_N))
coef2$Species <- factor(coef2$Species,levels=unique(coef2$Species))
coef2$Species <- get_species(coef2$Species)
#
coef_nadd <- coef2 %>% unite(nadd,c(Coef_N,Std.Err_N),sep="±") %>%
  select(Species,parm1,nadd,p.value_N) %>% 
  select(-parm1)
#
coef_neigh <- coef2 %>%
  unite(nadd,c(Coef_N_neigh,Std.Err_N_neigh),sep="±") %>%
  unite(neigh,c(Coef_neigh,Std.Err_neigh),sep="±") %>%
  select(Species,parm2,nadd,p.value_N_Neigh,parm3,neigh,p.value_Neigh)
#
coef_int <- coef2 %>% 
  select(Species,parm4,p.lrr,p.aic,p.int) %>% 
  arrange(p.int)
#
colnames(coef_nadd)[2:3] <- c("Effect Size ± SE","p-value")
#
coef_neigh$parm2 <- "ln[N+D]"
coef_neigh$parm3 <- "Neighborhood NumSp"
colnames(coef_neigh)[2:7] <- c("Parameter","Effect Size ± SE",
                               "p-value",
                               "Parameter",
                               "Effect Size ± SE",
                               "p-value")
#
colnames(coef_int)[2:5] <- c("Parameter","Log Response Ratio","AIC","p-value")
#
write.csv(x = coef_nadd,file = "Tables/SuppTable5_neighe001.csv",
          row.names=FALSE)
write.csv(x = coef_neigh,file = "Tables/SuppTable6_neighe001.csv",
          row.names=FALSE)
write.csv(x = coef_int,file = "Tables/SuppTable7_neighe001.csv",
          row.names=FALSE)

#
#######
#a function to test if any species have an edge effect 
edgefun <- function(x){
  # x <- "Agrorepe"
  # print(x)
  dat=e001.l
  response <- "Biomass"
  dat <- dat %>% filter(Year<2005)
  neigh <- paste(x,"neighbor",sep=".")
  dat$NAdd.fac <- as.factor(dat$NAdd)
  dat <- dat %>% filter(Specid==x)
  dat1 <- dat %>% filter(Specid==x) %>% 
    group_by(Plot,NAdd,NTrt,ln.NAtm.NAdd,NAdd.fac,edgeeffect) %>% 
    summarise(Biomass=mean(Biomass),
              neighbor=mean(!!sym(neigh)))
  tmp <- paste(response,"~ ",sep=" ")
  int <- paste("ln.NAtm.NAdd","neighbor",sep="*")
  parms1 <- c(int,"edgeeffect")
  form1 =as.formula(paste(tmp, paste(int, collapse= "+")))
  form2 =as.formula(paste(tmp, paste(parms1, collapse= "+")))
  mod1 <- do.call("gls", args = list(form1,
                                     method="ML",
                                     weights=varIdent(form=~1|NTrt),data=dat1))
  mod2 <- do.call("gls", args = list(form2,
                                     method="ML",
                                     weights=varIdent(form=~1|NTrt),data=dat1))
  
  #
  pval= anova(mod1,mod2)$`p-value`[2]
  lrr = anova(mod1,mod2)$L.Ratio[2]
  aic = round(diff(anova(mod1,mod2)$AIC),1)
  out1 <- data.frame(aic=aic,lrr=lrr,pval=pval)
  return(out1)
}
spp <- splist
names(spp) <- set_names(spp)
#
edge <- map_df(.x = spp,.f = edgefun,.id = "Species")
edge$`p-value` <- p.adjust(edge$pval,method="fdr")
#
coef2 <- coef2 %>% arrange(p.int)
########
edge$`p-value` <- round(edge$`p-value`,2)
edge$aic <- round(edge$aic,2)
edge$lrr <- round(edge$lrr,2)
edge$pval <- round(edge$pval,2)
edge$Species <- get_species(edge$Species)
out <- edge %>% arrange(`p-value`)
# write.csv(file = "Figures/Table2_neighe001.csv",coef2,row.names = FALSE)
write.csv(file = "Tables/SuppTable8_neighe001.csv",out,row.names = FALSE)
#
e001.l$func <- get_funcgroup(e001.l$Specid)# get functional groups and remove woody plants
ranks <- e001.l %>%
  filter(func!="W") %>%
  group_by(Specid)  %>% 
  filter(Year<2005) %>% 
  summarise(Biomass=mean(Biomass)) %>% filter(Biomass>1) %>% 
  filter(Specid %ni% c("Misclitt","Mosslich","Calacana","Heligiga","Rumeacet")) %>% 
  arrange(desc(Biomass))
ranks$rank <- seq(from=1,to=NROW(ranks))
ranksshort <- ranks 
abundance <- ggplot(ranks,aes(x=rank,y=Biomass,label=Specid))+geom_text_repel()+geom_point()+
  theme_bw()+
  scale_x_continuous(breaks=seq(from=0,to=20))
abundance
splist <- ranks$Specid[1:11]
splist <- splist[-9]#remove Apoccann as it is mostly zeros with a few large 
#create a function to take a given species
#then run a regression with the neighborhood effects
#plot with only low and middle levels of nitrogen
plots <- function(species){
  # species="Agrorepe"
  dat=e001.l
  response <- "Biomass"
  neigh <- paste(species,"neighbor",sep=".") # make neighborhood param
  #use mean of whole times series
  dat1 <- dat %>% filter(Specid==species) %>% 
    group_by(Plot,NAdd,NTrt,ln.NAtm.NAdd,Specid,edgeeffect) %>% 
    summarise(Biomass=mean(Biomass),
              neighbor=mean(!!sym(neigh)))
  #create interaction 
  int <- paste("neighbor","ln.NAtm.NAdd",sep="*")
  #get parms
  parms <- c(int)
  tmp <- paste(response,"~ ",sep=" ")
  #formula
  form =as.formula(paste(tmp, paste(parms, collapse= "+")))
  #use eval and subs due to error with gls
  #add in variance function
  mod1 <- eval(substitute(gls(form,
                              weights=varIdent(form=~1|NTrt),data=dat1), list(form=form)))
  #get standard error 
  fits <-  AICcmodavg::predictSE(mod1,newdata=dat1,se=TRUE)
  dat1$fit <- fits$fit
  dat1$se <- fits$se.fit
  #use only low and middle nitrogen 
  # unique(dat$NAdd)
  dat1 <- dat1 %>% filter(NAdd %in% c(0,9.52))
  title1 <- get_species(unique(dat1$Specid))
  dat1$NAdd_chr <- paste(dat1$NAdd,"Added N")
  #
  dat.mean <- dat1 %>% group_by(ln.NAtm.NAdd,NAdd_chr) %>% 
    summarise(Biomass=mean(Biomass),
              neighbor=mean(neighbor))
  print(dat.mean)
  #
  p1 <- ggplot(dat1,
               aes(x=neighbor,y=Biomass,fill=as.factor(ln.NAtm.NAdd)))+
    # geom_point(data = dat.mean,shape=3,size=5)+
    geom_line(aes(y=fit,col=as.factor(ln.NAtm.NAdd)),size=0.5)+
    geom_line(aes(y=fit+se,col=as.factor(ln.NAtm.NAdd)),size=0.5,linetype=2,alpha=0.3)+
    geom_line(aes(y=fit-se,col=as.factor(ln.NAtm.NAdd)),size=0.5,linetype=2,alpha=0.3)+
    theme_bw()+
    # facet_wrap(~NAdd_chr)+
    geom_point(shape=21,size=2)+
    ylab("")+xlab("")+
    scale_color_brewer(palette = "Dark2",name="",labels=
                         c("= 0.0 Added N"," = 9.52 Added N"))+
    scale_fill_brewer(palette = "Dark2",name="",labels=
                        c("= 0.0 Added N"," = 9.52 Added N"))+
    # scale_color_brewer(palette = "Dark2",name=expression(paste("N(g"%.%"m"^-2,")")),labels=c(0,9.52))+
    # scale_fill_brewer(palette = "Dark2",name=expression(paste("N(g"%.%"m"^-2,")")),labels=c(0,9.52))+
    guides(col=FALSE)+
    # guides(fill=FALSE)
    ggtitle(title1)+
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
    theme(plot.title = element_text(size=10,margin=margin(0,0,0,0)),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.spacing.x = unit(0.02, 'cm'),
          # axis.text.x = element_text(size=5),
          legend.position = "top",
          plot.margin=margin(t = 0, r = 7, b = 0, l = 0,),
          legend.margin=margin(c(0,0,0,0)),
          legend.box.margin=margin(0,0,0,0),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1)))
  p1
  return(p1)
} 
p1 <- plots(species=splist[1])
p1 <- p1+
  geom_text(aes(x=133,y=50,label="No N"),col="#1b9e77",size=3)+
  geom_text(aes(x=135,y=160,label="High N"),col="#d95f02",size=3)
#create a function to get the legend # I don't remember the stackoverflow for this code
#possibly
#https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
legendfunc<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend.1<-legendfunc(p1)
#
#
splist
p2 <- plots(species=splist[2])+ggtitle("Elymus repens")
p2 <- p2+  
  geom_text(aes(x=165,y=0,label="No N"),col="#1b9e77",size=3)+
  geom_text(aes(x=111,y=110,label="High N"),col="#d95f02",size=3)
#
p3 <- plots(species=splist[6])+
  ggtitle("Symphyotrichum oolentangiense")+
  geom_text(aes(x=19,y=19,label="No N"),col="#1b9e77",size=3)+
  geom_text(aes(x=11.5,y=4,label="High N"),col="#d95f02",size=3)
p3  
p4 <- plots(species=splist[9])
p4 <- p4+
  geom_text(aes(x=16.5,y=14,label="No N"),col="#1b9e77",size=3)+
  geom_text(aes(x=13,y=4,label="High N"),col="#d95f02",size=3 )
##########
bioindlab <- expression(paste("Each Species' Focal Plot Biomass", " ","(","g"%.%"m"^-2, ")"))
neighlab<- expression(paste("Each Species' Neighborhood Biomass", " ","(","g"%.%"m"^-2, ")"))
t <- textGrob(bioindlab,rot=90,x=unit(3,"mm"),y=unit(80,"mm"),
              gp=gpar(fontsize=10))
x <- textGrob(neighlab,x=unit(37,unit="mm"),y=unit(3,"mm"),
              gp=gpar(fontsize=10))

ggsave(filename = "Figures/Figure4_neighe001.eps",
       device = cairo_ps,
       dpi=600,
       plot = grid.arrange(mylegend.1,
                           p1 + theme(legend.position="none",
                                      plot.margin=margin(t = 0,
                                                         r = 6,
                                                         b = 0, 
                                                         l = 0.5,)),
                           p2 + theme(legend.position="none",
                                      plot.margin=margin(t = 0,
                                                         r = 6,
                                                         b = 0, 
                                                         l = 0.5,)),
                           p3 + theme(legend.position="none",
                                      plot.margin=margin(t = 0,
                                                         r = 2,
                                                         b = 0, 
                                                         l = 5,)),
                           p4 + theme(legend.position="none",
                                      plot.margin=margin(t = 0,
                                                         r = 2,
                                                         b = 0, 
                                                         l = 5,)),
                           left=t,bottom=x,
                           heights=c(0.2,1,1,1,1),
                           widths=c(1),
                           nrow=5,
                           ncol=1),
       height=160,width=82,unit="mm")

#################
unique(e001.l$Year)
agro <- e001.l %>% filter(Specid =="Agrorepe") %>% 
  group_by(Plot,NAdd,NTrt) %>% 
  summarise(Biomass=mean(Biomass),
            Agrorepe.neighbor=mean(Agrorepe.neighbor)) %>% 
  filter(NAdd==9.52)
one <- agro %>% filter(Plot %in% c(1,13,31)) %>% 
  group_by(NAdd) %>% summarise(Biomass=mean(Biomass),
                               Agrorepe.neighbor=mean(Agrorepe.neighbor))
two <- agro %>% filter(Plot %in% c(9,24,52)) %>% 
  group_by(NAdd) %>% summarise(Biomass=mean(Biomass),
                               Agrorepe.neighbor=mean(Agrorepe.neighbor))
one#lines 330-344
two#lines 330-344
e001 %>% filter(NAdd==9.52) %>% 
  filter(Plot %in% c(9,24,52)) %>% 
  summarise(Agrorepe=mean(Agrorepe),
            Agrorepe.neighbor=mean(Agrorepe.neighbor))

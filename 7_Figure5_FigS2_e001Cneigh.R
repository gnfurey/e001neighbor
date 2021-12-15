source("helperfunctions_e001Cneigh.R")
library(tidyverse)
library(nlme)
library(viridis)
#read in data
e001.l <- read_csv("2021-06-25-e001.l.csv")
#reorder factor
e001.l$NTrt <- as.factor(e001.l$NTrt)
e001.l$NTrt <- factor(e001.l$NTrt,levels(e001.l$NTrt)[c(9,1,2:8)])
#
#get three dominant species 
splist <- c("Schiscop","Agrorepe","Poaprate")
e001.l <- e001.l %>%
  filter(Specid %in% splist)
# write.csv(x = e001.l,file = "2021-06-25-e001.l.csv")
unique(e001.l$NTrt)
#########
e001.l$Year <- e001.l$Year-1982
#plotting function and model fits
plot_1 <- function(species,title1,c1,c2,ntrt){
  #get species 
  # species= "Poaprate"
  dat=e001.l
  #filter to species
  neigh <- paste(species,"neighbor",sep=".") # make neighborhood param
  dat <- dat %>% filter(Specid==species)
  #
  response <- "sqrt(Biomass)"
  # response <- "Biomass"
  #fit three models with nitrogen,neighborhood,and their interaction
  tmp <- paste(response,"~ ",sep=" ")
  #get parameters
  parms2 <- c("NTrt","Year",neigh)
  #get formula
  form2 =as.formula(paste(tmp, paste(parms2, collapse= "*")))
  #arrange properly for temporal autocorrelation function
  dat <- dat %>% arrange(Plot,Year)
  #########
  # mod2 <- lme(form2,
  #             random=~1|Plot,
  #             method="REML",
  #             correlation=corLin(form=~1|Plot,nugget = TRUE),
  #             control=list(maxIter = 10000,msMaxIter=10000,niterEM=10000,msMaxEval = 10000),
  #             weights=varIdent(form=~1|NTrt),
  #             data=dat)
  #https://stackoverflow.com/questions/7666807/anova-test-fails-on-lme-fits-created-with-pasted-formula
  mod2 <- do.call("lme", args = list(form2,
                                     random=~1|Plot,
                                     method="REML",
                                     correlation=corLin(form=~1|Plot,nugget = TRUE),
                                     control=list(maxIter = 100,
                                                  msMaxIter=100,
                                                  niterEM=100,
                                                  msMaxEval = 100),
                                     weights=varIdent(form=~1|NTrt),
                                     data=dat))
  #get anova with non-sequential results
  out1 <- as.data.frame(anova(mod2,type="marginal"))
  #add in row names
  out1$coef <- rownames(out1)
  #########
  #get N treatments
  NTrt_tmp <- unique(dat$NTrt)
  #get mean neigborhood
  m1 <- mean(pull(dat,neigh))
  #get sd
  sd <- sd(pull(dat,neigh))
  #get high and low neighborhood biomass
  val1 <- 0
  val2 <- m1+sd
  ###########
  library(emmeans)
  #Agrorepe.neighbor 
  #Schiscop.neighbor 
  pfun <- function(ntrt,title1){
    # ntrt="I-0.00"
  parmy <- as.formula(paste("~",paste(parms2, collapse= "*"),sep=""))
  print(parmy)
  ##
  list1 <- list(NTrt= c(ntrt),
                neigh=c(val1,val2),
                Year=seq(from=0,to=22))
  names(list1)[2] <- neigh
  list1
  
  ##
  out <- emmeans(object = mod2, 
           # specs = ~NTrt*Year*Schiscop.neighbor,
           specs = parmy,
           var="Year",
           df=5,#very important n=6 at each NTrt
          type="response",
          sigmaAdjust = TRUE,
          at=list1)
  out
  out1 <- as.data.frame(out)
  colnames(out1) <- str_replace_all(string = colnames(out1),
                                    pattern = neigh,
                                    replacement = "neighvar")
  # out1$neighvar <- c("0","mu+sd")
  out1$NTrt <- as.character(out1$NTrt)
    out1$NTrt <- ifelse(out1$NTrt=="I-0.00","0.00 g N",out1$NTrt)
    out1$NTrt <- ifelse(out1$NTrt=="D-3.40","3.40 g N",out1$NTrt)
    out1$NTrt <- ifelse(out1$NTrt=="F-9.52","9.52 g N",out1$NTrt)

  p1 <- ggplot(out1,aes(x=Year,y=response,
                  col=as.factor(neighvar),
                  fill=as.factor(neighvar)))+
    geom_line(size=1)+
    geom_ribbon(aes(ymin=lower.CL,ymax=upper.CL),
                size=0.5,
                col=NA,
                alpha=0.4)+
    scale_fill_viridis(begin = 0.2,end = 0.8,discrete = TRUE,
                       # labels=c("mu*0.25","mu*2"),
                       labels=c("0","Mean+1SD"),
                        option = "A", name="Neighborhood Abundance",
    )+
    scale_color_viridis(begin = 0.2,end = 0.8,discrete = TRUE,
                        # labels=c("mu*0.25","mu*2"),
                        labels=c("0","Mean+1SD"),
                        option = "A", name="Neighborhood Abundance")+
    guides(
           fill="none")+
    facet_wrap(~NTrt)+
    ggtitle(title1)+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
              strip.text.x = element_text(margin = margin(t =0.07,r = 0,b = 0.07,l = 0, "cm")),
                    axis.text = element_text(size=7),
                    plot.margin = margin(t = 0,r = 3,b = 0,l = 3),
                    plot.title = element_text(size=7,family = "Helvetica",
                                              face="italic",
                                              margin = margin(0,0,0,0)))
  p1
  return(p1)
  }
  #control to get plot
  if(c2=="plot"){
    # ntrt="I-0.00"
    p1 <- pfun(ntrt=ntrt,title1=title1)
    return(p1)
  }
  #control to get table
  if(c2=="tab"){
    return(out1)
  }
  if(c2=="coef"){
    out2 <- as.data.frame(summary(mod2)$tTable)
    out2$coef <- rownames(out2)
    return(out2)
  }
}
#
#get species list
splist <- set_names(splist)
#################
#get anova table for each species see ?map_dfr
tabs_3 <- map_dfr(.x = splist,
                  .f = plot_1,
                  .id = "Species",
                  c1="",c2="coef",
                  ntrt="",title1="")
#get table 
outtab_3 <- tabs_3
#round
colnames(outtab_3)
outtab_3$Value <- 
  ifelse(outtab_3$Species=="Agrorepe",
         round(outtab_3$Value,4),
         round(outtab_3$Value,3))#
outtab_3$Std.Error <- 
  ifelse(outtab_3$Species=="Agrorepe",
        round(outtab_3$Std.Error,4),
        round(outtab_3$Std.Error,3))
        
outtab_3$`t-value` <- 
  ifelse(outtab_3$Species=="Agrorepe",
        round(outtab_3$`t-value`,4),
        round(outtab_3$`t-value`,3))#
outtab_3$`p-value` <- ifelse(outtab_3$`p-value`==0,
                        "<0.001",
                        ifelse(outtab_3$`p-value`<0.01,
                               "<0.01",
                               round(outtab_3$`p-value`,3)))
#
colnames(outtab_3)
outtab_3 <- outtab_3 %>%
  select(coef,Value:`p-value`,Species)
#subset each table
schis <- outtab_3 %>% filter(Species=="Schiscop")
agro <- outtab_3 %>% filter(Species=="Agrorepe")
poa <- outtab_3 %>% filter(Species=="Poaprate")
library(flextable)
library(officer)
#format nice word tables
#
schis$Species <- NULL
schis1 <- flextable(schis) %>% fontsize(size=8)
set_table_properties(schis1, width = 1, layout = "autofit")
schis1 <- font(schis1,fontname = "Times")
doc <- read_docx() %>% 
  body_add_par(value = "Supplemental Table 5", style = "Normal") %>%
  body_add_flextable(value = schis1)
print(doc, target = "Tables/SupplementalTable5_coef.docx")
#
agro$Species <- NULL
#
agro$coef <- str_replace_all(string =agro$coef,pattern = "Agrorepe",
                             "Elymrepe")
#
agro1 <- flextable(agro) %>% fontsize(size=8)
set_table_properties(agro1, width = 1, layout = "autofit")
agro1 <- font(agro1,fontname = "Times")
doc1 <- read_docx() %>% 
  body_add_par(value = "Supplemental Table 7", style = "Normal") %>%
  body_add_flextable(value = agro1)
#
print(doc1, target = "Tables/SupplementalTable7_coef.docx")
#
poa$Species <- NULL
poa1 <- flextable(poa) %>% fontsize(size=8)
set_table_properties(poa1, width = 1, layout = "autofit")
poa1 <- font(poa1,fontname = "Times")
doc1 <- read_docx() %>% 
  body_add_par(value = "Supplemental Table 9", style = "Normal") %>%
  body_add_flextable(value = poa1)
#
print(doc1, target = "Tables/SupplementalTable9_coef.docx")
#
#
tabs_4 <- map_dfr(.x = splist,
                  .f = plot_1,
                  .id = "Species",
                  c1="",c2="tab",
                  ntrt="",title1="")
#get table 
outtab_4 <- tabs_4
#get table 
#round
colnames(outtab_3)
outtab_4$`F-value` <- round(outtab_4$`F-value`,2)
outtab_4$pval <- p.adjust(outtab_4$`p-value`,"fdr")
outtab_4$`p-value` <- NULL
outtab_4$pval <- ifelse(outtab_4$pval==0,
                             "<0.001",
                             ifelse(outtab_4$pval<0.01,
                                    "<0.01",
                                    round(outtab_4$pval,3)))
#
colnames(outtab_4)
outtab_4 <- outtab_4 %>%
  select(coef,numDF:pval,Species)
#subset each table
schis <- outtab_4 %>% filter(Species=="Schiscop")
agro <- outtab_4 %>% filter(Species=="Agrorepe")
poa <- outtab_4 %>% filter(Species=="Poaprate")
library(flextable)
library(officer)
#format nice word tables
#
schis$Species <- NULL
schis1 <- flextable(schis) %>% fontsize(size=8)
set_table_properties(schis1, width = 1, layout = "autofit")
schis1 <- font(schis1,fontname = "Times")
doc <- read_docx() %>% 
  body_add_par(value = "Supplemental Table 4", style = "Normal") %>%
  body_add_flextable(value = schis1)
print(doc, target = "Tables/SupplementalTable4.docx")
#
agro$Species <- NULL
#
agro$coef <- str_replace_all(string =agro$coef,pattern = "Agrorepe",
                             "Elymrepe")
#
agro1 <- flextable(agro) %>% fontsize(size=8)
set_table_properties(agro1, width = 1, layout = "autofit")
agro1 <- font(agro1,fontname = "Times")
doc1 <- read_docx() %>% 
  body_add_par(value = "Supplemental Table 6", style = "Normal") %>%
  body_add_flextable(value = agro1)
#
print(doc1, target = "Tables/SupplementalTable6.docx")
#
poa$Species <- NULL
poa1 <- flextable(poa) %>% fontsize(size=8)
set_table_properties(poa1, width = 1, layout = "autofit")
poa1 <- font(poa1,fontname = "Times")
doc1 <- read_docx() %>% 
  body_add_par(value = "Supplemental Table 8", style = "Normal") %>%
  body_add_flextable(value = poa1)
#
print(doc1, target = "Tables/SupplementalTable8.docx")
#
###########
#make plots
library(gridExtra)
#run function for elymus repends at NTrt 9.52
#I run this function 9 times for each panel in Figure 5
o1 <- plot_1(species="Agrorepe",
             c1="short",c2="plot",
             ntrt="F-9.52",title1="i: Elyus repens")
o1
#
o2 <- plot_1(species="Agrorepe",
             c1="short",c2="plot",
             ntrt="D-3.40",title1="h: Elyus repens")
o2
o3 <- plot_1(species="Agrorepe",
             c1="short",c2="plot",
             ntrt="I-0.00",title1="g: Elyus repens")
o3
#
#get legend as a grob 
#https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
legendfunc<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
#
o1 <- o1+theme(legend.position = "bottom")
o1
#get legend
leg <- legendfunc(o1)
#get scales
ggplot_build(o1)$layout$panel_scales_y
#set proper y scale 
o1 <- o1+
  scale_y_continuous(limits=c(0,325))+
  theme(axis.title = element_blank())+
  guides(color="none")
o2 <- o2+
  scale_y_continuous(limits=c(0,325))+
  theme(axis.title = element_blank())+
  guides(color="none")
o3 <- o3+
  scale_y_continuous(limits=c(0,325))+
  theme(axis.title = element_blank())+
  guides(color="none")
#
o4 <-  plot_1(species="Schiscop",
              c1="short",c2="plot",
              ntrt="I-0.00",title1="a: S. scorparium")
ggplot_build(o4)$layout$panel_scales_y
o5 <-  plot_1(species="Schiscop",
              c1="short",c2="plot",
              ntrt="D-3.40",title1="b: S. scorparium")
o6 <-  plot_1(species="Schiscop",
                c1="short",c2="plot",
                ntrt="F-9.52",title1="c: S. scorparium")
#
o4 <- o4+
  scale_y_continuous(limits=c(0,225))+
  theme(axis.title = element_blank())+
  guides(color="none")
o5 <- o5+
  scale_y_continuous(limits=c(0,225))+
  theme(axis.title = element_blank())+
  guides(color="none")
o6 <- o6+
  scale_y_continuous(limits=c(0,225))+
  theme(axis.title = element_blank())+
  guides(color="none")
o6
#####
o9 <- plot_1(species="Poaprate",
             c1="short",c2="plot",
             ntrt="F-9.52",title1="f: Poa pratensis")
o9
o7 <- plot_1(species="Poaprate",
             c1="short",c2="plot",
             ntrt="I-0.00",title1="d: Poa pratensis")
o7
#
o8 <- plot_1(species="Poaprate",
             c1="short",c2="plot",
             ntrt="D-3.40",title1="e: Poa pratensis")
ggplot_build(o8)$layout$panel_scales_y
o8
#####
o7 <- o7+
  scale_y_continuous(limits=c(0,375))+
  theme(axis.title = element_blank())+
  guides(color="none")
o8 <- o8+
  scale_y_continuous(limits=c(0,375))+
  theme(axis.title = element_blank())+
  guides(color="none")
o9 <- o9+
  scale_y_continuous(limits=c(0,375))+
  theme(axis.title = element_blank())+
  guides(color="none")
#####
#make label 
bioindlan1 <- expression(paste("Species' Fitted Focal Plot Biomass", " ", "(", 
                 "g" %.% "m"^-2, ")"))
library(grid)
#generate figure
#this uses some obscure grid and gridextra code
#why is making paneled figures so hard ? 
#1) arrange all plots with grid.arrange
#2) add the legend as a grob
#3) arrange each species in a grob with three panels
#4) add in a common x and y axis
#5) adjust the heights
#6) save as a eps file 
library(gridExtra)
ggsave(filename = "Figures/Figure5_neighe001.eps",
       device = cairo_ps,
       dpi=600,
       plot=grid.arrange(
         leg,
         arrangeGrob(o4,o5,o6,ncol=3),
         arrangeGrob(o7,o8,o9,ncol=3),
         arrangeGrob(o3,o2,o1,ncol=3),
         left=
           textGrob(bioindlan1,rot=90,
                    gp=gpar(fontsize=8,fontface="plain")),
         bottom=
           textGrob("Year",
                    gp=gpar(fontsize=8,fontface="plain")),
         heights=list(unit(8,"mm"),
                      unit(37,"mm"),
                      unit(37,"mm"),
                      unit(37,"mm")),
         padding=unit(1,"mm")),
       
       height=100,width=150,unit="mm")
#get species means for the last 10 years
meandat <- e001.l %>% 
  filter(Specid %in% splist) %>% 
  filter(Year>1994-1982) %>% 
  group_by(Specid,NTrt) %>% 
  summarise(
    se=my.stand(Biomass),
    Biomass=mean(Biomass))
#
# length(unique(meandat$Year))
#make a factor 
match <- data.frame(Specid=c("Agrorepe","Schiscop","Poaprate"),
                    Species=c("E. repens","S. scoparium","P. pratensis"))
#match in the factor
meandat <- left_join(meandat,match)
#make the plot
b1 <- ggplot(meandat,aes(x=NTrt,y=Biomass,col=Species,fill=Species))+
  geom_line(aes(group=Species))+
  ylab(bioindlab)+
  geom_errorbar(aes(ymax=Biomass+se,ymin=Biomass-se,width=0.1))+
  theme_classic()+
  scale_fill_brewer(palette = "Dark2",name="Species")+
  scale_color_brewer(palette = "Dark2",name="Species")+
  guides(color=FALSE)+
  xlab("Nutrient Addition Treatment")+
  theme(axis.title = element_text(size=8),
        axis.text  = element_text(size=6),
        legend.position = "bottom",
        legend.text = element_text(size=6),
        legend.title = element_text(size=7))+
  geom_point(shape=21,col="Black",size=2)
b1  
#
ggsave(filename = "Figures/SuppFigure2_neighe001.eps",
       device = cairo_ps,
       dpi=600,
       plot=b1,
       height=4,width=5,unit="in")
#
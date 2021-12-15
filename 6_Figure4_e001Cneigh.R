source("helperfunctions_e001Cneigh.R")
library(gridExtra)
library(tidyverse)
library(nlme)
library(spdep)
library(viridis)
library(gstat)#see package details for class that is required
#ln.NAtm.NAdd = log(Added N +1)
e001 <- read_csv("2021-06-25-e001.csv")
#Figure3
justI <- e001 %>% filter(NTrt=="I-0.00")#only control
ggplot(justI,aes(x=Year,y=sr))+
  geom_point()+
  facet_wrap(~Plot)
#run a regression for each plot through time
library(broom)
#
mods <- e001 %>% filter(NTrt=="I-0.00") %>%
  group_by(Plot) %>% 
  nest() %>%
  mutate(mod=map(data,~lm(sr~Year,data=.)),
         mods.t = map(mod,glance),
         fits = map(mod,augment),
         tid=map(mod,tidy))
#get model statistics
mods.t <- mods %>% 
  unnest(mods.t) %>%
  ungroup() %>% 
  arrange(desc(r.squared)) %>%
  mutate(p.value.fdr=p.adjust(p.value,"fdr")) %>%
  select(-c(data,mod))
##############
#get pvalues 
pvals <- mods.t %>% select(Plot,p.value.fdr) 
#
plotfun <- function(x,dat){
  # x="18";dat=pvals
  tmp <- justI %>% 
    filter(Plot==x)
  #run reg
  mod1 <- lm(sr~Year,data=tmp)
  #get stats
  pval <- round(summary(mod1)$coef[2,4],4)
  #get fits
  fit <- predict(mod1,se.fit=TRUE)
  #merge in
  tmp$fitsr <- fit$fit
  tmp$se <- fit$se.fit
  #get data
  dat <- dat %>% filter(Plot==x)
  #get pvalue corrected for six tests
  pval.fdr <- round(dat$p.value.fdr,3)
  #plot
  #
  means1 <- round(mean(tmp$means.Sr.4),1)
  means1
  #
  p0 <- ggplot(tmp,aes(x=Year,y=sr))+
    geom_point(size=2,shape=21,fill="grey",col="Black")+
    theme_bw()+
    geom_line(aes(y=fitsr))+
    geom_line(aes(y=fitsr+se),linetype=2)+
    geom_line(aes(y=fitsr-se),linetype=2)+
    scale_y_continuous(limits=c(0,24),breaks = c(0,5,10,15,20))+
    ylab("")+
    xlab("")+
    annotate(x=2000,y=23,geom="text",size=3,
             label=paste("p",pval.fdr,sep=" = "))+
    annotate(x=1985,y=23,geom="text",size=3,
             label=paste("Plot",unique(tmp$Plot),sep=" = "))+
    annotate(x=1989,y=2,geom="text",size=3,
             label=paste("NeighborSR",means1,sep=" = "))+
    scale_x_continuous(limits=c(1982,2004),breaks=c(1982,1992,2002))+
    theme(plot.margin = margin(t = 0.5,r = 0.5,b = 0.5,l = 0.5),
          panel.grid = element_blank(),
          axis.title = element_blank())
  p0
  return(p0)
}
#
plots <- unique(justI$Plot)
#get plots
out <- map(.x = plots,.f = plotfun,dat=pvals)
#
labfun <- function(x,lab){
  x <- x + ggtitle(lab)+
    theme(axis.text = element_text(size=7),
          axis.title = element_text(size=8),
          plot.title = element_text(size=10,family = "Helvetica",
                                    margin = margin(0,0,0,0)))
   return(x) 
}
labs <- c("a","b","c","d","e","f")
out <- map2(.x  = out,.f = labfun,.y = labs)
#
library(grid)
o1 <- exec(.fn = grid.arrange,
           grobs=out,
           ncol=3,
           bottom=textGrob("Year",
                           gp=gpar(cex=1,fontsize=9, 
                                   family="Helvitica")),
           left=
             textGrob(
               "Species Richness",rot=90, 
               gp=gpar(cex=1,fontsize=9,
                       family="Helvitica")),
           padding=unit(1,"mm"))
o1
#reviewer 2 move to main figure
# ggsave(filename = "Figures/SuppFig1_neighe001.eps",
#        device = cairo_ps,
#        dpi=600,
#        plot = o1,
#        height=4,width=6,unit="in")#############
############
#get mean across time
justImean <- justI %>% 
  group_by(Plot) %>% 
  summarise(means.Sr.4=mean(means.Sr.4))
#get coefs
coef <- mods %>% 
  select(Plot,tid) %>% 
  unnest(tid) %>% 
  filter(term=="Year")
#merge together
justImean <- left_join(justImean,coef)
##invert slopes to improve understanding
justImean$estimate <- justImean$estimate *-1
#means
justImean %>% 
  filter(Plot %in% c(18,28,30)) %>% 
  summarise(estimate=mean(estimate))
#major axis reg
sma1 <- smatr::sma(estimate~means.Sr.4,data=justImean)
sma1
#reported in text
sma1$coef
#get slope and intercept
m1 <- sma1$coef[[1]][2,1]
b1 <- sma1$coef[[1]][1,1]
#
o2 <- ggplot(justImean,aes(x=means.Sr.4,y=estimate,
                           label=Plot,
                           col=p.value>0.05))+
  geom_text()+
  scale_color_viridis(begin = 0.2,end = 0.8,discrete = TRUE,
                      option = "A",
                     name=expression(beta[Year]==0))+
  geom_errorbar(aes(ymax=estimate+std.error,
                ymin=estimate-std.error,width=0.1))+
  geom_hline(yintercept = 0,linetype=2)+
  theme_bw(base_size = 10)+
  ggtitle("g")+
  ylab("Rate of species loss per year")+
  xlab("Average Neighborhood Species Richness")+
  geom_abline(slope=m1,intercept = b1)+
  theme(legend.position = "bottom",
        axis.text = element_text(size=7),
        axis.title = element_text(size=9),
        plot.title = element_text(size=10,family = "Helvetica",
                                  margin = margin(0,0,0,0)))
o2
#
# write_rds(x = o2,file = "Fig4b_colorblind.rds")
#
#
ggsave(filename = "Figures/Figure4_neighe001.pdf",
       # device = cairo_ps,
       dpi=600,
       plot = grid.arrange(o1,o2,
                           heights=c(unit(3,"in"),unit(2,"in")),
                           ncol=1),
       height=6,width=6,unit="in")#############
                

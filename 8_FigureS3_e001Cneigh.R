library(tidyverse)
source("helperfunctions_e001Cneigh.R")
e001 <- read_csv("2021-06-25-e001.all.csv")
srNXdata <- read.csv("FINAL-Speciesliste001FieldC.csv",stringsAsFactors=FALSE)
srNXdata$Specid <- as.character(get_specid(srNXdata$Species))

e001 <- left_join(e001,srNXdata)

gam <- e001 %>% 
  filter(is.Species==1) %>% 
  filter(Biomass>0) %>% 
  group_by(Year) %>% 
  summarise(gamma=sum(NROW(unique(Specid))))
unique(e001$Specid)
ggplot(gam,aes(x=Year,y=gamma))+
  geom_point()
mod1 <- lm(gamma~Year,data=gam)
summary(mod1)
out <- predict(mod1,data=gam,se=TRUE)
gam$fit <- out$fit
gam$se <- out$se.fit
#
p1 <- ggplot(gam,aes(x=Year,y=gamma))+
  geom_line(alpha=0.3)+
  geom_point(size=3,shape=21,fill="grey",col="Black")+
  geom_line(aes(y=fit))+
  geom_line(aes(y=fit+se),linetype=2)+
  geom_line(aes(y=fit-se),linetype=2)+
  theme_bw()+
  ylab("Total Number of Species Observed")+
  scale_x_continuous(breaks=c(1982,1992,2002))
p1
ggsave(plot = p1,filename = "Figures/SuppFigure3_neighe001.pdf",
       height=4,width=6,unit="in")

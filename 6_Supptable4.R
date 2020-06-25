#
library(broom)
library(tidyverse)
# Package ID: knb-lter-cdr.7.8 Cataloging System:https://pasta.edirepository.org.
# Data set title: Soil nitrate and ammonium: Long-Term Nitrogen Deposition: Population, Community, and Ecosystem Consequences.
# Data set creator:  David Tilman -  
# Metadata Provider:    - Cedar Creek LTER 
# Contact:  Dan Bahauddin - Information Manager Cedar Creek Ecosystem Science Reserve  - webmaster@cedarcreek.umn.edu
# Stylesheet v2.7 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/7/8/c632431649614dd9e8dbd4cadbb0eb89" 
infile1 <- tempfile()
download.file(inUrl1,infile1,method="curl")
soilN <-read.csv(infile1,header=F 
                 ,skip=1
                 ,sep="\t"  
                 , col.names=c(
                   "Year",     
                   "Field",     
                   "Plot",     
                   "NTrt",     
                   "NAdd",     
                   "NitrAdd",     
                   "NAtm.plus.NAdd",     
                   "NO3Soil",     
                   "NH4Soil"    ), check.names=TRUE)
###############
e001 <- read_csv("2020-02-03-e001.csv")
#
soilN$NTrt <- as.factor(soilN$NTrt)
soilN$NTrt <- factor(soilN$NTrt,levels(soilN$NTrt)[c(9,1,2:8)])
#recode factor to match true treatments
soilN<- soilN  %>%  
  mutate(NTrt=recode(NTrt,
                     "1"="A-0.00",
                     "2"="B-1.02",
                     "3"="C-2.04",
                     "4"="D-3.40",
                     "5"="E-5.44",
                     "6"="F-9.52",
                     "7"="G-17.0",
                     "8"="H-27.2",
                     "9"="I-0.00"))
#
soilN <- soilN %>% filter(Field=="C")
justIN <- soilN %>% filter(NTrt=="I-0.00")
############
tmp <- e001 %>% filter(NTrt=="I-0.00") %>% 
  select(Plot,edgeeffect,Year,sr,Agrorepe,means.Sr.4)
justIN <- left_join(justIN,tmp)
#
mod1 <- lm(sr~Year+means.Sr.4+Agrorepe+NO3Soil+edgeeffect,
           data=justIN)
#
car::vif(mod1)
car::Anova(mod1,type="III")
write.csv(car::Anova(mod1,type="III"))
write.csv(x=car::Anova(mod1,type="III"),file="Tables/SuppTable4_neighe001.csv")

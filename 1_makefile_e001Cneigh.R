#
library(tidyverse)
library(vegan)
library(purrr)
library(furrr)
library(spdep)
library(lme4)
library(nlme)
library(viridis)
future::plan(multisession)
source("helperfunctions_e001Cneigh.R")
##########
# Package ID: knb-lter-cdr.14.8 Cataloging System:https://pasta.edirepository.org.
# Data set title: Plant aboveground biomass data: Long-Term Nitrogen Deposition: Population, Community, and Ecosystem Consequences.
# Data set creator:  David Tilman -  
# Metadata Provider:    - Cedar Creek LTER 
# Contact:  Dan Bahauddin - Information Manager Cedar Creek Ecosystem Science Reserve  - webmaster@cedarcreek.umn.edu
# Stylesheet v2.7 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/14/8/057e39850bd748d364df8a5ef60bb08d"
download.file(inUrl1,"ple001.csv",method="curl")
# #
dat <-read.csv("ple001.csv",header=F 
               ,skip=1
               ,sep="\t"  
               , col.names=c(
                 "Exp",     
                 "Year",     
                 "Field",     
                 "Plot",     
                 "NTrt",     
                 "NAdd",     
                 "NitrAdd",     
                 "NAtm.plus.NAdd",     
                 "Species",     
                 "Biomass"    ), check.names=TRUE)
############
#cedar creek species list
taxdat <- read.csv("ccesrPlantSpeciesData.csv")
#data table with species list
srNXdata <- read.csv("FINAL-Speciesliste001FieldC.csv",stringsAsFactors=FALSE)
#import row and column associations
rowcol <- read.csv("rowcol.csv",stringsAsFactors = FALSE)
####################
#filter only field C pre burning and fencing removal to get stable experimental conditions (1982-2004)
dat.c <- dat %>% filter(Field=="C") %>% filter(Year<2005)#field C only
#change type
dat.c$Year <- as.character(dat.c$Year)
#change type
dat.c$Species <- as.character(dat.c$Species)
#use custom function to get shorter species name 
dat.c$Specid <- as.character(get_specid(dat.c$Species))
#get functional groups
dat.c$funcgroup <- as.character(get_funcgroup(dat.c$Specid))#use function to get functional group, but add as character string to get rid of species not present
#arrange dataframe
dat.c <- dat.c %>% arrange(Plot,Year,Specid)
#change type
dat.c$NTrt <- as.factor(dat.c$NTrt)
#order properly note 9 is true control (Treatment I)
dat.c$NTrt <- factor(dat.c$NTrt,levels(dat.c$NTrt)[c(9,1,2:8)])
#recode factor to match true treatments
dat.c<- dat.c  %>%  
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
#join with species data
dat.c <- left_join(dat.c,srNXdata)
#join together row and columns
##########
#calculate species richness using dplyr
e001 <- dat.c %>%   filter(is.Species==1) %>% 
  group_by(Plot,Year,NAdd,NTrt) %>% 
  summarise(sr=sum(NROW(unique(Species)))) 
#merge in row column
e001 <- left_join(e001,rowcol)
#create a few useful variables
e001$NAtm.NAdd <- e001$NAdd+1
e001$ln.NAtm.NAdd <- log(e001$NAtm.NAdd)
#########
#get biomass dataframe
SpeciesBioData <- dat.c %>% select(Exp,Field,Year,Plot,Species,Specid,funcgroup,is.Species,Biomass)
###########
#calculate a few useful community properties using the Vegan package
#get matrix of biomass
SpeciesBioData.w <-  SpeciesBioData %>% 
  arrange(Year,Plot) %>%
  filter(is.Species==1) %>% 
  select(Year,Plot,Species,Biomass) %>% 
  spread(Species,Biomass,fill=0) %>% 
  select(-c(Year,Plot))
##arrange data so values line up
e001 <- e001 %>% arrange(Year,Plot)
#get counts with vegan
e001$Sr.vegan <- specnumber(SpeciesBioData.w)
#check dplyr and vegan SR
tmp <-  e001 %>% ungroup() %>% select(sr,Sr.vegan,Plot,Year)
all(as.integer(tmp$sr)==as.integer(tmp$Sr.vegan))#they are equal
#remove 
e001$Sr.vegan <- NULL
###########
#get spatial arrangment of plots using spdep
#create coordinate matrix
coor <- e001 %>% ungroup() %>% select(row,col) %>% distinct()
coor <- as.data.frame(coor)
#get neighborhood list
neighbor <- cell2nb(nrow=9,ncol = 6,type = "rook",torus = FALSE)#change to queen if you want eight instead of four
plot(neighbor,coords=coor)#check if you want to see the arrangement
#############
#create a vector of plots
plot_vec <- unique(e001$Plot)
#give the neighborhood matrix names from this vector
names(neighbor) <- c(plot_vec)
nfun <- function(x,dat,neighmat,var){
  # x=2;dat=e001[e001$Year=="1982",];neighmat=neighbor;var="sr"
  #create a vector where we take the neighbors of a given plot x 
  #and then calculate the mean values of a given variable
  tmp <- dat %>% filter(Plot %in% neighmat[[x]]) %>% ungroup() %>% 
    summarise(means.Sr.4=mean(!! var))
  return(tmp$means.Sr.4)
}
yrnfun <- function(yr,dat,vector,neighmat,var,yvar){
  # yr="1982";dat=e001;vector=plot_vec;neighmat=neighmat;var=quo(sr);yvar="means.Sr.4"
  #across years apply the previous function
  var <- enquo(var)
  neighmat <- neighmat
  tmp <- dat %>% filter(Year==yr)
  out <- future_map_dbl(.x = vector,.f = nfun,dat=tmp,neighmat=neighmat,var=var) %>% 
    enframe(name="Plot",value=yvar)
  return(out)
  # names(out) <- yr
}
yrlist <- unique(e001$Year)#year list
yrlist <- set_names(yrlist)#it needs names
#get neighbor SR
mat1 <- future_map_dfr(.x = yrlist,.f = yrnfun,
                       dat=e001,vector=plot_vec,
                       neighmat=neighbor,var=quo(sr),.id="Year",
                       yvar="means.Sr.4")
######checkfuns
#check using a corner, edge, and center plot
#first manually with a corner
check <- e001 %>% filter(Year==1983) %>% filter(Plot %in% c(2,7)) %>% 
  ungroup() %>% 
  select(Plot,sr) %>% 
  summarise(sr=mean(sr))
tmp <- mat1 %>% filter(Year==1983) %>% filter(Plot==1) %>% 
  select(means.Sr.4)
tmp==check$sr
#with functions across years
checkfun.corner <- function(yr,dat,var){
  # dat=e001;Year="1984"
  # myvar <- quo(var)
  check <- dat %>% filter(Year==yr) %>% filter(Plot %in% c(2,7)) %>% ungroup() %>% 
    summarise(sr=mean(!!var))
  a1 <- mat1 %>% filter(Year==yr & Plot==1) %>% select(means.Sr.4)
  a1==check$sr
}
#
checkfun.center <- function(yr,dat,var){
  check <- dat %>% filter(Year==yr) %>% filter(Plot %in% c(2,7,9,14)) %>% ungroup() %>% 
    summarise(sr=mean(!!var))
  a1 <- mat1 %>% filter(Year==yr & Plot==8) %>% select(means.Sr.4)
  a1==check$sr
}
#
checkfun.edge<- function(yr,dat,var){
  check <- dat %>% filter(Year==yr) %>% filter(Plot %in% c(7,14,19)) %>% ungroup() %>% 
    summarise(sr=mean(!!var))
  a1 <- mat1 %>% filter(Year==yr & Plot==13) %>% select(means.Sr.4)
  a1==check$sr
}
#run tests
all(future_map_lgl(.x = yrlist,.f = checkfun.center,dat=e001,var=quo(sr)))==TRUE
all(future_map_lgl(.x = yrlist,.f = checkfun.corner,dat=e001,var=quo(sr)))==TRUE
all(future_map_lgl(.x = yrlist,.f = checkfun.edge,dat=e001,var=quo(sr)))==TRUE
#join it up
e001$Year <- as.character(e001$Year)
e001 <- left_join(e001,mat1)
##############
#include zeros for neighbor matrix for species
SpeciesBioData.wfull <- SpeciesBioData %>%
  select(Plot,Year,Specid,Biomass) %>% 
  spread(Specid,Biomass,fill=0)
colnames(SpeciesBioData.wfull)
SpeciesBioData.l <- SpeciesBioData.wfull %>% gather(Specid,Biomass,Achimill:Violsagi)  
###########
#function to calculate neighborhood for all species
#currently VERY slow
#the current code creates a new dataframe for each interation that is useful for checking, but is slow 
iteratefun <- function(x,dat){
  # dat=SpeciesBioData.l;x="Arteludo"
  lab <- paste(x,".neighbor",sep="")
  dat <- dat %>% filter(Specid==x)
  out <- future_map_dfr(.x = yrlist,.f = yrnfun,
                        dat=dat,vector=plot_vec,
                        neighmat=neighbor,var=quo(Biomass),.id="Year",
                        yvar=lab)
  return(out)
}
#
#specieslist
splist <- c("Schiscop","Agrorepe","Poaprate")
#run function 
out <- future_map(.x = splist,.f =iteratefun,dat=SpeciesBioData.l)
#turn list of list into tibble
specs <- reduce(.x = out,.f = left_join,by=c("Plot","Year"))
#join together
e001 <- left_join(e001,SpeciesBioData.wfull)#biomass
#
#specieslist 
e001 <- left_join(e001,specs)#neighborbiomass
########
#check funs 
checkfun.edge_spec<- function(yr,dat,x){
  # dat=e001;x="Arteludo";yr=1982;
  lab <- paste(x,".neighbor",sep="")
  check <- dat %>% filter(Year==yr) %>% filter(Plot %in% c(7,14,19)) %>% ungroup() %>% 
    summarise(var=mean(!!sym(x)))
  check
  a1 <- e001 %>% filter(Year==yr & Plot==13)  %>% ungroup()%>% select(!!sym(lab))
  a1==check$var
}
#with functions across years
checkfun.corner_spec<- function(yr,dat,x){
  # dat=e001;x="Arteludo";yr=1982;
  lab <- paste(x,".neighbor",sep="")
  check <- dat %>% filter(Year==yr) %>% filter(Plot %in% c(2,7)) %>% ungroup() %>% 
    summarise(var=mean(!!sym(x)))
  check
  a1 <- e001 %>% filter(Year==yr & Plot==1)  %>% ungroup()%>% select(!!sym(lab))
  a1==check$var
}
#
checkfun.center_spec<- function(yr,dat,x){
  # dat=e001;x="Schiscop";yr=1982;
  lab <- paste(x,".neighbor",sep="")
  check <- dat %>% filter(Year==yr) %>% filter(Plot %in% c(2,7,9,14)) %>% ungroup() %>% 
    summarise(var=mean(!!sym(x)))
  a1 <- e001 %>% filter(Year==yr & Plot==8)  %>% ungroup()%>% select(!!sym(lab))
  a1==check$var
}
#
all(future_map_lgl(.x = yrlist,.f = checkfun.center_spec,dat=e001,x="Schiscop"))==TRUE
all(future_map_lgl(.x = yrlist,.f = checkfun.corner_spec,dat=e001,x="Schiscop"))==TRUE
all(future_map_lgl(.x = yrlist,.f = checkfun.edge_spec,dat=e001,x="Schiscop"))==TRUE
#
########
#get long data frame
colnames(e001)
e001.l <- e001 %>% gather(Specid,Biomass,Achimill:Violsagi)
#and functional groups
e001.l$func <- get_funcgroup(e001.l$Specid)
#
dat <- e001 %>%
  arrange(Plot,Year)
#caculate mean neighborhood across time series 
dat.m <- e001 %>% group_by(Plot) %>% 
  summarise(means.Sr.4_mean=mean(means.Sr.4))
dat <- left_join(dat,dat.m)
#
write_csv(file = "2021-06-25-e001.csv",x = dat)
# write_csv(file = "2021-06-25-e001.l.csv",x = e001.l)
#file is too big with all species present
#shorten to just those we use
e001.l <- e001.l %>%
  filter(Specid %in% splist)
write_csv(file = "2021-06-25-e001.l.csv",x = e001.l)
#clean it up
rm(list=setdiff(ls(), c("e001","e001.l","srNXdata","dat")))
####
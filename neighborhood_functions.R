#negate in as in NOT IN
"%ni%" <- Negate("%in%")
#get path
path <- "ccesrPlantSpeciesData.csv"
#get scaled values
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}
#matching function 
get_specid<-function(Species) {
  taxdat <-read.csv(path)
  taxdat[match(Species, taxdat$Species),"Specid"]
} 
#matching function
get_funcgroup<-function(Specid) {
  taxdat <-read.csv(path)
  fg <- taxdat[match(Specid, taxdat$Specid),"FunctionalGroup"]
  fg
}
#get AIC diff
aicfun <- function(mod1,mod2){
  print(round(diff(anova(mod1,mod2)$AIC),1))
  print(round(anova(mod1,mod2)$L.Ratio[2],1))
}
#get standard error 
#https://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean
my.stand <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
#matching function 
get_species<-function(Specid) {
  taxdat <-read.csv(path)
  taxdat[match(Specid, taxdat$Specid),"Species"]
} 
#
biolab <- expression(paste("Focal Plot Biomass", " ", "(", 
                          "g" %.% "m"^-2, ")"))
#
naddlab <- expression(paste("Added Nitrogen ","(g" %.% "m"^-2 %.% "yr"^-1,")"))
#
bioindlab <- expression(paste("Each Species' Focal Plot Biomass", " ", "(", 
                              "g" %.% "m"^-2, ")"))
#
neighlab<- expression(paste("Each Species' Neighborhood Biomass", " ","(","g"%.%"m"^-2, ")"))
#
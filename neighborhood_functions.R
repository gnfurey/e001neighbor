path <- "ccesrPlantSpeciesData.csv"
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}
get_specid<-function(Species) {
  taxdat <-read.csv(path)
  taxdat[match(Species, taxdat$Species),"Specid"]
} 
get_funcgroup<-function(Specid) {
  taxdat <-read.csv(path)
  fg <- taxdat[match(Specid, taxdat$Specid),"FunctionalGroup"]
  fg
}
scale_this <-function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}
aicfun <- function(mod1,mod2){
  print(round(diff(anova(mod1,mod2)$AIC),1))
  print(round(anova(mod1,mod2)$L.Ratio[2],1))
}
#https://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean
my.stand <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
get_species<-function(Specid) {
  taxdat <-read.csv(path)
  taxdat[match(Specid, taxdat$Specid),"Species"]
} 
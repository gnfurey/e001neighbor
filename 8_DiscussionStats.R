###########
e001 <- read_csv("2020-02-03-e001.csv")
source("/Users/furey034/Documents/Functions/funcs.R")
pchange <- function(new,old){
  x <- ((new-old)/abs(old))*100
  return(x)
}
means <- e001 %>%
  filter(Year>1994)%>%
  group_by(Plot,NTrt,NAdd) %>% 
  summarise(sr=mean(sr)) %>% 
  filter(NAdd<4)
unique(means$NAdd)

mod <- lm(sr~NAdd,
          data=means) 
summary(mod)
confint(mod)
ggplot(means,aes(x=NAdd,y=sr))+geom_point()+
  geom_smooth(method="lm")
pchange(5.6,4.6)#line 362
a <- 12.1*0.22
a#line 368
b <- 2.7*0.3
b
#
dat <- read.csv("2020-02-03-e001.csv")
Table1 <- lme(sr~
                ln.NAtm.NAdd+Year+
                means.Sr.4_mean,
              random=~1|Plot,
              weights=varIdent(form=~1|NTrt),
              correlation=corARMA(form=~Year|Plot,p=2,q=1),
              control=list(maxIter = 10000,msMaxIter=10000,niterEM=10000,msMaxEval = 10000),
              method="REML",
              data=dat)
summary(Table1)
intervals(Table1,which="fixed")
2.7*0.3
0.07*2.7
0.52*2.7
```
source("neighborhood_functions.R")
# library(animation)
library(tidyverse)
library(gganimate)
library(viridis)
e001.l <- read_csv("2021-06-25-e001.l.csv")
e001.l$NTrt <- as.factor(e001.l$NTrt)
e001.l$NTrt <- factor(e001.l$NTrt,levels(e001.l$NTrt)[c(9,1,2:8)])
tmp <- e001.l %>% select(Plot,NAdd,NTrt,row,col,Year)
plot.gif <- function(specid){
  # specid <- "Agrorepe"
  dat <- e001.l
  ###
  cleanspec <- get_species(specid)
  cleanspec
  cleanspec <- ifelse(cleanspec=="Agropyron repens",
                      "Elymus repens",cleanspec)
  cleanspec <- ifelse(cleanspec=="Aster azureus",
                      "S. oolentangiense",
                      cleanspec)
  
  ###
  dat1 <- e001.l %>% filter(Specid==specid) %>% filter(Biomass>0)
  plot <- ggplot(dat1, aes(x=as.factor(row),
                           y=as.factor(col),
                           frame=Year,fill=NTrt)) +
    geom_tile(col="black",data=tmp)+
    geom_point(aes(size=Biomass),
               shape=21,fill="grey")+
    scale_size(range = c(0, 6))+
    scale_fill_viridis(discrete=TRUE,option="A",direction=-1,begin=0.3,end = 1)+
    ylab("Rows")+xlab("Columns")+
    scale_x_discrete(labels=c(6:1))+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.key.height =  unit(0.4, "cm"),
          legend.key.width = unit(0.3,"cm"),
          legend.text = element_text(size=8),
          panel.border = element_rect(linetype = 1, fill = NA))+
    transition_manual(Year,cumulative = FALSE)+
    labs(title = paste(cleanspec,"Year: {current_frame}"))+
    scale_size_continuous(range = c(3,17),breaks=c(1,100,500))+
    theme(plot.title = element_text(size=15))
  plot
  # name <- paste(as.character(substitute(specid)),"plot.gif",sep="")
  return(plot)
}
#
splist <- c("Schiscop","Agrorepe","Poaprate")
animate(plot.gif("Agrorepe"),duration = 32,rewind=FALSE,
        renderer=ffmpeg_renderer(
          format = "mov",
          ffmpeg = NULL,
          options = list(pix_fmt = "yuv420p")),
        height = 6, width = 6, units = "in", res = 300)
animate(plot.gif("Poaprate"),duration = 32,rewind=FALSE,
        renderer=ffmpeg_renderer(
          format = "mov",
          ffmpeg = NULL,
          options = list(pix_fmt = "yuv420p")),
        height = 6, width = 6, units = "in", res = 300)
animate(plot.gif("Schiscop"),duration = 32,rewind=FALSE,
        renderer=ffmpeg_renderer(
          format = "mov",
          ffmpeg = NULL,
          options = list(pix_fmt = "yuv420p")),
        height = 6, width = 6, units = "in", res = 300)
# savefun <- function(x){
#   x="Agrorepe"
#   animate(plot.gif(specid=x),duration = 32,rewind=FALSE,
#           height = 6, width = 6, units = "in", res = 300)
#   anim_save(filename = paste("Figures/Species_gif/",x,".gif",sep=""))
# }
# walk(.x = splist,.f = savefun)
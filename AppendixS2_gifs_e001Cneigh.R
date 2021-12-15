source("helperfunctions_e001Cneigh.R")
# library(animation)
library(tidyverse)
library(gganimate)
library(viridis)
e001.l <- read_csv("2021-06-25-e001.l.csv")
e001.l$NTrt <- as.factor(e001.l$NTrt)
e001.l$NTrt <- factor(e001.l$NTrt,levels(e001.l$NTrt)[c(9,1,2:8)])
#
col2 <- data.frame(col2=c(1,2,3,4,5,6,7,8,9),col=c(9,8,7,6,5,4,3,2,1))
e001.l <- left_join(e001.l,col2)
row2 <- data.frame(row2=c(1,2,3,4,5,6),row=c(6,5,4,3,2,1))
e001.l <- left_join(e001.l,row2)
e001.l$row2 <- as.factor(e001.l$row2)
e001.l$col2 <- as.factor(e001.l$col2)
#
keydata <- e001.l %>% select(Plot,NAdd,NTrt,row,col,row2,col2) %>% 
  distinct()
#
plot.gif <- function(specid){
  # specid <- "Agrorepe"
  # dat <- e001.l
  ###
  cleanspec <- get_species(specid)
  cleanspec
  cleanspec <- ifelse(cleanspec=="Agropyron repens",
                      "Elymus repens",cleanspec)
  ###
  dat1 <- e001.l %>% filter(Specid==specid) %>% filter(Biomass>0)
  ###############
  cols <- viridis(n=8,option="A",direction=-1,begin=0.3,end = 1)
  tmpdat <- data.frame(hex=c("#FFFFFF",cols),NTrt=c("I-0.00",levels(keydata$NTrt)[-1]))
  tmpdat$NTrt <- as.factor(tmpdat$NTrt)
  tmpdat$NTrt <- factor(tmpdat$NTrt,levels(tmpdat$NTrt)[c(9,1,2:8)])
  #
  keydata <- left_join(keydata,tmpdat)
  #
  keydata <- keydata %>% arrange(NTrt)
  keydata$hex <- factor(keydata$hex,levels = unique(keydata$hex))
  levels(keydata$hex)
  #
  dat_col <- left_join(dat1,keydata)
  #
  plot <- ggplot(dat1, aes(x=row2,
                           y=col2,
                           frame=Year)) +
    scale_size(range = c(0, 6))+
    scale_fill_manual(values=levels(keydata$hex),labels=
                        c("Control",
                          "0 g N",
                          "1.02 g N",
                          "2.04 g N",
                          "3.40 g N",
                          "5.44 g N",
                          "9.52 g N",
                          "17.0 g N",
                          "27.2 g N"))+
    ylab("Rows")+xlab("Columns")+
    scale_x_discrete(labels=c(6:1))+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background = element_rect(fill = "black",
                                          colour = "black",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.key.height =  unit(0.4, "cm"),
          legend.key.width = unit(0.3,"cm"),
          legend.text = element_text(size=8),
          panel.border = element_rect(linetype = 1, fill = NA))+
    transition_manual(Year,cumulative = FALSE)+
    labs(title = paste(cleanspec,"Year: {current_frame}"))+
    scale_size_continuous(range = c(3,17),
                          breaks=c(0,10,100,500))+
    theme(plot.title = element_text(size=15))
  plot
  if(specid=="Schiscop"){
    plot <- plot+geom_tile(data=dat_col,col="Black",
                               aes(x=row2,y=col2,
                                   fill=NTrt),
                               inherit.aes = FALSE)+
      geom_point(aes(size=Biomass),
                 shape=21,fill="grey")
  }else{
    plot <- plot+geom_tile(data=keydata,col="Black",
              aes(x=row2,y=col2,
                  fill=NTrt),
              inherit.aes = FALSE)+
      geom_point(aes(size=Biomass),
                 shape=21,fill="grey")
  }
  plot
  return(plot)
}
#
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
####
plot.nocolor <- function(specid){
  # specid <- "Agrorepe"
  ###
  cleanspec <- get_species(specid)
  cleanspec
  cleanspec <- ifelse(cleanspec=="Agropyron repens",
                      "Elymus repens",cleanspec)
  ###
  dat1 <- e001.l %>% 
    filter(Specid==specid) %>% 
    filter(Biomass>0)
  #
  tmp1 <- e001.l %>% 
    select(Plot,row2,col2,Year)
  #
  plot <- ggplot(dat1, aes(x=row2,
                           y=col2,
                           frame=Year)) +
    geom_tile(fill="White",col="Black",data=tmp1)+
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
          panel.background = element_rect(fill = "black",
                                          colour = "black",
                                          size = 0.5, linetype = "solid"), 
          panel.grid.major=element_blank(),
          # panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.key.height =  unit(0.4, "cm"),
          legend.key.width = unit(0.3,"cm"),
          legend.text = element_text(size=8),
          panel.border = element_rect(linetype = 1, fill = NA))+
    transition_manual(Year,cumulative = FALSE)+
    labs(title = paste(cleanspec,"Year: {current_frame}"))+
    scale_size_continuous(range = c(3,17),breaks=c(0,1,100,500))+
    theme(plot.title = element_text(size=15))
  plot
  return(plot)
}
#
animate(plot.nocolor("Agrorepe"),duration = 32,rewind=FALSE,
        renderer=ffmpeg_renderer(
          format = "mov",
          ffmpeg = NULL,
          options = list(pix_fmt = "yuv420p")),
        height = 6, width = 6, units = "in", res = 300)
animate(plot.nocolor("Poaprate"),duration = 32,rewind=FALSE,
        renderer=ffmpeg_renderer(
          format = "mov",
          ffmpeg = NULL,
          options = list(pix_fmt = "yuv420p")),
        height = 6, width = 6, units = "in", res = 300)
animate(plot.nocolor("Schiscop"),duration = 32,rewind=FALSE,
        renderer=ffmpeg_renderer(
          format = "mov",
          ffmpeg = NULL,
          options = list(pix_fmt = "yuv420p")),
        height = 6, width = 6, units = "in", res = 300)
writeLines(capture.output(sessionInfo()), "Figures/VIDEOS1S6_sessionInfo.txt")
citation("viridis")
citation("ggplot2")
citation("gganimate")

# savefun <- function(x){
#   x="Agrorepe"
#   animate(plot.gif(specid=x),duration = 32,rewind=FALSE,
#           height = 6, width = 6, units = "in", res = 300)
#   anim_save(filename = paste("Figures/Species_gif/",x,".gif",sep=""))
# }
# walk(.x = splist,.f = savefun)
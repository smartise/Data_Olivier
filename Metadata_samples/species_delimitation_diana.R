library(readr)
library(ggplot2)
library(ggfortify)
library(dplyr)
library(stats)
library(tidyr)
library(readtext)
library(lubridate)
library(stringr)
library(rstatix)
library(purrr)
library(readxl)

setwd("C:\\Users\\olivi\\Desktop\\Blast")
HPLC<- read_excel("Blast_diana.xlsx") 
HPLC <- HPLC[,c(1, 2, 12)]

plot_style <-theme(axis.title.y = element_text(colour="Red", size=15),
                   axis.title.x = element_text(colour="Red", size=15),
                   axis.text.x = element_text(size=15),
                   axis.text.y = element_text(size=15),
                   strip.text = element_text(size=10),
                   legend.title=element_text(size=18),
                   legend.title.align = 0.5,
                   legend.text=element_text(size=18),
                   plot.title = element_text(colour = "DarkBlue", 
                                             size=15),
                   panel.border = element_rect(linetype = "dashed", fill=NA)) 


HPLC <- HPLC%>% group_by(pecies, zone)%>% mutate(number=nrow())






ggplot(HPLC, aes(x=zone, fill =as.factor(pecies)))+
  geom_bar()+
  ylab("Number")+
  xlab("Zone")+
  ggtitle(paste("Number of Corals depending on the species and the zone", sep=""))+
  plot_style
ggsave(filename="setup2.JPG", width=30, height=20, dpi=300, units="cm")

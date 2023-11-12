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
setwd("C:\\Users\\olivi\\OneDrive - Universite de Liege\\OBSIDIAN\\Project Panama\\")


samples<- read_excel("stress thermique.xlsx") 

###################################################plot###############################

p <- ggplot(samples, aes(x=day, y=temp, color=stress)) + 
  geom_point(size=3)+
  geom_line(aes(group=stress)) + 
  labs(title="Temperature over Days", x="Day", y="Temperature (°C)") + 
  scale_color_manual(values=c("cold stress"="blue", "control"="green", "heat stress"="red"))+
  theme_minimal()
  ggsave(filename="R_schema/temperature stress.JPG",plot=p,  width=20, height=10, dpi=300, units="cm")

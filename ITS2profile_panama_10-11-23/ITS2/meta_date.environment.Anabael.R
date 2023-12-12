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
library(vegan)
library(gplots)
library(fossil)
library(openxlsx)

setwd("C:\\Users\\olivi\\OneDrive\\OneDrive liège\\OBSIDIAN\\Project Panama\\Data\\")
Info <- read.csv("Metadata_samples\\data_matt.csv")

Info <- Info[c(2555:3094),]
Info <- Info[,c(3, 6,7, 8, 9)]

Info <- unique(Info)

write.xlsx(Info, "Metadata_samples\\environmental_data.xlsx")

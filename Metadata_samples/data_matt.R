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
library(writexl)
library(gridExtra)
library(jpeg)
library(grid) 
library(googlesheets4)
library(googledrive)
library(DT)
setwd("C:\\Users\\olivi\\OneDrive\\OneDrive liège\\OBSIDIAN\\Project Panama\\Data\\Metadata_samples")

#sheet <- drive_get(path = "https://docs.google.com/spreadsheets/d/1PQhLMpu4syPWe8bu-Oq2eGqD4wQr94ExNmKny_h8UpQ/edit?pli=1#gid=0")
#drive_download(sheet, path = "data_matt.csv", overwrite = TRUE, type = "csv")

########################## traitement csv matt#############################################
data <- read.csv("data_matt.csv")
data <- data[,c(2,3,4,5,6,10,15,29)]
colnames(data)[colnames(data) == "DNA.ID"] <- "ID"
colnames(data)[colnames(data) == "Sample.number...ARMS.Unit"] <- "colony"
colnames(data)[colnames(data) == "Station.number"] <- "replicate"
data$colony <- substr(data$colony, start=10, stop=12)

############################ ajoux des espèces #############################################
species <- read.csv("Diana_Blast_231023.csv")
species <- species[,c(2,6)]
colnames(species)[colnames(species) == "Colony"] <- "colony"
species$colony <- as.character(species$colony)

temp <- which(species$colony == 158) #correction diana 
species$colony[species$colony == 159] <- 158
species$colony[temp] <- 159
############################ fusion########################################################

info<- merge(data, species, by = "colony")
library(DT)

datatable(info) %>%
  formatStyle(columns = 'species', 
              target = 'cell', 
              fontStyle = 'italic')


write.csv(info, "info.csv", row.names = FALSE)

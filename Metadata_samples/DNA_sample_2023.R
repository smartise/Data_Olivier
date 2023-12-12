sheet <- drive_get(path = "https://docs.google.com/spreadsheets/d/1PQhLMpu4syPWe8bu-Oq2eGqD4wQr94ExNmKny_h8UpQ/edit?pli=1#gid=0")
drive_download(sheet, path = "data_matt.csv", overwrite = TRUE, type = "csv")

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

setwd("C:\\Users\\olivi\\OneDrive - Universite de Liege\\OBSIDIAN\\Project Panama\\")

samples<- read_excel("DNA_sample_2023.xlsx") 
Info<- read_excel("SymPortal_datasheet.xlsx") 
species <- read.csv("Diana_Blast_231023.csv")
data <- read.csv("data_matt.csv")

################################donnée info 2#############################################
colnames(data)[colnames(data) == "DNA.ID"] <- "ID"
colnames(data)[colnames(data) == "Sample.number...ARMS.Unit"] <- "colony"
colnames(data)[colnames(data) == "Qubit..ng.uL."] <- "cons"
data <- data[,c(4,5,6,10,15,29)]

data$colony <- substr(data$colony, start=9, stop=12)
data$colony <- as.integer(data$colony)
################################traitement donnée info#####################################

#colnames(Info)[colnames(Info) == "sample_name"] <- "ID"
#colnames(Info)[colnames(Info) == "Season"] <- "season"
#colnames(Info)[colnames(Info) == "Sample.number...ARMS.Unit"] <- "colony"
#Info <- Info[,c(1, 15, 16, 17, 18, 19)]

################################traitement donnée species#####################################

species <- species[,c(5,6)]
colnames(species)[colnames(species) == "Colony"] <- "colony"

################################Ther merging################################################

merged_df <- merge(samples, data, by = "ID")
merged_df <- merge(merged_df, species, by = "colony")
merged_df <- merged_df[!is.na(merged_df$Species), ] 

#########################################save as a CSV#######################################

write.csv(merged_df, "Diana_colony_identity.csv", row.names = FALSE)





##############keep the data that has upwelling and non upwelling value#########################

merged_df <- merged_df %>% group_by(colony) %>% filter(length(unique(season)) == 2)

#################################choose sample to use#######################################

selected_colonies <- merged_df %>%
  group_by(reagion, Species) %>%
  sample_n(3) %>%
  ungroup() %>%
  select(colony)

filtered_df <- merged_df %>%
  semi_join(selected_colonies, by = c("colony"))

#################################save the table as jpg#####################################

table_image <- tableGrob(filtered_df)

# Export the table grob as a JPG image
jpeg("R_schema/output.jpg", width = 800, height = 800, units = "px", res = 100)
grid.draw(table_image)
dev.off()

###############################info of the merged data###########################
num_rows <- length(merged_df$ID)
print(num_rows)



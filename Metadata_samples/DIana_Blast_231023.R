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

blasttype1 <- read.csv("Alignment_type1mitochondriallineage_231023_Diana_pocillopora.csv") 
blasttype3 <- read.csv("Alignment_type3mitochondriallineage_231023_Diana_pocillopora.csv")
datadiana <- read_excel("PCR_colony_ID.xlsx")

################################################process data 

blasttype1 <- blasttype1[,c(1,7)]
colnames(blasttype1)[colnames(blasttype1) == "Per..ident"] <- "type1"

blasttype3 <- blasttype3[,c(1,7)]
colnames(blasttype3)[colnames(blasttype3) == "Per..ident"] <- "type3"

mergedtipe <- merge(blasttype1, blasttype3, by = "Description")

#####################################################decide the tipe#############################

mergedtipe <- mergedtipe %>% mutate(mitochondrial = ifelse(type1 > type3, "type1", "type3"))
mergedtipe <- mergedtipe %>% group_by(mitochondrial) %>% mutate(species = ifelse(mitochondrial == "type1", "P.meandrina", "P.verrucosa") )

###################################################match the colony######################"
datadiana <- datadiana [,c(1,2)]
colnames(datadiana)[colnames(datadiana) == "RRR_ID"] <- "Description"
datadiana$Description <- as.integer(datadiana$Description)

mergedtipe <- merge(datadiana, mergedtipe, by = "Description", all.x = TRUE)

######################################################export data##########################

write.csv(mergedtipe, "Diana_Blast_231023.csv", row.names = FALSE)


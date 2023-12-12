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

setwd("C:\\Users\\olivi\\OneDrive\\OneDrive liège\\OBSIDIAN\\Project Panama\\Data\\")
relative_abund<- read_delim("ITS2profile_panama_10-11-23\\ITS2\\ITS2profile_data_profiles.absolute.abund_08.10.23.txt", delim="\t")  # assuming tab-delimited; adjust accordingly
absolute_abund<- read_delim("ITS2profile_panama_10-11-23\\ITS2\\ITS2profile_data_profiles.relative.abund_08.10.23.txt", delim="\t")
Info <- read.csv("Metadata_samples\\info.csv")
Corr_DL<- read_excel("ITS2profile_panama_10-11-23\\ITS2\\ITS2profile_meta_CD_11.11.23.xlsx")
#############################################traitement des donnés##############################################################
relative_abund <- relative_abund[-(1:5), -1 ]
relative_abund[1,1] <- "ID"
colnames(relative_abund) <- as.character(unlist(relative_abund[1, ]))
relative_abund <- relative_abund[-1,]
relative_abund <- relative_abund[-c(568, 569), ]

Relative<- relative_abund %>%
  pivot_longer(c(2:93), names_to = "symbtypes", values_to = "relative")%>%
  mutate(relative= as.numeric(relative)) %>%
  filter(relative != 0)  


absolute_abund <- absolute_abund[-(1:5), -1 ]
absolute_abund[1,1] <- "ID"
colnames(absolute_abund) <- as.character(unlist(absolute_abund[1, ]))
absolute_abund <- absolute_abund[-1,]
absolute_abund <- absolute_abund[-c(568, 569), ]

Absolute<- absolute_abund %>%
  pivot_longer(c(2:93), names_to = "symbtypes", values_to = "absolute")%>%
  mutate(absolute= as.numeric(absolute)) %>%
  filter(absolute != 0)  
#################################################fusion des données#########################################################

dataITS2 <- inner_join(Relative,Absolute,by=c("ID","symbtypes"))%>%
  left_join(.,Info, by=c("ID"))%>%
  mutate(Locality=as.factor(Locality),ID=as.factor(ID),Site=as.factor(Site), Season=as.factor(Season))%>%
  mutate(clade = case_when(
    str_detect(symbtypes, "A") ~ "A",
    str_detect(symbtypes, "B") ~ "B",
    str_detect(symbtypes, "C") ~ "C",
    str_detect(symbtypes, "D") ~ "D",
    str_detect(symbtypes, "E") ~ "E",
    str_detect(symbtypes, "F") ~ "F", 
    TRUE ~ NA_character_))%>%
  filter(clade=="C"| clade=="D")

###############################################correction diana ###########################################################

Corr_DL <- subset(Corr_DL, Corr_DL$`photo ID check` == "mixed types")
colony_extract <- Corr_DL$`Colony no.`
dataITS2<- dataITS2[!(dataITS2$colony %in% colony_extract), ]

################################################other adding######################################################

dataITS2$symbtypes_s <- sub("(/|-).*", "", dataITS2$symbtypes) #keep the first clade
dataITS2$species_colony <- paste(dataITS2$species, dataITS2$colony) #fuse species and colony
dataITS2 <- dataITS2[!grepl("BLANK-.", dataITS2$ID), ] #delete the blank

################################################save the whole data#######################################################

write.csv(dataITS2, file = "ITS2profile_panama_10-11-23\\ITS2\\ITS2profile_data_whole.csv", row.names = FALSE)

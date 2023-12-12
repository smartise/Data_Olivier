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
library(vegan)
library(fossil)
library(RColorBrewer)

setwd("C:\\Users\\olivi\\OneDrive - Universite de Liege\\OBSIDIAN\\Project Panama\\")

#abundance data tunred into community matrix
info <- read.csv("info.csv")
data1<-read.csv("bc05_rel-abundance.tsv", sep = "\t")
data1$ID <- "ML5825"
data2 <-read.csv("bc06_rel-abundance.tsv", sep = "\t")
data2$ID <- "ML6219"
data3<-read.csv("bc11_rel-abundance.tsv", sep = "\t")
data3$ID <- "ML5976"
data4<-read.csv("bc12_rel-abundance.tsv", sep = "\t")
data4$ID <- "ML6279"
emu2 <- rbind(data1, data2, data3, data4)

####### put the info on the dataframe####################################

emu2 <- merge(emu2, info, by = "ID")
emu2$other <- paste(emu2$species.y, emu2$colony, emu2$Locality)

#######################

emu2[,15] <- as.integer(emu2[,15]) #count column. # colonne 11 ? 
commu <- create.matrix(emu2, tax.name="species.x", locality="ID", abund = TRUE, abund.col="estimated.counts")
commu <- t(commu)

#rarefaction plots


rarefy(commu, 20,  se = TRUE, MARGIN = 1)
rarecurve(commu, 20, xlab = "Number of reads", ylab = "Species", label = TRUE)
rareslope(commu, 20)


S <- specnumber(commu) # observed number of species
raremax <- min(rowSums(commu))
Srare <- rarefy(commu, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(commu, step = 10, sample = raremax, col = "blue", cex = 0.6)

########################## je retire les taxa trop minime#######################
emu2 <- emu2[emu2$abundance >= 0.005, ]
emu2$species.x <- sapply(strsplit(emu2$species.x, " "), `[`, 1)
emu2 <- emu2 %>% group_by(species.x, ID) %>% mutate(relative_abundance = sum(abundance))

emu2$species_over_0_2 <- ifelse(emu2$relative_abundance > 0.03, emu2$species.x, "other")

##########################################plot des espèces 

# Create the plot

p <- ggplot(emu2, aes(x= Season, y = abundance, fill = species_over_0_2)) +
  geom_bar(stat = "identity", width = 0.95) +
  labs(x = "Season",
       y = "Relative abundances",
       fill = "Genus") +
  theme(legend.position = "bottom", legend.title = element_text(vjust = 1)) +
  facet_grid(cols = vars(other), scales ="free")+
  guides(fill = guide_legend(ncol = 6))+
  scale_fill_manual(values = colors)+
  scale_x_discrete(expand = c(0.5, 0))+
  #scale_y_discrete(expand = c(0.01, 0))+
  plot_style

# Save the plot to a file
ggsave(filename = paste("R_schema/plot_bact.jpg"), width=29, height=27, dpi=300, units="cm",plot = p)


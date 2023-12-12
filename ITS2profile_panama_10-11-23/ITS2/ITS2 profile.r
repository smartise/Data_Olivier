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

setwd("/home/olivier/OneDrive/OBSIDIAN/Project Panama/Data")
relative_abund<- read_delim("ITS2profile_panama_10-11-23\\profiles.relative.abund_and_meta.txt", delim="\t")  # assuming tab-delimited; adjust accordingly
absolute_abund<- read_delim("ITS2profile_panama_10-11-23\\profiles.absolute.abund_and_meta.txt", delim="\t")
Info <- read.csv("Metadata_samples\\info.csv")
species <- read_csv("Metadata_samples\\Diana_Blast_231023.csv")
################################traitement donnée de symportal################################

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



#################################merge of the profile #####################################
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

######################################### Checking number of reads ###############################

absol_sum <-absolute_abund%>% #data set with absolute read number values per sample ID
  mutate(Totalrow=rowSums (.[2:93] )) %>%
  filter(Totalrow>=250)%>% #standardizing number of reads to 2 std above and below mean
  filter(Totalrow<=16560) %>%
  select (ID) 



data <- data[!grepl("BLANK-.", data$ID), ]

#################################lower the tipe of clade#########################################

data$ITS2 <- sub("(/|-).*", "", data$ITS2)

##################################add the species###########################################
species <- species[,c(2,6)]
colnames(species)[colnames(species) == "Colony"] <- "colony"
species$colony <- as.character(species$colony)

data <- merge(species, data, by = "colony", all.x = TRUE)
data$species_colony <- paste(data$species, data$colony)
##############################plot##################################

# ggplot(data, aes(x=ID, fill=conc, color=ITS2))+
#   geom_bar()+
#   ylab("colonies")+
#   xlab("realtive concentation")+
#   ggtitle(paste("Measurement of rETR Day ",temp, " stress ",temp_stress, sep=""))+
# ggsave(filename=paste("R_schema/ITS2_profile.jpg", sep=""), width=100, height=100, dpi=300, units="cm")
# 
# 
# ggplot(merged_df, aes(x = ID, y = conc, fill = ITS2)) +
#   geom_bar(stat = "identity") +
#   labs(title = "Bar plot with ID, conc and ITS2",
#        x = "ID",
#        y = "Concentration",
#        fill = "ITS2 Type") +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   facet_wrap(~ site + season, scales = "free_x")
#   ggsave(filename=paste("R_schema/ITS2_profile.jpg"), width=100, height=100, dpi=300, units="cm", limitsize = FALSE)
  
  # Loop through each site

######## fro the label 

data <- data %>%
  arrange(ID, ITS2) %>%
  group_by(ID) %>%
  mutate(cumulative_conc = cumsum(conc) - (conc / 2))

sites <- unique(data$Site)

# Loop through each site
  for (s in sites) {
    
    # Filter data for the current site
    temp_data <- data[data$Site == s, ]
    
    # Create the plot
    p <- ggplot(temp_data, aes(x = replicate, y = conc, fill = ITS2)) +
      geom_bar(stat = "identity") +
      labs(title = paste("Bar plot for site", s),
           x = "ID",
           y = "Concentration",
           fill = "ITS2 Type") +
      theme_minimal() +
      theme(legend.position = "bottom", legend.title = element_text(vjust = 1)) +
      guides(fill = guide_legend(ncol = 2))+
      facet_grid(cols = vars(Season), rows = vars(species), scales ="free")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    
    
    # Save the plot to a file
    ggsave(filename = paste("R_schema/plot_for_sit_", s, ".jpg"), width=18, height=45, dpi=300, units="cm",plot = p)
  }

for (s in sites) {
  
  # Filter data for the current site
  temp_data <- data[data$Site == s, ]
  
  # Create the plot
  p <- ggplot(temp_data, aes(x = replicate, y = conc, fill = ITS2)) +
    geom_bar(stat = "identity") +
    labs(title = paste("Bar plot for site", s),
         x = "ID",
         y = "Concentration (%)",
         fill = "clade") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_text(vjust = 1)) +
    guides(fill = guide_legend(ncol = 2))+
    facet_grid(cols = vars(species_colony), rows = vars(Season), scales ="free")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  
  # Save the plot to a file
  ggsave(filename = paste("R_schema/clade_for_sit_", s, ".jpg"), width=45, height=30, dpi=300, units="cm",plot = p)
}

#######################################statistical analysis #########################################

matrixITS2 <- create.matrix(data, tax.name="ITS2", locality="ID", abund = TRUE, abund.col="conc")


########################################distance based redudancy####################################

########################################Non-Metric Multidimensional Scaling (NMDS)######################

nmds_result <- metaMDS(matrixITS2, distance = "bray")
plot(nmds_result) #there is 4 points that are hugly 


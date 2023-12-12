library(readr)
library(ggplot2)
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
relative_abund<- read_delim("ITS2profile_panama_10-11-23/ITS2/ITS2profile_data_profiles.absolute.abund_08.10.23.txt", delim="\t")  # assuming tab-delimited; adjust accordingly
absolute_abund<- read_delim("ITS2profile_panama_10-11-23/ITS2/ITS2profile_data_profiles.relative.abund_08.10.23.txt", delim="\t")
Info <- read.csv("Metadata_samples/info.csv")
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
  left_join(.,Info, by=c("ID"))

##################################delete the Blank###########################################

dataITS2 <- dataITS2[!grepl("BLANK-.", dataITS2$ID), ]

#################################lower the tipe of clade#########################################

dataITS2$clade <- sub("(/|-).*", "", dataITS2$symbtypes)

################################make a column with all the info fused######################

dataITS2$combined <- paste(dataITS2$colony, dataITS2$species, sep = " ")

################################defining it as factor for the colors on the plot ###########

dataITS2 <- dataITS2%>%mutate(Locality=as.factor(Locality),ID=as.factor(ID),Site=as.factor(Site), Season=as.factor(Season), clade=as.factor(clade))


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

dataITS2 <- dataITS2 %>%
  arrange(ID, ITS2) %>%
  group_by(ID) %>%
  mutate(cumulative_conc = cumsum(conc) - (conc / 2))

sites <- unique(dataITS2$Site)

# Loop through each site
  for (s in sites) {
    
    # Filter data for the current site
    temp_data <- dataITS2[dataITS2$Site == s, ]
    
    # Create the plot
    p <- ggplot(temp_data, aes(x = replicate, y = absolute, fill = clade)) +
      geom_bar(stat = "identity") +
      labs(title = paste("Bar plot for site", s),
           x = "ID",
           y = "Concentration",
           fill = "ITS2 Type") +
      scale_fill_manual(values = color_palette) +
      theme(legend.position = "bottom", legend.title = element_text(vjust = 1)) +
      guides(fill = guide_legend(ncol = 2))+
      facet_grid(cols = vars(Season), rows = vars(combined), scales ="free")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
      plot_style
    
    
    # Save the plot to a file
    ggsave(filename = paste("ITS2profile_panama_10-11-23/Graph/plot_for_sit_", s, ".pdf"), width=18, height=45, dpi=300, units="cm",plot = p)
  }

for (s in sites) {
  
  # Filter data for the current site
  temp_data <- dataITS2[dataITS2$Site == s, ]
  
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
  ggsave(filename = paste("ITS2profile_panama_10-11-23/Graph/clade_for_sit_", s, ".jpg"), width=45, height=30, dpi=300, units="cm",plot = p)
}



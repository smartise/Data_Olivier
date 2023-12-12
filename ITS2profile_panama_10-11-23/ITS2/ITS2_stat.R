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

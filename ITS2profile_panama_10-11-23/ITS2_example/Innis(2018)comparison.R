library(data.table)
library(effects)
library(lsmeans)
library(devtools)
library(plyr)
library(reshape2)
library(RgoogleMaps)
library(plotrix)
library(zoo)
library(rgdal)
library(car)
library(scales)
library(png)
library(ecodist)
library(nnet)

# Import Sample Metadata
setwd("~/Desktop/TInnis comparison")

Coral_Data <- read.csv("Data/Coral_Collection.csv")
Coral_Data$Depth..m. <- as.numeric(as.character(Coral_Data$Depth..m.))

# Import qPCR Data

###qpcr analysis

library(tidyverse) #install.packages("tidyverse")
library(janitor)   #install.packages("janitor")
library(readxl)    #install.packages("readxl")

#Cunning et al. 2012 - Excess algal symbionts increase the susceptibility of reef corals to bleaching
#Cunning et al. 2017 - Patterns of bleaching and recovery of Montipora capitata in Kaneohe Bay, Hawaii, USA

################################### plate 1 ######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_1_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("\\+",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

a<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 2 ######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_2_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

b<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 3######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_3_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

c<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 4######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_4_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

d<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 5######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_5_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

e<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 6######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_6_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

f<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 7######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_7_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

g<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 8######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_8_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

h<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 9######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_9_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

i<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 10######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_10_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

j<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 11######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_11_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

k<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 12######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_12_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

l<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 13######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_13_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

m<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 14######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_14_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

n<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)


###################################plate 15######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_15_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

o<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)


###################################plate 16######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_16_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

p<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 17######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_17_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

q<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)


###################################plate 18######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_18_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

r<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

###################################plate 19######################################################
# upload, cleanup 
####set your working directory
rawdata<-read_excel("Data/qPCR_data/Symcap_19_data.xlsx")                                      ####import file
rawdata<-rawdata[-(1:6),]                                                       ####deleted first six rows of header information
colnames(rawdata) <- rawdata[1,]                                                ####set first row as column names
rawdata <- rawdata[-1, ]                                                        ####delete first row
format_data<-clean_names(rawdata); colnames(format_data)[7] <- "ct"; colnames(format_data)[8] <- "ct_mean"; colnames(format_data)[9] <- "ct_sd"  ### this renames columns manually, there are some weird capital subscripts in column names, so this fixes them

##set arguments 
copy.n.C<-33                                             #from Cunning et al. 2017 Supplemental Materials
copy.n.D<-3                                              #from Cunning et al. 2017 Supplemental Materials
#copy.n.Mcap<-1                                           #from Cunning et al. 2017 Supplemental Materials
copy.n.ratioCD<-(copy.n.C/copy.n.D)
#copy.n.ratioCM<-(copy.n.C/copy.n.Mcap)
#copy.n.ratioDM<-(copy.n.D/copy.n.Mcap)
fluo.C<-2.26827                                          #from Cunning et al. 2017 Supplemental Materials
fluo.D<-0                                                #from Cunning et al. 2017 Supplemental Materials
#fluo.Mcap<-0.84815                                       #from Cunning et al. 2017 Supplemental Materials

# process data 

data<-format_data%>%  
  select(sample_name,target_name,ct_mean, ct_sd)%>%                                                                                         ###take only the columns we care about
  mutate_at(vars(starts_with("ct")),funs(as.numeric))%>%                                                                                    ###set all the ct columns to numeric values (they are originally stored as characters)
  group_by(sample_name,target_name)%>%                                                                                                      ###group technical replicates
  summarise_all(funs(mean(., na.omit = TRUE)))%>%                                                                                           ###take mean of technical replicates
  filter(!grepl("Control",sample_name))%>%
  filter(sample_name!="NA")

C<-filter(data,target_name=="C")
D<-filter(data,target_name=="D")
#Mcap<-filter(data,target_name=="Mcap")

s<-left_join(C,D,by="sample_name")%>%
  select(-target_name.x,-target_name.y)%>%                                                                                     ###final formatting, get rid of redundant sample name column from previous step
  dplyr::rename("d_mean" = "ct_mean.y")%>%                                                                                                               ###rename columns from rbind
  dplyr::rename(d_sd=ct_sd.y)%>%
  dplyr::rename(c_mean=ct_mean.x)%>%
  dplyr::rename(c_sd=ct_sd.x)%>%
  #rename(mcap_mean=ct_mean)%>%
  #rename(mcap_sd=ct_sd)%>%  
  mutate(c_mean=c_mean-fluo.C)%>%                                                          #from Cunning et al. 2017
  mutate(d_mean=d_mean-fluo.D)%>%                                                          #from Cunning et al. 2017
  mutate(presence=case_when(d_mean>0&is.na(c_mean)~"D",
                            is.na(d_mean)&c_mean>0~"C",
                            d_mean>0&c_mean>0~"CD"))%>%
  mutate(cd_ratio=(c_mean/d_mean))%>%                           #from Cunning et al. 2017                         ###calculate ratio, log of 2^(difference in ct values)
  #mutate(sh_ratio=sum(cm_ratio,dm_ratio,na.rm=TRUE))%>%                                    #from Cunning et al. 2012                         ###add C:host ratio + D:host ratio for total symbiont:host ratio
  #mutate(colony=as.numeric(str_extract(sample_name,"[0-9]{1,3}")))%>%                                                                        ###extract numbers from sample name, set as colony number
  #mutate(phenotype=case_when(colony%%2==0~"nonbleached",colony%%2==1~"bleached"))%>%                                                         ###set even numbers as nonbleached, odd as bleached. this works by using the remainder command (%%), if 0 it must be an even number
  mutate(prop_c=cd_ratio/(cd_ratio+1))%>%                                                  #from Cunning et al. 2012
  mutate(prop_c=case_when(presence=="C"~1,presence=="D"~0,TRUE ~ as.numeric(prop_c)))%>%             #from Cunning et al. 2012
  mutate(prop_d=1-prop_c)

################### MERGE ALL PLATE RESULTS #################

output<-dplyr::bind_rows(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s)%>%
  filter(sample_name!="-")%>%
  filter(sample_name!="+")%>%
  filter(sample_name!= "NTC")%>%
  select(sample_name,  prop_d)

output[is.na(output)] <- 0

outputclean<- output%>%
  dplyr::rename(Colony=sample_name)%>%
  mutate(Colony = as.numeric(Colony))

dataallqpcr<-left_join(outputclean, Coral_Data, by="Colony")%>%
  select(-Date, -Time, -Reef.Type, -Reef.Area, -Bay.Area, -Reef.ID)


############################### Making map ####################################################

setwd("~/Desktop/Map/")
library(sf);library(ggrepel);library(ggsn);library(cowplot);library(tidyverse);library(patchwork);library(jpeg);library(ggnewscale);library(ggcharts);library(spdep)
library(tidyverse);library(janitor)  

oahu<-st_transform(st_read("coast_n83.shp"),crs=4326)
cropped<-st_crop(oahu,xmin=-158.8,xmax=-157.4,ymin=20.8,ymax=21.8)
fringe<-st_read("Fringing Reef.shp")%>%select(id,geometry)%>%mutate(zone='Fringing Reef')
habitat<-st_read("haw_benthic_habitat.shp")%>%clean_names()%>%filter(zone=='Backreef'|zone=="Reef Flat")%>%select(id,zone,geometry)
patch<-st_read("Patches2.shp")%>%mutate(type="Reef")%>%clean_names()%>%select(id,zone,block,geometry)%>%mutate(zone='Patch Reef')


out<-bind_rows(st_transform(fringe,crs=4326),st_transform(habitat,crs=4326))%>%
  mutate(zone=case_when(zone=="Reef Flat"~ "Fringing Reef",
                        zone=="Reef Crest"~"Fringing Reef",
                        zone=="Backreef"|zone=="Fringing Reef"~"Back/Fringing Reef",
                        TRUE~as.character(zone)))

#lines here to 32 are all dealing with rotating the map, you can exclude them if you dont want the tall/skinny version

any(is.na(dataallqpcr$Longitude))

dataallqpcr[is.na(dataallqpcr)] <- 0

dataallclean<- dataallqpcr%>%
  filter(!(Longitude==0 | Latitude==0))

rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
tcrop<-st_geometry(st_transform(cropped%>%filter(COASTLIN_1==17),crs=32604))
tpoints<-st_geometry(st_transform(st_as_sf(dataallclean, coords = c("Longitude", "Latitude"),crs=4326),crs=32604))
tout<-st_geometry(st_transform(out,crs=32604))
tpatch<-st_geometry(st_transform(patch,crs=32604))
cntrdc = st_centroid(tcrop)

d2 = st_as_sf(((tcrop - cntrdc) *rot(0.75) *2 + cntrdc))
p2 = st_as_sf(((tpoints - cntrdc) *rot(0.75) *2 + cntrdc))%>%bind_cols(.,dataallclean)
t2 = st_as_sf(((tout - cntrdc) *rot(0.75) *2 + cntrdc))
pat2= st_as_sf(((tpatch - cntrdc) *rot(0.75) *2 + cntrdc))%>%bind_cols(.,patch)


###PLOT


plot_Innis<-ggplot()+
  geom_sf(data=t2,fill="lightgray",color="lightgray")+
  geom_sf(data=d2,fill="darkgray")+
  geom_sf(fill="lightgray",dat=pat2,color=NA)+
  geom_sf(aes(fill=prop_d),data=p2,size=2,pch=21,color="black")+
  theme_classic(base_size=8)+
  coord_sf(xlim=c(627000,636000),ylim=c(2338000,2365000))+
  theme(legend.position=c(0.16,0.9),
        legend.key.size=unit(0.2,"cm"),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank())+
  scale_fill_gradient(low="yellow",high="blue",name="D Innis")+
  annotate("text",x=634700,y=2364000,label="N",size=3)+
  geom_segment(aes(x=635000,xend=636200,y=2364500,yend=2366000),arrow = arrow(length = unit(0.2,"cm")))
# theme(plot.margin=margin(l=-0.8,unit="cm"))






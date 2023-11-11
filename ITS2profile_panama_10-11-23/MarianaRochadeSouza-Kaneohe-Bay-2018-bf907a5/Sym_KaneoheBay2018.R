# Clonality 2018 samples ########

########################################## Loading packages####################################################################################
dependencies <- c("tidyverse",
                  "plotly",
                  "here",
                  "readxl",
                  "janitor",
                  "vegan", 
                  "ggplot2", 
                  "dplyr"
)

#check if packages are installed - load if so, install+load if not)
for (i in dependencies) {
  if (i %in% row.names(installed.packages())){
    eval(bquote(library(.(i))))
    message(paste("loaded package",i))
  } else {
    install.packages(i)
    eval(bquote(library(.(i))))
    message(paste("installed and loaded package",i))
  }
}

setwd("C:\\Users\\olivi\\OneDrive\\OneDrive liège\\OBSIDIAN\\Project Panama\\Data\\ITS2profile_panama_10-11-23\\MarianaRochadeSouza-Kaneohe-Bay-2018-bf907a5")


######################################### Adding individual datasets ##############
## Symbiodiniaceae datasets
Profile_absolute<-readRDS("Symportal_proabsolute") #Profile absolute read numbers

Type_absolute<-readRDS("Symportal_abso") #Type absolute read numbers

Type_relative<-readRDS("Symportal_relat")#Type relative read numbers


## Metadata

meta<-readRDS("depthmeta")


#__making Symbiodiniaceae datasets tidy __

Relative<-Type_relative %>%
  pivot_longer(c("B1": "F3.1"), names_to = "symbtypes", values_to = "relative")%>%
  mutate(relative= as.numeric(relative)) %>%
  filter(relative != 0)  

Absolute<-Type_absolute %>%
  pivot_longer(c("B1": "F3.1"), names_to = "symbtypes", values_to = "absolute")%>%
  mutate(absolute= as.numeric(absolute)) %>%
  filter(absolute != 0)  

ProfileAb<-Profile_absolute%>%
  pivot_longer(c(2:32), names_to = "Profiles", values_to = "absoluteepro")%>%
  mutate(absoluteepro= as.numeric(absoluteepro)) %>%
  filter(absoluteepro != 0)  


#Merging all tables into a single table, with all information and saving it as RDS

Processeddata2018<-inner_join(Relative,Absolute,by=c("ID","symbtypes"))%>%
  left_join(.,ProfileAb,by=c("ID"))%>%
  left_join(.,meta,by=c("ID"))%>%
  mutate(block=as.factor(block),ID=as.factor(ID),site=as.factor(site))%>%
  mutate(clade = case_when(
    str_detect(symbtypes, "A") ~ "A",
    str_detect(symbtypes, "B") ~ "B",
    str_detect(symbtypes, "C") ~ "C",
    str_detect(symbtypes, "D") ~ "D",
    str_detect(symbtypes, "E") ~ "E",
    str_detect(symbtypes, "F") ~ "F", 
    TRUE ~ NA_character_))%>%
  filter(clade=="C"| clade=="D")

str(PProcesseddata2018)
saveRDS(Processeddata2018,"Processeddata2018")
 
######################################### Checking number of reads ###############################

absol_sum<-Type_absolute%>% #data set with absolute read number values per sample ID
  mutate(Totalrow=rowSums (.[2:284] )) %>%
  filter(Totalrow>=250)%>% #standardizing number of reads to 2 std above and below mean
  filter(Totalrow<=16560) %>%
  select (ID) 

Processeddata2018_clean<-left_join(absol_sum, Processeddata2018, by= "ID") #only keeping samples that follow the requirement for number of reads

######################################### nMDS per block ####################################
wide2018<- Type_relative%>%
  right_join(meta,., by="ID")%>%
  left_join(absol_sum, .)

numerical<-wide2018%>%
  select(5:278)%>%
  select_if(colSums(.)!=0)%>%mutate(sum=rowSums(.)) #selects only numerical (clade) abundance data

meta_nmds<-wide2018%>%
  select(1:4) %>%
  rename(Block=block)

set.seed(216)
NMDS=metaMDS(numerical, distance='bray') 
nmdsall<- bind_cols(meta_nmds, NMDS$points[,1],NMDS$points[,2]) %>%
  rename(NMDS1=5, NMDS2=6)

a<-ggplot(nmdsall, aes(x=NMDS1, y=NMDS2, color=Block)) + 
  stat_ellipse(aes(fill = Block,color=Block), geom = "polygon", alpha = 0.2) +
  scale_fill_manual(values=c("#4B0055", "#00588B", "#009B95" ,"#53CC67", "#FDE333"))+
  scale_color_manual(values=c("#4B0055", "#00588B", "#009B95" ,"#53CC67", "#FDE333"))+
  geom_point(alpha = 0.4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  geom_hline(yintercept = 0, linetype="dotted", color="grey") + 
  geom_vline(xintercept = 0, linetype="dotted", color="grey") +
  labs(color = "Blocks", fill = "Blocks")  +
  theme(plot.title = element_text(hjust = 0.5))+
  annotate(geom="text", x= -0.65, y=-0.50, label="Stress=0.062",
           color="black"); a

######################################### Permanova per block ######################################################

install.packages("devtools")
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

siteblock<-wide2018%>%
  select(1:4)%>%#p.env contains only the environmental (independant) variables 
  mutate(site=as.factor(site))%>%
  mutate(block=as.factor(block))

#running the permanova test

#for all sites in block (general)
Permanova1 <- adonis(numerical ~ site:block+block,data=siteblock, permutations=999, method="bray"); Permanova1 

#for pairwise permanova 
library(pairwiseAdonis)
permanova<-pairwise.adonis(numerical,factors=siteblock$block, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni',reduce=NULL,perm=999)

######################################### Heatmap #######################################

heatmap<-read_excel("PermanovaR2.xlsx")%>%
  select(c(2:6))

pheatmap(heatmap, color=colorRampPalette(c("Orange" ,"white"))(50))

######################################### Barplot Symbiondiniaceae community per block ###########################

####_________________________Symbiont types ________________________________________
Data2018<-Processeddata2018_clean %>% 
  group_by(ID)%>% 
  mutate(totalreads=sum(absolute)) %>% 
  group_by(ID, clade) %>% 
  mutate(cladereads=sum(absolute)) %>% 
  mutate(proportion=cladereads/totalreads) %>% 
  spread(clade,proportion)%>%
  dplyr::rename(prop_c=C,prop_d=D)%>%
  mutate(prop_d = replace_na(prop_d, 0)) %>%
  mutate(prop_c = 1-prop_d) 


Majorsymbtypes18 <-Data2018%>% 
  mutate(typessummary = case_when(
    str_detect(symbtypes, "D1") ~ "D1",
    str_detect(symbtypes, "D4") ~ "D4",
    str_detect(symbtypes, "D6") ~ "D6",
    str_detect(symbtypes, "C17") ~ "C17",
    str_detect(symbtypes, "D3") ~ "D3", 
    str_detect(symbtypes, "_D") ~ "D", 
    str_detect(symbtypes, "C31") ~ "C31", 
    str_detect(symbtypes, "C21") ~ "C21", 
    str_detect(symbtypes, "C15") ~ "C15", 
    str_detect(symbtypes, "_C") ~ "C", 
    str_detect(symbtypes, "C3") ~ "C3", 
    str_detect(symbtypes, "C1") ~ "C1", 
    str_detect(symbtypes, "D2") ~ "D2" ))%>% 
  group_by(site)%>% 
  mutate(totalreadstype=sum(absolute)) %>% 
  group_by(site, typessummary) %>% 
  mutate(summarytypes=sum(absolute)) %>% 
  mutate(proportionsymb=summarytypes/totalreadstype) %>% 
  spread(typessummary,proportionsymb)%>%
  mutate(D1 = replace_na(D1, 0)) %>%
  mutate(D4 = replace_na(D4, 0))   %>%  
  mutate(D6 = replace_na(D6, 0))   %>%
  mutate(C17 = replace_na(C17, 0))   %>%
  mutate(D3 = replace_na(D3, 0))   %>%
  mutate(D = replace_na(D, 0))%>%
  mutate(C31 = replace_na(C31, 0))%>%
  mutate(C21 = replace_na(C21, 0))%>%
  mutate(C15 = replace_na(C15, 0))%>%
  mutate(C = replace_na(C, 0))%>%
  mutate(C3 = replace_na(C3, 0))%>%
  mutate(C1 = replace_na(C1, 0))%>%
  mutate(D2 = replace_na(D2, 0))%>%
  dplyr::select(D1,D4,D6,C17,D3,D,C31,C21,C15,C,C3,C1,D2, site, block) %>% 
  distinct() 

#making it tidy

tidysym<-Majorsymbtypes18 %>%
  pivot_longer(c("D1": "D2"), names_to = "majortypes", values_to = "proportion")%>%
  mutate(proportion= as.numeric(proportion)) %>%
  filter(proportion != 0)  


colours1=c("#14505C" ,"#226B70" ,"#348781", "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           "#B13F63", "#CC5762", "#DE7461", "#E8956D", "#F0B384", "#F5D1A8")

quartz(height = 4)
mx4 = ggplot(tidysym, aes(x = site, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") + 
  facet_wrap (~block, ncol=5, scales="free") +
  theme(axis.text.x = element_text(angle = 90, size = 7, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 9),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm'), legend.margin=margin(-0.4,0,0,0,unit="cm"))  + 
  scale_y_continuous(expand = c(0,0)) + 
  #ggtitle(label = "A") +
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours1) +
  guides(fill = guide_legend(nrow = 2)); mx4

###_____________________________________Symbiont profiles

Data18Pro<-Processeddata2018_clean %>% 
  drop_na(Profiles)%>%
  group_by(site)%>% 
  mutate(totalreadstype=sum(absoluteepro))%>%
  group_by(site, Profiles) %>% 
  mutate(summarytypes=sum(absoluteepro)) %>% 
  mutate(proportionsymb=summarytypes/totalreadstype)%>%
  select(site,Profiles,proportionsymb, block)%>%
  distinct()

colours2<-c( "#14505C", "#175762" ,"#1B5F67", "#1F666C", "#236E71", "#287576", "#2D7D7B", "#328580", "#388C84", "#3E9488", "#459B8C", "#4CA38F", "#54AA92", "#5CB295", "#65B998", "#6EC09A", "#7AC79D", "#88CCA1", "#96D2A6", "#A3D8AB",
             "#AFDDB1", "#BBE2B7", "#C7E5BE", "#D6765D" ,"#F8DFC1", "#772C4B")

quartz(height = 6)
mx5 = ggplot(Data18Pro, aes(x = site, fill = Profiles, y = proportionsymb)) + 
  geom_bar(stat = "identity") + 
  facet_wrap (~block, ncol=5, scales="free") +
  theme(axis.text.x = element_text(angle = 90, size = 7, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 9), legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 7),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm'), legend.margin=margin(-0.4,0,0,0,unit="cm")) +
  scale_y_continuous(expand = c(0,0)) + 
  #ggtitle(label = "B") +
  labs(x = "", y = "Relative proportion", fill = "Profiles") +
  scale_fill_manual(values = colours2) +
  guides(fill = guide_legend(nrow = 10)); mx5

#___making a panel

install.packages(cowplot)
quartz()
plot_grid(mx4, NULL,mx5, labels = c('A','', 'B') ,align = "v" ,label_size = 12, ncol=1, rel_heights = c(1,0,1.3))

######################################### dbRDA############################################
#_____adding sedimentation and temperature data _______________________________________

mastersedimentation<-read_excel("mastersedimentation.xlsx")%>%
  select(c("site","S_mean","S_min", "S_max", "S_range", "Std"))  #sedimentation

mastertemperature<-read_excel("tempsummary.xlsx")#temperature stats
  
environmentalmeta<-left_join(mastertemperature, mastersedimentation, by="site") #merged

Symbmeta<-left_join(wide2018, environmentalmeta, by= "site") #adding symbiont data

CD<-Processeddata2018_clean%>% 
  group_by(ID)%>% 
  mutate(totalreads=sum(absolute)) %>% 
  group_by(ID, clade) %>% 
  mutate(cladereads=sum(absolute)) %>% 
  mutate(proportion=cladereads/totalreads)%>%
  select(ID,clade,proportion)%>%distinct()%>%
  spread(clade,proportion)%>%
  mutate(D = replace_na(D, 0))%>%
  mutate(C=replace_na(C, 0))%>%
  mutate(Majorclade=case_when(
    C>0.8~"C",
    D>0.8~"D",
    TRUE~"CD"))%>%
  select(ID, Majorclade)%>%
  distinct() #calculatuing majority CD proportion

SymbmetaCD<-left_join(Symbmeta, CD, by="ID")

Symbmeta2018 <- SymbmetaCD[c(1,2,4, 300,3, 288:299, 5:287)] #changing order


species<-Symbmeta2018[,18:300] 
environment= Symbmeta2018[,5:17] 
scaled_clim = scale(environment)
species001= (species + 0.001)

str(species001)
str(environment)
rankindex(environment, species001, indices = c("euc", "man", "gow","bra", "kul"), stepacross= FALSE, method = "spearman") #Bra is 0.07

dbRDA = capscale(species001 ~ depth_m + T_DHW + T_mean + T_max + T_min + T_range + T_drange + T_dsd +  S_mean+ S_min +S_max  + Std ,environment, dist="bray")


#Permutation tests to access significant of constraints
set.seed(999)
anova(dbRDA) # overall test of the significant of the analysis 0.001
anova(dbRDA, by="terms", permu=999) # test for sign. environmental. variables 

#making the graph in ggplot
x <- as.data.frame(scores(dbRDA, display = "sites"))
Symbmeta2018$CAP1 <- x$CAP1
Symbmeta2018$CAP2 <- x$CAP2 
d <- data.frame(Variables = rownames(dbRDA$CCA$biplot), dbRDA$CCA$biplot)%>%
  filter( Variables=="T_max"| Variables=="depth_m"| Variables=="T_DHW" | Variables== "T_dsd")


quartz()
i<-ggplot(Symbmeta2018, aes(x= CAP1, y= CAP2, color=Majorclade)) + 
  stat_ellipse(aes(fill = Majorclade,color=Majorclade), geom = "polygon", alpha = 0.2) +
  scale_fill_manual(values=c("#14505C", "orange1", "#772C4B"))+
  scale_color_manual(values=c("#14505C", "orange1", "#772C4B"))+
  geom_point(alpha = 0.4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  coord_cartesian(xlim=c(-4, 5), ylim=c(-7, 8))  +
  geom_hline(yintercept = 0, linetype="dotted", color="grey") + 
  geom_vline(xintercept = 0, linetype="dotted", color="grey") +
  labs(color = "Major Clades", fill = "Major Clades")  +
  theme(plot.title = element_text(hjust = 0.5))+ 
  geom_segment(data = d, aes(x = 0, y = 0, xend = (CAP1*5),yend = (CAP2*5)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +
  annotate("text", x = (d$CAP1*5), y = (d$CAP2*5),label = d$Variables);i

######################################### Correlation env factors #########################
depth = Symbmeta2018$depth_m.x
maxtemp = Symbmeta2018$T_max     
mintemp = Symbmeta2018$T_min
DHW = Symbmeta2018$T_DHW
mintemp = Symbmeta2018$T_min
depth = Symbmeta2018$depth_m
maxtemp = Symbmeta2018$T_max 
mean = Symbmeta2018$T_mean
dailyrange = Symbmeta2018$T_drange
T_dsd = Symbmeta2018$T_dsd
S_mean= Symbmeta2018$S_mean
S_max= Symbmeta2018$S_max
S_min= Symbmeta2018S_min
Std= Symbmeta2018$Std
S= Symbmeta2018$Std


res_mintemp <- cor.test(depth, mintemp, 
                method = "pearson")

res_maxtemp <- cor.test(depth, maxtemp, 
                method = "pearson")

res_DHW <- cor.test(depth, DHW, 
                method = "pearson")

res_mean <- cor.test(depth, mean, 
                method = "pearson")

res_dailyrange <- cor.test(depth, dailyrange, 
                method = "pearson")

res_T_dsd <- cor.test(depth, T_dsd, 
                method = "pearson")

res_S_mean <- cor.test(depth, S_mean, 
                method = "pearson")

res_S_max <- cor.test(depth, S_max, 
                method = "pearson")

res_S_min <- cor.test(depth, S_min, 
                method = "pearson")

res_Std <- cor.test(depth, Std, 
                method = "pearson")

######################################### Barplot Symbiodiniaceae in each colony per site #############

Majorsymbtypes_all <-Data2018%>% 
  mutate(typessummary = case_when(
    str_detect(symbtypes, "D1") ~ "D1",
    str_detect(symbtypes, "D4") ~ "D4",
    str_detect(symbtypes, "D6") ~ "D6",
    str_detect(symbtypes, "C17") ~ "C17",
    str_detect(symbtypes, "D3") ~ "D3", 
    str_detect(symbtypes, "_D") ~ "D", 
    str_detect(symbtypes, "C31") ~ "C31", 
    str_detect(symbtypes, "C21") ~ "C21", 
    str_detect(symbtypes, "C15") ~ "C15", 
    str_detect(symbtypes, "_C") ~ "C", 
    str_detect(symbtypes, "C3") ~ "C3", 
    str_detect(symbtypes, "C1") ~ "C1", 
    str_detect(symbtypes, "D2") ~ "D2" ))%>% 
  group_by(ID)%>% 
  mutate(totalreadstype=sum(absolute)) %>% 
  group_by( ID, typessummary) %>% 
  mutate(summarytypes=sum(absolute)) %>% 
  mutate(proportionsymb=summarytypes/totalreadstype) %>% 
  spread(typessummary,proportionsymb)%>%
  mutate(D1 = replace_na(D1, 0)) %>%
  mutate(D4 = replace_na(D4, 0))   %>%  
  mutate(D6 = replace_na(D6, 0))   %>%
  mutate(C17 = replace_na(C17, 0))   %>%
  mutate(D3 = replace_na(D3, 0))   %>%
  mutate(D = replace_na(D, 0))%>%
  mutate(C31 = replace_na(C31, 0))%>%
  mutate(C21 = replace_na(C21, 0))%>%
  mutate(C15 = replace_na(C15, 0))%>%
  mutate(C = replace_na(C, 0))%>%
  mutate(C3 = replace_na(C3, 0))%>%
  mutate(C1 = replace_na(C1, 0))%>%
  mutate(D2 = replace_na(D2, 0))%>%
  dplyr::select(D1,D4,D6,C17,D3,D,C31,C21,C15,C,C3,C1,D2, site, ID, block) %>% 
  distinct()  


#_________________________________Block 1 _____________________________________

tidysym1<-Majorsymbtypes_all %>%
  pivot_longer(c("D1": "D2"), names_to = "majortypes", values_to = "proportion")%>%
  mutate(proportion= as.numeric(proportion)) %>%
  filter(proportion != 0) %>%
  filter (block==1)

tidysym1$site_f = factor(tidysym1$site, levels=c('1_2','1_3','1_4','1_6', '1_9', '1_10'))


colours1=c("#14505C" ,"#226B70" ,"#348781", "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762", "#DE7461","#E8956D", "#F0B384", "#F5D1A8")



quartz(height = 4)
p1 = ggplot(tidysym1, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site_f, ncol=3, scales="free") +
  theme(axis.title.y = element_text(size = 9),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 9),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  ggtitle(label = "A")+
  guides(fill = guide_legend(nrow = 1)) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours1); p1

###_____________________________Block 2 _______________________________________

tidysym2<-Majorsymbtypes_all %>%
  pivot_longer(c("D1": "D2"), names_to = "majortypes", values_to = "proportion")%>%
  mutate(proportion= as.numeric(proportion)) %>%
  filter(proportion != 0) %>%
  filter (block==2)

colours2=c("#14505C" ,"#226B70" , "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762", "#E8956D", "#F0B384", "#F5D1A8")

quartz(height = 4)
p2 = ggplot(tidysym2, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site, ncol=3, scales="free") +
  theme(axis.title.y = element_text(size = 9),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 9),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  ggtitle(label = "B") + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours2) +
  guides(fill = guide_legend(nrow = 1)) ; p2

###_____________________________Block 3 _______________________________________

tidysym3<-Majorsymbtypes_all %>%
  pivot_longer(c("D1": "D2"), names_to = "majortypes", values_to = "proportion")%>%
  mutate(proportion= as.numeric(proportion)) %>%
  filter(proportion != 0) %>%
  filter (block==3)



colours3=c("#14505C" ,"#226B70" ,"#348781", "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762", "#DE7461","#E8956D", "#F0B384", "#F5D1A8")


quartz(height = 4)
p3 = ggplot(tidysym3, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site, ncol=3, scales="free") +
  theme(axis.title.y = element_text(size = 9),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 9),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  ggtitle(label = "C") + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours3) +
  guides(fill = guide_legend(nrow = 1)) ; p3

###_____________________________Block 4 _______________________________________

tidysym4<-Majorsymbtypes_all  %>%
  pivot_longer(c("D1": "D2"), names_to = "majortypes", values_to = "proportion")%>%
  mutate(proportion= as.numeric(proportion)) %>%
  filter(proportion != 0) %>%
  filter (block==4)



colours4=c("#14505C" ,"#226B70" , "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762","#E8956D", "#F0B384", "#F5D1A8")


quartz(height = 4)
p4 = ggplot(tidysym4, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site, ncol=3, scales="free") +
  theme(axis.title.y = element_text(size = 9),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 9),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  ggtitle(label = "D") + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours4) +
  guides(fill = guide_legend(nrow = 1)) ; p4

###_____________________________Block 5 _______________________________________

tidysym5<-Majorsymbtypes_all  %>%
  pivot_longer(c("D1": "D2"), names_to = "majortypes", values_to = "proportion")%>%
  mutate(proportion= as.numeric(proportion)) %>%
  filter(proportion != 0) %>%
  filter (block==5)



colours5=c("#14505C" ,"#226B70", "#4CA38F", "#6BBE99", "#9AD4A8", "#C7E5BE", 
           
           "#B13F63", "#CC5762","#E8956D", "#F0B384", "#F5D1A8")


quartz(height = 4)
p5 = ggplot(tidysym5, aes(x = ID, fill = majortypes, y =proportion)) + 
  geom_bar(stat = "identity") +
  facet_wrap (~site, ncol=3, scales="free") +
  theme(axis.title.y = element_text(size = 9),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size = 7), 
        axis.text.y = element_text(colour = "black", size = 9),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position ="bottom",
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, 'cm')) +
  ggtitle(label = "E") + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative proportion", fill = "ITS2 types") +
  scale_fill_manual(values = colours5) +
  guides(fill = guide_legend(nrow = 1)) ; p5






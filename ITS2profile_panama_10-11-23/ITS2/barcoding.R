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


#######################################load the data########################################
setwd("C:\\Users\\olivi\\OneDrive\\OneDrive liège\\OBSIDIAN\\Project Panama\\Data\\")
barcode123 <- read_excel("ITS2profile_panama_10-11-23\\ITS2\\ITS2profile_meta_index.barcode.P123_30.07.22.xlsx")
barcode456 <- read_excel("ITS2profile_panama_10-11-23\\ITS2\\ITS2profile_meta_index.barcode.P456_22.08.2022.xlsx")
fasta_doc <- read_excel("ITS2profile_panama_10-11-23\\ITS2\\ITS2profile_meta_fastq.file.name_23.10.23.xlsx")

colnames(barcode123) <- as.character(unlist(barcode123[10, ]))
barcode123 <- barcode123[-c(1:10),c(1:5)]
barcode123$Plate <- "Plate123"

colnames(barcode456) <- as.character(unlist(barcode456[10, ]))
barcode456<- barcode456[-c(1:10),c(1:5)]
barcode456$Plate <- "Plate456"

plate_marta <- rbind(barcode123, barcode456)

colnames(plate_marta)[1] <- "ID"
colnames(plate_marta)[2] <- "V4"
colnames(plate_marta)[4] <- "V3"
plate_marta$ID <- paste0("ML0", plate_marta$ID)

######################################matrix of the PCR plate####################################

PCRplate <- as.matrix(read_excel("General_meta\\Meta_PCRplate.xlsx"))
PCRplate <- PCRplate [,-1]

#####################################matrix for the sample#######################################
excel_file <- "ITS2profile_panama_10-11-23\\ITS2\\ITS2profile_meta_Platebarcoding_14.11.23.xlsx"
sheet_names <- excel_sheets("ITS2profile_panama_10-11-23\\ITS2\\ITS2profile_meta_Platebarcoding_14.11.23.xlsx")
Samples <- list()

for (sheet_name in sheet_names) {
  sheet_data <- read_excel(excel_file, sheet = sheet_name)
  Samples[[sheet_name]] <- as.matrix(sheet_data)
}
#####################################matrix for the barcode + PCR plate#########################
excel_file <- "General_meta\\Meta_PCRplate.xlsx"
sheet_names <- excel_sheets("General_meta\\Meta_PCRplate.xlsx")
barcode <- list()

for (sheet_name in sheet_names) {
  sheet_data <- read_excel(excel_file, sheet = sheet_name)
  barcode[[sheet_name]] <- as.matrix(sheet_data)
}
####################################superposition of the sample/barcode###########################

Plate1 <- list(Samples[["Plate1"]], barcode[["PCRplate"]], barcode[["barcodeBF"]], barcode[["barcodeBR"]])
Plate1 <- do.call(cbind, mapply(function(mat) as.vector(mat), Plate1, SIMPLIFY = FALSE))
Plate2 <- list(Samples[["Plate2"]], barcode[["PCRplate"]], barcode[["barcodeCF"]], barcode[["barcodeCR"]])
Plate2 <- do.call(cbind, mapply(function(mat) as.vector(mat), Plate2, SIMPLIFY = FALSE))
Plate3 <- list(Samples[["Plate3"]], barcode[["PCRplate"]], barcode[["barcodeDF"]], barcode[["barcodeDR"]])
Plate3 <- do.call(cbind, mapply(function(mat) as.vector(mat), Plate3, SIMPLIFY = FALSE))
Plate123 <- as.data.frame(rbind(Plate1, Plate2, Plate3))
Plate123$Plate <- "Plate123"

Plate4 <- list(Samples[["Plate4"]], barcode[["PCRplate"]], barcode[["barcodeBF"]], barcode[["barcodeBR"]])
Plate4 <- do.call(cbind, mapply(function(mat) as.vector(mat), Plate4, SIMPLIFY = FALSE))
Plate5 <- list(Samples[["Plate5"]], barcode[["PCRplate"]], barcode[["barcodeCF"]], barcode[["barcodeCR"]])
Plate5 <- do.call(cbind, mapply(function(mat) as.vector(mat), Plate5, SIMPLIFY = FALSE))
Plate6 <- list(Samples[["Plate6"]], barcode[["PCRplate"]], barcode[["barcodeDF"]], barcode[["barcodeDR"]])
Plate6 <- do.call(cbind, mapply(function(mat) as.vector(mat), Plate6, SIMPLIFY = FALSE))
Plate456 <- as.data.frame(rbind(Plate4, Plate5, Plate6))
Plate456$Plate <- "Plate456"

meta_barcoding <- rbind(Plate123, Plate456)
meta_barcoding

rm(list = setdiff(ls(), c("meta_barcoding", "plate_marta"))) # i keep what i need
####################################the one that i gave marta bouhouhou###########################

merged_barcoding <- meta_barcoding %>%
  inner_join(plate_marta, by = c("V3", "V4", "Plate")) #merge de dataframe 
merged_barcoding$barcodeI7 <- toupper(merged_barcoding$barcodeI7)

different_rows <- which(merged_barcoding$V1 != merged_barcoding$ID)
rows_with_differences <- merged_barcoding[different_rows, ]

####################################correction using the Tags#####################################

correction <- read_csv("ITS2profile_panama_10-11-23\\ITS2\\SampleSheetUsed.csv")
correction <- correction[,c(1,4,6)]
colnames(correction)[2] <- "barcodeI7"
colnames(correction)[3] <- "barcodeI5"
correction$Sample_ID <- paste0("ML0", correction$Sample_ID)

####################################for plate 456########################################

merged_barcoding_P456 <- merged_barcoding[merged_barcoding$Plate != "Plate123", ]
merged_corr <- correction %>% inner_join(merged_barcoding_P456, by = c("barcodeI7", "barcodeI5"))

different_rows <- which(merged_corr$Sample_ID != merged_corr$V1)
rows_with_differences <- merged_corr[different_rows, ]

different_rows <- which(merged_corr$Sample_ID != merged_corr$V1)


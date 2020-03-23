#Data pre-processing for Guppy data
#Network inference
library(openxlsx)
## Data loading
## Load normalized counts for genes used in DE analysis

###############
#Aripo
counts_genome_Aripo <-read.csv(file = "~/Dropbox/Project_Kim/data/OLD_counts_matrix.csv", 
                        row.names=1, stringsAsFactors = FALSE) %>% ceiling()
# Old behaviroal data
behave_Aripo <- read.xlsx(xlsxFile = "~/Dropbox/Project_Kim/data/behavioral_data_updated.xlsx",1)

#Quare
counts_genome_Quare  <- read.csv(file = "~/Dropbox/Project_Kim/data/NEW_counts_matrix.csv", 
                          row.names=1, stringsAsFactors = FALSE) %>% ceiling()
# New behaviroal data
behave_Quare <- read.xlsx(xlsxFile = "~/Dropbox/Project_Kim/data/behavioral_data_2.0.xlsx",1)

###
#Behavioral data cleaning
behave_Aripo$fish[1:9] <- paste("X00",behave_Aripo$fish[1:9],sep="")
behave_Aripo$fish[-(1:9)] <- paste("X0",behave_Aripo$fish[-(1:9)],sep="")
rownames(behave_Aripo)<- behave_Aripo$fish
behave_Aripo$pop <- factor(behave_Aripo$pop,levels = c("LP","HP"))
behave_Aripo$rear <- factor(behave_Aripo$rear, levels  = c("P","NP"))
behave_Aripo$group = factor(paste0(behave_Aripo$pop,behave_Aripo$rear))
behave_Aripo <- behave_Aripo[colnames(counts_genome_Aripo),]

# Add cross effect in behaviroal data
behave_Aripo$cross <- 1:40
behave_Aripo$cross[behave_Aripo$family=="HP216"] <- 31
behave_Aripo$cross[behave_Aripo$family=="LP205"] <- 16
behave_Aripo$cross <- factor(behave_Aripo$cross)
behave_Aripo$cross
behave_Aripo$group = paste(behave_Aripo$pop,behave_Aripo$rear,sep="")

behave_Quare<- behave_Quare[match(colnames(counts_genome_Quare),behave_Quare$fish),]

#Remove problematic samples - M17, M15
counts_genome_Quare <- counts_genome_Quare[,-c(14,16)]
behave_Quare<-behave_Quare[-c(14,16),]

behave_Quare$cross = 1:58
behave_Quare$cross[behave_Quare$family %in% c("207B")] = 38
behave_Quare$pop = ifelse(behave_Quare$pop == "CM","LP","HP") %>% factor(levels=c("LP","HP"))
behave_Quare$rear = factor(behave_Quare$rear,levels=c("P","NP"))
behave_Quare$group = factor(paste(behave_Quare$pop,behave_Quare$rear,sep=""))

behave_Quare["42",]$family = "207A" # Relabel wrong familiy
behave_Quare$family_new = substr(behave_Quare$family,1,3)
behave_Quare$family_new[which(behave_Quare$family_new =="207")]=behave_Quare$family[which(behave_Quare$family_new =="207")]
behave_Quare$group = paste(behave_Quare$pop,behave_Quare$rear,sep="")


#Define group names
#In Aripo
LPP_Aripo = which(behave_Aripo$group=="LPP")
LPNP_Aripo = which(behave_Aripo$group=="LPNP")
HPP_Aripo = which(behave_Aripo$group=="HPP")
HPNP_Aripo = which(behave_Aripo$group=="HPNP")
#In Quare
LPP_Quare = which(behave_Quare$group=="LPP")
LPNP_Quare = which(behave_Quare$group=="LPNP")
HPP_Quare = which(behave_Quare$group=="HPP")
HPNP_Quare = which(behave_Quare$group=="HPNP")
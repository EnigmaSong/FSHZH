#Chi-square tests

#read xlsx

# Aripo_DE_res_genome = read.csv("output/Aripo_glmer_effects_genome.csv")
# Quare_DE_res_genome = read.csv("output/Quare_glmer_genome.csv")

library(openxlsx)

Aripo_DE_res=read.xlsx("../../outputs/DE_GLMM_genome.xlsx",sheet=2)
Quare_DE_res=read.xlsx("../../outputs/DE_GLMM_genome.xlsx",sheet=3)

#If we exclude significant interaction genes, 
merged_DE_res = merge(Aripo_DE_res[!Aripo_DE_res$inter_DE_flag,],
                      Quare_DE_res[!Quare_DE_res$inter_DE_flag,],
                      # by = "gene",
                      by = "Id",
                      suffixes = c("_Aripo","_Quare"))
##Re-obtain direction and DE/NDE information 
Aripo_pop_ud2 = ifelse(merged_DE_res$pop_DE_flag_Aripo & unlist(merged_DE_res$pop_estim_Aripo >0), "DEU", 
                       ifelse(merged_DE_res$pop_DE_flag_Aripo & unlist(merged_DE_res$pop_estim_Aripo <0), "DED",
                              ifelse(!merged_DE_res$pop_DE_flag_Aripo & unlist(merged_DE_res$pop_estim_Aripo >0), 
                                     "NDEU","NDED")))
Aripo_rear_ud2 = ifelse(merged_DE_res$rear_DE_flag_Aripo & unlist(merged_DE_res$rear_estim_Aripo >0), "DEU", 
                        ifelse(merged_DE_res$rear_DE_flag_Aripo & unlist(merged_DE_res$rear_estim_Aripo <0), "DED",
                               ifelse(!merged_DE_res$rear_DE_flag_Aripo & unlist(merged_DE_res$rear_estim_Aripo >0), 
                                      "NDEU","NDED")))

Quare_pop_ud2 = ifelse(merged_DE_res$pop_DE_flag_Quare & unlist(merged_DE_res$pop_estim_Quare >0), "DEU", 
                       ifelse(merged_DE_res$pop_DE_flag_Quare & unlist(merged_DE_res$pop_estim_Quare <0), "DED",
                              ifelse(!merged_DE_res$pop_DE_flag_Quare & unlist(merged_DE_res$pop_estim_Quare >0), 
                                     "NDEU","NDED")))
Quare_rear_ud2 = ifelse(merged_DE_res$rear_DE_flag_Quare & unlist(merged_DE_res$rear_estim_Quare >0), "DEU", 
                        ifelse(merged_DE_res$rear_DE_flag_Quare & unlist(merged_DE_res$rear_estim_Quare <0), "DED",
                               ifelse(!merged_DE_res$rear_DE_flag_Quare & unlist(merged_DE_res$rear_estim_Quare >0), 
                                      "NDEU","NDED")))
###############
#Pop: Aripo vs Quare
table(Aripo_pop_ud2,Quare_pop_ud2)
chisq.test(table(Aripo_pop_ud2,Quare_pop_ud2),
           simulate.p.value = TRUE)#Up/Down regulated
#DE vs NDE in two data set
aa=table(Aripo_pop_ud2,Quare_pop_ud2)
chisq.test(aa,
           simulate.p.value = TRUE)#DE vs NDE in two data set
fisher.test(aa,
            simulate.p.value = TRUE)#DE vs NDE in two data set

aaa= matrix(c(sum(aa[1:2,1:2]),sum(aa[3:4,1:2]),sum(aa[1:2,3:4]),sum(aa[3:4,3:4])) ,2,2)
chisq.test(aaa,
           simulate.p.value = TRUE)#DE vs NDE in two data set
fisher.test(aaa,
            simulate.p.value = TRUE)#DE vs NDE in two data set

aaaa= matrix(c(sum(aa[c(1,3),c(1,3)]),sum(aa[c(2,4),c(1,3)]),sum(aa[c(1,3),c(2,4)]),sum(aa[c(2,4),c(2,4)])) ,2,2)
chisq.test(aaaa,
           simulate.p.value = TRUE)#DE vs NDE in two data set
fisher.test(aaaa,
            simulate.p.value = TRUE)#DE vs NDE in two data set

aaaaa = matrix(c(sum(aa[c(1,3),c(1,3)])-aa[3,3], sum(aa[c(2,4),c(1,3)])-aa[4,3],
                 sum(aa[c(1,3),c(2,4)])-aa[3,4],sum(aa[c(2,4),c(2,4)])-aa[4,4])
               ,2,2)
chisq.test(aaaaa, simulate.p.value = TRUE)#
fisher.test(aaaaa,
            simulate.p.value = TRUE)#DE vs NDE in two data set
chisq.test(aa[1:2,1:2],simulate.p.value=TRUE)
fisher.test(aa[1:2,1:2],simulate.p.value=TRUE)
###############
#Rear: Aripo vs Quare
aa= table(Aripo_rear_ud2,Quare_rear_ud2)
chisq.test(aa,
           simulate.p.value = TRUE)#Up/Down regulated
fisher.test(aa,
            simulate.p.value = TRUE)#Up/Down regulated

aaa= matrix(c(sum(aa[1:2,1:2]),sum(aa[3:4,1:2]),sum(aa[1:2,3:4]),sum(aa[3:4,3:4])) ,2,2)
chisq.test(aaa,
           simulate.p.value = TRUE)#DE vs NDE in two data set
fisher.test(aaa,
            simulate.p.value = TRUE)#DE vs NDE in two data set

aaaa= matrix(c(sum(aa[c(1,3),c(1,3)]),sum(aa[c(2,4),c(1,3)]),sum(aa[c(1,3),c(2,4)]),sum(aa[c(2,4),c(2,4)])) ,2,2)
chisq.test(aaaa,
           simulate.p.value = TRUE)#DE vs NDE in two data set
fisher.test(aaaa,
            simulate.p.value = TRUE)#DE vs NDE in two data set

aaaaa = matrix(c(sum(aa[c(1,3),c(1,3)])-aa[3,3], sum(aa[c(2,4),c(1,3)])-aa[4,3],
                 sum(aa[c(1,3),c(2,4)])-aa[3,4],sum(aa[c(2,4),c(2,4)])-aa[4,4])
               ,2,2)
chisq.test(aaaaa, simulate.p.value = TRUE)#
fisher.test(aaaaa,
            simulate.p.value = TRUE)#DE vs NDE in two data set

chisq.test(aa[1:2,1:2],simulate.p.value=TRUE)
fisher.test(aa[1:2,1:2],simulate.p.value=TRUE)


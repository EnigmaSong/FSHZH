# Chi-square test for independence (excluding significant interaction genes)
# Pop vs rearing in the same data set
# Pop in two data set

######################
# Pop vs rearing in the same data set
#Aripo: Excluding significant interaction gene
Aripo_pop_DE = unlist(Aripo_DE_res$pop_DE_flag)
Aripo_rear_DE = unlist(Aripo_DE_res$rear_DE_flag)

#Quare: Excluding significant interaction gene
Quare_pop_DE = unlist(Quare_DE_res$pop_DE_flag)
Quare_rear_DE = unlist(Quare_DE_res$rear_DE_flag)

##Including up&down regulated
#Aripo
Aripo_pop_ud = ifelse(Aripo_pop_DE & unlist(Aripo_DE_res$pop_estim >0), "DEU", 
                      ifelse(Aripo_pop_DE & unlist(Aripo_DE_res$pop_estim <0), "DED",
                             ifelse(!Aripo_pop_DE & unlist(Aripo_DE_res$pop_estim >0), "NDEU","NDED")))
Aripo_rear_ud = ifelse(Aripo_rear_DE & unlist(Aripo_DE_res$rear_estim >0), "DEU", 
                       ifelse(Aripo_rear_DE & unlist(Aripo_DE_res$rear_estim <0), "DED",
                              ifelse(unlist(Aripo_DE_res$rear_estim >0), "NDEU","NDED")))

#Excluding interaction genes
table(Aripo_pop_ud[!Aripo_DE_res$inter_DE_flag],Aripo_rear_ud[!Aripo_DE_res$inter_DE_flag])
chisq.test(table(Aripo_pop_ud[!Aripo_DE_res$inter_DE_flag],
                 Aripo_rear_ud[!Aripo_DE_res$inter_DE_flag]),
           simulate.p.value = TRUE)
# table(Aripo_pop_ud,Aripo_rear_ud)
# chisq.test(table(Aripo_pop_ud,Aripo_rear_ud),simulate.p.value = TRUE)#Recap : using all genes (chisq = 491.96 pval = 0.00049)
#Quare
Quare_pop_ud = ifelse(Quare_pop_DE & unlist(Quare_DE_res$pop_estim >0), "DEU", 
                      ifelse(Quare_pop_DE & unlist(Quare_DE_res$pop_estim <0), "DED",
                             ifelse(unlist(Quare_DE_res$pop_estim >0), "NDEU","NDED")))
Quare_rear_ud = ifelse(Quare_rear_DE & unlist(Quare_DE_res$rear_estim >0), "DEU", 
                       ifelse(Quare_rear_DE & unlist(Quare_DE_res$rear_estim <0), "DED",
                              ifelse(unlist(Quare_DE_res$rear_estim >0), "NDEU","NDED")))

#Excluding interaction genes
table(Quare_pop_ud[!Quare_DE_res$inter_DE_flag],Quare_rear_ud[!Quare_DE_res$inter_DE_flag])
chisq.test(table(Quare_pop_ud[!Quare_DE_res$inter_DE_flag],Quare_rear_ud[!Quare_DE_res$inter_DE_flag]),
           simulate.p.value = TRUE)

#########################
# Pop/rear in two data set
#If we include all genes
# merged_DE_res = merge(Aripo_DE_res,Quare_DE_res,by = "gene",suffixes = c("_Aripo","_Quare"))

#If we exclude significant interaction genes, 19218
merged_DE_res = merge(Aripo_DE_res[!Aripo_DE_res$inter_DE_flag,],
                      Quare_DE_res[!Quare_DE_res$inter_DE_flag,],by = "gene",
                      suffixes = c("_Aripo","_Quare"))
##Re-obtain direction and DE/NDE information for 19218 genes
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

chisq.test(table(Aripo_pop_ud2,Aripo_rear_ud2),simulate.p.value = TRUE)
chisq.test(table(Quare_pop_ud2,Quare_rear_ud2),simulate.p.value = TRUE)

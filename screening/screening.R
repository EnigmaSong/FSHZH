# library(bigmemory)

##S_b : Bootstrap samples
##S_hat : testing statistics
##B : # bootstrap sample
##alpha : type-I error
#######################
## Under H0, the test statistic can be drawed from bootstrap distribution.
## That is, # total false rejection ratio from bootstrap with type-I error alpha = 0.05 should be around 
## G*theta_hat = G*(1/(G*B)*sum(S_b >= t_check,na.rm=TRUE)) = sum(ratio of rejection among B bootstrap samples for each tests) = approximated # of False rejections.
## # rejection is sum(S_hat >= t_check,na.rm=TRUE).

get_t_hat<-function(t = c(0,30), S_b, S_hat, B, alpha=0.05, na.rm=TRUE,track=FALSE){
  #Bisection method
  # t = seq(t_min,t_max,t_int)
  G = length(S_hat)
  
  #searching infimum
  while(t[2]-t[1]> 10^(-6)){
    t_check = mean(t)
    # #FP/ #P
    #one sided
    theta_hat = 1/(G*B) * sum(S_b >= t_check,na.rm=TRUE)
    
    ##############################
    ##To debug code
    if(track){
      message(cat("t:", t_check))
      message(cat("G*theta_hat:", G*theta_hat))
      message(cat("Numerator: ",max(sum(S_hat >= t_check,na.rm=TRUE),1)))
      message(cat("FDP: ",theta_hat*G/max(sum(S_hat >= t_check,na.rm=TRUE),1)))
    }
    ##############################
    
    cond = (theta_hat*G/max(sum(S_hat >= t_check,na.rm=TRUE),1) <= alpha)
    
    #If t_check satisfies the condition above, we try smaller one. 
    if(cond){
      t[2]=t_check
    }
    #If t_check does not satisfy the condition above, we try larger one. 
    else{
      t[1]=t_check
    }
  }
  return(t_check)
}


norm_Ctm_counts <- norm_counts_Aripo[ctm_gene_name,] ## Targets
norm_counts_Aripo <- norm_counts_Aripo[ctm_candidate,] ##Candidates

#Get Kendall's tau between targets and candidiates for each group
ctm_gene_cor_LPP = cor(t(norm_Ctm_counts[,LPP_Aripo]),
                       t(norm_counts_Aripo[,LPP_Aripo]), method = "kendall")
ctm_gene_cor_LPNP = cor(t(norm_Ctm_counts[,LPNP_Aripo]),
                        t(norm_counts_Aripo[,LPNP_Aripo]), method = "kendall")
ctm_gene_cor_HPP = cor( t(norm_Ctm_counts[,HPP_Aripo]),
                        t(norm_counts_Aripo[,HPP_Aripo]), method = "kendall")
ctm_gene_cor_HPNP = cor(t(norm_Ctm_counts[,HPNP_Aripo]),
                        t(norm_counts_Aripo[,HPNP_Aripo]), method = "kendall")
### Get Testing statistics
# For biological reason, we only consider LPNP,HPP, and HPNP.
#For one-sided
S_hat = colSums(ctm_gene_cor_LPNP,na.rm=TRUE) + 
  colSums(ctm_gene_cor_HPP,na.rm=TRUE) + 
  colSums(ctm_gene_cor_HPNP,na.rm=TRUE)
S_hat <- S_hat/sqrt(12)
head(S_hat)

#####################
### Bootstrap
B=2000 #Number of bootstrap samples
set.seed(2017)
p = dim(norm_counts_Aripo)[1]
S_b = matrix(0,B,p)
t_norm_Ctm_counts = t(norm_Ctm_counts)
t_norm_counts = t(norm_counts_Aripo)
for(i in 1:B){
  if(i%%50==0) print(paste(round(i/B*100,2),"%")) #print progress
  #To check independence between groups in heterogeneity world.
  #Re-Sampling in targets
  b_group12_ind_target = sample(LPNP_Aripo,replace=TRUE)
  b_group21_ind_target = sample(HPP_Aripo,replace=TRUE)
  b_group22_ind_target = sample(HPNP_Aripo,replace=TRUE)
  #Re-Sampling in candidates
  b_group12_ind = sample(LPNP_Aripo,replace=TRUE)
  b_group21_ind = sample(HPP_Aripo,replace=TRUE)
  b_group22_ind = sample(HPNP_Aripo,replace=TRUE)
  
  #Get cross-correlations to get Bootstrap samples
  B_cor_LPNP = cor(t_norm_Ctm_counts[b_group12_ind_target,], 
                   t_norm_counts[b_group12_ind,], method = "kendall")
  B_cor_HPP = cor( t_norm_Ctm_counts[b_group21_ind_target,], 
                   t_norm_counts[b_group21_ind,], method = "kendall")
  B_cor_HPNP = cor(t_norm_Ctm_counts[b_group22_ind_target,],
                   t_norm_counts[b_group22_ind,], method = "kendall")
  #Get bootstrap samples
  S_b[i,] = colSums(B_cor_LPNP,na.rm=TRUE)+
               colSums(B_cor_HPP,na.rm=TRUE)+
               colSums(B_cor_HPNP,na.rm=TRUE)
}
#Studentize.
S_b <- S_b/sqrt(12)

########
# Apply FDR
# Reference: Cai et. al (2015)
#########
t_hat2=get_t_hat(t=c(0,6),as.matrix(S_b),S_hat,B=2000,alpha = 0.2)
cont2 <- (S_hat>t_hat2)
cont2[is.na(cont2)] = FALSE
sum(cont2) #Check the number of screened genes

##
# Construct normalized count matrices for DE analysis
norm_counts_Aripo=norm_counts_Aripo[cont_tc_ind,] #Aripo
norm_counts_Quare=norm_counts_Quare[cont_tc_ind,]#Quare

############################
# DE Analysis
############################

############################
# Load packages
library(DESeq2)
library(dplyr)
library(lme4)
library(MASS)
library(parallel)
library(doSNOW)
library(foreach)
library(optimx)
library(lsmeans)

############################
# Load user-defined functions
source("Data preprocessing/normalize.R")

#############################
# Step 0: data pre-processing 
source("Data preprocessing/pre_processing.R")

###################################################
## Step 1
## Filtering trivially non-contaminated genes
## No variation or low counts

### Aripo
dd_Aripo <- DESeqDataSetFromMatrix(countData = counts_Aripo,colData = behave_Aripo, design= ~ pop + rear)
dd_Aripo <- dd_Aripo[ rowSums(counts(dd_Aripo)) > 50, ]

ctm_candidate <- rownames(counts(dd_Aripo))
ctm_gene_name = c("TRINITY_DN97289_c0_g1","TRINITY_DN230859_c0_g1",
                  "TRINITY_DN19667_c0_g2","TRINITY_DN99305_c0_g1")
ctm_gene_idx = match(ctm_gene_name,rownames(counts))
### We use 4 transcriptom as seeding genes.
Ctm_counts = counts_Aripo[ctm_gene_name,]
rownames(Ctm_counts)<-1:4

### Quare
dd_Quare <- DESeqDataSetFromMatrix(countData = counts_Quare,colData = behave_Quare, design= ~ pop + rear)
dd_Quare <- dd_Quare[ rowSums(counts(dd_Quare)) > 75, ]
ctm_candidate_Quare <- rownames(counts(dd_Quare))
length(ctm_candidate_Quare)

###################################################
## Step 2
## Filtering contaminated genes with Kendall's tau
###################################################
# Reference for normalization methods:
# Dillies, Rau, Aubert et al. (2013)
###################################################
# Use Aripo data to determine filtered genes
normalizing_gene_list <- c(ctm_candidate, ctm_gene_name)

#Use DESeq normalization
norm_method="DESeq"
#Aripo
norm_counts_Aripo <- normalize(counts_Aripo[normalizing_gene_list,],method=norm_method)
rownames(norm_counts_Aripo)<- normalizing_gene_list
colnames(norm_counts_Aripo)<- colnames(counts_Aripo)

#Quare
normalizing_contig_list_Quare <- c(ctm_candidate_Quare, ctm_gene_name)
#Use Q method for guppy_kallisto_genes.counts.matrix
norm_counts_Quare <- normalize(counts_Quare[normalizing_contig_list_Quare,],method=norm_method)
rownames(norm_counts_Quare)<- normalizing_contig_list_Quare
colnames(norm_counts_Quare)<- colnames(counts_Quare)
norm_counts_Quare <- norm_counts_Quare[ctm_candidate_Quare,] ##Candidates

#It takes more than 8 hours in Macbook 2013 Pro.
source("screening/screening.R")
#Output: 
## cont2: the screend gene indicator in Aripo 
# List of non-contaminated genes
cont_tc_ind <- !cont2 #rep(TRUE,nrow(norm_counts_Aripo)) #List of genes

#Construct normalized count matrices
#Aripo
norm_counts_Aripo <- norm_counts_Aripo[cont_tc_ind,] #non-contimated genes (Normalized)

#Quare

cont_tc_ind_Quare <- rep(TRUE,nrow(norm_counts_Quare))
names(cont_tc_ind_Quare)<- rownames(counts_Quare[ctm_candidate_Quare,])

#We decide to use the screened genes from old data
cont_tc_ind_Quare[intersect(names(cont_tc_ind[!cont_tc_ind]),names(cont_tc_ind_Quare))]<-FALSE

norm_counts_Quare <- norm_counts_Quare[cont_tc_ind_Quare,]

###################################
## Part 3: fitting the model
##############
n_cluster=2
##############
## interaction model, Aripo
message("Part1: Interaction model, Aripo")
# Get non-converged index from previous simulation
p = dim(norm_counts_Aripo)[1]
#Start sim
begin_t=Sys.time()
cl<- makeCluster(n_cluster,type="SOCK")
registerDoSNOW(cl)
#LRT results + DE results
Aripo_DE_res = foreach(i=1:p, .combine = rbind,.packages = c("lme4","MASS","lsmeans"),
                .errorhandling = "remove",.verbose=TRUE) %dopar%{
  #If fail to converge by Negative binomial, use poisson.
  #If fail to converge by Negative binomial, change tolPwrss. use poisson.
  warn_flag= FALSE
  #If any convergence warning, set warn_flag = TRUE
  #If any errors occur, use different optimization method for the second stage optimization.
  #Fit by default algorithm, if fails, set error_flag=TRUE
  err_flag=FALSE
  tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(norm_counts_Aripo[i,]~pop*rear + (1|week)+(1|family),
                                                    data = behave_Aripo)},
                               warning = function(w){
                                 warn_flag<<- TRUE
                                 invokeRestart("muffleWarning")
                               })
           ,error = function(e) err_flag<<-TRUE,fianlly=(if(!err_flag)alg_optim <<- c("bobyqa", "Nelder_Mead")))
  #Try lme4 optimizers
  optim_alg = c("bobyqa","Nelder_Mead","nloptwrap","nlminbwrap")
  j=1
  while(err_flag & (j<=length(optim_alg))){
    err_flag = FALSE
    tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(norm_counts_Aripo[i,]~pop*rear + (1|week)+(1|family),
                                                      data = behave_Aripo, control=glmerControl(optimizer=optim_alg[j]))},
                                 warning = function(w){
                                   warn_flag<<- TRUE
                                   invokeRestart("muffleWarning")
                                 })
             ,error = function(e) err_flag<<-TRUE, finally = {if(!err_flag) alg_optim<<-optim_alg[j]})
    j=j+1
  }
  #If no error
  if(err_flag) optimx_flag=TRUE
  else optimx_flag=FALSE
  #Try optimx optimizers if still errors
  j=1
  optimx_alg = c('Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm',
                 'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa',
                 'nmkb', 'hjkb', 'Rcgmin', 'Rvmmin')
  while(err_flag & (j<=length(optimx_alg))){
    err_flag = FALSE
    tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(norm_counts_Aripo[i,]~pop*rear + (1|week)+(1|family),
                                                      data = behave_Aripo,
                                                      control=glmerControl(optimizer="optimx",
                                                                           optCtrl=list(method=optimx_alg[j])))},
                                 warning = function(w){
                                   warn_flag<<- TRUE
                                   invokeRestart("muffleWarning")
                                 })
             ,error = function(e) err_flag<<-TRUE, finally = {alg_optim<<-optimx_alg[j]})
    j=j+1
  }
  #Update to remove warnings: up to 10 times
  #If error occures in updates, ignore it.
  j=0
  while(warn_flag & j<10){
    warn_flag = FALSE
    if(!optimx_flag){
      tryCatch(withCallingHandlers(fit_glmer<-update(fit_glmer,start=getME(fit_glmer,c("theta","fixef")),
                                                     control=glmerControl(optimizer=alg_optim)),
                                   warning=function(w){warn_flag<<-TRUE;invokeRestart("muffleWarning")})
               ,error = function(e) warn_flag=FALSE)
    }else{
      tryCatch(withCallingHandlers(fit_glmer<-update(fit_glmer,start=getME(fit_glmer,c("theta","fixef")),
                                                     control=glmerControl(optimizer="optimx",
                                                                          optCtrl=list(method=optimx_alg[j]))),
                                   warning=function(w){warn_flag<<-TRUE;invokeRestart("muffleWarning")})
               ,error = function(e) warn_flag=FALSE)
    }
    j=j+1
  }
  
  #Fit
  sum_glmer = summary(fit_glmer)
  #test contrast to check population and rearing effect
  # main_eff=test(lsmeans(fit_glmer,~group,contr=list(pop=c(0.5,0.5,-0.5,-0.5),rear=c(0.5,-0.5,0.5,-0.5),
  #                                                   interaction=c(1,-1,-1,1)))$contrasts)
  main_eff=test(lsmeans(fit_glmer,~pop*rear,contr=list(pop=c(0.5,-0.5,0.5,-0.5),rear=c(-0.5,-0.5,0.5,0.5),
                                                       inter=c(1,-1,-1,1),
                                                       simple_pop_fix_rearNP= c(1,-1,0,0),
                                                       simple_pop_fix_rearP= c(0,0,1,-1),
                                                       simple_rear_fix_popHP= c(-1,0,1,0),
                                                       simple_rear_fix_popLP= c(0,-1,0,1)
  ))$contrasts)
  
  return(c(i, id = rownames(norm_counts_Aripo)[i],
           coef_glmer = getME(fit_glmer,"beta"), theta_glmer = getME(fit_glmer,"glmer.nb.theta"),
           var_re = getME(fit_glmer,"theta")^2,
           estim_pop = main_eff[1,2], se_pop = main_eff[1,3], stat_pop_wald=main_eff[1,5], pval_pop_wald = main_eff[1,6],
           estim_rear = main_eff[2,2], se_rear = main_eff[2,3], stat_rear_wald=main_eff[2,5], pval_rear_wald = main_eff[2,6],
           estim_int = main_eff[3,2], se_int = main_eff[3,3], stat_int_wald=main_eff[3,5], pval_int_wald = main_eff[3,6],
           estim_simple_pop_fix_rearNP = main_eff[4,2], se_simple_pop_fix_rearNP = main_eff[4,3],
           stat_simple_pop_fix_rearNP_wald=main_eff[4,5], pval_simple_pop_fix_rearNP_wald = main_eff[4,6],
           estim_simple_pop_fix_rearP = main_eff[5,2], se_simple_pop_fix_rearP= main_eff[5,3],
           stat_simple_pop_fix_rearP_wald=main_eff[5,5], pval_simple_pop_fix_rearP_wald = main_eff[5,6],
           estim_simple_rear_fix_popHP  = main_eff[6,2], se_simple_rear_fix_popHP = main_eff[6,3],
           stat_simple_rear_fix_popHP_wald=main_eff[6,5], pval_simple_rear_fix_popHP_wald = main_eff[6,6],
           estim_simple_rear_fix_popLP= main_eff[7,2], se_simple_rear_fix_popLP = main_eff[7,3],
           stat_simple_rear_fix_popLP_wald=main_eff[7,5], pval_simple_rear_fix_popLP_wald = main_eff[7,6],
           conv_opt = sum_glmer$optinfo$conv$opt,
           conv_lme4 = ifelse(is.null(sum_glmer$optinfo$conv$lme4$code),0,-1),
           warn_flag=warn_flag,
           resid = resid(fit_glmer)
  ))
}
stopCluster(cl)
end_t = Sys.time()
end_t-begin_t
write.csv(Aripo_DE_res,"output/Aripo_glmer_effects_u_int.csv")

message("End of Part1: Interaction model, Aripo")

# ## Interaction model, new data
message("Part2: Interaction model, Quare")
p = dim(norm_counts_Quare)[1]
begin_t=Sys.time()
cl<- makeCluster(n_cluster,type="SOCK")
registerDoSNOW(cl)
#New data DE result in additive
Quare_DE_res = foreach(i=1:p, .combine = rbind,.packages = c("lme4","MASS","lsmeans"),
                .errorhandling = "remove",.verbose=TRUE) %dopar%{
  #If fail to converge by Negative binomial, use poisson.
  #If fail to converge by Negative binomial, change tolPwrss. use poisson.
  warn_flag= FALSE
  #If any convergence warning, set warn_flag = TRUE
  #If any errors occur, use different optimization method for the second stage optimization.
  #Fit by default algorithm, if fails, set error_flag=TRUE
  err_flag=FALSE
  tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(norm_counts_Quare[i,]~pop*rear + (1|week)+(1|family_new),
                                                    data = behave_Quare)},
                               warning = function(w){
                                 warn_flag<<- TRUE
                                 invokeRestart("muffleWarning")
                               })
           ,error = function(e) err_flag<<-TRUE,
           fianlly=(if(!err_flag)alg_optim <<- c("bobyqa", "Nelder_Mead")))
  #Try lme4 optimizers
  optim_alg = c("bobyqa","Nelder_Mead","nloptwrap","nlminbwrap")
  j=1
  while(err_flag & (j<=length(optim_alg))){
    err_flag = FALSE
    tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(norm_counts_Quare[i,]~pop*rear + (1|week)+(1|family_new),
                                                      data = behave_Quare, control=glmerControl(optimizer=optim_alg[j]))},
                                 warning = function(w){
                                   warn_flag<<- TRUE
                                   invokeRestart("muffleWarning")
                                 })
             ,error = function(e) err_flag<<-TRUE, finally = {if(!err_flag) alg_optim<<-optim_alg[j]})
    j=j+1
  }
  #If no error 
  if(err_flag) optimx_flag=TRUE
  else optimx_flag=FALSE
  #Try optimx optimizers if still errors
  j=1
  optimx_alg = c('Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm', 
                 'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa', 
                 'nmkb', 'hjkb', 'Rcgmin', 'Rvmmin')
  while(err_flag & (j<=length(optimx_alg))){
    err_flag = FALSE
    tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(norm_counts_Quare[i,]~pop*rear + (1|week)+(1|family_new),
                                                      data = behave_Quare, 
                                                      control=glmerControl(optimizer="optimx",
                                                                           optCtrl=list(method=optimx_alg[j])))},
                                 warning = function(w){
                                   warn_flag<<- TRUE
                                   invokeRestart("muffleWarning")
                                 })
             ,error = function(e) err_flag<<-TRUE, finally = {alg_optim<<-optimx_alg[j]})
    j=j+1
  }
  #Update to remove warnings: up to 10 times
  #If error occures in updates, ignore it.
  j=0
  while(warn_flag & j<10){
    warn_flag = FALSE
    if(!optimx_flag){
      tryCatch(withCallingHandlers(fit_glmer<-update(fit_glmer,start=getME(fit_glmer,c("theta","fixef")),
                                                     control=glmerControl(optimizer=alg_optim)),
                                   warning=function(w){warn_flag<<-TRUE;invokeRestart("muffleWarning")})
               ,error = function(e) warn_flag=FALSE)
    }else{
      tryCatch(withCallingHandlers(fit_glmer<-update(fit_glmer,start=getME(fit_glmer,c("theta","fixef")),
                                                     control=glmerControl(optimizer="optimx",
                                                                          optCtrl=list(method=optimx_alg[j]))),
                                   warning=function(w){warn_flag<<-TRUE;invokeRestart("muffleWarning")})
               ,error = function(e) warn_flag=FALSE)
    }
    j=j+1
  }
  #Fit 
  sum_glmer = summary(fit_glmer)
  #test contrast to check population and rearing effect
  #pop: HP-LP
  #rear: P-NP
  
  main_eff=test(lsmeans(fit_glmer,~pop*rear,contr=list(pop=c(0.5,-0.5,0.5,-0.5),rear=c(-0.5,-0.5,0.5,0.5),
                                                       inter=c(1,-1,-1,1),
                                                       simple_pop_fix_rearNP= c(1,-1,0,0),
                                                       simple_pop_fix_rearP= c(0,0,1,-1),
                                                       simple_rear_fix_popHP= c(-1,0,1,0),
                                                       simple_rear_fix_popLP= c(0,-1,0,1)
  ))$contrasts)
  
  return(c(i, id = rownames(norm_counts_Quare)[i],
           coef_glmer = getME(fit_glmer,"beta"), theta_glmer = getME(fit_glmer,"glmer.nb.theta"),
           var_re = getME(fit_glmer,"theta")^2, 
           estim_pop = main_eff[1,2], se_pop = main_eff[1,3], stat_pop_wald=main_eff[1,5], pval_pop_wald = main_eff[1,6],
           estim_rear = main_eff[2,2], se_rear = main_eff[2,3], stat_rear_wald=main_eff[2,5], pval_rear_wald = main_eff[2,6],
           estim_int = main_eff[3,2], se_int = main_eff[3,3], stat_int_wald=main_eff[3,5], pval_int_wald = main_eff[3,6],
           estim_simple_pop_fix_rearNP = main_eff[4,2], se_simple_pop_fix_rearNP = main_eff[4,3], 
           stat_simple_pop_fix_rearNP_wald=main_eff[4,5], pval_simple_pop_fix_rearNP_wald = main_eff[4,6],
           estim_simple_pop_fix_rearP = main_eff[5,2], se_simple_pop_fix_rearP= main_eff[5,3],
           stat_simple_pop_fix_rearP_wald=main_eff[5,5], pval_simple_pop_fix_rearP_wald = main_eff[5,6],
           estim_simple_rear_fix_popHP  = main_eff[6,2], se_simple_rear_fix_popHP = main_eff[6,3],
           stat_simple_rear_fix_popHP_wald=main_eff[6,5], pval_simple_rear_fix_popHP_wald = main_eff[6,6],
           estim_simple_rear_fix_popLP= main_eff[7,2], se_simple_rear_fix_popLP = main_eff[7,3],
           stat_simple_rear_fix_popLP_wald=main_eff[7,5], pval_simple_rear_fix_popLP_wald = main_eff[7,6],
           conv_opt = sum_glmer$optinfo$conv$opt, 
           conv_lme4 = ifelse(is.null(sum_glmer$optinfo$conv$lme4$code),0,-1),
           warn_flag=warn_flag,
           resid = resid(fit_glmer)
  ))
}
stopCluster(cl)
end_t = Sys.time()
end_t-begin_t

write.csv(Quare_DE_res,"output/Quare_glmer_w_new_family_labels.csv")
message("End of Part2: Interaction model, Quare")



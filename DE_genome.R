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
# Load user-defined functions, change path accordingly
setwd("DE_github")
source("Data preprocessing/normalize.R")
source("Data preprocessing/miscellaneous.R")
#############################
# Step 0: data pre-processing 
# To change data file path, edit pre_processing_genome.R.
source("Data preprocessing/pre_processing_genome.R")

###
#Draw MDS plot
#Aripo
MDS_Aripo = log2(1+counts_genome_Aripo) %>% t() %>% dist() %>% cmdscale(k=2)
plot(MDS_Aripo,xlab="Dim 1", ylab = "Dim 2",type="n")
text(MDS_Aripo,labels=behave_Aripo$fish,col = c("blue","red")[as.numeric(behave_Aripo$pop)])
legend("topright",fill=c("blue","red"),legend=levels(behave_Aripo$pop))

#Quare
MDS_Quare = log2(1+counts_genome_Quare) %>% t() %>% dist() %>% cmdscale(k=2)
plot(MDS_Quare,xlab="Dim 1", ylab = "Dim 2",type="n")
text(MDS_Quare,labels=behave_Quare$fish,col = c("blue","red")[as.numeric(behave_Quare$pop)])
legend("topleft",fill=c("blue","red"),legend=levels(behave_Quare$pop))

###################################################
## Step 1
## Filtering trivially non-contaminated genes
## No variation or low counts,

### Aripo
dd_genome_Aripo <- DESeqDataSetFromMatrix(countData = counts_genome_Aripo,
                                   colData = behave_Aripo, design= ~ pop + rear)
dd_genome_Aripo <- dd_genome_Aripo[ rowSums(counts(dd_genome_Aripo)) > 50, ]

ctm_candidate_genome <- rownames(counts(dd_genome_Aripo))
#Load seed-gene names. Change the file path accordingly
ctm_gene_name_genome <- read.csv('../../data/seedGenes.csv', stringsAsFactor =FALSE)[,1]
# ctm_gene_genome_idx = match(ctm_gene_name_genome,ctm_candidate_genome)
### We use 7 genomes as seeding genes.
Ctm_counts_genome = counts_genome_Aripo[ctm_gene_name_genome,]

### Quare
dd_genome_Quare <- DESeqDataSetFromMatrix(countData = counts_genome_Quare,colData = behave_Quare, design= ~ pop + rear)
dd_genome_Quare <- dd_genome_Quare[ rowSums(counts(dd_genome_Quare)) > 75, ]

###################################################
## Step 2
## Filtering contaminated genes with Kendall's tau
###################################################
# Reference for normalization methods:
# Dillies, Rau, Aubert et al. (2013)
###################################################
# Use Aripo data to determine filtered genes
normalizing_gene_list_genome <- c(ctm_candidate_genome, ctm_gene_name_genome)

#Use DESeq normalization
norm_method="DESeq"
#Aripo
norm_counts_Aripo_genome <- normalize(counts_genome_Aripo[normalizing_gene_list_genome,],
                                      method=norm_method)
rownames(norm_counts_Aripo_genome)<- normalizing_gene_list_genome
colnames(norm_counts_Aripo_genome)<- colnames(counts_genome_Aripo)

#It takes more than 8 hours in Macbook 2013 Pro.
norm_Ctm_counts_genome <- norm_counts_Aripo_genome[ctm_gene_name_genome,] ## Targets
norm_counts_Aripo_genome <- norm_counts_Aripo_genome[ctm_candidate_genome,] ##Candidates

#Get Kendall's tau between targets and candidiates for each group
# ctm_gene_cor_LPP_genome = cor(t(norm_Ctm_counts_genome[,LPP_Aripo]),
#                        t(norm_counts_Aripo_genome[,LPP_Aripo]), method = "kendall")
ctm_gene_cor_LPNP_genome = cor(t(norm_Ctm_counts_genome[,LPNP_Aripo]),
                        t(norm_counts_Aripo_genome[,LPNP_Aripo]), method = "kendall")
ctm_gene_cor_HPP_genome = cor( t(norm_Ctm_counts_genome[,HPP_Aripo]),
                        t(norm_counts_Aripo_genome[,HPP_Aripo]), method = "kendall")
ctm_gene_cor_HPNP_genome = cor(t(norm_Ctm_counts_genome[,HPNP_Aripo]),
                        t(norm_counts_Aripo_genome[,HPNP_Aripo]), method = "kendall")
### Get Testing statistics
# For biological reason, we only consider LPNP,HPP, and HPNP.
#For one-sided
S_hat_genome = colSums(ctm_gene_cor_LPNP_genome,na.rm=TRUE) + 
  colSums(ctm_gene_cor_HPP_genome,na.rm=TRUE) + 
  colSums(ctm_gene_cor_HPNP_genome,na.rm=TRUE)
S_hat_genome <- S_hat_genome/sqrt(3*7)
head(S_hat_genome)
hist(S_hat_genome)
#####################
### Bootstrap
B=2000 #Number of bootstrap samples
set.seed(2017)
p = dim(norm_counts_Aripo_genome)[1]
S_b_genome = matrix(0,B,p)
t_norm_Ctm_counts_genome = t(norm_Ctm_counts_genome)
t_norm_counts_genome = t(norm_counts_Aripo_genome)
for(i in 1:B){
  if(i%%20==0) print(paste(round(i/B*100,2),"%")) #print progress
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
  B_cor_LPNP = cor(t_norm_Ctm_counts_genome[b_group12_ind_target,], 
                   t_norm_counts_genome[b_group12_ind,], method = "kendall")
  B_cor_HPP = cor( t_norm_Ctm_counts_genome[b_group21_ind_target,], 
                   t_norm_counts_genome[b_group21_ind,], method = "kendall")
  B_cor_HPNP = cor(t_norm_Ctm_counts_genome[b_group22_ind_target,],
                   t_norm_counts_genome[b_group22_ind,], method = "kendall")
  #Get bootstrap samples
  S_b_genome[i,] = colSums(B_cor_LPNP,na.rm=TRUE)+
    colSums(B_cor_HPP,na.rm=TRUE)+
    colSums(B_cor_HPNP,na.rm=TRUE)
}
#Studentize.
S_b_genome <- S_b_genome/sqrt(3*7)

########
# Apply FDR
# Reference: Cai et. al (2015)
#########
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
t_hat2_genome=get_t_hat(t=c(0,6),as.matrix(S_b_genome),S_hat_genome,B=2000,alpha = 0.2)
cont2_genome <- (S_hat_genome>t_hat2_genome)
cont2_genome[is.na(cont2_genome)] = FALSE
sum(cont2_genome) #Check the number of screened genes
rm(S_b_genome)#Free memory
#Output: 
## cont2: the screend gene indicator in Aripo 
# List of non-contaminated genes
cont_tc_ind_genome <- !cont2_genome #rep(TRUE,nrow(norm_counts_Aripo_genome)) #List of genes

##
#Save screened genes and the average correlation
screened_gene_genome = names(which(cont2_genome))
gene_annote = read.csv("~/Dropbox/Project_Kim/data/GenomeAnnotations.csv",
                       stringsAsFactors = FALSE)
screened = subset(gene_annote, subset = (GenestableID %in% screened_gene_genome))
screened$avg.cor = S_hat_genome[screened$GenestableID]/sqrt(21)
screened = screened[order(screened$avg.cor,decreasing=TRUE),]
rownames(screened)<-NULL
write.csv(screened, file = "~/Dropbox/Project_Kim/outputs/screened_genes_genome.csv",
          row.names = FALSE)

#Construct normalized count matrices
#Aripo
norm_counts_Aripo_genome <- norm_counts_Aripo_genome[cont_tc_ind_genome,] #non-contimated genes (Normalized)

#Quare
#Remove the genes from Aripo filtering
cont_tc_ind_Quare <- rep(TRUE,nrow(counts(dd_genome_Quare)))
names(cont_tc_ind_Quare)<- rownames(counts(dd_genome_Quare))

#We decide to use the screened genes from Aripo data
cont_tc_ind_Quare[intersect(names(cont_tc_ind_genome[!cont_tc_ind_genome]),
                            names(cont_tc_ind_Quare))]<-FALSE
#Normalize, then wipe out the screened genes
norm_counts_Quare_genome <- counts(dd_genome_Quare) %>% normalize(method="DESeq") %>% 
                              .[cont_tc_ind_Quare,]

rownames(norm_counts_Quare_genome)<- rownames(counts(dd_genome_Quare)[cont_tc_ind_Quare,])
colnames(norm_counts_Quare_genome)<- colnames(counts_genome_Quare)

###################################
## Part 3: fitting the model
##
# DESeq2: For comparison
## Aripo
DESeq2_Aripo_genome = dd_genome_Aripo[!cont2_genome,] %>% DESeq()
res_DESeq2_Aripo_genome_pop = lfcShrink(DESeq2_Aripo_genome,contrast= c("pop","HP","LP"))
which(res_DESeq2_Aripo_genome_pop$padj<0.05)
res_DESeq2_Aripo_genome_rear = lfcShrink(DESeq2_Aripo_genome,contrast = c("rear","P","NP"))
which(res_DESeq2_Aripo_genome_rear$padj<0.05)
## Quare
DESeq2_Quare_genome = dd_genome_Quare[cont_tc_ind_Quare,] %>% DESeq()
res_DESeq2_Quare_genome_pop = lfcShrink(DESeq2_Quare_genome,contrast= c("pop","HP","LP"))
sum(res_DESeq2_Quare_genome_pop$padj<0.05)
res_DESeq2_Quare_genome_rear = lfcShrink(DESeq2_Quare_genome,contrast = c("rear","P","NP"))
sum(res_DESeq2_Quare_genome_rear$padj<0.05)

#w/ interaction

##############
n_cluster=2
##############
## interaction model, Aripo
message("Part1: Interaction model, Aripo")
# Get non-converged index from previous simulation
p = dim(norm_counts_Aripo_genome)[1]
#Start sim
begin_t=Sys.time()
cl<- makeCluster(n_cluster,type="SOCK")
registerDoSNOW(cl)
#LRT results + DE results
Aripo_DE_res_genome = foreach(i=1:p, .combine = rbind,.packages = c("lme4","MASS","lsmeans"),
                       .errorhandling = "remove",.verbose=TRUE) %dopar%{
                         #If fail to converge by Negative binomial, use poisson.
                         #If fail to converge by Negative binomial, change tolPwrss. use poisson.
                         warn_flag= FALSE
                         #If any convergence warning, set warn_flag = TRUE
                         #If any errors occur, use different optimization method for the second stage optimization.
                         #Fit by default algorithm, if fails, set error_flag=TRUE
                         err_flag=FALSE
                         tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(norm_counts_Aripo_genome[i,]~pop*rear + (1|week)+(1|family),
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
                           tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(norm_counts_Aripo_genome[i,]~pop*rear + (1|week)+(1|family),
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
                           tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(norm_counts_Aripo_genome[i,]~pop*rear + (1|week)+(1|family),
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
                         
                         ## Model comparison test between glmer.nb and glm.nb to check random effects
                         fit_glm <- glm.nb(norm_counts_Aripo_genome[i,]~pop*rear,
                                             data = behave_Aripo)
                         
                         model_comp = anova(fit_glmer, fit_glm)[2,6:8]
                         #test contrast to check population and rearing effect
                         # main_eff=test(lsmeans(fit_glmer,~group,contr=list(pop=c(0.5,0.5,-0.5,-0.5),rear=c(0.5,-0.5,0.5,-0.5),
                         #                                                   interaction=c(1,-1,-1,1)))$contrasts)
                         main_eff=test(lsmeans(fit_glmer,~pop*rear,contr=list(pop=c(-0.5,0.5,-0.5,0.5),rear=c(0.5,0.5,-0.5,-0.5),
                                                                              inter=c(1,-1,-1,1),
                                                                              simple_HPP_LPP= c(-1,1,0,0),
                                                                              simple_HPNP_LPNP= c(0,0,-1,1),
                                                                              simple_LPP_LPNP= c(1,0,-1,0),
                                                                              simple_HPP_HPNP= c(0,1,0,-1)
                                                                              
                         ))$contrasts)
                         
                         return(c(i, id = rownames(norm_counts_Aripo_genome)[i],
                                  coef_glmer = getME(fit_glmer,"beta"), theta_glmer = getME(fit_glmer,"glmer.nb.theta"),
                                  var_re = getME(fit_glmer,"theta")^2,
                                  estim_pop = main_eff[1,2], se_pop = main_eff[1,3], stat_pop_wald=main_eff[1,5], pval_pop_wald = main_eff[1,6],
                                  estim_rear = main_eff[2,2], se_rear = main_eff[2,3], stat_rear_wald=main_eff[2,5], pval_rear_wald = main_eff[2,6],
                                  estim_int = main_eff[3,2], se_int = main_eff[3,3], stat_int_wald=main_eff[3,5], pval_int_wald = main_eff[3,6],
                                  estim_simple_HPP_LPP = main_eff[4,2], se_simple_HPP_LPP = main_eff[4,3], 
                                  stat_simple_HPP_LPP_wald=main_eff[4,5], pval_simple_HPP_LPP_wald = main_eff[4,6],
                                  estim_simple_HPNP_LPNP = main_eff[5,2], se_simple_HPNP_LPNP= main_eff[5,3],
                                  stat_simple_HPNP_LPNP_wald=main_eff[5,5], pval_simple_HPNP_LPNP_wald = main_eff[5,6],
                                  estim_simple_LPP_LPNP = main_eff[6,2], se_simple_LPP_LPNP = main_eff[6,3],
                                  stat_simple_LPP_LPNP_wald=main_eff[6,5], pval_simple_LPP_LPNP_wald = main_eff[6,6],
                                  estim_simple_HPP_HPNP= main_eff[7,2], se_simple_HPP_HPNP = main_eff[7,3],
                                  stat_simple_HPP_HPNP_wald=main_eff[7,5], pval_simple_HPP_HPNP_wald = main_eff[7,6],
                                  model_comp_stat = model_comp[1],model_comp_df = model_comp[2],  
                                  model_comp_pval = model_comp[3],
                                  conv_opt = sum_glmer$optinfo$conv$opt,
                                  conv_lme4 = ifelse(is.null(sum_glmer$optinfo$conv$lme4$code),0,-1),
                                  warn_flag=warn_flag,
                                  resid = resid(fit_glmer)
                         ))
                       }
stopCluster(cl)
end_t = Sys.time()
end_t-begin_t
write.csv(Aripo_DE_res_genome,"output/Aripo_glmer_effects_genome.csv")

message("End of Part1: Interaction model, Aripo")

# ## Interaction model, new data
message("Part2: Interaction model, Quare")
p = dim(norm_counts_Quare_genome)[1]
begin_t=Sys.time()
n_cluster=2
cl<- makeCluster(n_cluster,type="SOCK")
registerDoSNOW(cl)
#New data DE result in additive
Quare_DE_res_genome1 = foreach(i=1:p, .combine = rbind,.packages = c("lme4","MASS","lsmeans"),
                       .errorhandling = "remove",.verbose=TRUE) %dopar%{
                         #If fail to converge by Negative binomial, use poisson.
                         #If fail to converge by Negative binomial, change tolPwrss. use poisson.
                         warn_flag= FALSE
                         #If any convergence warning, set warn_flag = TRUE
                         #If any errors occur, use different optimization method for the second stage optimization.
                         #Fit by default algorithm, if fails, set error_flag=TRUE
                         err_flag=FALSE
                         tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(norm_counts_Quare_genome[i,]~pop*rear + (1|week)+(1|family_new),
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
                           tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(norm_counts_Quare_genome[i,]~pop*rear + (1|week)+(1|family_new),
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
                           tryCatch(withCallingHandlers({fit_glmer<-glmer.nb(norm_counts_Quare_genome[i,]~pop*rear + (1|week)+(1|family_new),
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
                         
                         ## Model comparison test between glmer.nb and glm.nb to check random effects
                         fit_glm <- glm.nb(norm_counts_Quare_genome[i,]~pop*rear, 
                                           data = behave_Quare)
                        
                         model_comp = anova(fit_glmer, fit_glm)[2,6:8]
                         #test contrast to check population and rearing effect
                         #pop: HP-LP
                         #rear: P-NP
                         
                         main_eff=test(lsmeans(fit_glmer,~pop*rear,contr=list(pop=c(-0.5,0.5,-0.5,0.5),rear=c(0.5,0.5,-0.5,-0.5),
                                                                              inter=c(1,-1,-1,1),
                                                                              simple_HPP_LPP= c(-1,1,0,0),
                                                                              simple_HPNP_LPNP= c(0,0,-1,1),
                                                                              simple_LPP_LPNP= c(1,0,-1,0),
                                                                              simple_HPP_HPNP= c(0,1,0,-1)
                         ))$contrasts)
                         
                         return(c(i, id = rownames(norm_counts_Quare_genome)[i],
                                  coef_glmer = getME(fit_glmer,"beta"), theta_glmer = getME(fit_glmer,"glmer.nb.theta"),
                                  var_re = getME(fit_glmer,"theta")^2, 
                                  estim_pop = main_eff[1,2], se_pop = main_eff[1,3], stat_pop_wald=main_eff[1,5], pval_pop_wald = main_eff[1,6],
                                  estim_rear = main_eff[2,2], se_rear = main_eff[2,3], stat_rear_wald=main_eff[2,5], pval_rear_wald = main_eff[2,6],
                                  estim_int = main_eff[3,2], se_int = main_eff[3,3], stat_int_wald=main_eff[3,5], pval_int_wald = main_eff[3,6],
                                  estim_simple_HPP_LPP = main_eff[4,2], se_simple_HPP_LPP = main_eff[4,3], 
                                  stat_simple_HPP_LPP_wald=main_eff[4,5], pval_simple_HPP_LPP_wald = main_eff[4,6],
                                  estim_simple_HPNP_LPNP = main_eff[5,2], se_simple_HPNP_LPNP= main_eff[5,3],
                                  stat_simple_HPNP_LPNP_wald=main_eff[5,5], pval_simple_HPNP_LPNP_wald = main_eff[5,6],
                                  estim_simple_LPP_LPNP = main_eff[6,2], se_simple_LPP_LPNP = main_eff[6,3],
                                  stat_simple_LPP_LPNP_wald=main_eff[6,5], pval_simple_LPP_LPNP_wald = main_eff[6,6],
                                  estim_simple_HPP_HPNP= main_eff[7,2], se_simple_HPP_HPNP = main_eff[7,3],
                                  stat_simple_HPP_HPNP_wald=main_eff[7,5], pval_simple_HPP_HPNP_wald = main_eff[7,6],
                                  model_comp_stat = model_comp[1],model_comp_df = model_comp[2],  
                                  model_comp_pval = model_comp[3],
                                  conv_opt = sum_glmer$optinfo$conv$opt, 
                                  conv_lme4 = ifelse(is.null(sum_glmer$optinfo$conv$lme4$code),0,-1),
                                  warn_flag=warn_flag,
                                  resid = resid(fit_glmer)
                         ))
                       }
stopCluster(cl)
end_t = Sys.time()
end_t-begin_t

write.csv(Quare_DE_res_genome,"output/Quare_glmer_genome.csv")
message("End of Part2: Interaction model, Quare")

# Count DE genes
# Aripo_DE_res_genome
# Quare_DE_res_genome
# Aripo_DE_res_genome = read.csv("output/Aripo_glmer_effects_genome.csv")
# Quare_DE_res_genome = read.csv("output/Quare_glmer_genome.csv")



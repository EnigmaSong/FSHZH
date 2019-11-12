#Mean/Median comparison between DE and NDE (population effect)

"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y

# DE list by Benjamini & Hochberg (1995) with FDR level 0.05
Aripo_pop_DE = Aripo_DE_res$id[which(p.adjust(Aripo_DE_res$pval_pop_wald,method = "BH")<0.05)]
Quare_pop_DE = Quare_DE_res$id[which(p.adjust(Quare_DE_res$pval_pop_wald,method = "BH")<0.05)]

##Two-sided
#two-sample t
Aripo_mean_t = t.test(norm_counts_Aripo[Aripo_pop_DE,], 
                    norm_counts_Aripo[Aripo_DE_res$id %w/o% Aripo_pop_DE,])
Quare_mean_t = t.test(norm_counts_Quare[Quare_pop_DE,], 
                    norm_counts_Quare[Quare_DE_res$id %w/o% Quare_pop_DE,])
# Wilcoxon-rank sum
Aripo_ws = wilcox.test(norm_counts_Aripo[Aripo_pop_DE,], 
                     norm_counts_Aripo[Aripo_DE_res$id %w/o% Aripo_pop_DE,])
Quare_ws = wilcox.test(norm_counts_Quare[Quare_pop_DE,], 
                     norm_counts_Quare[Quare_DE_res$id %w/o% Quare_pop_DE,])

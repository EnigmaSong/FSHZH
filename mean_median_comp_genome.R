#two-sample t-test and wilcoxon test

# Aripo_DE_res_genome = read.csv("output/Aripo_glmer_effects_genome.csv")
# Quare_DE_res_genome = read.csv("output/Quare_glmer_genome.csv")


#Aripo pop
aripo_pop_DE_id = Aripo_DE_res$Id[Aripo_DE_res$pop_DE_flag==1]
aripo_pop_NDE_id = Aripo_DE_res$Id[Aripo_DE_res$pop_DE_flag==0]

#Mean (2-sample t)
t.test(norm_counts_Aripo_genome[aripo_pop_DE_id,],
       norm_counts_Aripo_genome[aripo_pop_NDE_id,])
t.test(norm_counts_Aripo_genome[aripo_pop_DE_id,],
       norm_counts_Aripo_genome[aripo_pop_NDE_id,],alternative = "less")

#Median (Wilcox)
wilcox.test(norm_counts_Aripo_genome[aripo_pop_DE_id,],
       norm_counts_Aripo_genome[aripo_pop_NDE_id,])
median(norm_counts_Aripo_genome[aripo_pop_DE_id,])#median
median(norm_counts_Aripo_genome[aripo_pop_NDE_id,])

#Quare pop
quare_pop_DE_id = Quare_DE_res$Id[Quare_DE_res$pop_DE_flag==1]
quare_pop_NDE_id = Quare_DE_res$Id[Quare_DE_res$pop_DE_flag==0]

#Mean (2-sample t)
t.test(norm_counts_Quare_genome[quare_pop_DE_id,],
       norm_counts_Quare_genome[quare_pop_NDE_id,])
#Median (Wilcox)
wilcox.test(norm_counts_Quare_genome[quare_pop_DE_id,],
            norm_counts_Quare_genome[quare_pop_NDE_id,])
median(norm_counts_Quare_genome[quare_pop_DE_id,])#median
median(norm_counts_Quare_genome[quare_pop_NDE_id,])

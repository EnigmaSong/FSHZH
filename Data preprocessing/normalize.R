normalize<-function(counts, method=c("TC","UQ","Med","DESeq","TMM","Q","RPKM")){
  method <- match(toupper(method), c("TC","UQ","MED","DESEQ","TMM","Q","RPKM"))
  if(is.na(method)){
    message("Use DESeq as the defalut")
    method = 4
  }
  
  ##TC
  if(method==1){
    # norm_counts[,i] = counts[normalizing_contig_list,i]/col_total[i]*mean(col_total)
    counts = apply(counts,2,function(counts){counts/sum(counts)})*mean(colSums(counts))
  }
  #UQ
  else if(method==2){
    meanUQ = mean(apply(counts,2, quantile, prob=0.75))
    counts = apply(counts,2,function(counts){counts/quantile(counts[counts>0],prob=0.75)})*meanUQ
    
  }
  #Med
  else if(method==3){
    meanMed = mean(apply(counts,2, median))
    counts = apply(counts,2,function(counts){counts/median(counts[counts>0])})*meanMed
  }
  #DESeq2
  else if(method==4){
    norm_factor=estimateSizeFactorsForMatrix(counts)
    counts <-mapply("/",as.data.frame(counts),norm_factor)
  }
  #TMM
  else if(method==5){
    norm_factor=calcNormFactors(counts)
    counts <-mapply("/",as.data.frame(counts),norm_factor)
  }
  #Q
  else if(method==6){
    counts <-as.matrix(normalizeQuantiles(counts))
  }
  #RPKM
  else if(method==7){
    counts <-rpkm(counts,gene.length=nrow(counts),lib.size=ncol(counts))
  }
  return(counts)
}

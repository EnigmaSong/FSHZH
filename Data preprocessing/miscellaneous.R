##############################
#colSum of positive elements
##############################
positive_colSums<-function(X1){
  X1<-X1*(X1>0)
  return(colSums(X1))
}


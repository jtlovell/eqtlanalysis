counts.summary<-function(counts){
  sums<-apply(counts,1,function(x){
    x<-x[!is.na(x)]
    n.na<-length(x[x==0])
    mean.out<-mean(x[x!=0])
    skew.out<-skewness(x[x!=0])
    skew.log<-skewness(log(x[x!=0]))
    return(cbind(n.na,mean.out,skew.out,skew.log))
  }
  )
  sums<-data.frame(t(sums))
  colnames(sums)<-c("n.na","mean","skewness","skew.log")
  sums<-sums[complete.cases(sums),]
  #have a look at the distribution of expression for each gene
  par(mfrow=c(2,2))
  hist(sums$n.na)
  hist(log(sums$mean))
  hist(sums$skewness)
  hist(sums$skew.log)
  return(sums)
}


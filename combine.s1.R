combine.s1<-function(s1s,n=100){
  sets<-split(1:length(s1s), ceiling(seq_along(1:length(s1s))/n))
  cat("combining scanones by sets of", n, "\n")
  out.byset<-list()
  for(i in 1:length(sets)){
    out.byset[[i]]<-s1.test[[sets[[i]][1]]]
    for(j in (min(sets[[i]])+1): max(sets[[i]])){
      out.byset[[i]]<-cbind(out.byset[[i]], s1.test[[j]])
    }
    cat(length(out.byset)*n,"...")
  }
  out.all<-out.byset[[1]]
  cat("\n combining sets (total =",length(s1s),") : \n")
  for(i in 2:length(out.byset))  {
    out.all<-cbind(out.all, out.byset[[i]])
    cat(length(colnames(out.all)),"...")
  }
  return(out.all)
}
#combine 20 scanones
combine.s1<-function(s1s,n){
  s1.out<-s1s[[1]]
  if(n>length(s1s)){
    for(j in 2:length(s1s)){
      s1.out<-c(s1.out,s1s[[j]])
    }
  }else{
    for(j in 2:n){
      s1.out<-c(s1.out,s1s[[j]])
    }
  }
  return(s1.out)
}

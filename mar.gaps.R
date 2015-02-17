mar.gaps<-function(map){
  dout<-data.frame()
  for(i in 1:(length(map[,1])-1)){
    pos.1<-map[i,]
    pos.2<-map[i+1,]
    if(pos.2$pos<.001){
      d<-0
    }else{
      d<-pos.2$pos-pos.1$pos
    }
    d.out<-data.frame(pos.1$chr, pos.1$pos, pos.2$pos,d)
    colnames(d.out)<-c("chr","pos1","pos2","dist")
    dout<-rbind(dout,d.out)
  }
  return(dout)
}
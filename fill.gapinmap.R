fill.gapinmap<-function(cross, chr, pos.low, pos.high, genos, genepos){
  gp.chr<-genepos[genepos$lg==chr & genepos$pos>pos.low & genepos$pos < pos.high,]
  if(length(gp.chr[,1])<1){
    return(cross)
  }else{
    wind<-seq(from=min(gp.chr$pos), to=max(gp.chr$pos), by=.2)
    if(length(wind)==1){
      wind<-c(min(gp.chr$pos),max(gp.chr$pos))
    }
    for(i in 1:(length(wind)-1)){
      gp<-gp.chr[with(gp.chr, pos>wind[i] & pos <wind[i+1]),]
      if(length(gp[,1])<1){
        next()
      }else{
        gp.mar<-gp[gp$lod==max(gp$lod),]
        if(length(gp.mar$lod)>1){
          gp.mar<-gp.mar[1,]
        }
      }
      mar<-as.character(gp.mar$gene)
      als<-as.character(genos[,mar])
      als[als=="fil"]<-1
      als[als=="hal"]<-3
      als[als=="het"]<-2
      als<-as.numeric(as.character(als))
      cross<-addmarker(cross, genotypes=als, markername=mar, chr=gp.mar$lg, pos=gp.mar$pos)
    }
    return(cross)
  }
}

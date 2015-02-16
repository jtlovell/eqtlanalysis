circos.mapcompare<-function(temp, toplot, tolines){
  colnames(temp)<-c("gene","chr","start")
  tp<-rbind(toplot,temp)
  tp$chr<-as.character(tp$chr)
  tp$start<-as.numeric(tp$start)
  tp<-tp[with(tp, order(chr)),]
  tp<-tp[-grep("super",tp$chr),]
  
  xlim = cbind(rep(0,length(unique(tp$chr))),
               aggregate(start~chr,data=tp,max)$start)
  chromosome = unique(tp[,"chr"])
  
  par(mar = c(1, 1, 1, 1), lwd = 0.5)
  circos.par(cell.padding = c(0, 0, 0, 0))
  circos.initialize(factors = factor(chromosome, levels = chromosome), xlim = xlim)
  circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.1)
  col=c(rainbow(length(chromosome[grep("Chr",chromosome)]), alpha=.2),
        rainbow(length(chromosome[grep("lg",chromosome)]), alpha=1))
  for(i in 1:length(unique(tp$chr))) {
    # data in current `chr`
    chr=as.character(unique(tp$chr)[order(unique(tp$chr))][i])
    d2 = tp[tp$chr == chr, ]
    n = nrow(d2)
    circos.rect(min(d2[,"start"]), 0, max(d2[,"start"]), 0.4, sector.index = chr, border = "black", col=col[i])
    circos.text((max(d2$start)-min(d2$start))/2, 1, labels = gsub("chr", "", chr),
                sector.index = chr, cex = 0.8)
  }
  
  for(i in 1:9) {
    # data in current `chr`
    chr=paste("lg_0",c(1:9), sep="")[i]
    d2 = tolines[tolines$lg == chr, ]
    d2<-d2[,]
    for(j in d2$gene){
      cm.pos<-d2[d2$gene==j,"pos.cor"]
      cm.lg<-d2[d2$gene==j,"lg"]
      bp.pos<-d2[d2$gene==j,"start"]
      bp.chr<-d2[d2$gene==j,"chr"]
      circos.link(cm.lg, cm.pos, bp.chr, bp.pos, col=col[i+9])
    }
  }
  
  circos.clear()
}

change.lgnames<-function(cross){
  pm<-pull.map(cross)
  gene.lg<-data.frame(names(unlist(pm)))
  gene.lg$gene.id<-sapply(names(unlist(pm)), function(x) strsplit(x, "[.]")[[1]][2])
  gene.lg$ps<-sapply(as.character(gene.lg[,1]), function(x) strsplit(x, "[.]")[[1]][3])
  gene.lg$lg<-sapply(as.character(gene.lg[,1]), function(x) strsplit(x, "[.]")[[1]][1])
  gene.lg$gene<-paste(gene.lg$gene.id,gene.lg$ps, sep=".")
  
  lg.info<-merge(annot,gene.lg,by="gene")
  for(i in c("Chr_01","Chr_02","Chr_03","Chr_04","Chr_05","Chr_06","Chr_07","Chr_08","Chr_09")){
    print(i)
    print(table(lg.info$lg[lg.info$chr==i]))
  }
  
  lg.call<-data.frame()
  out.all<-data.frame()
  for (i in unique(lg.info$lg)){
    tab<-table(lg.info$chr[lg.info$lg==i])
    out<-data.frame(i,names(tab),as.numeric(tab))
    colnames(out)<-c("lg.new","chr.annot","count")
    out.all<-rbind(out.all,out)
    out2<-data.frame(i,as.character(out$chr.annot[out$count==max(out$count)]))
    colnames(out2)<-c("lg.new","chr.annot")
    lg.call<-rbind(lg.call,out2)
  }
  #make lg names reflext the chrs that are most represented in the lg
  nam <- names(cross$geno)
  lg.call$chr.annot.num<-gsub("chr0", "", lg.call$chr.annot)
  nam[1]<-lg.call$chr.annot.num[lg.call$lg.new==1]
  nam[2]<-lg.call$chr.annot.num[lg.call$lg.new==2]
  nam[3]<-lg.call$chr.annot.num[lg.call$lg.new==3]
  nam[4]<-lg.call$chr.annot.num[lg.call$lg.new==4]
  nam[5]<-lg.call$chr.annot.num[lg.call$lg.new==5]
  nam[6]<-lg.call$chr.annot.num[lg.call$lg.new==6]
  nam[7]<-lg.call$chr.annot.num[lg.call$lg.new==7]
  nam[8]<-lg.call$chr.annot.num[lg.call$lg.new==8]
  nam[9]<-lg.call$chr.annot.num[lg.call$lg.new==9]
  names(cross$geno) <- nam
  return(cross)
}


Rqtl2_ped.genABEL<-function(cross, phys, 
                            ped.filename="out.ped", 
                            map.filename="out.map",
                            ga.phefilename="phe.txt",
                            return.ga=T){
  g<-pull.geno(cross)
  m<-data.frame(pull.map(cross, as.table=T))
  m$gene<-rownames(m)
  phys<-phys[,c("gene","start")]
  map<-merge(phys,m, by="gene")
  map<-map[,c("chr","gene","pos","start")]
  write.table(map, file=map.filename, sep="\t", row.names=F, col.names=F)
  
  i<-data.frame(getid(cross)); colnames(i)<-"sample_id"
  ped<-data.frame(cbind("FAM1",i$sample_id))
  colnames(ped)<-c("family_id","sample_id")
  ped$paternal_id<-0
  ped$maternal_id<-0
  ped$sex<-0
  ped$affection=0
  
  g[is.na(g)]<-"0/0"
  g[g==1]<-"A/A"
  g[g==2]<-"A/T"
  g[g==3]<-"T/T"
  ped<-cbind(ped,g)
  write.table(ped, file=ped.filename, sep="\t", row.names=F, col.names=F, quote=F)
  
  phe<-pull.pheno(cross)
  phe<-phe[,c("id","Treatment",m$gene)]
  phe$sex<-0
  write.table(phe, file=ga.phefilename, sep="\t", quote=F, row.names=F)
  if(return.ga){
    convert.snp.ped(ped=ped.filename, mapfile=map.filename, outfile="genos.raw", wslash=T)
    ga<-load.gwaa.data(phenofile = "phe.txt", genofile = "genos.raw", 
                       force = TRUE, makemap = FALSE, sort = TRUE, id = "id")
    return(ga)
  }
}


#function to call alleles from ASE RNA-seq data
#requires the name of the gene and two files with the ASE data with colnames as genes
genobyMclust<-function(id, gene.name, a1se, a2se,  
                       p1, p2,
                       a1.name="fil", a2.name="hal",
                       n.clust=3,p.tol=1e-10, plot=FALSE, logT=T, uncertaintyOnly=T){
  geno.in<-data.frame(a1se, a2se,id, stringsAsFactors = F)
  parents<-data.frame(rbind(cbind(p1,0),cbind(0,p2)))
  parents$id<-c(rep("p1", length(p1)), rep("p2", length(p2)))
  colnames(parents)<-colnames(geno.in)
  dat<-rbind(parents,geno.in)
  dat$id[grep("F1",dat$id)]<-"F1"
  
  if(logT){
    dat$a1se<-log10(dat$a1se+1)
    dat$a2se<-log10(dat$a2se+1)
  }
  clust1<-Mclust(dat[,"a1se"],G=2, prior = priorControl())
  clust2<-Mclust(dat[,"a2se"],G=2, prior = priorControl())
  dat$allele1<-clust1$classification
  dat$allele2<-clust2$classification
  
  dat$uncertainty1<-clust1$uncertainty
  dat$uncertainty2<-clust2$uncertainty
  
  p1.group1<-unique(dat$allele1[dat$id=="p1"])
  p1.group2<-unique(dat$allele2[dat$id=="p1"])
  
  p2.group1<-unique(dat$allele1[dat$id=="p2"])
  p2.group2<-unique(dat$allele2[dat$id=="p2"])
  
  if(length(p1.group1)>1 | length(p1.group2)>1 | length(p2.group1)>1 | length(p2.group2)>1){
    genotype<-rep(NA, nrow(dat))
  }else{
    dat$allele1[dat$uncertainty1>p.tol]<-0
    dat$allele2[dat$uncertainty2>p.tol]<-0
    genotype<-with(dat, ifelse(allele1==p1.group1 & allele2==p1.group2, a1.name,
                               ifelse(allele1==p2.group1 & allele2==p2.group2, a2.name,
                                      ifelse(allele1==p2.group1 & allele2==p1.group2 |
                                               allele1==p1.group1 & allele2==p2.group2, 
                                             paste(a1.name,a2.name,sep=""), NA))))
    genotype[dat$a1se==0 & dat$a2se>0]<-a2.name
    genotype[dat$a2se==0 & dat$a1se>0]<-a1.name
    dat$group2<-genotype
    if(plot){
      
      print(ggplot(dat, aes(x=a2se, y=a1se,col=as.factor(group2)))+
              geom_point()+
              scale_y_continuous(name="FIL-allele-specific expression (log10)")+
              scale_x_continuous(name="HAL-allele-specific expression (log10)")+
              theme_bw()+
              ggtitle("expression categorization by  genotype"))
    }
  }
  return(dat$group2)
}
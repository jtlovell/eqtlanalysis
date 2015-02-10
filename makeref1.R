#function to make a qtl object at a specified location, but allow it to "wiggle" n cM away from that position
make.1<-function(cross, chr, pos, phe, cov, wiggle=10){
  m<-makeqtl(cross=cross, chr=chr, pos=pos, what="prob")
  r<-refineqtl(cross=cross, qtl=m, pheno.col=phe, 
               covar=cov, method="hk", model="normal", verbose=F)
  if(abs(m$pos-r$pos)>wiggle){
    return(m)
  }else{return(r)}
}


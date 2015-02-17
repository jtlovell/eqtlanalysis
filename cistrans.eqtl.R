#Author JT Lovell
#Date 9-Feb 2014
#Version 1.2

#This function fits an exhaustive set of cis/trans QTL models
#Here, the cis QTL term is not penalized, but all other terms are
#Penalties must be calculated prior and must be a vector of 9 (numeric)
#chromosome and postion are the locations of the cis eQTL to fit without a penalty
#cov is a dataframe with a covariate (only one possible) that can be fit without penalty

# #this returns a list of 4 objects
#   the best model (integer 1:9)
#   the LOD score of the non-covariate terms in the best model
#   the best QTL object
#   the best fitted QTL model with statistics
cistrans.eqtl<-function(cross, chromosome, position, phe, pens=NULL, forms.in=NULL,cov=covar){
  
  #for a standard cis-trans eqtl, we search exhaustively through mode space
  #   the cis eQTL ("Q1") and treatment covariate ("trt") are always included. 
  if(is.null(forms.in)){
    forms<-c("y ~ Q1 + trt",
             "y ~ Q1 + Q1*trt + trt",
             "y ~ Q1 + Q2 + trt",
             "y ~ Q1 + Q2 + Q1*trt + trt",
             "y ~ Q1 + Q2 + Q2*trt + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + Q1*trt + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + Q2*trt + trt",
             "y ~ Q1 + Q2 + Q1*Q2 + Q1*trt + Q2*trt + trt")
  }
  #get name for phe, so that we index by character not number
  if(is.numeric(phe)){
    phe<-phenames(cross)[phe]
  }
  mods<-list()
  fits<-list()
  # the first two models need to be built manually because of a single QTL
  cis.eqtl<-makeqtl(cross=cross, chr=chromosome, pos=position, what="prob")
  mods[[1]]<-cis.eqtl
  mods[[2]]<-cis.eqtl
  fit.1<-fitqtl(cross, qtl=cis.eqtl, 
                formula=forms[1], pheno.col=phe, 
                covar=cov, method="hk", dropone=T)
  
  fit.2<-fitqtl(cross, qtl=cis.eqtl, 
                formula=forms[2], pheno.col=phe, 
                covar=cov, method="hk", dropone=T)
  
  lod.1<-data.frame(fit.1$result.drop)
  fits[[1]]<-lod.1
  lod.1<-sum(lod.1[which(rownames(lod.1)!=colnames(cov)),"LOD"])
  lod.2<-data.frame(fit.2$result.drop)
  fits[[2]]<-lod.2
  lod.2<-sum(lod.2[which(rownames(lod.2)!=colnames(cov)),"LOD"])
  lod.all<-data.frame(lod.1,lod.2)
  # fit the remaining 7 models
  for (i in 3:length(forms)){
    form.in<-forms[i]
    scan <- addqtl(cross, qtl=cis.eqtl, formula=form.in, 
                   method="hk", covar=cov, pheno.col=phe)
    mod <- addtoqtl(cross, cis.eqtl, max(scan)$chr, max(scan)$pos)
    fit <- fitqtl(cross, qtl=mod, 
                  formula=form.in, pheno.col=phe, 
                  covar=cov, method="hk", dropone=T, get.ests=T)
    lod.out<-data.frame(fit$result.drop)
    lod.all<-cbind(lod.all,lod.out)
    mods[[i]]<-mod
  }
  #deterime which mode fit is best by taking the best pLOD score from the QTL and QTL*trt interactions
  best<-which(lod.all-pens == max(lod.all-pens))
  best.lod<-lod.all[best]
  best.mod<-mods[[best]]
  best.fit<-fits[[best]]
  
  #extract component and return table of stats
  pos.out<-sapply(rownames(lod.out), function(x) strsplit(x,"@")[[1]][2])
  pos.out[grep(":",pos.out)]<-NA
  lod.out$pos<-pos.out
  
  chr.out<-sapply(rownames(lod.out), function(x) strsplit(x,"@")[[1]][1])
  chr.out[is.na(pos.out)]<-NA
  lod.out$chr<-chr.out
  
  cat1<-rownames(lod.out)
  category<-rownames(lod.out)
  category[intersect(grep("trt", cat1), grep(":", cat1))]<-"trt.int"
  category[intersect(grep("trt", cat1), grep(":", cat1,invert=TRUE))]<-"trt"
  category[intersect(grep("trt", cat1, invert=TRUE), grep(":", cat1))]<-"epi"
  cis.index<-which(as.numeric(cis.eqtl$pos)==pos & as.numeric(cis.eqtl$chr)==chr)
  trans.index<-which(as.numeric(pos.out)==max(scan)$pos & as.numeric(chr.out)==max(scan)$chr)
  category[cis.index]<-"cis"
  category[trans.index]<-"trans"
  
  lod.out$category<-category
  
  lod.out$phenotype<-phe
  
  ests.all<-summary(fit)$ests
  rows.dom<-intersect(grep("d",rownames(ests.all)), grep(":",rownames(ests.all),invert=TRUE))
  rows.add<-intersect(intersect(grep("a",rownames(ests.all)), grep(":",rownames(ests.all),invert=TRUE)), 
                      grep(colnames(covar),rownames(ests.all),invert=TRUE))
  rows.cov<-which(rownames(ests.all) %in% colnames (covar))
  if(nqtls==1){
    covar.ests_out<-data.frame(t(c(ests.all[rows.cov,],NA,NA,NA)))
    colnames(covar.ests_out)<-c("est.add","SE.add","t.add","est.dom","SE.dom","t.dom")
    dom.ests_out<-data.frame(t(ests.all[rows.dom,1:3])); colnames(dom.ests_out)<-c("est.dom","SE.dom","t.dom")
    add.ests_out<-data.frame(t(ests.all[rows.add,1:3])); colnames(add.ests_out)<-c("est.add","SE.add","t.add")
    qtl.ests_out<-cbind(add.ests_out,dom.ests_out)
    ests.out<-rbind(covar.ests_out,qtl.ests_out)
  }else{
    covar.ests_out<-data.frame(t(c(ests.all[rows.cov,],NA,NA,NA)))
    colnames(covar.ests_out)<-c("est.add","SE.add","t.add","est.dom","SE.dom","t.dom")
    rownames(covar.ests_out)<-colnames(covar)
    dom.ests_out<-data.frame(ests.all[rows.dom,1:3]); colnames(dom.ests_out)<-c("est.dom","SE.dom","t.dom")
    add.ests_out<-data.frame(ests.all[rows.add,1:3]); colnames(add.ests_out)<-c("est.add","SE.add","t.add")
    qtl.ests_out<-cbind(add.ests_out,dom.ests_out)
    ests.out<-rbind(covar.ests_out,qtl.ests_out)
  }
  #add nas if epistasis
  if(n.epi>0){
    for (i in 1:n.epi) ests.out<-rbind(ests.out,NA)     
  }
  
  fits[[i]]<-lod.out
  lod.out<-sum(lod.out[which(rownames(lod.out)!=colnames(cov)),"LOD"])
  lod.out<-data.frame(lod.out)
  colnames(lod.out)<-paste("lod",i, sep=".")
  rownames(lod.out)<-phe
  
  return(list(best,lod.out,best.mod,best.fit))
}
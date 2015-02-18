#Author JT Lovell
#Date 9-Feb 2014
#Version 2.2

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
cistrans.eqtl<-function(cross, chromosome, position, phe, pens=NULL, forms.in=NULL,trt=covar, wiggle=1){
  
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
  }else{
    forms<-forms.in
  }
  #get name for phe, so that we index by character not number
  if(is.numeric(phe)){
    phe<-phenames(cross)[phe]
  }
  mods<-list()
  fits<-list()
  
  #if our position is not 100%, allow the position to move
  if(wiggle>0){
    s1.1<-as.data.frame(scanone(cross, pheno.col=phe, method="hk",  addcovar=trt))
    s1.2<-as.data.frame(scanone(cross, pheno.col=phe, method="hk",  addcovar=trt, intcovar=trt))
    s1.wiggle<-cbind(s1.1[s1.1$chr==chromosome & s1.1$pos<position+wiggle & s1.1$pos>position-wiggle,],
                     s1.2$lod[s1.1$chr==chromosome & s1.1$pos<position+wiggle & s1.1$pos>position-wiggle])
    wig.sum<-s1.wiggle[,3]+s1.wiggle[,4]
    position.new<-s1.wiggle$pos[which(wig.sum==max(wig.sum))[1]]
    wiggle.move<-abs(position.new-position)
    position<-position.new
  }
  # the first two models need to be built manually because of a single QTL
  cis.eqtl<-makeqtl(cross=cross, chr=chromosome, pos=position, what="prob")
  mods[[1]]<-cis.eqtl
  mods[[2]]<-cis.eqtl
  fit.1<-fitqtl(cross, qtl=cis.eqtl, 
                formula=forms[1], pheno.col=phe, 
                covar=trt, method="hk", dropone=T)
  
  fit.2<-fitqtl(cross, qtl=cis.eqtl, 
                formula=forms[2], pheno.col=phe, 
                covar=trt, method="hk", dropone=T)
  
  lod.1<-data.frame(fit.1$result.drop)
  fits[[1]]<-lod.1
  lod.nullmod<-lod.1[cis.eqtl$name,"LOD"]
  lod.1<-sum(lod.1[which(rownames(lod.1)!=colnames(trt)),"LOD"])-lod.nullmod
  lod.2<-data.frame(fit.2$result.drop)
  fits[[2]]<-lod.2
  lod.2<-sum(lod.2[which(rownames(lod.2)!=colnames(trt)),"LOD"])-lod.nullmod
  lod.all<-data.frame(lod.1,lod.2)
  cis.eqtl$name
  # fit the remaining 7 models
  for (i in 3:length(forms)){
    form.in<-forms[i]
    scan <- addqtl(cross, qtl=cis.eqtl, formula=form.in, 
                   method="hk", covar=trt, pheno.col=phe)
    mod <- addtoqtl(cross, cis.eqtl, max(scan)$chr, max(scan)$pos)
    if(length(unique(mod$chr))==1 & abs(diff(mod$pos)) < 35){
      scan <- addqtl(cross, qtl=cis.eqtl, formula=form.in, 
                     chr = chrnames(cross)[-which(chrnames(cross)==cis.eqtl$chr)],
                     method="hk", covar=trt, pheno.col=phe)
      mod <- addtoqtl(cross, cis.eqtl, max(scan)$chr, max(scan)$pos)
    }
    fit <- fitqtl(cross, qtl=mod, 
                  formula=form.in, pheno.col=phe, 
                  covar=trt, method="hk", dropone=T, get.ests=T)
    lod.out<-data.frame(fit$result.drop)
    fits[[i]]<-lod.out
    lodi<-sum(lod.out[which(rownames(lod.out)!=colnames(trt)),"LOD"])-lod.nullmod
    lod.all<-cbind(lod.all,lodi)
    mods[[i]]<-mod
  }
  #deterime which mode fit is best by taking the best pLOD score from the QTL and QTL*trt interactions
  plods<-lod.all-pens
  best<-which(plods == max(plods))
  best.form<-forms[best]
  best.lod<-as.numeric(lod.all[best])
  best.mod<-mods[[best]]
  best.fit<-fits[[best]]
  #make the output object
  all.out<-data.frame(rownames(best.fit)); colnames(all.out)[1]<-"term.id"
  
  #add chr and pos to output
  pos.out<-sapply(rownames(best.fit), function(x) strsplit(x,"@")[[1]][2])
  pos.out[grep(":",pos.out)]<-NA
  all.out$pos<-pos.out
  chr.out<-sapply(rownames(best.fit), function(x) strsplit(x,"@")[[1]][1])
  chr.out[is.na(pos.out)]<-NA
  all.out$chr<-chr.out
  #add category information
  cat1<-rownames(best.fit)
  category<-rownames(best.fit)
  category[intersect(grep("trt", cat1), grep(":", cat1,invert=TRUE))]<-"trt"
  category[intersect(grep("trt", cat1, invert=TRUE), grep(":", cat1))]<-"epi"
  
  category[category==cis.eqtl$name]<-"cis"
  cis.trt.int<-category[intersect(grep("trt", cat1), grep(cis.eqtl$name, cat1))]
  category[category==cis.trt.int]<-"cis.trt.int"
  
  trans.id<-category[intersect(grep("@", category), grep(":", category,invert=TRUE))]
  category[category==trans.id]<-"trans"
  trans.trt.int<-category[intersect(grep("@", category), grep(":", category))]
  category[category==trans.trt.int]<-"trans.trt.int"
  all.out$category<-category
  
  #add phenotype name
  all.out$phenotype<-phe
  
  #add estimates
  ests.all<-summary(fit)$ests
  rows.dom<-intersect(grep("d",rownames(ests.all)), grep(":",rownames(ests.all),invert=TRUE))
  rows.add<-intersect(intersect(grep("a",rownames(ests.all)), grep(":",rownames(ests.all),invert=TRUE)), 
                      grep(colnames(covar),rownames(ests.all),invert=TRUE))
  rows.cov<-which(rownames(ests.all) == "trt")
  covar.ests_out<-data.frame(t(c(ests.all[rows.cov,],NA,NA,NA)))
  colnames(covar.ests_out)<-c("est.add","SE.add","t.add","est.dom","SE.dom","t.dom")
  rownames(covar.ests_out)<-"trt"
  dom.ests_out<-data.frame(ests.all[rows.dom,1:3]); colnames(dom.ests_out)<-c("est.dom","SE.dom","t.dom")
  add.ests_out<-data.frame(ests.all[rows.add,1:3]); colnames(add.ests_out)<-c("est.add","SE.add","t.add")
  qtl.ests_out<-cbind(add.ests_out,dom.ests_out)
  ests.out<-rbind(covar.ests_out,qtl.ests_out)
  ests.out$term.id<-gsub("a","",rownames(ests.out))
  
  all.out<-merge(all.out, ests.out, by="term.id", all.x=T)
  
  lod.df<-data.frame(best.fit)
  lod.df$term.id<-rownames(lod.df)
  all.out<-merge(all.out, lod.df, by="term.id", all.x=T)
  
  move<-wiggle.move
  return(list(formula=best.form,plods=plods,ciseqtl.lod=lod.nullmod,model=best.mod,stats=all.out,cis.position.move=move))
}
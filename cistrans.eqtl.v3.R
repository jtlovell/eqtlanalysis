#Author JT Lovell
#Date 23-Feb 2014
#Version 2.3

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
cistrans.eqtl<-function(cross, chromosome, position, phe, pens=NULL, forms.in=NULL,trt=covar, wiggle=1, calc.polygenic=FALSE, fit3=FALSE){
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
  scans<-list()
  colnames(trt)<-"trt"
  #if our position is not 100%, allow the position to move
  if(wiggle>0){
    s1.1a<-scanone(cross, pheno.col=phe, method="hk",  addcovar=trt)
    scans[[1]]<-s1.1a
    s1.2a<-scanone(cross, pheno.col=phe, method="hk",  addcovar=trt, intcovar=trt)
    scans[[2]]<-s1.2a
    s1.1<-as.data.frame(s1.1a)
    s1.2<-as.data.frame(s1.2a)
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
  # fit the remaining 7 models
  for (i in 3:length(forms)){
    form.in<-forms[i]
    scan <- addqtl(cross, qtl=cis.eqtl, formula=form.in, 
                   method="hk", covar=trt, pheno.col=phe)
    scans[[i]]<-scan
    mod <- addtoqtl(cross, cis.eqtl, max(scan)$chr, max(scan)$pos)
    if(length(unique(mod$chr))==1 & abs(diff(mod$pos)) < 35){
      scan <- addqtl(cross, qtl=cis.eqtl, formula=form.in, 
                     chr = chrnames(cross)[-which(chrnames(cross)==cis.eqtl$chr)],
                     method="hk", covar=trt, pheno.col=phe)
      scans[[i]]<-scan
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
  best.scan<-scans[[best]]
  
#   if(refine.mod2){
#     r<-refineqtl(cross,)
#   }
  ### 
  # if we want to fit a 3rd QTL look here - 
  if(fit3){
    q3forms<-c("+ Q3", "+ Q3 + Q3*trt", "+ Q3 + Q1*Q3", "+ Q3 + Q3*trt + Q3*Q1")
    for (i in 1:length(q3forms)){
      form.in3<-paste(best.form,q3forms[i])
      scan3 <- addqtl(cross, qtl=best.mod, formula=form.in3, 
                     method="hk", covar=trt, pheno.col=phe)
      scans[[i+9]]<-scan
      mod3 <- addtoqtl(cross, best.mod, max(scan3)$chr, max(scan3)$pos)
      if(length(unique(mod3$chr))<3){
        chr.dup<-mod3$chr[(duplicated(mod3$chr))]
        diff.pos<-abs(diff(mod3$pos[mod3$chr==chr.dup]))
      }else{
        mod3a<-mod3
        diff.pos<-NA
      }
      if(diff.pos < 50 & !is.na(diff.pos)){
          scan3.2 <- addqtl(cross, qtl=best.mod, formula=form.in3,
                          chr = chrnames(cross)[-which(chrnames(cross)==chr.dup)],
                          method="hk", covar=trt, pheno.col=phe)
          scans[[i+9]]<-scan3.2
          mod3a <- addtoqtl(cross, best.mod, max(scan3.2)$chr, max(scan3.2)$pos)
        }else{
          mod3a<-mod3
        }
      fit3 <- fitqtl(cross, qtl=mod3a, 
                    formula=form.in3, pheno.col=phe, 
                    covar=trt, method="hk", dropone=T, get.ests=T)
      lod.out3<-data.frame(fit3$result.drop)
      fits[[i+9]]<-lod.out3
      lodi<-sum(lod.out3[which(rownames(lod.out3)!=colnames(trt)),"LOD"])-lod.nullmod
      lod.all<-cbind(lod.all,lodi)
      mods[[i+9]]<-mod3
    }
    pens3<-c((pens[best] + pens[3]),  
             (pens[best] + pens[4]),
             (pens[best] + pens[6]),
             (pens[best] + pens[7]))
    plods3<-lod.all[10:13]-pens3
    if(max(plods3)>max(plods)){
      best3<-which(plods3 == max(plods3))
      best.form<-paste(forms[best], q3forms[best3])
      best.lod<-as.numeric(lod.all[9+best3])
      best.mod<-mods[[9+best3]]
      best.fit<-fits[[9+best3]]
      best.scan<-scans[[9+best3]]
    }
  }
  #make the output object
  all.out<-data.frame(rownames(best.fit)); colnames(all.out)[1]<-"term.id"
  mgsub <- function(pattern, replacement, x, ...) {
    if (length(pattern)!=length(replacement)) {
      stop("pattern and replacement do not have the same length.")
    }
    result <- x
    for (i in 1:length(pattern)) {
      result <- gsub(pattern[i], replacement[i], result, ...)
    }
    result
  }
  cis.name<-cis.eqtl$name
  cis.ints<-c(paste(":",cis.name, sep=""),paste(cis.name,":", sep=""))
  qnames<-rownames(best.fit)
  qnames<-mgsub(c(":trt","trt:",cis.ints),rep("",4),qnames)
  chr.out<-sapply(qnames, function (x) strsplit(x,"@")[[1]][1]) 
  chr.out[chr.out=="trt"]<-NA
  pos.out<-sapply(qnames, function (x) strsplit(x,"@")[[1]][2]) 
  pos.out[pos.out=="trt"]<-NA

  all.out$chr<-chr.out
  all.out$pos<-pos.out
  
  #add category information
  trans.name<-best.mod$name[!best.mod$name==cis.name]
  category<-rownames(best.fit)
  category[intersect(grep("trt", category), grep(paste(trans.name,collapse="|"), category))]<-"trans.trt.int"
  category[intersect(grep("trt", category), grep(":", category))]<-"cis.trt.int"
  category[grep(":", category)]<-"epi"
  category[grep(paste(trans.name,collapse="|"), category)]<-"trans"
  category[category==cis.name]<-"cis"
  all.out$category<-category
  #add phenotype name
  all.out$phenotype<-phe
  
  #add estimates
  fit<-fitqtl(cross, qtl=best.mod, 
              formula=best.form, pheno.col=phe, 
              covar=trt, method="hk", dropone=T, get.ests=T)
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
  
  if(calc.polygenic){
    g<-pull.genoprob(cross, include.pos.info=T, rotate=T, omit.first.prob=F)
    chrs<-best.mod$chr
    poss<-best.mod$pos
    mars<-find.marker(cross,chrs, poss)
    to.vc<-data.frame(cbind(pull.pheno(cross,pheno.col=phe),trt))
    colnames(to.vc)[1]<-c("y")
    for(i in mars){
      gs<-g[g$marker ==i,]
      gs.out<-as.data.frame(as.numeric(unlist(apply(gs[,-c(1:4)],2,function(x) which(x==max(x))[1]))))
      colnames(gs.out)<-i
      to.vc<-cbind(to.vc,gs.out)
    }
    colnames(to.vc)[3:length(colnames(to.vc))]<-paste("Q",seq(from=1, to=length(colnames(to.vc))-2), sep="")
    snps<-as.data.frame(pull.geno(cross))
    print(head(to.vc))
    to.vc$SNP<-snps

    v1<-varComp(as.formula(best.form) , to.vc,  random= ~ibs(SNP))
    pg<-as.numeric(v1$varComps)
    err<-as.numeric(v1$sigma2)
    prop.pg<-pg/(pg+err)
    return(list(formula=best.form,
                plods=plods,
                ciseqtl.lod=lod.nullmod,
                scans=best.scan,
                stats=all.out,
                cis.position.move=move,
                polygenic=c(pg,err,prop.pg)))
  }else{
    return(list(formula=best.form,
                plods=plods,
                ciseqtl.lod=lod.nullmod,
                scans=best.scan,
                stats=all.out,
                cis.position.move=move))
  }
}
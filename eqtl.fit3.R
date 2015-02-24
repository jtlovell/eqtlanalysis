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
#final cis-trans eQTL
#author: jtlovell
#data: 18-Feb 2014

#######################################
# Part 1: Set up environment
#######################################
rm(list=ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/ph2015_eqtl")
pkg <- c("RCurl","e1071","qtl","snow","rlecuyer","plyr","ggplot2","reshape","reshape2","lme4","nlme","car","doBy","scales","foreach","doMC","doSNOW","plyr","qtlhot")
invisible(lapply(pkg, function(x) {cat(x,"..."); library(x, character.only=T, verbose=F, warn.conflicts=F,quietly=T)} ))
sessionInfo()
# import necessary functions from github
function.names<-c("counts.summary.R","quantnorm.no0.R","cistrans.eqtl.v3.R", "calc.polygenic.R")
us<-paste("https://raw.githubusercontent.com/jtlovell/eqtlanalysis/master/",function.names, sep="")
for(i in 1:length(function.names)){
  script <- getURL(us[i], ssl.verifypeer = FALSE)
  eval(parse(text = script))
}
#read in the cross object with empty phenotypes
cross<-read.cross("csvs",dir="", genotypes=c("filfil","filhal","halhal"),
                  genfile="ph2015_finalF2map_gen.csv", 
                  phefile= "ph2015_finalF2map_phe.csv",
                  alleles=c("fil","hal"),
                  na.strings="-")


#read in counts data
ge<-read.delim("Phallii_eQTL_Normalized.counts") #normalized counts from Jerry

#######################################
# Part 2: Make phenotype matrices w/ counts data
#######################################
#look at the distribution of counts for each gene- cull to a set of genes with mean non-0 exp of >5
sums<-counts.summary(ge) #get summary of data
genes.tokeep<-rownames(sums[sums$mean>=5,]) #pull out genes with expression that satisfies 
ge.cull<-data.frame(t(ge[genes.tokeep,])) #covert to dataframe
ge.cull$full.name<-rownames(ge.cull)
phes<-pull.pheno(cross)
eqtl.phes<-merge(phes, ge.cull, by="full.name")
dim(eqtl.phes)

info<-read.csv("JGI_Samples_3_6_14.csv") #experimental design information from Lowry/TEJ, via Scott
info$jgi.id<-paste(info$id,info$N,info$Plant.Number, sep="_")
idtrt<-info[,c("id","Treatment","jgi.id")]
colnames(idtrt)[c(1,3)]<-c("FH.id","full.name")
#change some typos so that it merges with our gene expression data
idtrt$full.name[grep("FH147",idtrt$full.name)]<-"FH147_G973_87"
idtrt$full.name[grep("FH078_H079.2_265",idtrt$full.name)]<-"FH078_H079_265.2"
idtrt$full.name[grep("FIL2_H328",idtrt$full.name)]<-"FIL2_H328_83"
eqtl.phes<-merge(idtrt, eqtl.phes, by="full.name")

#################################
#transform all genes by quantiles
genes.tonorm<-rownames(sums)[sums$n.na<=180 &  sums$mean>=5] #remove genes w/ >50% nas and mean <5
ge.qn<-data.frame(apply(eqtl.phes[,genes.tonorm], 2, quantnorm.no0))
ge.qn<-cbind(eqtl.phes[,1:4], ge.qn)
rownames(ge.qn)<-ge.qn$id


#################################
#transform to binary 
#gene must have a median of >=2 normalized counts and >15% 0s and <85% 0s
genes.tobinary<-rownames(sums)[sums$n.na>=54 & sums$n.na<=306 &  sums$mean>=5]
ge.binary<-eqtl.phes[,genes.tobinary]
set.seed(42)
prop.0<-apply(ge.binary,2,function(x) length(x[x==0])/length(x))
prop.0<-data.frame(colnames(ge.binary),prop.0)
colnames(prop.0)[1]<-"gene"
ge.binary<-data.frame(apply(ge.binary, 2,function(x) {
  x[x!=0]<-rnorm(length(x[x!=0]), mean=1, sd=0.0001)
  x[x==0]<-rnorm(length(x[x==0]), mean=0, sd=0.0001)
  x}))
ge.binary<-cbind(eqtl.phes[,1:4], ge.binary)
rownames(ge.binary)<-ge.binary$id

#######################################
# Part 3: Make QTL-ready cross objects
#######################################
#quantile normalized (by gene) normalized (by id) counts
cross.norm<-add.phenos(cross, newdata = ge.qn)
#binary-ized counts
cross.bin<-add.phenos(cross, newdata = ge.binary)
save(cross.bin, prop.0, file="ph2015_binary.eqtl_input.RData")
#make a cross object with a random phenotype to run permutations

#make covariates (for treatment and cis-eqtl)
covar<-data.frame(as.numeric(cross.norm$phe$Treatment)); colnames(covar)<-"covar"
cis.trt<-cbind(covar,sample(c(0,1,2),nind(cross.norm), replace=T)); colnames(cis.trt)<-c("trt","cis")
#penalty for cis-trt interaction:
cis.trt.pen<-qchisq(0.95,df=2)/2/log(10)
epi.pen<-qchisq(0.95,df=4)/2/log(10)
trt<-covar; colnames(trt)<-"trt"
save(cross, cross.norm, cross.bin,covar,cis.trt,cis.trt.pen,trt,epi.pen, file="ph2015_inputfor10kperms.RData")
#######################################
# Part 4: Run permutation (do on cluster)
#######################################
load("ph2015_inputfor10kperms.RData")
#for normal traits
# #need four sets of s1 penalties
perm.phe<-data.frame(sample(c(0,1,2),nind(cross.norm), replace=T))
colnames(perm.phe)<-"perm.phe"
cross.perm<-add.phenos(cross,perm.phe)
cross.perm<-calc.genoprob(cross.perm, step=1,error.prob=0.001, map.function="kosambi")
pens<-fourperms(cross=cross.perm,phename="perm.phe", np=10, alpha=.05, nc=1, fit3=T, covar=covar)
# 
# #for binary, need to produce a scanone curve for each penalty
cross.norm<-calc.genoprob(cross.norm, step=1,error.prob=0.001, map.function="kosambi")
series<-c(seq(from=.15, to=.25, by=.01),seq(from=.26, to=.35, by=.03),.4,.45,.5)
binperms.out<-data.frame()
for ( i in series){
  print(i)
  prop<-as.character(i)
  zeros<-rep(0,i*nind(cross))
  ones<-rep(1,nind(cross)-length(zeros))
  perm.bin<-data.frame(sample(c(zeros,ones), replace=F)); colnames(perm.bin)<-"perm.bin"
  cross.binperm<-add.phenos(cross,perm.bin)
  pens<-fourperms(cross=cross.binperm,phename="perm.bin", np=10000, alpha=.05, nc=20)
  pens.out<-data.frame(rep(i,4),paste("pen",1:4,sep=""),pens); colnames(pens.out)<-c("prop.0","penalty.type","penalty")
  binperms.out<-rbind(binperms.out,pens.out)
}
# save(perms.out,pens file="phe2015_allperm10k.RData")
# load("phe2015_binaryperm10ksim.RData")
# ggplot(binperms.out, aes(x=prop.0, y=penalty, col=penalty.type))+
#   geom_point()+geom_line()+theme_bw()+facet_grid(penalty.type~.)
load("phe2015_allperm10k.RData")
#make a vector of 9 penalties for normal
pens<-vector()
pens[1]<-3.967144
pens[2]<-6.144188-cis.trt.pen
pens[3]<-6.922426-epi.pen
pens[4]<-9.576449-(epi.pen+cis.trt.pen)
pens.all<-c(0,  # plod = 0
        cis.trt.pen, #plod = lod (Q1 *trt)- chi2.2
        pens[1], #plod = lod (Q2) - pen1
        pens[1] + cis.trt.pen, #plod= lod(Q2) +lod(Q1*trt) - (pen1 + chi2.2)
        pens[2] + cis.trt.pen, #plod= lod(Q2) +lod(Q2*trt) - (pen2 + chi2.2)
        pens[3] + epi.pen, #plod= lod(Q2) +lod(Q2*trt) - (pen3 + chi2.4) 
        pens[3] + epi.pen + cis.trt.pen, #plod= lod(Q2) +lod(Q2*trt) - (pen4 + chi2.4 + chi2.2) 
        pens[4] + epi.pen + cis.trt.pen, #plod= lod(Q2) +lod(Q2*trt) - (pen4 + chi2.4 + chi2.2) 
        pens[4] + epi.pen + cis.trt.pen + cis.trt.pen) #plod= lod(Q2) +lod(Q2*trt) - (pen4 + chi2.4 + chi2.2 + chi2.2) 

pens.all<-c(0.000000,  1.301030,  3.967144,  5.268174,  6.144188,  6.922426,  8.223456,  9.576449, 10.877479)
cross.norm<-calc.genoprob(cross.norm, step=1,error.prob=0.001, map.function="kosambi")
which(genpos$gene == "Pahalv11b044716m.g")
i=7369

trt.num<-as.numeric(trt[,1])
for(i in c(124,3215,6234,1437,3457,7235,9465)){
  test<-cistrans.eqtl(cross=cross.norm, phe=i,
                      chromosome=genpos$lg[i],
                      position=genpos$pos[i], 
                      pens=pens.all, trt=trt, wiggle=3)
  test2<-calc.polygenic(cross=cross.norm, chrs=test$chr, poss=test$pos, best.form=test$formula, trt=trt.num)
  
  print(i)
  test["formula"]; test["polygenic"]; test["stats"]
}

i=3457

pens.all3<-c(pens.all, pens3)

#for binary - generate a model for each penalty. 
binpen1<-
bin.1<-data.frame(predict(interpSpline( penalty ~ prop.0,  binperms.out[binperms.out$penalty.type=="pen1",] ), nseg = 200))
bin.2<-data.frame(predict(interpSpline( penalty ~ prop.0,  binperms.out[binperms.out$penalty.type=="pen2",] ), nseg = 200))
bin.3<-data.frame(predict(interpSpline( penalty ~ prop.0,  binperms.out[binperms.out$penalty.type=="pen3",] ), nseg = 200))
bin.4<-data.frame(predict(interpSpline( penalty ~ prop.0,  binperms.out[binperms.out$penalty.type=="pen4",] ), nseg = 200))
pens.bin<-cbind(0,  # plod = 0
        cis.trt.pen, #plod = lod (Q1 *trt)- chi2.2
        bin.1, #plod = lod (Q2) - pen1
        bin.1 + cis.trt.pen, #plod= lod(Q2) +lod(Q1*trt) - (pen1 + chi2.2)
        bin.2 + cis.trt.pen, #plod= lod(Q2) +lod(Q2*trt) - (pen2 + chi2.2)
        bin.3 + epi.pen, #plod= lod(Q2) +lod(Q2*trt) - (pen3 + chi2.4) 
        bin.3 + epi.pen + cis.trt.pen, #plod= lod(Q2) +lod(Q2*trt) - (pen4 + chi2.4 + chi2.2) 
        bin.4 + epi.pen + cis.trt.pen, #plod= lod(Q2) +lod(Q2*trt) - (pen4 + chi2.4 + chi2.2) 
        bin.4 + epi.pen + cis.trt.pen + cis.trt.pen) #plod= lod(Q2) +lod(Q2*trt) - (pen4 + chi2.4 + chi2.2 + chi2.2) 


#merge gene expression data and cross

#read in annotation
genpos<-read.csv("ph2015_inferredmappos1.csv")
length(phenames(cross.norm))
length(genpos$gene)
genes.norm<-phenames(cross.norm)[phenames(cross.norm) %in% genpos$gene]
genes.bin<-phenames(cross.bin)[phenames(cross.bin) %in% genpos$gene]
#read in mapping positions for genes with that information

# pipeline - inputs
gene.name<-i
v1.1chr<-annot
v1.1pos<-annot
v0.5scaf<-annot
map.chr<-annot
map.pos<-annot
marker<-yn
genoinfer<-yn



phys<-annot
map<-map.anot
#######################################
# Part 5 - Post processing and plotting
#######################################
#post processing:
#load the output from cistrans.eqtl (run on server)
load("ph2015_stats.failed.10keqtlnorm.RData")
genes.norm[which(genes.norm %in%failed)]
badphes<-which(genes.norm %in%failed)

colnames(plods.df)<-paste("form",1:9,sep="")
test<-apply(plods.df, 1, function(x) x/max(x))
plot(1:9,test[,1], type="n", ylim=c(0,1))

test1<-as.character(formulae)
cols<-c("red","orange","yellow","lawngreen",
        "darkturquoise","forestgreen","blue","purple","black")
for(i in 1:9) test1[test1==unique(test1)[i]]<-cols[i]
for(i in 1:length(formulae))lines(1:9,test[,i], col=test1[i])
test2<-fac2col(test1)
form.counts<-data.frame(table(formulae))
ggplot(form.counts, aes(y=Freq, x=formulae, fill=formulae))+
  geom_bar(stat="identity")+
  theme_bw()+
  scale_x_discrete(limits=c(levels(form.counts$formulae)[c(9,1,8,6,7,5,3,4,2)]))+
  scale_fill_manual(values=c("orange","lightgreen","lightgreen","lightgreen",
                              "blue","blue","blue","blue", "orange"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none")+
  ggtitle("number of phenotypes fit best by each potential model")
head(stats.df)
table(stats.df$category)
trans.qtl<-stats.df[stats.df$category=="trans",]
cis.qtl<-stats.df[stats.df$category=="cis",]
ggplot(trans.qtl, aes(x=pos, y=LOD, col=chr))+
  geom_point()+
  theme_bw()+
  facet_wrap(~chr)
ggplot(trans.qtl, aes(x=pos, col=chr, fill=chr))+
  geom_histogram()+
  theme_bw()+
  facet_wrap(~chr)
ggplot(cis.qtl, aes(x=pos, col=chr, fill=chr))+
  geom_histogram()+
  theme_bw()+
  facet_wrap(~chr)

counts.trans<-data.frame(melt(table(trans.qtl$chr, trans.qtl$pos)))
counts.trans$type<-"trans"
counts.cis$type<-"cis"
cistrans.counts<-rbind(counts.trans,counts.cis)
colnames(cistrans.counts)<-c("chr","pos","count","type")
ggplot(cistrans.counts, aes(x=pos, y=count, col=type, fill=type))+
  #geom_line()+
  theme_bw()+
  stat_smooth(method="loess")+
  facet_wrap(~chr, scale="free")
ggplot(cistrans.counts, aes(x=pos, y=count, col=type, fill=type))+
  geom_point()+
  theme_bw()+
  facet_wrap(~chr)
gp.toplot<-genpos
colnames(gp.toplot)<-c("phenotype","v1.1chr","v1.1start","chr.gp","pos.gene","lod.gp")
stats.wpos<-merge(gp.toplot, stats.df, by="phenotype")
stats.wpos<-stats.wpos[complete.cases(stats.wpos),]
dim(stats.wpos)
stats.wpos$pos<-as.numeric(as.character(stats.wpos$pos))
stats.wpos$chr.gp<-as.factor(stats.wpos$chr.gp)
stats.wpos$chr.gp <- factor(stats.wpos$chr.gp, levels=rev(levels(stats.wpos$chr.gp)) )
ggplot(stats.wpos, aes(x=pos.gene, y=pos, col=category))+
  geom_point()+
  facet_grid(chr~chr.gp, scales="free", space="free")+
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "ghostwhite")
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,strip.background = element_blank()
  ) +
  scale_x_continuous("mapping position of each gene")+
  scale_y_continuous("mapping position of each QTL")+
  theme(strip.text.y = element_text(angle = 0,hjust=0))

ggplot(stats.wpos, aes(x=pos, col=category, fill=category))+
  geom_histogram(binwidth=1)+
  facet_grid(category~chr, scales="free", space="free")+
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "ghostwhite")
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,strip.background = element_blank()
  ) +
  #scale_x_continuous("mapping position of each gene")+
  #scale_y_continuous("mapping position of each QTL")+
  theme(strip.text.y = element_text(angle = 0,hjust=0))

#extract each aspect of the output to its own list

form.list<-lapply(first500, function(x) x[[1]])
plods.list<-lapply(first500, function(x) x[[2]])
ciseqtl.lod.list<-lapply(first500, function(x) x[[3]])
model.list<-lapply(first500, function(x) x[[4]])
stats.list<-lapply(first500, function(x) x[[5]])
wiggle.list<-lapply(first500, function(x) x[[6]])
#extract stats, ciseqtllod and plods to df
stats.df<-ldply(stats.list,data.frame)
plods.df<-ldply(plods.list,data.frame)
ciseqtl<-unlist(ciseqtl.lod.list)
load("ph2015_inputfor10kperms.RData")
cross.null<-read.cross("csvs",dir="", genotypes=c("filfil","filhal","halhal"),
                  genfile="ph2015_finalF2map_gen.csv", 
                  phefile= "ph2015_finalF2map_pheNULL.csv",
                  alleles=c("fil","hal"),
                  na.strings=c("NA","-"))
save(cross.null, file="ph2015_cross.null.RData")
test<-add.phenos(cross.null, pull.pheno(cross.norm)[,5:25])
hotperm1 <- hotperm(cross = test,
                    n.quant = 300,
                    n.perm = 100,
                    lod.thrs = 3.95,
                    alpha.levels = c(0.0001,0.001,0.01,0.05,0.1),
                    drop.lod = 1.5,
                    verbose = T)

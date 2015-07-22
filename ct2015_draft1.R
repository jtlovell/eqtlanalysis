# Differential expression analysis of Hal Fil F1 population
# Script #3
# Author- JT Lovell
# Data - 10-March 2015

####################
####################
#Part 1: Raw counts processing:
rm(list=ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/ph2015_eqtl")
pkg <- c("RCurl","plyr","mclust","qtl","DESeq2","GenomicRanges","car")
invisible(lapply(pkg, function(x) {cat(x,"..."); library(x, character.only=T, verbose=F, warn.conflicts=F,quietly=T)} ))
sessionInfo()
options(warn=-1)
# import necessary functions from github
function.names<-c("genobyMclust.R","multiplot.R")
us<-paste("https://raw.githubusercontent.com/jtlovell/eqtlanalysis/master/",function.names, sep="")
for(i in 1:length(function.names)){
  script <- getURL(us[i], ssl.verifypeer = FALSE)
  eval(parse(text = script))
}
#read in ASE raw count data
counts.fil<-read.delim("Phallii-FIL.counts")
counts.hal<-read.delim("Phallii-HAL.counts")
colnames(counts.fil)[grep("FIL2_H328_83",colnames(counts.fil))]<-"FIL2_H328_383"
colnames(counts.hal)[grep("FIL2_H328_83",colnames(counts.hal))]<-"FIL2_H328_383"
bad.lines<-vector()
for(i in c("HAL","FIL","F1")){
  dat.hal<-counts.hal[,grep(i,colnames(counts.hal))]
  dat.fil<-counts.fil[,grep(i,colnames(counts.fil))]
  mean.hal<-apply(dat.hal,2,mean) ;  mean.fil<-apply(dat.fil,2,mean)
  plot(mean.hal,mean.fil, main=i)
  mod<-lm(mean.fil~mean.hal)
  ol<-outlierTest(mod)
  if(ol$bonf.p<0.0001){
    out<-names(ol$bonf);    bad.lines<-c(bad.lines, out)
  }
}

#drop contaminated individuals w/ high Cooks D for outlier test
counts.hal<-counts.hal[,-which(colnames(counts.hal) %in% bad.lines)]
counts.fil<-counts.fil[,-which(colnames(counts.fil) %in% bad.lines)]

#get info
info<-read.csv("JGI_Samples_3_6_14.csv") #experimental design information from Lowry/TEJ, via Scott
info$jgi.id<-paste(info$id,info$N,info$Plant.Number, sep="_")
idtrt<-info[,c("id","Treatment","jgi.id","VWC_July","Day_July","JGI.Batch")]
idtrt$generation<-ifelse(grepl("F1",idtrt$jgi.id),"F1", 
                         ifelse(grepl("FIL",idtrt$jgi.id) | grepl("HAL",idtrt$jgi.id) ,"F0",
                                ifelse(grepl("FH",idtrt$jgi.id),"F2","NA")))


#cull to just F1/parents
counts.fil<-counts.fil[,-grep("FH",colnames(counts.fil))]
counts.hal<-counts.hal[,-grep("FH",colnames(counts.hal))]
ids<-colnames(counts.hal)
genes<-rownames(counts.hal)
colnames(counts.fil)<-paste(colnames(counts.fil),"_fil",sep="")
colnames(counts.hal)<-paste(colnames(counts.hal),"_hal",sep="")

#quick summary of bias in F1
f1.counts.hal<-counts.hal[,grep("F1", colnames(counts.hal))]
f1.counts.fil<-counts.fil[,grep("F1", colnames(counts.fil))]
f0.counts.hal<-counts.hal[,-grep("F1", colnames(counts.hal))]
f0.counts.fil<-counts.fil[,-grep("F1", colnames(counts.fil))]
f1.genesums.hal<-apply(f1.counts.hal,1,sum)
f1.genesums.fil<-apply(f1.counts.fil,1,sum)
f0.genesums.hal<-apply(f0.counts.hal,1,sum)
f0.genesums.fil<-apply(f0.counts.fil,1,sum)

plot(log10(f1.genesums.fil+1),log10(f1.genesums.hal+1), pch=".", 
     main="total counts / gene \n (black=1:1, red=F0, blue=F1)", bty='n',
     ylim=c(0,max(log10(f0.genesums.hal+1),log10(f1.genesums.hal+1))), xlim=c(0,max(log10(f1.genesums.fil+1),log10(f0.genesums.fil+1))))
points(log10(f0.genesums.fil+1),log10(f0.genesums.hal+1), pch=".", col="red")
abline(0,1, lty=2, col="black", lwd=3)
abline(lm(log10(f1.genesums.hal+1)~log10(f1.genesums.fil+1)), lty=1, lwd=2, col="blue")
abline(lm(log10(f0.genesums.hal+1)~log10(f0.genesums.fil+1)), lty=1, lwd=2, col="red")
summary(lm(log10(f1.genesums.hal+1)~log10(f1.genesums.fil+1)))
summary(lm(log10(f0.genesums.hal+1)~log10(f0.genesums.fil+1)))
#there is more read mapping bias in the parents than the F1
#merge datasets
all<-cbind(counts.fil,counts.hal)
tot.counts<-apply(all,1,sum)
med.counts<-apply(all,1,median)
mean.counts<-apply(all,1,mean)
hist(log10(tot.counts+1), breaks=100)
hist(log10(med.counts+1), breaks=100)
hist(log10(mean.counts+1), breaks=100)

#filter to only genes w/ mean expression >5
all<-all[mean.counts>=5,]
#make a new dataset with some highly sigificant genes
library(MASS)
aov.out<-data.frame()
for(i in 1:100){
  line.info$dat<-as.numeric(all[i,])
  line.info$cis<-line.info$allele
  line.info$trans1<-ifelse(line.info$generation == "F0" & line.info$allele =="fil","Fil.F0","not")
  line.info$trans2<-ifelse(line.info$generation == "F0" & line.info$allele =="hal","Hal.F0","not")
  with(line.info,by(dat, cis, mean))
  tapply(line.info$dat, c(line.info[,c("cis","trt")]), mean)
  ddply(line.info$dat, ~ line.info$cis, mean)
  aov1<-aov(dat~cis+trans1+trans2+cis*Treatment+trans1*Treatment+trans2*Treatment, data=line.info)
  boxplot(dat ~ cis, data=line.info)
  sa<-summary(aov1)
  ps<-data.frame(rownames(all)[i],t(data.frame(summary(aov1)[[1]][5][[1]][-8])))
  colnames(ps)<-c("gene.id","cis","trans1","trans2","Treatment","cis_trt","trans1_trt","trans2_trt")
  rownames(ps)<-i
  aov.out<-rbind(aov.out,ps)  
}


#check matchup
ids[!ids %in% idtrt$jgi.id] #genes in counts data, not in info file

#create two additional files - annotation for each gene, experimental design by sample
idtrt<-idtrt[idtrt$jgi.id %in% ids,]
halcols<-idtrt; halcols$jgi.id<-paste(halcols$jgi.id,"_hal",sep="")
filcols<-idtrt; filcols$jgi.id<-paste(filcols$jgi.id,"_fil",sep="")
halcols$allele<-"hal";filcols$allele<-"fil"

line.info<-rbind(halcols,filcols)
line.info$allele<-as.factor(line.info$allele)
line.info$generation<-as.factor(line.info$generation)

annot<-read.csv("ph2015_v1.1annot.edited.csv")
annot<-annot[,c("chr","start","strand","id", "end")]
rowData<-annot[annot$id %in% rownames(all),]


#files:
dim(all)
dim(rowData)
dim(line.info)
chrs<-data.frame(rowData[,"chr"]); colnames(chrs)<-"chrs"
gr<-GRanges(seqnames = rowData$id, 
            ranges = IRanges(
              start=rowData$start, 
              end=rowData$end, 
              names=rowData$id),
            strand = rowData$strand)

counts<-data.matrix(all)

line.info<-DataFrame(line.info)
se<-SummarizedExperiment(assays = counts,
                     rowData = gr, 
                     colData = line.info,
                     verbose=T)

dds <- DESeqDataSet(se = se, design = ~ Treatment + generation + allele + generation : allele + Treatment : allele + generation : Treatment)
dds$Treatment <- relevel(dds$Treatment, "Dry")
dds$generation <- relevel(dds$generation, "F0")
dds$allele <- relevel(dds$allele, "fil")

#####################
# Part 2: Cis-trans test... workflow from Cublillos et al 2014 (Loudet's group)
# 2.1 library-size adjust data

# 2.2 binom.test - summed read counts / allele in F1; q-value adjust results

# 2.3 Variance estimation
# Notes:
# Use variance stabilized expression values (from DES)
#get normalized counts

#test for ASE in F1

#variance component modelling - separately normalize hal/fil alleles separately (count median...)

cds = estimateDispersions(dds)
vsd = getVarianceStabilizedData(cds)

vst<-varianceStabilizingTransformation(dds, blind = F)
par(mfrow=c(1,2))
plot(rank(rowMeans(counts(dds))), genefilter::rowVars(log2(counts(dds)+1)),
     main="log2(x+1) transform")
plot(rank(rowMeans(assay(vst))), genefilter::rowVars(assay(vst)),
     main="VST")

test<-getVarianceStabilizedData(vst)

dds1 <- makeExampleDESeqDataSet(m=6)
vsd1 <- varianceStabilizingTransformation(dds1, blind=TRUE)
par(mfrow=c(1,2))
plot(rank(rowMeans(counts(dds1))), genefilter::rowVars(log2(counts(dds1)+1)),
     main="log2(x+1) transform")
plot(rank(rowMeans(assay(vsd1))), genefilter::rowVars(assay(vsd1)),
     main="VST")


dds <- DESeq(dds)
save(dds, file="ct2015_deseq1.RData")
res <- results(dds)
# treatment main effect
resMFType <- results(dds, contrast=c("Treatment","Wet","Dry"))
# allele main effect 
resMFType <- results(dds, name="allele_hal_vs_fil")
# generation by allele (Trans)
resMFType <- results(dds, name="generationF1.allelehal")

all.norm<-data.frame(counts(vst, normalized=F))
colnames(all.norm)<-colnames(all)

fil.ase<-all.norm[,grep("_fil",colnames(all.norm))]
hal.ase<-all.norm[,grep("_hal",colnames(all.norm))]
fil.means<-apply(fil.ase,2,mean)
hal.means<-apply(hal.ase,2,mean)
plot(fil.means, hal.means)
fil.means<-apply(counts.fil,2,mean)
hal.means<-apply(counts.hal,2,mean)
plot(fil.means, hal.means)


id-idtrt$jgi.id[!grep("FH",idtrt$jgi.id)]
table(id.f2 %in% colnames(counts.fil))



all<-cbind()
#merge by genes






#merge to do normalizations
counts.all<-rbind(counts.fil,counts.hal)

#QC
feature <- data.frame(gc=yeastGC,length=yeastLength)
data <- newSeqExpressionSet(counts=as.matrix(geneLevelData[common,]),
                            featureData=feature[common,],
                            phenoData=data.frame(
                              conditions=c(rep("mut",2),rep("wt",2)),
                              row.names=colnames(geneLevelData)))
data




lines(cnvrt.coords( c(0,0,.5), c(.5,0,0), input='plt')$usr)


hist

plot(fil.counts.f1[,"Pahalv11b000009m.g"],hal.counts.f1[,"Pahalv11b000009m.g"])

#less stringent cuttoffs - 40% increase in markers
genes.parents.lax<-genes[med.hal.hal>10 &
                           med.fil.fil>10 &
                           med.hal.fil<2 &
                           med.fil.hal<2]
length(genes.parents.lax)
#genes.parents<-genes.parents.lax

#for remaining analyses, use 1/50 cutoffs - these will most likely give the strongest clustering for later
####################
#################### 
#Part 3: of the genes that look good in parents, which have both alleles in F1s?
fil.counts.f1<-counts.fil[,grep("F1",colnames(counts.fil))]
hal.counts.f1<-counts.hal[,grep("F1",colnames(counts.hal))]
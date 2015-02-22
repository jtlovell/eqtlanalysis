fourperms<-function(cross, phename, alpha=0.05, np=100, nc=1){
  set.seed(42)
  cat("running perm set: 1..")
  s1perm.1<-scanone(cross, phe=phename, addcovar=cis.trt, method="hk",
                    n.perm=np, n.cluster=nc,verbose=F)
  set.seed(42)
  cat("2..")
  s1perm.2<-scanone(cross, phe=phename, addcovar=cis.trt, intcovar=cis.trt$trt,  method="hk",
                    n.perm=np, n.cluster=nc, verbose=F)
  set.seed(42)
  cat("3..")
  s1perm.3<-scanone(cross, phe=phename, addcovar=cis.trt, intcovar=cis.trt$cis,  method="hk",
                    n.perm=np, n.cluster=nc, verbose=F)
  set.seed(42)
  cat("4..\n")
  s1perm.4<-scanone(cross, phe=phename, addcovar=cis.trt, intcovar=cis.trt,  method="hk",
                    n.perm=np, n.cluster=nc, verbose=F)
  pens.1<-summary(s1perm.1, alpha=alpha)[1]
  pens.2<-summary(s1perm.2, alpha=alpha)[1]
  pens.3<-summary(s1perm.3, alpha=alpha)[1]
  pens.4<-summary(s1perm.4, alpha=alpha)[1]
  return(c(pens.1,pens.2,pens.3,pens.4))
}
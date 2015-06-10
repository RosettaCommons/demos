rmsd.a<-read.table("rmsd_allVsAll.txt", row.names=1)
# minimum number of elements for clustering is 2
if(dim(rmsd.a)[1]<2){
  rmsd.c<-rmsd.a
} else {
  rmsd.d<-as.dist(rmsd.a)
  rmsd.h=hclust(rmsd.d)
  rmsd.c<-cutree(rmsd.h, h=20)
}
write.table(rmsd.c, file="cluster.txt", quote=FALSE, col.names=FALSE, row.names=TRUE, sep="\t")

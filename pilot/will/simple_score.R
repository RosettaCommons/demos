tmp = as.matrix(read.table('tmp.dat',row.names=NULL,header=F)[,-1:-3])
dis = tmp[1,]
tmp = tmp[-1,]

# for(i in 1:7) assign(paste("ref",i,sep=''),tmp[i,])
# save(list=c("ref1","ref2","ref3","ref4","ref5","ref6","ref7"),file="refs.rdata")
load('refs.rdata')

par(mfrow=c(2,4))
par(mar=c(2,2,2,2))

# plot( dis, tmp[1,], type='l', xlim=c(3,6), ylim=c(-1,1)/3 )
# plot( dis, tmp[2,], type='l', xlim=c(3,6), ylim=c(-1,1)/3 )
# plot( dis, tmp[3,], type='l', xlim=c(3,6), ylim=c(-1,1)/3 )
# plot( dis, tmp[4,], type='l', xlim=c(3,6), ylim=c(-1,1)/3 )

plot( dis, tmp[ 8,], type='l', xlim=c(3,max(dis)), ylim=c(-1,0.1)/4 ); lines(dis,tmp[1,],col=4); lines(dis,c(ref1,rep(0,length(dis)-length(ref1))),col=2)
plot( dis, tmp[ 9,], type='l', xlim=c(3,max(dis)), ylim=c(-0.1,1)/4 ); lines(dis,tmp[2,],col=4); lines(dis,c(ref2,rep(0,length(dis)-length(ref2))),col=2)
plot( dis, tmp[10,], type='l', xlim=c(3,max(dis)), ylim=c(-0.1,1)/3 ); lines(dis,tmp[3,],col=4); lines(dis,c(ref3,rep(0,length(dis)-length(ref3))),col=2)
plot( dis, tmp[11,], type='l', xlim=c(3,max(dis)), ylim=c(-1,0.1)/3 ); lines(dis,tmp[4,],col=4); lines(dis,c(ref4,rep(0,length(dis)-length(ref4))),col=2)
plot( dis, tmp[12,], type='l', xlim=c(3,max(dis)), ylim=c(-1,1)/2 ); lines(dis,tmp[5,],col=4); lines(dis,c(ref5,rep(0,length(dis)-length(ref5))),col=2)
plot( dis, tmp[13,], type='l', xlim=c(3,max(dis)), ylim=c(-1,1)/2 ); lines(dis,tmp[6,],col=4); lines(dis,c(ref6,rep(0,length(dis)-length(ref6))),col=2)
plot( dis, tmp[14,], type='l', xlim=c(3,max(dis)), ylim=c(-1,1)/2 ); lines(dis,tmp[7,],col=4); lines(dis,c(ref7,rep(0,length(dis)-length(ref7))),col=2)



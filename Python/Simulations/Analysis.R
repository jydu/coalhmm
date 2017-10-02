datum<-read.table("ILS09_gamma.csv", header=T)
datum$e1<-datum$pf1-datum$rf1
datum$e2<-datum$pf2-datum$rf2
datum$e3<-datum$pf3-datum$rf3
datum$e4<-datum$pf4-datum$rf4

library(wvioplot)
with(datum, wvioplot(e1, e2, e3, e4, at=1:4, clip=F, adjust=3)); abline(h=0)
with(datum, wvioplot(rf1, pf1, rf2, pf2, rf3, pf3, rf4, pf4, at=1:8, clip=F, adjust=1)); abline(h=0)
 
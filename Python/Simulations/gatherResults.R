library(ape)

nan<-NA

results<-data.frame(
    coal.tau1=numeric(0),
    coal.tau2=numeric(0),
    coal.tau3=numeric(0),
    coal.theta1=numeric(0),
    coal.theta2=numeric(0),
    coal.rho=numeric(0),
    GTR.a=numeric(0),
    GTR.b=numeric(0),
    GTR.c=numeric(0),
    GTR.d=numeric(0),
    GTR.e=numeric(0),
    GTR.theta=numeric(0),
    GTR.theta1=numeric(0),
    GTR.theta2=numeric(0),
    Gamma.alpha=numeric(0),
    logL=numeric(0),
    time=numeric(0),
    nbIter=numeric(0)
)
for (i in 1:4) results[, paste("pf", i, sep="")] <- numeric(0)
for (i in 1:4) results[, paste("rf", i, sep="")] <- numeric(0)


for (i in 1:10)
{
  path<-paste("Results_ILS09/Sim", i, "/scaledtrees.dnd",sep="")
  cat(path,"\n")
  newick<-lapply(strsplit(readLines(path), " "), function(x) paste(x[-(1:2)], collapse=""))
  pos<-lapply(strsplit(readLines(path), " "), function(x) return(round(1e6*as.numeric(x[2]) - round(1e6*as.numeric(x[1])))))
  rf<-numeric(4)
  for (j in 1:length(newick)) {
    tree<-read.tree(text=newick[[j]])
    mat<-cophenetic(tree)
    dHC<-mat["'1'","'2'"]
    dHG<-mat["'1'","'0'"]
    dCG<-mat["'2'","'0'"]
    if (dHG < dHC) rf[3]<-rf[3]+pos[[j]]
    else if (dCG < dHC) rf[4]<-rf[4]+pos[[j]]
    else if (dHC/2 < 0.0066) rf[1]<-rf[1]+pos[[j]]
    else rf[2]<-rf[2]+pos[[j]]
  }

  path<-paste("Results_ILS09/Sim", i, "/seq.ras.user.txt",sep="")
  if (file.exists(path))
  {
    cat(path,"\n")
    source(path)
    probpath<-paste("Results_ILS09/Sim", i, "/seq.ras.hmm.postdecode.csv", sep="")
    n<-NA
    if (! is.na(logL) & file.exists(probpath))
    {
      postdecode<-read.table(probpath, header=TRUE, sep="\t")
      n<-nrow(postdecode)
      table(postdecode[,"State"])->tbl
      pf<-numeric(4)
      for (i in 1:4)
      {	      
        pf[i]<-tbl[as.character(i-1)]
        if (is.na(pf[i])) pf[i]<-0
      }
    }
    else
    {
      pf<-rep(NA,4)
    }
    results[nrow(results)+1,]<-c(tau1, tau2, tau3, theta1, theta2, rho, GTR.a, GTR.b, GTR.c, GTR.d, GTR.e, GTR.theta, GTR.theta1, GTR.theta2, Gamma.alpha, logL, time, nbIterations, pf, rf)
  }
  #write.table(results, file=paste("ILS09_gamma.tmp", chr, ".csv", sep=""), sep="\t", row.names=FALSE, quot=FALSE)
}

write.table(results, file=paste("ILS09_gamma.csv", sep=""), sep="\t", row.names=FALSE, quot=FALSE)


source("../GenerateSSCorr.R");

numrnds=10;
MethNames=c("LIONESS", "SSN", "SWEET", "BONOBO");
set.seed(1);
eps=runif(1)*10^-16;

# vary number of genes, hold the number of samples fixed
NumG=c(50,100,200, 400, 800); NumSamples=100;
runT1=matrix(0,length(MethNames), length(NumG));
for(gcnt in seq(1, length(NumG),1)){
  for(rndcnt in seq(1, numrnds,1)){
    ExpData=matrix(rnorm(NumG[gcnt]*NumSamples, mean = 0, sd = 1), nrow = NumG[gcnt], ncol = NumSamples);
    ExpData=ExpData-min(ExpData)+eps;
    
    start_time=proc.time()
    GenerateLIONESS(ExpData);
    runT1[1,gcnt]=runT1[1,gcnt]+(proc.time()-start_time)[3];
    
    start_time=proc.time()
    GenerateSSN(ExpData);
    runT1[2,gcnt]=runT1[2,gcnt]+(proc.time()-start_time)[3];   
  
    start_time=proc.time()
    GenerateSWEET(ExpData);
    runT1[3,gcnt]=runT1[3,gcnt]+(proc.time()-start_time)[3];   
 
    start_time=proc.time()
    GenerateBONOBO(ExpData);
    runT1[4,gcnt]=runT1[4,gcnt]+(proc.time()-start_time)[3];  
  }
  runT1[,gcnt]=runT1[,gcnt]/numrnds;
  print(gcnt);
}

plot(NumG,runT1[1,], type="b", col="blue", lwd=2, lty=1, ylim=c(0,max(runT1)+0.5), xlim=c(0,NumG[length(NumG)]+1), xlab = "Number of Genes", ylab="Average runtime (s)", main = "Number of Samples = 100")
lines(NumG, runT1[2,], type="b", col="red", lwd=2, lty=1)
lines(NumG, runT1[3,], type="b", col="orange", lwd=2, lty=1)
lines(NumG, runT1[4,], type="b", col="purple", lwd=2, lty=1)
legend("topleft", legend = MethNames, col = c("blue", "red", "orange", "purple"), lty = c(1, 1, 1, 1), pch = c(1, 1, 1, 1), lwd=c(2,2,2,2))


# vary number of samples, hold the number of genes fixed
NumS=c(50,100,200,400,800); NumGenes=200;
runT2=matrix(0,length(MethNames), length(NumS));
for(scnt in seq(1, length(NumS),1)){
  for(rndcnt in seq(1, numrnds,1)){
    ExpData=matrix(rnorm(NumGenes*NumS[scnt], mean = 0, sd = 1), nrow = NumGenes, ncol = NumS[scnt]);
    ExpData=ExpData-min(ExpData)+eps;
    
    start_time=proc.time()
    GenerateLIONESS(ExpData);
    runT2[1,scnt]=runT2[1,scnt]+(proc.time()-start_time)[3];
    
    start_time=proc.time()
    GenerateSSN(ExpData);
    runT2[2,scnt]=runT2[2,scnt]+(proc.time()-start_time)[3];   
    
    start_time=proc.time()
    GenerateSWEET(ExpData);
    runT2[3,scnt]=runT2[3,scnt]+(proc.time()-start_time)[3];   
    
    start_time=proc.time()
    GenerateBONOBO(ExpData);
    runT2[4,scnt]=runT2[4,scnt]+(proc.time()-start_time)[3];  
  }
  runT2[,scnt]=runT2[,scnt]/numrnds;
  print(scnt);
}

plot(NumS,runT2[1,], type="b", col="blue", lwd=2, lty=1, ylim=c(0,max(runT2)+0.5), xlim=c(0,NumS[length(NumS)]+1), xlab = "Number of Samples", ylab="Average runtime (s)", main = "Number of Genes = 200")
lines(NumS, runT2[2,], type="b", col="red", lwd=2, lty=1)
lines(NumS, runT2[3,], type="b", col="orange", lwd=2, lty=1)
lines(NumS, runT2[4,], type="b", col="purple", lwd=2, lty=1)
legend("topleft", legend = MethNames, col = c("blue", "red", "orange", "purple"), lty = c(1, 1, 1, 1), pch = c(1, 1, 1, 1), lwd=c(2,2,2,2))












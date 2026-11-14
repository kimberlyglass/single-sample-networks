GenerateSSCorr <- function(LocExp){
  
  # determine size of input data
  NumNode=dim(LocExp)[1];
  NumSamp=dim(LocExp)[2];
  
  # create filter to extract the upper diagonal
  uF=upper.tri(matrix(0, NumNode, NumNode));
  
  # pre-computation step for LIONESS / SSN / SWEET
  AgNet=cor(t(LocExp)); AgNet=AgNet[uF];
  
  # pre-computation step for BONOBO (simplified for large N)
  AgCov=cov(t(LocExp));
  
  # additional pre-computation step for SWEET
  PCCmat=cor(LocExp); muPCC=(rowSums(PCCmat)-1)/(NumSamp-1);
  minPCC=min(muPCC); maxPCC=max(muPCC);
  Sp=(muPCC-minPCC+0.01)/(maxPCC-minPCC+0.01); # x=0.01
  
  # matrices to store the computed networks
  LION=matrix(0, nrow = length(AgNet), ncol = NumSamp);
  SSN=matrix(0, nrow = length(AgNet), ncol = NumSamp);
  SWEET=matrix(0, nrow = length(AgNet), ncol = NumSamp);
  BONO=matrix(0, nrow = length(AgNet), ncol = NumSamp);
  
  idxvec=c(1:NumSamp);
  diag_idx=seq(1,NumNode*NumNode,NumNode+1);
  for(idx in seq(1,NumSamp, 1)){
    # compute Pearson correlation with all but sample idx;  used for LIONESS / SSN / SWEET
    LocNet=cor(t(LocExp[,idxvec[-idx]])); LocNet=LocNet[uF];
    
    # LIONESS
    LION[,idx]=NumSamp*(AgNet-LocNet)+LocNet;
    
    # SSN
    SSN[,idx]=(NumSamp-2)*(AgNet-LocNet)/(1-LocNet^2);
    
    # SWEET
    SWEET[,idx]=Sp[idx]*0.1*(NumSamp-1)*(AgNet-LocNet)+LocNet;
    
    # BONOBO
    LocCov=cov(t(LocExp[,idxvec[-idx]]));
    cov_diag=LocCov[diag_idx];
    delta=1/(3+2*mean(sqrt(cov_diag))/var(cov_diag));
    sscov=delta*NumSamp*(AgCov-LocCov)+LocCov;
    sds=matrix(0,NumNode,NumNode); sds[diag_idx]=sscov[diag_idx]^(-1/2);
    BON=sds %*% sscov %*% sds;
    BONO[,idx]=BON[uF];
  }
  AllSSN=list(LION,SSN,SWEET,BONO);
  return(AllSSN);
}


GenerateLIONESS <-function(LocExp){
  
  # determine size of input data
  NumNode=dim(LocExp)[1];
  NumSamp=dim(LocExp)[2];
  
  # create filter to extract the upper diagonal
  uF=upper.tri(matrix(0, NumNode, NumNode));
  
  # pre-computation step for LIONESS / SSN / SWEET
  AgNet=cor(t(LocExp)); AgNet=AgNet[uF];
  
  # matrices to store the computed networks
  LION=matrix(0, nrow = length(AgNet), ncol = NumSamp); 
  
  idxvec=c(1:NumSamp);
  for(idx in seq(1,NumSamp, 1)){
    # compute Pearson correlation with all but sample idx;  used for LIONESS / SSN / SWEET
    LocNet=cor(t(LocExp[,idxvec[-idx]])); LocNet=LocNet[uF];
    
    # LIONESS
    LION[,idx]=NumSamp*(AgNet-LocNet)+LocNet;
  }
  return(LION);
}


GenerateSSN <- function(LocExp){
  
  # determine size of input data
  NumNode=dim(LocExp)[1];
  NumSamp=dim(LocExp)[2];
  
  # create filter to extract the upper diagonal
  uF=upper.tri(matrix(0, NumNode, NumNode));
  
  # pre-computation step for LIONESS / SSN / SWEET
  AgNet=cor(t(LocExp)); AgNet=AgNet[uF];
  
  # matrices to store the computed networks
  SSN=matrix(0, nrow = length(AgNet), ncol = NumSamp);

  idxvec=c(1:NumSamp);
  for(idx in seq(1,NumSamp, 1)){
    # compute Pearson correlation with all but sample idx;  used for LIONESS / SSN / SWEET
    LocNet=cor(t(LocExp[,idxvec[-idx]])); LocNet=LocNet[uF];
    
    # SSN
    SSN[,idx]=(NumSamp-2)*(AgNet-LocNet)/(1-LocNet^2);
  }
  return(SSN);
}


GenerateSWEET <- function(LocExp){
  
  # determine size of input data
  NumNode=dim(LocExp)[1];
  NumSamp=dim(LocExp)[2];
  
  # create filter to extract the upper diagonal
  uF=upper.tri(matrix(0, NumNode, NumNode));
  
  # pre-computation step for LIONESS / SSN / SWEET
  AgNet=cor(t(LocExp)); AgNet=AgNet[uF];

  # additional pre-computation step for SWEET
  PCCmat=cor(LocExp); muPCC=(rowSums(PCCmat)-1)/(NumSamp-1);
  minPCC=min(muPCC); maxPCC=max(muPCC);
  Sp=(muPCC-minPCC+0.01)/(maxPCC-minPCC+0.01); # x=0.01
  
  # matrices to store the computed networks
  SWEET=matrix(0, nrow = length(AgNet), ncol = NumSamp);

  idxvec=c(1:NumSamp);
  for(idx in seq(1,NumSamp, 1)){
    # compute Pearson correlation with all but sample idx;  used for LIONESS / SSN / SWEET
    LocNet=cor(t(LocExp[,idxvec[-idx]])); LocNet=LocNet[uF];

    # SWEET
    SWEET[,idx]=Sp[idx]*0.1*(NumSamp-1)*(AgNet-LocNet)+LocNet;
  }
  return(SWEET);
}


GenerateBONOBO <- function(LocExp){
  
  # determine size of input data
  NumNode=dim(LocExp)[1];
  NumSamp=dim(LocExp)[2];
  
  # create filter to extract the upper diagonal
  uF=upper.tri(matrix(0, NumNode, NumNode));
  
  # pre-computation step for BONOBO (simplified for large N)
  AgCov=cov(t(LocExp));

  # matrices to store the computed networks
  BONO=matrix(0, nrow = sum(uF), ncol = NumSamp);
  
  idxvec=c(1:NumSamp);
  diag_idx=seq(1,NumNode*NumNode,NumNode+1);
  for(idx in seq(1,NumSamp, 1)){
    # BONOBO
    LocCov=cov(t(LocExp[,idxvec[-idx]]));
    cov_diag=LocCov[diag_idx];
    delta=1/(3+2*mean(sqrt(cov_diag))/var(cov_diag));
    sscov=delta*NumSamp*(AgCov-LocCov)+LocCov;
    sds=matrix(0,NumNode,NumNode); sds[diag_idx]=sscov[diag_idx]^(-1/2);
    BON=sds %*% sscov %*% sds;
    BONO[,idx]=BON[uF];
  }
  return(BONO);
}
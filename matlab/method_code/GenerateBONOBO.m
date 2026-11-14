function [AllSSN]=GenerateBONOBO(LocExp)

% determine size of input data
[NumNode,NumSamp]=size(LocExp); 

% create filter to extract upper diagonal since these are all symmetric matrices / saves memory; NodeIdx are the indices of the elements returned
uF=triu(ones(NumNode),1);
uF=uF(:)==1;

% pre-computation step for BONOBO (simplified for large N)
AgCov=cov(LocExp');

% matrices to store the computed networks
BONO=zeros(sum(uF),NumSamp);

for(idx=1:NumSamp)
	% vector to extract all but sample idx
	idxvec=[1:(idx-1),(idx+1):NumSamp];

	% BONOBO
	LocCov=cov(LocExp(:,idxvec)'); % BONOBO: covariance_matrix
	cov_diag=LocCov(1:(NumNode+1):end);
	delta=1/(3+2*mean(sqrt(cov_diag))/var(cov_diag)); % BONOBO: delta
	sscov=delta*NumSamp*(AgCov-LocCov)+LocCov;
	sds=zeros(NumNode); sds(1:(NumNode+1):end)=(sscov(1:(NumNode+1):end).^(-1/2)); % sample-specific standard deviations
    BON=sds*sscov*sds;
	BONO(:,idx)=BON(uF);
end
% merge matrices into one structure to return
AllSSN=BONO;


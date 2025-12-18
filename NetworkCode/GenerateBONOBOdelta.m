function [BONO]=GenerateBONOBOdelta(LocExp,delta);

% determine size of input data
[NumNode,NumSamp]=size(LocExp); 

% create filter to extract upper diagonal since these are all symmetric matrices / saves memory; NodeIdx are the indices of the elements returned
uF=triu(ones(NumNode),1); [i,j]=ind2sub([NumNode,NumNode], find(uF)); NodeIdx=[i,j];
uF=uF(:)==1;

% pre-computation step for BONOBO (simplified for large N)
AgCov=cov(LocExp');


BONO=zeros(sum(uF),NumSamp,length(delta));

for(idx=1:NumSamp);
	% vector to extract all but sample idx
	idxvec=[1:(idx-1),(idx+1):NumSamp];

	% BONOBO
	LocCov=cov(LocExp(:,idxvec)'); % BONOBO: covariance_matrix
	for(dcnt=1:length(delta))
		sscov=delta(dcnt)*NumSamp*(AgCov-LocCov)+LocCov;
		sds=zeros(NumNode); sds(1:(NumNode+1):end)=(sscov(1:(NumNode+1):end).^(-1/2)); % sample-specific standard deviations
		BON=sds*sscov*sds;
		BONO(:,idx,dcnt)=BON(uF);
	end
end

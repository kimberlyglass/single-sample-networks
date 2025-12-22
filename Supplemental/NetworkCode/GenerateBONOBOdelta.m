function [BONO]=GenerateBONOBOdelta(LocExp,delta);

% determine size of input data
[NumNode,NumSamp]=size(LocExp); 

% create filter to extract upper diagonal since these are all symmetric matrices / saves memory; NodeIdx are the indices of the elements returned
uF=triu(ones(NumNode),1); [i,j]=ind2sub([NumNode,NumNode], find(uF)); NodeIdx=[i,j];
uF=uF(:)==1;

BONO=zeros(sum(uF),NumSamp,length(delta));

for(idx=1:NumSamp);
	% vector to extract all but sample idx
	idxvec=[1:(idx-1),(idx+1):NumSamp];

	% BONOBO
	cov_mat=cov(LocExp(:,idxvec)'); % BONOBO: covariance_matrix
	cov_diag=cov_mat(1:(NumNode+1):end);
	A=LocExp(:,idx)-mean(LocExp,2);
	for(dcnt=1:length(delta))
		sscov=delta(dcnt)*A*A'+(1-delta(dcnt))*cov_mat; % BONOBO: sscov
		sscov_diag=sqrt(sscov(1:(NumNode+1):end)); sscov_diag(sscov_diag==0)=1; % pull out the diagonal of sscov, set any zeros equal to 1
		sscov_eye=eye(NumNode); sscov_eye(1:(NumNode+1):end)=sscov_diag; % turn this into a diagonal matrix; BONOBO: diag
		sds=inv(sscov_eye);
		BON=sds*sscov*sds;
		BONO(:,idx,dcnt)=BON(uF);
	end
end


function [AllSSN, AgNet, NodeIdx]=GenerateLinearSSN(LocExp);

% determine size of input data
[NumNode,NumSamp]=size(LocExp); 

% create filter to extract upper diagonal since these are all symmetric matrices / saves memory; NodeIdx are the indices of the elements returned
uF=triu(ones(NumNode),1); [i,j]=ind2sub([NumNode,NumNode], find(uF)); NodeIdx=[i,j];
uF=uF(:)==1;

% pre-computation step for LIONESS / SSN / SWEET
AgNet=corr(LocExp'); AgNet=AgNet(uF);

% additional pre-computation step for SWEET
PCCmat=corr(LocExp); muPCC=(sum(PCCmat)-1)/(NumSamp-1);
minPCC=min(muPCC); maxPCC=max(muPCC);

% matrices to store the computed networks
LION=zeros(length(AgNet),NumSamp);
SSN=zeros(length(AgNet),NumSamp);
SWEET=zeros(length(AgNet),NumSamp);
BONO=zeros(length(AgNet),NumSamp);

for(idx=1:NumSamp);
	% vector to extract all but sample idx
	idxvec=[1:(idx-1),(idx+1):NumSamp];

	% compute Pearson correlation with all but sample idx; used for LIONESS / SSN / SWEET
	LocNet=corr(LocExp(:,idxvec)'); LocNet=LocNet(uF);

	% LIONESS
	LION(:,idx)=(NumSamp)*(AgNet-LocNet)+LocNet; % wq=1/N

	% SSN
	SSN(:,idx)=(NumSamp-2)*(AgNet-LocNet)./(1-LocNet.*LocNet);

	% SWEET
	% compute Pearson correlation with all samples with sample idx duplicated; used for SWEET
	LocNet=corr([LocExp, LocExp(:,idx)]'); LocNet=LocNet(uF);
	Sp=(muPCC(idx)-minPCC+0.01)/(maxPCC-minPCC+0.01); % x=0.01
	SWEET(:,idx)=Sp*0.1*NumSamp*(LocNet-AgNet)+AgNet; % K=0.1

	% BONOBO
	cov_mat=cov(LocExp(:,idxvec)'); % BONOBO: covariance_matrix
	cov_diag=cov_mat(1:(NumNode+1):end);
	delta=1/(3+2*mean(sqrt(cov_diag))/var(cov_diag)); % BONOBO: delta
	A=LocExp(:,idx)-mean(LocExp,2);
	sscov=delta*A*A'+(1-delta)*cov_mat; % BONOBO: sscov
	sscov_diag=sqrt(sscov(1:(NumNode+1):end)); sscov_diag(sscov_diag==0)=1; % pull out the diagonal of sscov, set any zeros equal to 1
	sscov_eye=eye(NumNode); sscov_eye(1:(NumNode+1):end)=sscov_diag; % turn this into a diagonal matrix; BONOBO: diag
	sds=inv(sscov_eye);
	BON=sds*sscov*sds;
	BONO(:,idx)=BON(uF);
end
% merge matrices into one structure to return
AllSSN={LION, SSN, SWEET, BONO};


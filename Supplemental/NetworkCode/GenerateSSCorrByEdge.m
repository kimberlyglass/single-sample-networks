function [AllSSN, AgNet]=GenerateLinearSSN(LocExp,NodeIdx);

% determine size of input data
[NumNode,NumSamp]=size(LocExp); 
LocNode=length(NodeIdx); % NodeIdx is a list of nodes between which to calculate the SSCorr

% create filter to extract upper diagonal since these are all symmetric matrices / saves memory
uF=triu(ones(LocNode),1); 
uF=uF(:)==1;

% pre-computation step for LIONESS / SSN / SWEET
AgNet=corr(LocExp(NodeIdx,:)'); AgNet=AgNet(uF);

% additional pre-computation step for SWEET; requires all nodes in dataset
PCCmat=corr(LocExp); muPCC=(sum(PCCmat)-1)/(NumSamp-1);
minPCC=min(muPCC); maxPCC=max(muPCC);

% additional pre-computation step for BONOBO; requires all nodes in dataset
delta=zeros(NumSamp,1);
for(idx=1:NumSamp)
	idxvec=[1:(idx-1),(idx+1):NumSamp];
	cov_diag=var(LocExp(:,idxvec),[],2)';
	delta(idx)=1/(3+2*mean(sqrt(cov_diag))/var(cov_diag)); % BONOBO: delta
end

% matrices to store the computed networks
LION=zeros(length(AgNet),NumSamp);
SSN=zeros(length(AgNet),NumSamp);
SWEET=zeros(length(AgNet),NumSamp);
BONO=zeros(length(AgNet),NumSamp);

for(idx=1:NumSamp);
	% vector to extract all but sample idx
	idxvec=[1:(idx-1),(idx+1):NumSamp];

	% compute Pearson correlation with all but sample idx; used for LIONESS / SSN
	LocNet=corr(LocExp(NodeIdx,idxvec)'); LocNet=LocNet(uF);

	% LIONESS
	LION(:,idx)=(NumSamp)*(AgNet-LocNet)+LocNet; % wq=1/N

	% SSN
	SSN(:,idx)=(NumSamp-2)*(AgNet-LocNet)./(1-LocNet.*LocNet);

	% SWEET
	% compute Pearson correlation with all samples with sample idx duplicated; used for SWEET
	LocNet=corr([LocExp(NodeIdx,:), LocExp(NodeIdx,idx)]'); LocNet=LocNet(uF);
	Sp=(muPCC(idx)-minPCC+0.01)/(maxPCC-minPCC+0.01); % x=0.01
	SWEET(:,idx)=Sp*0.1*NumSamp*(LocNet-AgNet)+AgNet; % K=0.1

	% BONOBO
	cov_mat=cov(LocExp(NodeIdx,idxvec)'); % BONOBO: covariance_matrix
	A=LocExp(NodeIdx,idx)-mean(LocExp(NodeIdx,:),2);
	sscov=delta(idx)*A*A'+(1-delta(idx))*cov_mat; % BONOBO: sscov
	sscov_diag=sqrt(sscov(1:(LocNode+1):end)); sscov_diag(sscov_diag==0)=1; % pull out the diagonal of sscov, set any zeros equal to 1
	sscov_eye=eye(LocNode); sscov_eye(1:(LocNode+1):end)=sscov_diag; % turn this into a diagonal matrix; BONOBO: diag
	sds=inv(sscov_eye);
	BON=sds*sscov*sds;
	BONO(:,idx)=BON(uF);
end
% merge matrices into one structure to return
AllSSN={LION, SSN, SWEET, BONO};


function [AllSSN, BgNet]=GenerateSSCorrBG(AllExp,bgsamp);

bgsamp=ismember([1:size(AllExp,2)], bgsamp); % turn the index values of bgsamp into a boolean vector

% determine size of input data
BgExp=AllExp(:,bgsamp); % expression for the 'background' samples
BgSamp=size(BgExp,2); % number of 'background' samples
LocExp=AllExp(:,~bgsamp); % expression for the 'samples of interest'
[NumNode,NumSamp]=size(LocExp);

% create filter to extract upper diagonal since these are all symmetric matrices / saves memory; NodeIdx are the indices of the elements returned
uF=triu(ones(NumNode),1); [i,j]=ind2sub([NumNode,NumNode], find(uF)); NodeIdx=[i,j];
uF=uF(:)==1;

% pre-computation step for LIONESS / SSN / SWEET
BgNet=corr(BgExp'); BgNet=BgNet(uF);

% pre-computation step for BONOBO (simplified for large N)
BgCov=cov(BgExp');

% additional pre-computation step for SWEET
PCCmat=corr([LocExp, BgExp]);
muPCC=mean(PCCmat(1:NumSamp,(NumSamp+1):(NumSamp+BgSamp)),2); % how correlated, on average, each (non-bg) sample (i.e. a sample in rows 1:NumSamp) is with the bg samples (those in columns NumSamp+1 to the end)
minPCC=min(muPCC); maxPCC=max(muPCC);

% additional pre-computation step for BONOBO
% delta is equal to a constant since it is computed on all but the sample of interest (ie the bg samples)
cov_mat=cov(BgExp'); % BONOBO: covariance_matrix
cov_diag=cov_mat(1:(NumNode+1):end);
delta=1/(3+2*mean(sqrt(cov_diag))/var(cov_diag)); % BONOBO: delta

% matrices to store the computed networks
LION=zeros(length(BgNet),NumSamp);
SSN=zeros(length(BgNet),NumSamp);
SWEET=zeros(length(BgNet),NumSamp);
BONO=zeros(length(BgNet),NumSamp);

for(idx=1:NumSamp);
	% compute Pearson correlation with all bg samples plus idx; used for LIONESS / SSN / SWEET
	LocNet=corr([LocExp(:,idx), BgExp]'); LocNet=LocNet(uF);

	% LIONESS
	LION(:,idx)=(BgSamp+1)*(LocNet-BgNet)+BgNet; % wq=1/N

	% SSN
	SSN(:,idx)=(BgSamp-1)*(LocNet-BgNet)./(1-BgNet.*BgNet);

	% SWEET
	Sp=(muPCC(idx)-minPCC+0.01)/(maxPCC-minPCC+0.01); % x=0.01
	SWEET(:,idx)=Sp*0.1*BgSamp*(LocNet-BgNet)+BgNet; % K=0.1

	% BONOBO
	LocCov=cov([LocExp(:,idx), BgExp]'); % BONOBO: covariance_matrix with bg samples + idx
	sscov=delta*(BgSamp+1)*(LocCov-BgCov)+BgCov;
	sds=zeros(NumNode); sds(1:(NumNode+1):end)=(sscov(1:(NumNode+1):end).^(-1/2)); % sample-specific standard deviations
	BON=sds*sscov*sds;
	BONO(:,idx)=BON(uF);
end
% merge matrices into one structure to return
AllSSN={LION, SSN, SWEET, BONO};


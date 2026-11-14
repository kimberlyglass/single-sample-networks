function [AllSSN]=GenerateSWEET(LocExp)

% determine size of input data
[NumNode,NumSamp]=size(LocExp); 

% create filter to extract upper diagonal since these are all symmetric matrices / saves memory; NodeIdx are the indices of the elements returned
uF=triu(ones(NumNode),1);
uF=uF(:)==1;

% pre-computation step for LIONESS / SSN / SWEET
AgNet=corr(LocExp'); AgNet=AgNet(uF);

% additional pre-computation step for SWEET
PCCmat=corr(LocExp); muPCC=(sum(PCCmat)-1)/(NumSamp-1);
minPCC=min(muPCC); maxPCC=max(muPCC);

% matrices to store the computed networks
SWEET=zeros(sum(uF),NumSamp);

for(idx=1:NumSamp)
	% vector to extract all but sample idx
	idxvec=[1:(idx-1),(idx+1):NumSamp];

	% compute Pearson correlation with all but sample idx; used for LIONESS / SSN / SWEET
	LocNet=corr(LocExp(:,idxvec)'); LocNet=LocNet(uF);

	% SWEET
	Sp=(muPCC(idx)-minPCC+0.01)/(maxPCC-minPCC+0.01); % x=0.01
	SWEET(:,idx)=Sp*0.1*(NumSamp-1)*(AgNet-LocNet)+LocNet; % K=0.1
end
% merge matrices into one structure to return
AllSSN={SWEET};


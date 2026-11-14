function [AllSSN]=GenerateLIONESS(LocExp)

% determine size of input data
[NumNode,NumSamp]=size(LocExp); 

% create filter to extract upper diagonal since these are all symmetric matrices / saves memory; NodeIdx are the indices of the elements returned
uF=triu(ones(NumNode),1);
uF=uF(:)==1;

% pre-computation step for LIONESS / SSN / SWEET
AgNet=corr(LocExp'); AgNet=AgNet(uF);

% matrices to store the computed networks
LION=zeros(sum(uF),NumSamp);

for(idx=1:NumSamp)
	% vector to extract all but sample idx
	idxvec=[1:(idx-1),(idx+1):NumSamp];

	% compute Pearson correlation with all but sample idx; used for LIONESS / SSN / SWEET
	LocNet=corr(LocExp(:,idxvec)'); LocNet=LocNet(uF);

	% LIONESS
	LION(:,idx)=(NumSamp)*(AgNet-LocNet)+LocNet; % wq=1/N
end
% merge matrices into one structure to return
AllSSN=LION;


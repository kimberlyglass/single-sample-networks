function [LION, AgNet]=GenerateLIONESS_WPCC(LocExp, weights)
% note, the values in weights must not be negative and must sum to one; this code skips data checks to save time.

numTFs=size(LocExp,1);
uF=triu(ones(numTFs),1); uF=uF(:)==1;

AgNet=weightedcorrs(LocExp', weights');
AgNet=AgNet(uF);

LION=zeros(length(AgNet),size(LocExp,2));
for(idx=1:size(LocExp,2));
	idxvec=[1:(idx-1),(idx+1):size(LocExp,2)];
	LocNetW=weightedcorrs(LocExp(:,idxvec)', weights(idxvec)'/sum(weights(idxvec)));
	LocNetW=LocNetW(uF);

	% LIONESS
	LION(:,idx)=(AgNet-LocNetW)/weights(idx)+LocNetW;
end

end

function R = weightedcorrs(Y, w)

[nobs, nvar] = size(Y); % nobs: number of observations; nvar: number of variables
wmean = w' * Y; % weighted means of X
temp = Y - repmat(w' * Y, nobs, 1); % center X by remove weighted means
temp = temp' * (temp .* repmat(w, 1, nvar)); % weighted covariance matrix normalized to nobs
temp = 0.5 * (temp + temp'); % Must be exactly symmetric
R = diag(temp); % weighted variances normalized to nobs
R = temp ./ sqrt(R * R'); % Matrix of Weighted Correlation Coefficients

end

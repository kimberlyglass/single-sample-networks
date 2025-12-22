function [AllSSN, AgMI, NodeIdx]=GenerateSSNL(LocExp);

[NumNode,NumSamp] = size(LocExp);
uF=triu(ones(NumNode),1); [i,j]=ind2sub([NumNode,NumNode], find(uF)); NodeIdx=[i,j];
uF=uF(:)==1;

% things calculated ahead of time
AgMI=CalcMIM(LocExp');
AgMI=AgMI(uF);

% precomputation step for CSN
% Define the neighborhood of each plot
boxsize=0.1; % default parameter value in CSN code; could be changed
upper=zeros(NumNode,NumSamp);
lower=zeros(NumNode,NumSamp);
for(i=1:NumNode)
	[s1,s2] = sort(LocExp(i,:));
	n3 = NumSamp-sum(sign(s1)); % note that sign(s1) assumes data is positive valued; could change to abs(sign(s1)) if data has negative values
	h = round(boxsize/2*sum(sign(s1))); % same comment as above regarding sign(s1)
	k = 1;
	while k <= NumSamp
		s = 0;
		while k+s+1 <= NumSamp && s1(k+s+1) == s1(k)
			s = s+1;
		end
		if s >= h
			upper(i,s2(k:k+s)) = LocExp(i,s2(k));
			lower(i,s2(k:k+s)) = LocExp(i,s2(k));
		else
			upper(i,s2(k:k+s)) = LocExp(i,s2(min(NumSamp,k+s+h)));
			lower(i,s2(k:k+s)) = LocExp(i,s2(max(n3*(n3>h)+1,k-h)));
		end
		k = k+s+1;
	end
end


LIONMI=zeros(length(AgMI),NumSamp);
CSN=zeros(length(AgMI),size(LocExp,2));
for(idx=1:NumSamp);
	idxvec=[1:(idx-1),(idx+1):NumSamp];

	% LIONESS::MI
	LocMI=CalcMIM(LocExp(:,idxvec)'); LocMI=LocMI(uF);
	LIONMI(:,idx)=(NumSamp)*(AgMI-LocMI)+LocMI;

	% CSN
	B = LocExp <= repmat(upper(:,idx),1,NumSamp) & LocExp >= repmat(lower(:,idx),1,NumSamp);
	a = sum(B,2);
	d = (B*B'*NumSamp-a*a')./sqrt((a*a').*((NumSamp-a)*(NumSamp-a)')/(NumSamp-1)+eps);
	CSN(:,idx)=d(uF); % CSN original code had a step removed elements of d<0: d=d.*(d>0);

end
AllSSN={LIONMI, CSN};

end


% the functions below calculate mutual information

function MIM=CalcMIM(data)
% calculates the mutual information between all columns in data
% returns upper diagonal of the MIM

p=0.05; % percentage of data points to use to calculate k
[ndat,n]=size(data);
MIM=nan(n);
for(n1=1:n)
	for(n2=(n1+1):n)
		MIM(n1,n2)=mi_cont_cont(data(:,n1),data(:,n2),ceil(p*ndat)); % sets k (k-nearest neighbors) to p% of the total number of data points
        end
end

end


% the functions below are based on code from https://github.com/otoolej/mutual_info_kNN

function mi = mi_cont_cont(x,y,k)

% define parameters
N = length(x);
xy = [x(:) y(:)];

% calculate the k-th nearest neighbour:
nn_mdlx = createns(xy, 'Distance','chebychev');
[~, dist_kth] = knnsearch(nn_mdlx, xy, 'k', k + 1);
dist_kth = dist_kth(:, k + 1);

% find all points within distance
[xs, idx] = sort(x); [ys, idy] = sort(y);
[~, idx] = sort(idx); [~, idy] = sort(idy);
nx_range = rangesearch_var_radius(xs, idx, dist_kth);
ny_range = rangesearch_var_radius(ys, idy, dist_kth);

% calculate MI, based on: Kraskov, A., Stögbauer, H., & Grassberger, P. (2004). Estimating mutual information. Physical Review E, 69(6), 16. https://doi.org/10.1103/PhysRevE.69.066138
mi = psi(k) + psi(N) - mean(psi(nx_range + 1) + psi(ny_range + 1));
mi(mi < 0) = 0;

end


function nx_range = rangesearch_var_radius(x, idx, d)

N = length(x);
nx_range = zeros(1, N);
for n = 1:N
	nx_range(n) = rangesearch_1point(x, idx(n), d(n), N);
end

end


function p = rangesearch_1point(x, ix0, d, N)

x0 = x(ix0);
p = 0;
for n = (ix0 + 1):N
	if(abs(x(n) - x0) >= d)
		break;
	else
		p = p + 1;
	end
end
for n = (ix0 - 1):-1:1
	if(abs(x(n) - x0) >= d)
		break;
	else
		p = p + 1;
	end
end

end

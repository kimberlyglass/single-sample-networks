function EvalMethodsMixtureNL(datafile, NumRGenes, seltype, RandIdx, dtag);

MethName={'LIONESS::MI', 'CSN'};
NumMeth=length(MethName);
[NumTiss, NumRnds, ~]=size(RandIdx);
NumSel=length(RandIdx{1,1,1});

rng(1);
% read in the data
MixData=load(datafile);
[uTissue,~,tLoc]=unique(MixData.Tissue(:));

% select genes from original data to use for the analysis
sigma=std(MixData.ExpData,[],2);
if(strcmp(seltype, 'rand'))
	gidx=randperm(length(MixData.GeneSymbols)); gidx=gidx(1:NumRGenes); % random genes
elseif(strcmp(seltype, 'high'))
	[~,gidx]=sort(sigma, 'descend'); gidx=gidx(1:NumRGenes); % most variable genes
elseif(strcmp(seltype, 'low'))
	[~,gidx]=sort(sigma, 'ascend'); gidx=gidx(1:NumRGenes); % least variable genes
else
	disp('sel variable not correctly defined, defaulting to selecting random genes');
	gidx=randperm(length(MixData.GeneSymbols)); gidx=gidx(1:NumRGenes); % random genes
end
AllExp=MixData.ExpData(gidx,:);

% assumes all samples from the same tissue have the same network
nfilter=triu(ones(NumRGenes),1)==1;
AllGS=zeros(sum(nfilter(:)),size(AllExp,2));
for(tcnt=1:length(uTissue))
       	locGS=CalcMIM(AllExp(:, tLoc==tcnt)');
        AllGS(:,tLoc==tcnt)=repmat(locGS(nfilter), 1, sum(tLoc==tcnt));
end

% calculate correlation of inferred sample-specific network with `correct' tissue network and with a `different' tissue's network
CorrVals=cell(NumTiss,1);
RandCorrVals=cell(NumTiss,1);
for(cnt=1:1:NumTiss)
	disp(cnt);
	CorrVals{cnt}=zeros(NumMeth, NumRnds, NumSel);
	RandCorrVals{cnt}=zeros(NumMeth, NumRnds, NumSel);
	for(rcnt=1:NumRnds)
		tfilter=RandIdx{cnt,rcnt,1};
		rfilter=RandIdx{cnt,rcnt,2};
		LocExp=AllExp(:,tfilter); % expression for these samples
		LocGS=AllGS(:,tfilter); % correct tissue network for each of these samples
		RandGS=AllGS(:,rfilter); % different tissue network for each of these samples

		% generate single-sample networks
		[SSNets,AgNet]=GenerateSSNL(LocExp);

		% evaluation
		for(mcnt=1:NumMeth)
			NetCorr=corr(SSNets{mcnt},LocGS);
			CorrVals{cnt}(mcnt,rcnt,:)=NetCorr(1:(NumSel+1):end);
			NetCorr=corr(SSNets{mcnt},RandGS);
			RandCorrVals{cnt}(mcnt,rcnt,:)=NetCorr(1:(NumSel+1):end);
		end
	end
end

save([dtag, '_', seltype, '.mat'], 'CorrVals', 'RandCorrVals');

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

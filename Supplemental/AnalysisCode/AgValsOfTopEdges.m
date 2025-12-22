function AgValsOfTopEdges(datafile, selTissue, NodeIdx, AgValFile);

MixData=load(datafile);

InSilicoExp1=MixData.ExpData(:,strcmp(MixData.Tissue, selTissue{1}));
InSilicoExp2=MixData.ExpData(:,strcmp(MixData.Tissue, selTissue{2}));
NumEdges=size(NodeIdx,1);

CorrVals=zeros(NumEdges,2);
MIVals=zeros(NumEdges,2);
for(ecnt=1:NumEdges)
	temp=corr(InSilicoExp1(NodeIdx(ecnt,:),:)');
	CorrVals(ecnt,1)=temp(1,2);
	temp=corr(InSilicoExp2(NodeIdx(ecnt,:),:)');	
	CorrVals(ecnt,2)=temp(1,2);
	temp=CalcMIM(InSilicoExp1(NodeIdx(ecnt,:),:)');
	MIVals(ecnt,1)=temp(1,2);
	temp=CalcMIM(InSilicoExp2(NodeIdx(ecnt,:),:)');
	MIVals(ecnt,2)=temp(1,2);
end

save(AgValFile, 'CorrVals', 'MIVals', 'selTissue', 'NodeIdx');

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


function MakeInSilicoData(numGenes,numsampOrg,numsampRnd,datafile);

rng(1);
% numsampOrg=250;
% numGenes=6;
sidx=1:2*numsampOrg;

% make expression data with three groups; all corr, two modules, random/no corr
vec=sort(normrnd(0,1,[1,numsampOrg]));
xy=[vec,vec;vec,-vec];
ExpData=repmat(xy,[numGenes/2,1])+normrnd(0,0.1,[numGenes,2*numsampOrg]);
ExpData=ExpData([1:2:numGenes,2:2:numGenes],:); % reorder the nodes for easier visualization
ExpData=[ExpData,normrnd(0,1, numGenes, numsampRnd)];
ExpData=ExpData-min(ExpData(:))+eps; % NOTE: THIS IS BECAUSE CSN DOES NOT WORK WITH NEGATIVE DATA VALUES!!
numsamp=2*numsampOrg+numsampRnd;

% calculate the sample-specific networks; these functions report the upper diagonal of the square matrix
[SSNetsC,~,NodeIdx]=GenerateSSCorr(ExpData);
[SSNetsNL,~]=GenerateSSNL(ExpData);
SSNets=[SSNetsC,SSNetsNL];

save(datafile, 'ExpData', 'SSNets', 'numsampOrg', 'numsampRnd', 'NodeIdx');

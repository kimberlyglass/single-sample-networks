function CalcBGNets(datafile, netfile)

load(datafile);
[numGenes, numsamp]=size(ExpData);

% {sample set 1, sample set 2, sample set 3}
bgidx={[1:numsampOrg], [(numsampOrg+1):2*numsampOrg], [(2*numsampOrg+1):numsamp]};

SSNets=cell(length(bgidx),1);
bgsamp=cell(length(bgidx),1);
for(bgcnt=1:length(bgidx))
	[SSNets{bgcnt},~]=GenerateSSCorrBG(ExpData,bgidx{bgcnt});
	bgsamp{bgcnt}=ismember([1:size(ExpData,2)], bgidx{bgcnt}); % turn the index values of bgidx into a boolean vector

	% set the values of the background samples to nan
	for(mcnt=1:length(SSNets{bgcnt}))
		temp=nan(size(SSNets{bgcnt}{mcnt},1), numsamp);
        	temp(:,~bgsamp{bgcnt})=SSNets{bgcnt}{mcnt};
		SSNets{bgcnt}{mcnt}=temp;
	end
end

save(netfile, 'ExpData', 'SSNets', 'numsampOrg', 'numsampRnd', 'NodeIdx', 'bgsamp');

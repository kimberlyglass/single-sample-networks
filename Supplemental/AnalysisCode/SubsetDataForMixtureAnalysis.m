function RandIdx=SubsetDataForMixtureAnalysis(datafile, NumTissues, NumRnds, NumSel);

% read in the data
MixData=load(datafile);
[uTissue,~,tLoc]=unique(MixData.Tissue(:));

RandIdx=cell(length(NumTissues), NumRnds, 2);
for(cnt=1:1:length(NumTissues))
	disp(cnt);
	for(rcnt=1:NumRnds)
		randTiss=randperm(length(uTissue));
		randTiss=randTiss(1:NumTissues(cnt));
		tfilter=find(ismember(tLoc, randTiss)); % location of samples associated with the(se) tissue(s)
		tfilter=tfilter(randperm(length(tfilter))); % randomized location
		tfilter=tfilter(1:NumSel); % randomized selection of sample associated with the(se) tissue(s)
		% LocGS=AllGS(:,tfilter); % networks for these samples

		% for the random GS makes sure the assigned 'random' GS is from the other (or another for N>2) tissue in this subsample
		locT=tLoc(tfilter);
		[uT,~,uTidx]=unique(locT);
		rfilter=zeros(NumSel,1);
		for(tcnt=1:length(uT))
			selOpt=find(ismember(tLoc, randTiss) & ~ismember(tLoc,uT(tcnt))); % samples that are *not* in this tissue
			selOpt=selOpt(randperm(length(selOpt))); % randomize the order
			% RandGS(:,uTidx==tcnt)=AllGS(:,selOpt(1:sum(uTidx==tcnt))); % select a random set to be the 'random' GS
			rfilter(uTidx==tcnt)=selOpt(1:sum(uTidx==tcnt));
		end
		RandIdx{cnt,rcnt,1}=tfilter;
		RandIdx{cnt,rcnt,2}=rfilter;
	end
end

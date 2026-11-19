function MethodsMixture2_Summary(dtag, seltype, RandIdxFile)

load(RandIdxFile);

MethName={'LIONESS::PCC', 'SSN', 'SWEET', 'BONOBO'};
NumMeth=length(MethName);
[NumTiss, NumRnds, ~]=size(RandIdx);

allsel=seltype;

for(scnt=1:length(allsel))
	seltype=allsel{scnt};
	load([dtag, '_', seltype, '.mat'])

	CorrData=zeros(NumMeth, NumRnds, NumTiss);
	DiffData=zeros(NumMeth, NumRnds, NumTiss);
	for(cnt=1:NumTiss)
		CorrData(:,:,cnt)=mean(CorrVals{cnt},3);
		DiffData(:,:,cnt)=mean((CorrVals{cnt}-RandCorrVals{cnt})./(1-BGCorrVals{cnt}),3);
	end

    disp(['Median across sample-set with highest number of tissues, ', seltype, ' variance genes']);
    disp('Method: median accuracy; median specificity');
    for(mcnt=1:NumMeth)
        cc=median(squeeze(CorrData(mcnt,:,NumTiss)));
        dd=median(squeeze(DiffData(mcnt,:,NumTiss)));
        disp([MethName{mcnt}, ':', num2str(cc), '; ', num2str(dd)]);
    end
end

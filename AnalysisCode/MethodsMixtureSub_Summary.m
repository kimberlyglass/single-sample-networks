function MethodsMixtureSub_Summary(dtag, seltype, RandIdxFile, datafile)

load(RandIdxFile);

MethName={'LIONESS::PCC', 'SSN', 'SWEET', 'BONOBO'};
NumMeth=length(MethName);
[~, NumRnds, ~]=size(RandIdx);

SubtypeName=load(datafile, 'Tissue');
SubtypeName=SubtypeName.Tissue;
[uSubtypes, ~, Sloc]=unique(SubtypeName);
uSubtypes=uSubtypes([3,4,1,2]);
NumSub=length(uSubtypes);

allsel=seltype;
for(scnt=1:length(allsel))
	seltype=allsel{scnt};
	load([dtag, '_', seltype, '.mat']);

	CorrData=zeros(NumMeth,NumRnds,NumSub);
	DiffData=zeros(NumMeth,NumRnds,NumSub);

	for(rcnt=1:NumRnds)
		orgSub=SubtypeName(RandIdx{1,rcnt,1});
		for(tcnt=1:NumSub)
			locF=strcmp(orgSub, uSubtypes{tcnt});
			for(mcnt=1:NumMeth)
        		        CorrData(mcnt,rcnt,tcnt)=mean(CorrVals{1}(mcnt,rcnt,locF),3);
				DiffData(mcnt,rcnt,tcnt)=mean((CorrVals{1}(mcnt,rcnt,locF)-RandCorrVals{1}(mcnt,rcnt,locF))./(1-BGCorrVals{1}(mcnt,rcnt,locF)),3);
			end
		end
    	end
    	disp(['Median across samples in smallest subtype, ', seltype, ' variance genes']);
    	disp('Method: median accuracy; median specificity');
    	for(mcnt=1:NumMeth)
		cc=median(squeeze(CorrData(mcnt,:,NumSub)));
	        dd=median(squeeze(DiffData(mcnt,:,NumSub)));
        	disp([MethName{mcnt}, ':', num2str(cc), '; ', num2str(dd)]);
	end
end


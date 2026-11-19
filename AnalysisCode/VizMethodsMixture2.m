function VizMethodsMixture2(dtag, seltype, RandIdxFile, TissueLegend, otag);

if(strcmp(otag, ''))
	savefigs=0;
else
	savefigs=1;
end

load(RandIdxFile);

MethName={'LIONESS::PCC', 'SSN', 'SWEET', 'BONOBO'};
NumMeth=length(MethName);
[NumTiss, NumRnds, ~]=size(RandIdx);

% common visualization parameters
col=winter(NumTiss+1);
col=col(end:-1:1,:);

methloc=[1:1.5:1.5*NumMeth];
methidx=repmat(methloc', 1, NumRnds, NumTiss);
tissidx=zeros(NumMeth, NumRnds, NumTiss);
for(cnt=1:NumTiss)
	tissidx(:,:,cnt)=cnt;
end

allsel=seltype;
allCorr=zeros(NumMeth*NumRnds*NumTiss,length(allsel));
allDiff=zeros(NumMeth*NumRnds*NumTiss,length(allsel));
allMidx=zeros(NumMeth*NumRnds*NumTiss,length(allsel));
allTidx=zeros(NumMeth*NumRnds*NumTiss,length(allsel));
allMloc=zeros(NumMeth, length(allsel));
for(scnt=1:length(allsel))

	seltype=allsel{scnt};
	load([dtag, '_', seltype, '.mat'])

	CorrData=zeros(NumMeth, NumRnds, NumTiss);
	DiffData=zeros(NumMeth, NumRnds, NumTiss);
	for(cnt=1:NumTiss)
		CorrData(:,:,cnt)=mean(CorrVals{cnt},3);
		DiffData(:,:,cnt)=mean((CorrVals{cnt}-RandCorrVals{cnt})./(1-BGCorrVals{cnt}),3);
	end

	allCorr(:,scnt)=CorrData(:);
	allDiff(:,scnt)=DiffData(:);
	allMidx(:,scnt)=methidx(:)+(scnt-1)*7;
	allTidx(:,scnt)=tissidx(:);
	allMloc(:,scnt)=methloc(:)+(scnt-1)*7;
end

figure(1), bc=boxchart(allMidx(:), allCorr(:), 'groupbycolor', allTidx(:));
set(gca, 'xtick', allMloc(:), 'xticklabel', repmat(MethName, 1, length(allsel)), 'FontSize', 10)
xtickangle(45)
hold on, plot([0,max(allMidx(:))+1],[0,0], 'k--')
hold off
set(gca, 'ylim', [-0.1, 1.01], 'xlim', [0, max(allMidx(:))+1])
ylabel({'Correlation', '(Same Tissue Network)'}, 'FontSize', 10)
for(cnt=1:NumTiss)
	bc(cnt).BoxFaceColor=col(cnt,:);
	bc(cnt).MarkerColor=col(cnt,:);
end
legend(TissueLegend, 'location', 'northwest', 'FontSize', 10);

if(savefigs)
	set(gcf, 'PaperSize', [12,3], 'PaperPosition', [0,0,12,3]);
	print(gcf, [otag, '_OverallCorr.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_OverallCorr.png'], '-dpng', '-painters');
end

figure(2), bc=boxchart(allMidx(:), allDiff(:), 'groupbycolor', allTidx(:));
set(gca, 'xtick', allMloc(:), 'xticklabel', repmat(MethName, 1, length(allsel)), 'FontSize', 10)
xtickangle(45)
hold on, plot([0,max(allMidx(:))+1],[0,0], 'k--')
hold off
set(gca, 'ylim', [-0.101, 0.351], 'xlim', [0, max(allMidx(:))+1], 'ytick', [-0.1:.1:0.3])
ylabel({'Difference in Correlation', '(Same vs Other Tissue)'}, 'FontSize', 10)
for(cnt=1:NumTiss)
	bc(cnt).BoxFaceColor=col(cnt,:);
	bc(cnt).MarkerColor=col(cnt,:);
end
% legend(TissueLegend, 'location', 'northeast', 'FontSize', 12);

if(savefigs)
	set(gcf, 'PaperSize', [12,3], 'PaperPosition', [0,0,12,3]);
	print(gcf, [otag, '_DiffCorr.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_DiffCorr.png'], '-dpng', '-painters');
end

function VizMethodsMixture(dtag, seltype, RandIdxFile, TissueLegend, otag, datafile);

if(strcmp(otag, ''))
	savefigs=0;
else
	savefigs=1;
end

load(RandIdxFile);

MethName={'LIONESS::PCC', 'SSN', 'SWEET', 'BONOBO'};
NumMeth=length(MethName);
[~, NumRnds, ~]=size(RandIdx);

SubtypeName=load(datafile, 'Tissue');
SubtypeName=SubtypeName.Tissue;
[uSubtypes, ~, Sloc]=unique(SubtypeName);
uSubtypes=uSubtypes([3,4,1,2]);
[~,Sloc]=ismember(SubtypeName, uSubtypes);
NumSub=length(uSubtypes);
SubPer=100*hist(Sloc, 1:1:length(uSubtypes))/length(Sloc);

% common visualization parameters
% col=winter(NumSub+1);
col=lines(NumSub+2);
% col=col(end:-1:1,:);
col=col([1,2,4:end],:);

methloc=[1:1.5:1.5*NumMeth];
methidx=repmat(methloc', 1, NumRnds, NumSub);
subidx=zeros(NumMeth, NumRnds, NumSub);
SubLabs=cell(NumSub,1);
for(cnt=1:NumSub)
	subidx(:,:,cnt)=cnt;
    SubLabs{cnt}=[uSubtypes{cnt}, ' (', num2str(round(SubPer(cnt))), '%)'];
end


allsel=seltype;
allCorr=zeros(NumMeth*NumSub*NumRnds,length(allsel));
allRand=zeros(NumMeth*NumSub*NumRnds,length(allsel));
allDiff=zeros(NumMeth*NumSub*NumRnds,length(allsel));
allAg=zeros(NumMeth*NumSub*NumRnds,length(allsel));

allMidx=zeros(NumMeth*NumSub*NumRnds,length(allsel));
allSidx=zeros(NumMeth*NumSub*NumRnds,length(allsel));
allMloc=zeros(NumMeth,length(allsel));

for(scnt=1:length(allsel))

	seltype=allsel{scnt};
	load([dtag, '_', seltype, '.mat']);

	CorrData=zeros(NumMeth,NumRnds,NumSub);
	RandData=zeros(NumMeth,NumRnds,NumSub);
	DiffData=zeros(NumMeth,NumRnds,NumSub);
	AgData=zeros(NumMeth,NumRnds,NumSub);

	for(rcnt=1:NumRnds)
		orgSub=SubtypeName(RandIdx{1,rcnt,1});

		for(tcnt=1:NumSub)
			locF=strcmp(orgSub, uSubtypes{tcnt});
			for(mcnt=1:NumMeth)
				CorrData(mcnt,rcnt,tcnt)=mean(CorrVals{1}(mcnt,rcnt,locF),3);
				RandData(mcnt,rcnt,tcnt)=mean(RandCorrVals{1}(mcnt,rcnt,locF),3);
				DiffData(mcnt,rcnt,tcnt)=mean((CorrVals{1}(mcnt,rcnt,locF)-RandCorrVals{1}(mcnt,rcnt,locF))./(1-BGCorrVals{1}(mcnt,rcnt,locF)),3);
				AgData(mcnt,rcnt,tcnt)=mean(AgCorrVals{1}(mcnt,rcnt,locF),3);
			end
		end
	end

	allCorr(:,scnt)=CorrData(:);
	allRand(:,scnt)=RandData(:);
	allDiff(:,scnt)=DiffData(:);
	allAg(:,scnt)=AgData(:);
	allMidx(:,scnt)=methidx(:)+(scnt-1)*7;
	allSidx(:,scnt)=subidx(:);
	allMloc(:,scnt)=methloc(:)+(scnt-1)*7;
end

figure(1), bc=boxchart(allMidx(:), allDiff(:), 'groupbycolor', allSidx(:));
set(gca, 'xtick', allMloc(:), 'xticklabel', repmat(MethName, 1, length(allsel)), 'FontSize', 10)
xtickangle(45)
hold on, plot([0,max(allMidx(:))+1],[0,0], 'k--')
hold off
set(gca, 'ylim', [-0.601, 0.601], 'ytick', -0.6:0.2:0.6)
set(gca, 'xlim', [0, max(allMidx(:))+1])
ylabel({'Difference in Correlation', '(Same vs Other Subtype)'}, 'FontSize', 10)
for(cnt=1:NumSub)
    bc(cnt).BoxFaceColor=col(cnt,:);
    bc(cnt).MarkerColor=col(cnt,:);
end
% legend(uSubtypes, 'location', 'southwest', 'FontSize', 12);
legend(SubLabs, 'location', 'southwest', 'FontSize', 10);

if(savefigs)
    set(gcf, 'PaperSize', [12,3], 'PaperPosition', [0,0,12,3]);
    print(gcf, [otag, '_DiffCorr.pdf'], '-dpdf', '-painters');
    print(gcf, [otag, '_DiffCorr.png'], '-dpng', '-painters');
end


figure(2), bc=boxchart(allMidx(:), allCorr(:), 'groupbycolor', allSidx(:));
set(gca, 'xtick', allMloc(:), 'xticklabel', repmat(MethName, 1, length(allsel)), 'FontSize', 10)
xtickangle(45)
hold on, plot([0,max(allMidx(:))+1],[0,0], 'k--')
hold off
set(gca, 'ylim', [-0.151, 1.01], 'xlim', [0, max(allMidx(:))+1])
set(gca, 'xlim', [0, max(allMidx(:))+1])
ylabel({'Correlation', '(Same Subtype Network)'}, 'FontSize', 10)
for(cnt=1:NumSub)
	bc(cnt).BoxFaceColor=col(cnt,:);
	bc(cnt).MarkerColor=col(cnt,:);
end
% legend(uSubtypes, 'location', 'northwest', 'FontSize', 12);
legend(SubLabs, 'location', 'northwest', 'FontSize', 10);

if(savefigs)
    set(gcf, 'PaperSize', [12,3], 'PaperPosition', [0,0,12,3]);
    print(gcf, [otag, '_SubCorr.pdf'], '-dpdf', '-painters');
    print(gcf, [otag, '_SubCorr.png'], '-dpng', '-painters');
end

figure(3), bc=boxchart(allMidx(:), allRand(:), 'groupbycolor', allSidx(:));
set(gca, 'xtick', allMloc(:), 'xticklabel', repmat(MethName, 1, length(allsel)), 'FontSize', 10)
xtickangle(45)
hold on, plot([0,max(allMidx(:))+1],[0,0], 'k--')
hold off
set(gca, 'ylim', [-0.151, 1.01], 'xlim', [0, max(allMidx(:))+1])
ylabel({'Correlation', '(Other Subtype Network)'}, 'FontSize', 10)
for(cnt=1:NumSub)
	bc(cnt).BoxFaceColor=col(cnt,:);
	bc(cnt).MarkerColor=col(cnt,:);
end
% legend(uSubtypes, 'location', 'northeast', 'FontSize', 12);

if(savefigs)
    set(gcf, 'PaperSize', [12,3], 'PaperPosition', [0,0,12,3]);
    print(gcf, [otag, '_OthCorr.pdf'], '-dpdf', '-painters');
    print(gcf, [otag, '_OthCorr.png'], '-dpng', '-painters');
end

figure(4), bc=boxchart(allMidx(:), allAg(:), 'groupbycolor', allSidx(:));
set(gca, 'xtick', allMloc(:), 'xticklabel', repmat(MethName, 1, length(allsel)), 'FontSize', 10)
xtickangle(45)
hold on, plot([0,max(allMidx(:))+1],[0,0], 'k--')
hold off
set(gca, 'ylim', [-0.151, 1.01], 'xlim', [0, max(allMidx(:))+1])
ylabel({'Correlation', '(Aggregate Network)'}, 'FontSize', 10)
for(cnt=1:NumSub)
	bc(cnt).BoxFaceColor=col(cnt,:);
	bc(cnt).MarkerColor=col(cnt,:);
end
% legend(uSubtypes, 'location', 'northeast', 'FontSize', 12);

if(savefigs)
	set(gcf, 'PaperSize', [12,3], 'PaperPosition', [0,0,12,3]);
	print(gcf, [otag, '_AggCorr.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_AggCorr.png'], '-dpng', '-painters');
end

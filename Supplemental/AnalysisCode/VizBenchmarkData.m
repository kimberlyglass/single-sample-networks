function EvalMethodsMixture(datafile, NumRGenes, seltype, RandIdx, otag);

MethName={'LIONESS-Pearson', 'SSN', 'SWEET', 'BONOBO'};
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
       	locGS=corr(AllExp(:, tLoc==tcnt)');
        AllGS(:,tLoc==tcnt)=repmat(locGS(nfilter), 1, sum(tLoc==tcnt));
end

ecmap=hsv2rgb([repmat([39/360, 0.9], 64,1), [0.9:-(0.9/63):0]']); % orange to black
ecmap=ecmap(end:-1:1,:);

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
		[SSNets,AgNet]=GenerateSSCorr(LocExp);

		% evaluation
		for(mcnt=1:NumMeth)
			NetCorr=corr(SSNets{mcnt},LocGS);
			CorrVals{cnt}(mcnt,rcnt,:)=NetCorr(1:(NumSel+1):end);
			NetCorr=corr(SSNets{mcnt},RandGS);
			RandCorrVals{cnt}(mcnt,rcnt,:)=NetCorr(1:(NumSel+1):end);
		end


		[~,ridx]=sort(RandIdx{cnt,rcnt,1});
		[~,eidx]=sort(mean(LocGS,2), 'descend');

		figure(1),
		imagesc(LocExp(:,ridx), [0,15])
		colormap(ecmap);
		set(gca, 'xticklabel', '', 'yticklabel', '');
		xlabel('Samples'); ylabel('Genes');
		title('Gene Expression');
		colorbar

		if(length(otag)>0)
			set(gcf, 'PaperSize', [4,3], 'PaperPosition', [0,0,4,3]);
			print(gcf, [otag, '_SelData.pdf'], '-dpdf', '-painters');
			print(gcf, [otag, '_SelData.png'], '-dpng', '-painters');
		end

		figure(2),
		subplot(2,1,1), imagesc(LocGS(eidx,ridx),[-1,1]);
		colormap jet;
		set(gca, 'xticklabel', '', 'yticklabel', '');
		xlabel('Samples'); ylabel('Edges');
		title('Correct Associated Network');
		colorbar

		subplot(2,1,2), imagesc(RandGS(eidx,ridx),[-1,1]);
		colormap jet;
		set(gca, 'xticklabel', '', 'yticklabel', '');
		xlabel('Samples'); ylabel('Edges');
		title('Other Tissue Network');
		colorbar

		if(length(otag)>0)
			set(gcf, 'PaperSize', [4,4], 'PaperPosition', [0,0,4,4]);
			print(gcf, [otag, '_SelNets.pdf'], '-dpdf', '-painters');
			print(gcf, [otag, '_SelNets.png'], '-dpng', '-painters');
		end

		figure(3),
		for(mcnt=1:NumMeth)
			subplot(NumMeth,1,mcnt), imagesc(SSNets{mcnt}(eidx,ridx), [-1,1]);
			colormap jet;
			set(gca, 'xticklabel', '', 'yticklabel', '');
			xlabel('Samples'); ylabel('Edges');
			title(MethName{mcnt})
			colorbar
		end

		if(length(otag)>0)
			set(gcf, 'PaperSize', [4,8], 'PaperPosition', [0,0,4,8]);
			print(gcf, [otag, '_PredNets.pdf'], '-dpdf', '-painters');
			print(gcf, [otag, '_PredNets.png'], '-dpng', '-painters');
		end
	end
end

function VizSSNPatternsBG(netfile, otag, NodePairs, bgcnt);

MethNames={'LIONESS::Pearson', 'SSN', 'SWEET', 'BONOBO', 'Aggregate'};
NumMeth=length(MethNames)-1;

if(strcmp(otag, ''))
	savefigs=0;
else
	otag=[otag, num2str(bgcnt)];
	savefigs=1;
end

load(netfile);
[numGenes, numsamp]=size(ExpData);
sidx=1:2*numsampOrg;
SSNets=SSNets{bgcnt}; bgsamp=bgsamp{bgcnt};

% here are the two edges / node pairs for visualization
PairIdx=zeros(length(NodePairs),1);
for(cnt=1:length(NodePairs))
	PairIdx(cnt)=find(ismember(NodeIdx, NodePairs{cnt}, 'rows'));
end

% compare the edge weights computed by each algorithm for these two node pairs
figure(1)
lcolor=[118,42,131]/255; % [153,112,171]/255;
nlcolor=[27,120,55]/255;
for(mcnt1=1:NumMeth)
	for(mcnt2=1:NumMeth)
		if(mcnt1<mcnt2)
			subplot(NumMeth,NumMeth,sub2ind([NumMeth,NumMeth],mcnt1,mcnt2)), plot(SSNets{mcnt1}(PairIdx(1),sidx), SSNets{mcnt2}(PairIdx(1),sidx), '.', 'MarkerSize', 10, 'color', lcolor);
			set(gca, 'fontsize', 12)
		elseif(mcnt1>mcnt2)
			subplot(NumMeth,NumMeth,sub2ind([NumMeth,NumMeth],mcnt1,mcnt2)), plot(SSNets{mcnt1}(PairIdx(2),sidx), SSNets{mcnt2}(PairIdx(2),sidx), '.', 'MarkerSize', 10, 'color', nlcolor);
			set(gca, 'fontsize', 12)
		else
			locmax=max(max(SSNets{mcnt1}(:,sidx)));
			locmin=min(min(SSNets{mcnt1}(:,sidx)));
			locbins=[locmin:(locmax-locmin)/24:locmax];
			E1=hist(SSNets{mcnt1}(PairIdx(1),sidx), locbins);
			E2=hist(SSNets{mcnt1}(PairIdx(2),sidx), locbins);
			subplot(NumMeth,NumMeth,sub2ind([NumMeth,NumMeth],mcnt1,mcnt2)), plot(locbins, E1, '-', 'linewidth', 2, 'color', lcolor);
			hold on
			plot(locbins, E2, '-', 'linewidth', 2, 'color', nlcolor);
			hold off
			set(gca, 'fontsize', 12)
			% legend({'Nonlinear', 'Linear'}, 'location', 'northeast', 'FontSize', 10)
		end
		if(mcnt1==1)
			ylabel(MethNames{mcnt2}, 'FontSize', 14)
		end
		if(mcnt2==NumMeth)
			xlabel(MethNames{mcnt1}, 'FontSize', 14)
		end
%		if(mcnt2==1)
%			title(MethNames{mcnt1}, 'FontSize', 14)
%		end
	end
end
if(savefigs)
	set(gcf, 'PaperSize', [18,15], 'PaperPosition', [0,0,18,15]);
	print(gcf, [otag, '_EdgeComparisons.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_EdgeComparisons.png'], '-dpng', '-painters');
end


% visualize the edge weights compared to the original input data
figure(2)
nbins=64;
limvals=[min(min(ExpData(1:2,:))), max(max(ExpData(1:2,:)))];
binedges=-2:(2+2)/(nbins-1):2;
cols=jet(nbins);
binbyval=1;
% linear / set-scale
for(mcnt=1:NumMeth)
	if(binbyval==1)
		[~,locedges,bin]=histcounts(SSNets{mcnt}(PairIdx(2),:), binedges);
		bin(SSNets{mcnt}(PairIdx(2),:)<binedges(1))=1;
		bin(SSNets{mcnt}(PairIdx(2),:)>binedges(end))=nbins;
	else
		[~,locedges,bin]=histcounts(SSNets{mcnt}(PairIdx(2),:), nbins);
	end
	bin(bgsamp)=nan;

	subplot(NumMeth,4,4*mcnt-1), plot(ExpData(NodePairs{2}(1), isnan(bin)),ExpData(NodePairs{2}(2), isnan(bin)),'.', 'color', [0.7,0.7,0.7], 'MarkerSize', 20);
	hold on
	subplot(NumMeth,4,4*mcnt-1), plot(ExpData(NodePairs{2}(1), bin==1),ExpData(NodePairs{2}(2), bin==1),'.', 'color', cols(1,:), 'MarkerSize', 20);
	for(cnt=2:nbins)
		subplot(NumMeth,4,4*mcnt-1), plot(ExpData(NodePairs{2}(1), bin==cnt),ExpData(NodePairs{2}(2), bin==cnt),'.', 'color', cols(cnt,:), 'MarkerSize', 20);
	end
	set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', 14, 'xlim', limvals, 'ylim', limvals);
	hold off
	collims=[min(locedges):(max(locedges)-min(locedges))/4:max(locedges)];
	colorbar('ytick', 0:.25:1, 'yticklabel', round(collims,3), 'location', 'eastoutside');
	colormap jet;
	% title(MethNames{mcnt}, 'FontSize', 16);
	xlabel('Gene 1');
	ylabel('Gene 2');
end

% X-pattern / set-scale
for(mcnt=1:NumMeth)
	if(binbyval==1)
		[~,locedges,bin]=histcounts(SSNets{mcnt}(PairIdx(1),:), binedges);
		bin(SSNets{mcnt}(PairIdx(1),:)<binedges(1))=1;
		bin(SSNets{mcnt}(PairIdx(1),:)>binedges(end))=nbins;
	else
		[~,locedges,bin]=histcounts(SSNets{mcnt}(PairIdx(1),:), nbins);
	end
	bin(bgsamp)=nan;

	subplot(NumMeth,4,4*mcnt), plot(ExpData(NodePairs{1}(1), isnan(bin)),ExpData(NodePairs{1}(2), isnan(bin)),'.', 'color', [0.7,0.7,0.7], 'MarkerSize', 20);
	hold on
	subplot(NumMeth,4,4*mcnt), plot(ExpData(NodePairs{1}(1), bin==1),ExpData(NodePairs{1}(2), bin==1),'.', 'color', cols(1,:), 'MarkerSize', 20);
	for(cnt=2:nbins)
		subplot(NumMeth,4,4*mcnt), plot(ExpData(NodePairs{1}(1), bin==cnt),ExpData(NodePairs{1}(2), bin==cnt),'.', 'color', cols(cnt,:), 'MarkerSize', 20);
	end
	set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', 14, 'xlim', limvals, 'ylim', limvals);
	hold off
	collims=[min(locedges):(max(locedges)-min(locedges))/4:max(locedges)];
	colorbar('ytick', 0:.25:1, 'yticklabel', round(collims,3), 'location', 'eastoutside');
	colormap jet;
	% title(MethNames{mcnt}, 'FontSize', 16);
	xlabel('Gene 1');
	ylabel('Gene 4');
end

binbyval=0;
% linear / free scale
for(mcnt=1:NumMeth)
	if(binbyval==1)
		[~,locedges,bin]=histcounts(SSNets{mcnt}(PairIdx(2),:), binedges);
		bin(SSNets{mcnt}(PairIdx(2),:)<binedges(1))=1;
		bin(SSNets{mcnt}(PairIdx(2),:)>binedges(end))=nbins;
	else
		[~,locedges,bin]=histcounts(SSNets{mcnt}(PairIdx(2),:), nbins);
	end
	bin(bgsamp)=nan;

	subplot(NumMeth,4,4*mcnt-3), plot(ExpData(NodePairs{2}(1), isnan(bin)),ExpData(NodePairs{2}(2), isnan(bin)),'.', 'color', [0.7,0.7,0.7], 'MarkerSize', 20);
	hold on
	subplot(NumMeth,4,4*mcnt-3), plot(ExpData(NodePairs{2}(1), bin==1),ExpData(NodePairs{2}(1), bin==1),'.', 'color', cols(1,:), 'MarkerSize', 20);
	for(cnt=2:nbins)
		subplot(NumMeth,4,4*mcnt-3), plot(ExpData(NodePairs{2}(1), bin==cnt),ExpData(NodePairs{2}(2), bin==cnt),'.', 'color', cols(cnt,:), 'MarkerSize', 20);
	end
	set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', 14, 'xlim', limvals, 'ylim', limvals);
	hold off
	collims=[min(locedges):(max(locedges)-min(locedges))/4:max(locedges)];
	colorbar('ytick', 0:.25:1, 'yticklabel', round(collims,3), 'location', 'eastoutside');
	colormap jet;
	% title(MethNames{mcnt}, 'FontSize', 16);
	xlabel('Gene 1');
	ylabel('Gene 2');
end

% X pattern / free scale
for(mcnt=1:NumMeth)
	if(binbyval==1)
		[~,locedges,bin]=histcounts(SSNets{mcnt}(PairIdx(1),:), binedges);
		bin(SSNets{mcnt}(PairIdx(1),:)<binedges(1))=1;
		bin(SSNets{mcnt}(PairIdx(1),:)>binedges(end))=nbins;
	else
		[~,locedges,bin]=histcounts(SSNets{mcnt}(PairIdx(1),:), nbins);
	end
	bin(bgsamp)=nan;

	subplot(NumMeth,4,4*mcnt-2), plot(ExpData(NodePairs{1}(1), isnan(bin)),ExpData(NodePairs{1}(2), isnan(bin)),'.', 'color', [0.7,0.7,0.7], 'MarkerSize', 20);
	hold on
	subplot(NumMeth,4,4*mcnt-2), plot(ExpData(NodePairs{1}(1), bin==1),ExpData(NodePairs{1}(2), bin==1),'.', 'color', cols(1,:), 'MarkerSize', 20);
	for(cnt=2:nbins)
		subplot(NumMeth,4,4*mcnt-2), plot(ExpData(NodePairs{1}(1), bin==cnt),ExpData(NodePairs{1}(2), bin==cnt),'.', 'color', cols(cnt,:), 'MarkerSize', 20);
	end
	set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', 14, 'xlim', limvals, 'ylim', limvals);
	hold off
	collims=[min(locedges):(max(locedges)-min(locedges))/4:max(locedges)];
	colorbar('ytick', 0:.25:1, 'yticklabel', round(collims,3), 'location', 'eastoutside');
	colormap jet;
	% title(MethNames{mcnt}, 'FontSize', 16);
	xlabel('Gene 1');
	ylabel('Gene 4');
end
if(savefigs)
	set(gcf, 'PaperSize', [24,20], 'PaperPosition', [0,0,24,20]);
	print(gcf, [otag, '_EdgeValues.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_EdgeValues.png'], '-dpng', '-painters');
end

% look at these values across the full dataset
filter=triu(ones(numGenes),1);
OverallCorr=corr(ExpData');
LeftCorr=corr(ExpData(:,1:numsampOrg)');
RightCorr=corr(ExpData(:,(numsampOrg+1):2*numsampOrg)');
NullCorr=corr(ExpData(:,(2*numsampOrg+1):end)');
GSCorr=[repmat(LeftCorr(filter==1), [1, numsampOrg]), repmat(RightCorr(filter==1), [1, numsampOrg]), repmat(NullCorr(filter==1), [1,100])];
eridx=[1,2,3,10,14,15,4,5,6,7,8,9,11,12,13];

figure(3)
subplot(NumMeth+1,1,1), imagesc(GSCorr(eridx,:), [-2,2]);
set(gca, 'xticklabel', '', 'yticklabel', '');
title('Replicated Pearson Networks', 'FontSize', 20);
xlabel('Samples', 'FontSize', 20);
ylabel('Edges', 'FontSize', 20);
colorbar('FontSize', 14)
colormap(jet);

for(mcnt=1:NumMeth)
	subplot(NumMeth+1,1,mcnt+1), imagesc(SSNets{mcnt}(eridx,:),[-2,2]);
	set(gca, 'xticklabel', '', 'yticklabel', '');
	xlabel('Samples', 'FontSize', 20);
	ylabel('Edges', 'FontSize', 20);
	% title(MethNames{mcnt}, 'FontSize', 20);
	colorbar('FontSize', 14)
	colormap(jet);

	locData=SSNets{mcnt}(eridx,:);
	maskmatrix=ones(size(locData));
	maskmatrix(isnan(locData))=nan;
	maskim=repmat(maskmatrix~=1, 1, 1, 3)*0.7;
	hold on
	h=imagesc(maskim);
	h.AlphaData = isnan(maskmatrix);
	hold off
end
if(savefigs)
	set(gcf, 'PaperSize', [8,17.1429], 'PaperPosition', [0,0,8,17.1429]);
	print(gcf, [otag, '_PredictedNetworks.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_PredictedNetworks.png'], '-dpng', '-painters');
end

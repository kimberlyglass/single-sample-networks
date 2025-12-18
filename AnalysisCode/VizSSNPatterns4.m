function VizSSNPatterns4(datafile, otag, NodePairs); %, PairIdx);

% to improve:
% the 'Gene 1', 'Gene 4', etc labels are hard-wired. make them automatic based on the information in 'Node Pairs'
% MethNames would be better as a variable attached to the input data, right now its declared in the code

MethNames={'LIONESS::PCC', 'SSN', 'SWEET', 'BONOBO', 'LIONESS::MI', 'CSN', 'Aggregate'};
NumMeth=length(MethNames)-1;

if(strcmp(otag, ''))
	savefigs=0;
else
	savefigs=1;
end

load(datafile);
[numGenes, numsamp]=size(ExpData);
sidx=1:2*numsampOrg;

% here are the two edges / node pairs for visualization
PairIdx=zeros(length(NodePairs),1);
for(cnt=1:length(NodePairs))
	PairIdx(cnt)=find(ismember(NodeIdx, NodePairs{cnt}, 'rows'));
end

% plot the expression data for each group for these two edges
limvals=[min(min(ExpData(1:2,:))), max(max(ExpData(1:2,:)))];
idx1=1:numsampOrg;
idx2=[(numsampOrg+1):2*numsampOrg];
idx3=[(2*numsampOrg+1):numsamp];
figure(1)
subplot(3,2,1),plot(ExpData(NodePairs{2}(1), idx1),ExpData(NodePairs{2}(2), idx1), 'k.', 'MarkerSize', 10);
set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', 16, 'xlim', limvals, 'ylim', limvals);
xlabel('Gene 1', 'FontSize', 20); ylabel('Gene 2', 'FontSize', 20);
subplot(3,2,2),plot(ExpData(NodePairs{1}(1), idx1),ExpData(NodePairs{1}(2), idx1), 'k.', 'MarkerSize', 10);
set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', 16, 'xlim', limvals, 'ylim', limvals);
xlabel('Gene 1', 'FontSize', 20); ylabel('Gene 4', 'FontSize', 20);
subplot(3,2,3),plot(ExpData(NodePairs{2}(1), idx2),ExpData(NodePairs{2}(2), idx2), 'k.', 'MarkerSize', 10);
set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', 16, 'xlim', limvals, 'ylim', limvals);
xlabel('Gene 1', 'FontSize', 20); ylabel('Gene 2', 'FontSize', 20);
subplot(3,2,4),plot(ExpData(NodePairs{1}(1), idx2),ExpData(NodePairs{1}(2), idx2), 'k.', 'MarkerSize', 10);
set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', 16, 'xlim', limvals, 'ylim', limvals);
xlabel('Gene 1', 'FontSize', 20); ylabel('Gene 4', 'FontSize', 20);
subplot(3,2,5),plot(ExpData(NodePairs{2}(1), idx3),ExpData(NodePairs{2}(2), idx3), 'k.', 'MarkerSize', 10);
set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', 16, 'xlim', limvals, 'ylim', limvals);
xlabel('Gene 1', 'FontSize', 20); ylabel('Gene 2', 'FontSize', 20);
subplot(3,2,6),plot(ExpData(NodePairs{1}(1), idx3),ExpData(NodePairs{1}(2), idx3), 'k.', 'MarkerSize', 10);
set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', 16, 'xlim', limvals, 'ylim', limvals);
xlabel('Gene 1', 'FontSize', 20); ylabel('Gene 4', 'FontSize', 20);
if(savefigs)
	set(gcf, 'PaperSize', [8,10], 'PaperPosition', [0,0,8,10]);
	print(gcf, [otag, '_EdgeExamples.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_EdgeExamples.png'], '-dpng', '-painters');
end

% compare the edge weights computed by each algorithm for these two node pairs
figure(2)
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
figure(3)
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

	subplot(NumMeth,4,4*mcnt-1), plot(ExpData(NodePairs{2}(1), bin==1),ExpData(NodePairs{2}(2), bin==1),'.', 'color', cols(1,:), 'MarkerSize', 20);
	hold on
	for(cnt=2:nbins)
		subplot(NumMeth,4,4*mcnt-1), plot(ExpData(NodePairs{2}(1), bin==cnt),ExpData(NodePairs{2}(2), bin==cnt),'.', 'color', cols(cnt,:), 'MarkerSize', 20);
	end
	set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', 16, 'xlim', limvals, 'ylim', limvals);
	hold off
	collims=[min(locedges):(max(locedges)-min(locedges))/4:max(locedges)];
	colorbar('ytick', 0:.25:1, 'yticklabel', round(collims,2), 'location', 'eastoutside');
	colormap jet;
	% title(MethNames{mcnt}, 'FontSize', 16);
	xlabel('Gene 1', 'FontSize', 20);
	ylabel('Gene 2', 'FontSize', 20);
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

	subplot(NumMeth,4,4*mcnt), plot(ExpData(NodePairs{1}(1), bin==1),ExpData(NodePairs{1}(2), bin==1),'.', 'color', cols(1,:), 'MarkerSize', 20);
	hold on
	for(cnt=2:nbins)
		subplot(NumMeth,4,4*mcnt), plot(ExpData(NodePairs{1}(1), bin==cnt),ExpData(NodePairs{1}(2), bin==cnt),'.', 'color', cols(cnt,:), 'MarkerSize', 20);
	end
	set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', 16, 'xlim', limvals, 'ylim', limvals);
	hold off
	collims=[min(locedges):(max(locedges)-min(locedges))/4:max(locedges)];
	colorbar('ytick', 0:.25:1, 'yticklabel', round(collims,2), 'location', 'eastoutside');
	colormap jet;
	% title(MethNames{mcnt}, 'FontSize', 16);
	xlabel('Gene 1', 'FontSize', 20);
	ylabel('Gene 4', 'FontSize', 20);
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

	subplot(NumMeth,4,4*mcnt-3), plot(ExpData(NodePairs{2}(1), bin==1),ExpData(NodePairs{2}(1), bin==1),'.', 'color', cols(1,:), 'MarkerSize', 20);
	hold on
	for(cnt=2:nbins)
		subplot(NumMeth,4,4*mcnt-3), plot(ExpData(NodePairs{2}(1), bin==cnt),ExpData(NodePairs{2}(2), bin==cnt),'.', 'color', cols(cnt,:), 'MarkerSize', 20);
	end
	set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', 16, 'xlim', limvals, 'ylim', limvals);
	hold off
	collims=[min(locedges):(max(locedges)-min(locedges))/4:max(locedges)];
	colorbar('ytick', 0:.25:1, 'yticklabel', round(collims,2), 'location', 'eastoutside');
	colormap jet;
	% title(MethNames{mcnt}, 'FontSize', 16);
	xlabel('Gene 1', 'FontSize', 20);
	ylabel('Gene 2', 'FontSize', 20);
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

	subplot(NumMeth,4,4*mcnt-2), plot(ExpData(NodePairs{1}(1), bin==1),ExpData(NodePairs{1}(2), bin==1),'.', 'color', cols(1,:), 'MarkerSize', 20);
	hold on
	for(cnt=2:nbins)
		subplot(NumMeth,4,4*mcnt-2), plot(ExpData(NodePairs{1}(1), bin==cnt),ExpData(NodePairs{1}(2), bin==cnt),'.', 'color', cols(cnt,:), 'MarkerSize', 20);
	end
	set(gca, 'box', 'off', 'linewidth', 1, 'fontsize', 16, 'xlim', limvals, 'ylim', limvals);
	hold off
	collims=[min(locedges):(max(locedges)-min(locedges))/4:max(locedges)];
	colorbar('ytick', 0:.25:1, 'yticklabel', round(collims,2), 'location', 'eastoutside');
	colormap jet;
	% title(MethNames{mcnt}, 'FontSize', 16);
	xlabel('Gene 1', 'FontSize', 20);
	ylabel('Gene 4', 'FontSize', 20);
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
GSCorr=[repmat(LeftCorr(filter==1), [1, numsampOrg]), repmat(RightCorr(filter==1), [1, numsampOrg]), repmat(NullCorr(filter==1), [1,numsampRnd])];
eridx=[1,2,3,10,14,15,4,5,6,7,8,9,11,12,13];

figure(4)
% ecmap=hsv2rgb([repmat([277/360, 0.5], 64,1), [0.94:-(0.94/63):0]']); % purple to black
ecmap=hsv2rgb([repmat([39/360, 0.9], 64,1), [0.9:-(0.9/63):0]']); % orange to black
ecmap=ecmap(end:-1:1,:);
imagesc(ExpData)%,[-3,3])
colorbar('FontSize', 16)
colormap(ecmap);
title('Expression Data', 'FontSize', 20)
set(gca, 'xticklabel', '', 'yticklabel', 1:numGenes, 'FontSize', 20);
xlabel('Samples', 'FontSize', 20);
% ylabel('Genes', 'FontSize', 20);
if(savefigs)
	set(gcf, 'PaperSize', [8,4], 'PaperPosition', [0,0,8,4]);
	print(gcf, [otag, '_ExpData.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_ExpData.png'], '-dpng', '-painters');
end

figure(5)
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
end
if(savefigs)
	set(gcf, 'PaperSize', [8,24], 'PaperPosition', [0,0,8,24]);
	print(gcf, [otag, '_PredictedNetworks.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_PredictedNetworks.png'], '-dpng', '-painters');
end

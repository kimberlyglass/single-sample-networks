function VizDataMixture(datafile, NumRGenes, seltype, otag);

if(strcmp(otag, ''))
	savefigs=0;
else
	savefigs=1;
end

rng(1);
% read in the data
MixData=load(datafile);
[uTissue,~,tLoc]=unique(MixData.Tissue(:));
gidx=randperm(length(MixData.GeneSymbols)); RandIdx=gidx(1:NumRGenes); % index of random genes

for(scnt=1:length(seltype))
	% select genes from original data to use for the analysis
	sigma=std(MixData.ExpData,[],2);
	if(strcmp(seltype{scnt}, 'rand'))
		gidx=RandIdx;
	elseif(strcmp(seltype{scnt}, 'high'))
		[~,gidx]=sort(sigma, 'descend'); gidx=gidx(1:NumRGenes); % most variable genes
	elseif(strcmp(seltype{scnt}, 'low'))
		[~,gidx]=sort(sigma, 'ascend'); gidx=gidx(1:NumRGenes); % least variable genes
	end

	AllExp=MixData.ExpData(gidx,:);
	% assumes all samples from the same tissue have the same network
	nfilter=triu(ones(NumRGenes),1)==1;
	AllGS=zeros(sum(nfilter(:)),size(AllExp,2));
	for(tcnt=1:length(uTissue))
        	locGS=corr(AllExp(:, tLoc==tcnt)');
	        AllGS(:,tLoc==tcnt)=repmat(locGS(nfilter), 1, sum(tLoc==tcnt));
	end

	gsMu=mean(AllGS,2);
	[~,gsidx]=sort(gsMu, 'descend');

	figure(scnt), imagesc(AllGS(gsidx,:), [-1,1])
	hold on
	for(cnt=1:length(uTissue)-1)
		plot([250*cnt, 250*cnt], [0.5, nchoosek(NumRGenes,2)+0.5], 'k-', 'linewidth', 0.5);
	end
	hold off
	set(gca, 'linewidth', 1, 'xtick', [0,length(uTissue)*250], 'xticklabel', '', 'ytick', [0,NumRGenes+1], 'yticklabel', '');
	% xlabel('Samples', 'FontSize', 24);
	% ylabel('Edges', 'FontSize', 24);
	colormap(jet);
	colorbar('ytick', -1:0.5:1, 'FontSize', 16)
	if(savefigs)
		set(gcf, 'PaperSize', [6,4], 'PaperPosition', [0,0,6,4]);
		print(gcf, [otag, '_', seltype{scnt}, '_DerivedNet.pdf'], '-dpdf', '-painters');
		print(gcf, [otag, '_', seltype{scnt}, '_DerivedNet.png'], '-dpng', '-painters');
	end
end

ecmap=hsv2rgb([repmat([39/360, 0.9], 64,1), [0.9:-(0.9/63):0]']); % orange to black
ecmap=ecmap(end:-1:1,:);

sigma=std(MixData.ExpData,[],2);
[~,ridx]=sort(sigma, 'ascend');
figure(length(seltype)+1), subplot(1,3,1:2), imagesc(MixData.ExpData(ridx,:), [0,15]);
hold on
for(cnt=1:length(uTissue)-1)
        plot([250*cnt, 250*cnt], [0.5, length(ridx)+0.5], 'k-', 'linewidth', 0.5);
end
hold off
set(gca, 'linewidth', 1, 'xtick', [0,length(uTissue)*250], 'xticklabel', '', 'ytick', [0,length(ridx)+1], 'yticklabel', '');
% xlabel('Samples (grouped by tissue)')
% ylabel({'Genes'; '(\leftarrow increasing variance)'})
colormap(ecmap);

ridx=ridx(end:-1:1);
f=ismember(ridx, RandIdx); % the location of the randomly selected genes in the sorted data
subplot(1,3,3), plot(sigma(ridx), 1:1:length(ridx), 'k.');
hold on
% plot(sigma(ridx(f)), find(f), 'r.');
plot(repmat([5,5.25], sum(f), 1)', [find(f), find(f)]', 'b-', 'linewidth', 0.1)
plot([0,5.25],[NumRGenes, NumRGenes]+0.5, 'r--', 'linewidth', 1);
plot([0,5.25],length(ridx)-[NumRGenes, NumRGenes]+0.5, 'r--', 'linewidth', 1);
hold off
% xlabel('Expression Variance')
set(gca, 'ylim', [0.5,length(ridx)+0.5], 'xlim', [0,5.25], 'ytick', [0,length(ridx)+1], 'yticklabel', '', 'xtick', [0, 2.5,5], 'FontSize', 16);

if(savefigs)
	set(gcf, 'PaperSize', [10,8], 'PaperPosition', [0,0,10,8]);
	print(gcf, [otag, '_OrgData.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_OrgData.png'], '-dpng', '-painters');
end

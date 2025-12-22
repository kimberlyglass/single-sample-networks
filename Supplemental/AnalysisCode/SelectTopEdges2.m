function NodeIdx=SelectTopEdges(datafile, otag, selTissue, numRedges);

if(strcmp(otag, ''))
	savefigs=0;
else
	savefigs=1;
end

rng(2);
MixData=load(datafile);

InSilicoExp1=MixData.ExpData(:,strcmp(MixData.Tissue, selTissue{1}));
InSilicoExp2=MixData.ExpData(:,strcmp(MixData.Tissue, selTissue{2}));
[NumGenes, Num1]=size(InSilicoExp1); Num2=size(InSilicoExp2,2);

% identify exemplar edges
i=zeros(4,1); j=zeros(4,1);
CorrA=corr(InSilicoExp1'); CorrB=corr(InSilicoExp2');
CorrA(1:(NumGenes+1):end)=nan; CorrB(1:(NumGenes+1):end)=nan;
[~,~,~,stat]=ttest2(InSilicoExp1', InSilicoExp2');
DEup=repmat(abs(stat.tstat), NumGenes, 1)>25; DEup=DEup.*DEup'; % edges where both genes are DE
DEdn=repmat(abs(stat.tstat), NumGenes, 1)<4; DEdn=DEdn.*DEdn'; % eges where neither gene is DE
minE=repmat(min([InSilicoExp1,InSilicoExp2],[],2)>1, 1, NumGenes); minE=minE.*minE'; % edges where both genes have a min expression >0

% both highly positively correlated & not DE.
[i1,j1]=find(CorrA>0.85 & CorrB>0.85 & DEdn & minE);
ridx=randperm(length(i1)); i1=i1(ridx); j1=j1(ridx);
i(1)=i1(1); j(1)=j1(1);

% biggest diff in corr (and not DE)
CorrDiff=abs(CorrA-CorrB).*double(DEdn).*double(minE); maxDiff=max(CorrDiff(:));
[i2,j2]=find(CorrDiff==maxDiff); i(2)=i2(1); j(2)=j2(1);

% both low corr separately and when combined (both not DE)
[i3,j3]=find(abs(CorrA)<0.05 & abs(CorrB)<0.05 & DEdn & minE);
ridx=randperm(length(i3)); i3=i3(ridx); j3=j3(ridx);
i(3)=i3(1); j(3)=j3(1);

% both low corr originally, but high corr when combined (both strongly DE)
[i4,j4]=find(abs(CorrA)<0.05 & abs(CorrB)<0.05 & DEup & minE);
ridx=randperm(length(i4)); i4=i4(ridx); j4=j4(ridx);
i(4)=i4(1); j(4)=j4(1);

i=i([1,3,2,4]); j=j([1,3,2,4]);
NodeIdx=[i,j];


% Visualization Code
ecmap=hsv2rgb([repmat([39/360, 0.9], 64,1), [0.9:-(0.9/63):0]']); % orange to black
ecmap=ecmap(end:-1:1,:);
selTissue=regexprep(selTissue, 'esophagus_', '');

redges=ceil(numel(CorrA)*rand(numRedges,1));
figure(1), plot(CorrA(redges), CorrB(redges), 'k.');
hold on
for(ecnt=1:length(i))
        plot(CorrA(i(ecnt),j(ecnt)), CorrB(i(ecnt),j(ecnt)), 'r.', 'MarkerSize', 30);
end
hold off
set(gca, 'fontsize', 12)
set(gca, 'xtick', -1:0.5:1, 'ytick', -1:0.5:1);
xlabel(['Correlation: ', regexprep(selTissue{1}, '_', ' ')], 'FontSize', 20)
ylabel(['Correlation: ', regexprep(selTissue{2}, '_', ' ')], 'FontSize', 20)
if(savefigs)
	set(gcf, 'PaperSize', [6,5], 'PaperPosition', [0,0,6,5]);
	print(gcf, [otag, '_SelectedEdgesScatter.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_SelectedEdgesScatter.png'], '-dpng', '-painters');
end


figure(2),
for(ecnt=1:length(i))
	subplot(length(i)+1,1, ecnt),
	plot(InSilicoExp1(i(ecnt),:), InSilicoExp1(j(ecnt),:), 'b.');
	hold on
	plot(InSilicoExp2(i(ecnt),:), InSilicoExp2(j(ecnt),:), 'r.');
	hold off
	set(gca, 'fontsize', 13)
	xvals=[InSilicoExp1(i(ecnt),:), InSilicoExp2(i(ecnt),:)];
	yvals=[InSilicoExp1(j(ecnt),:), InSilicoExp2(j(ecnt),:)];
	xvals=[min(xvals), max(xvals)]; xpad=0.05*(xvals(2)-xvals(1)); xvals=[xvals(1)-xpad, xvals(2)+xpad];
	yvals=[min(yvals), max(yvals)]; ypad=0.05*(yvals(2)-yvals(1)); yvals=[yvals(1)-ypad, yvals(2)+ypad];
	set(gca, 'xlim', xvals, 'ylim', yvals);
	xlabel(MixData.GeneSymbols(i(ecnt)), 'FontSize', 18);
	ylabel(MixData.GeneSymbols(j(ecnt)), 'FontSize', 18);
end;
subplot(length(i)+1, 1, length(i)+1),
plot(InSilicoExp1(i(ecnt),:), InSilicoExp1(j(ecnt),:), 'b.');
hold on
plot(InSilicoExp2(i(ecnt),:), InSilicoExp2(j(ecnt),:), 'r.');
hold off
legend(regexprep(selTissue, '_', ' '), 'location', 'northoutside', 'FontSize', 16);
ymaxval=max(yvals);
set(gca, 'ylim', [ymaxval+1, ymaxval+2]);
axis off
if(savefigs)
	set(gcf, 'PaperSize', [4,15], 'PaperPosition', [0,0,4,15]);
	print(gcf, [otag, '_SelectedEdges.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_SelectedEdges.png'], '-dpng', '-painters');
end

[~,eridx]=sort([mean(InSilicoExp1,2)-mean(InSilicoExp2,2)].*[mean(InSilicoExp1,2)+mean(InSilicoExp2,2)], 'descend');
figure(3), imagesc([InSilicoExp1(eridx,:), InSilicoExp2(eridx,:)],[0,15]);
colormap(ecmap);
colorbar
set(gca, 'xtick', [Num1/2,Num2/2+Num1], 'xticklabel', regexprep(selTissue, '_', ' '), 'yticklabel', '', 'ytick', [0,NumGenes+1], 'FontSize', 16);
ylabel('Genes', 'FontSize', 20)
colorbar
if(savefigs)
	set(gcf, 'PaperSize', [6,7], 'PaperPosition', [0,0,6,7]);
	print(gcf, [otag, '_OrgExp.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_OrgExp.png'], '-dpng', '-painters');
end

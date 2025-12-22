function VizSubpopByVarOverTopEdges(vizdatafile, otag);

if(strcmp(otag, ''))
	savefigs=0;
else
	savefigs=1;
end

MethNames={'LIONESS::PCC', 'SSN', 'SWEET', 'BONOBO', 'LIONESS::MI', 'CSN'};
load(vizdatafile); % PerLabs, Gvar1, Gvar2
[~,NumEdges]=size(Gvar1);
[~,NumRnds]=size(Gvar1{1});

% Visualization
methidx=repmat([1:1.5:1.5*length(MethNames)]',1,length(PerLabs), NumRnds);
peridx=repmat(1:length(PerLabs),length(MethNames),1,NumRnds);

for(ecnt=1:NumEdges)
	ratdata=zeros(length(MethNames), length(PerLabs), NumRnds);
	for(mcnt=1:length(MethNames))
		ratdata(mcnt,:,:)=log2(Gvar1{mcnt,ecnt}./Gvar2{mcnt,ecnt});
	end
	xtickvals=[1:1.5:1.5*length(MethNames)];
	figure(1), subplot(NumEdges+1,1,ecnt),
	bc=boxchart(methidx(:), ratdata(:), 'groupbycolor', peridx(:));
	hold on
	plot([0,1.5*length(MethNames)+1],[0,0], 'k--');
	hold off
	set(gca, 'linewidth', 2, 'FontSize', 14, 'xlim', [min(xtickvals)-0.75,max(xtickvals)+0.75]);
	set(gca, 'xtick', xtickvals, 'xticklabel', MethNames);
	set(gca, 'ylim', [-20,20]);
	% legend(PerLabs, 'location', 'eastoutside', 'FontSize', 10)
	ylabel('log_2(\sigma^2_{1} / \sigma^2_{2})')

	col=turbo(length(PerLabs));
	for(pcnt=1:length(PerLabs))
		bc(pcnt).BoxFaceColor=col(pcnt,:);
		bc(pcnt).MarkerColor=col(pcnt,:);
	end
end
% this makes a replicate of the final plot in the last subplot space; this is only done to get the legend printed and positioned correctly; the plot is then "erased"
figure(1), subplot(NumEdges+1,1,NumEdges+1),
bc=boxchart(methidx(:), ratdata(:), 'groupbycolor', peridx(:));
legend(PerLabs, 'location', 'northoutside', 'FontSize', 10, 'Orientation', 'horizontal')
col=turbo(length(PerLabs));
for(pcnt=1:length(PerLabs))
	bc(pcnt).BoxFaceColor=col(pcnt,:);
	bc(pcnt).MarkerColor=col(pcnt,:);
end
% these next two steps "hide" the plot information, by shifting the y-limits and then turning off the axis.
set(gca, 'ylim', [20,40])
axis off

if(savefigs)
	set(gcf, 'PaperSize', [10,12], 'PaperPosition', [0,0,10,12]);
	print(gcf, [otag, '_ByVar_OverRnds.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_ByVar_OverRnds.png'], '-dpng', '-painters');
end

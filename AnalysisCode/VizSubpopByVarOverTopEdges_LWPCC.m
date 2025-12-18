function VizSubpopByVarOverTopEdges_LWPCC(vizdatafile, otag);

if(strcmp(otag, ''))
	savefigs=0;
else
	savefigs=1;
end

MethNames={'LIONESS::WPCC'};
EdgeNames={'Edge#1', 'Edge#2', 'Edge#3', 'Edge#4'};
load(vizdatafile); % PerLabs, Gvar1, Gvar2
[~,NumEdges]=size(Gvar1);
[~,NumRnds]=size(Gvar1{1});

% Visualization
edgeidx=repmat([1:1.5:1.5*NumEdges]',1,length(PerLabs), NumRnds);
peridx=repmat(1:length(PerLabs),NumEdges,1,NumRnds);

ratdata=zeros(NumEdges, length(PerLabs), NumRnds);
for(ecnt=1:NumEdges)
	ratdata(ecnt,:,:)=log2(Gvar1{1,ecnt}./Gvar2{1,ecnt});
end
xtickvals=[1:1.5:1.5*NumEdges];

figure(1), 
bc=boxchart(edgeidx(:), ratdata(:), 'groupbycolor', peridx(:));
hold on
plot([0,1.5*NumEdges+1],[0,0], 'k--');
hold off
set(gca, 'linewidth', 2, 'FontSize', 14, 'xlim', [min(xtickvals)-0.75,max(xtickvals)+0.75]);
set(gca, 'xtick', xtickvals, 'xticklabel', EdgeNames);
set(gca, 'ylim', [-20,20]);
legend(PerLabs, 'location', 'eastoutside', 'FontSize', 10)
ylabel('log_2(\sigma^2_{1} / \sigma^2_{2})')
title(MethNames{1}, 'FontSize', 16)

col=turbo(length(PerLabs));
for(pcnt=1:length(PerLabs))
	bc(pcnt).BoxFaceColor=col(pcnt,:);
	bc(pcnt).MarkerColor=col(pcnt,:);
end


if(savefigs)
	set(gcf, 'PaperSize', [6,3], 'PaperPosition', [0,0,6,3]);
	print(gcf, [otag, '_ByVar_OverRnds.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_ByVar_OverRnds.png'], '-dpng', '-painters');
end

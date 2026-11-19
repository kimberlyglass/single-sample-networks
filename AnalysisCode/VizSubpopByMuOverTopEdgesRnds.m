function VizSubpopByMuOverTopEdgesRnds(vizdatafile, otag, AgValFile)

if(strcmp(otag, ''))
	savefigs=0;
else
	savefigs=1;
end

MethNames={'LIONESS::PCC', 'SSN', 'SWEET', 'BONOBO', 'LIONESS::MI', 'CSN'};
load(vizdatafile); % PerLabs, Gmu1, Gmu2
[~,NumEdges]=size(Gvar1);
X=load(AgValFile);

% Visualization
figure(1),
pcnt=1;
for(mcnt=1:length(MethNames))
	for(ecnt=1:NumEdges)
		subplot(length(MethNames), NumEdges, pcnt),
		imagesc([Gmu1{mcnt,ecnt}'; Gmu2{mcnt,ecnt}']', [-2,2])
		hold on
	        plot([100.5,100.5],[0.5,length(PerLabs)+0.5], 'k-')
		hold off
		% ylabel('Percent Samples From Mucosa', 'FontSize', 8)
		colormap jet;
		set(gca, 'ytick', [1:1:length(PerLabs)], 'yticklabel', PerLabs, 'FontSize', 10);
        	set(gca, 'xtick', [50.5,150.5], 'xticklabel', {'mucosa', 'muscularis'})
		% title([MethNames{mcnt}, ' - Edge#', num2str(ecnt)])
	        title([MethNames{mcnt}], 'FontSize', 12)
		pcnt=pcnt+1;
    end
end

if(savefigs)
	set(gcf, 'PaperSize', [12,12], 'PaperPosition', [0,0,12,12]);
	print(gcf, [otag, '_ByMu_OverRndsHeat.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_ByMu_OverRndsHeat.png'], '-dpng', '-painters');
end

% make heat maps colored by the 'correct' Pearson / MI vals
figure(2),
for(ecnt=1:NumEdges)
    subplot(1,NumEdges,ecnt),
    imagesc([X.CorrVals(ecnt,1),X.CorrVals(ecnt,2);X.MIVals(ecnt,1),X.MIVals(ecnt,2)],[-2,2])
    text(1,1,num2str(round(X.CorrVals(ecnt,1),2)), 'HorizontalAlignment','center', 'FontSize',10)
    text(2,1,num2str(round(X.CorrVals(ecnt,2),2)), 'HorizontalAlignment','center', 'FontSize',10)
    text(1,2,num2str(round(X.MIVals(ecnt,1),2)), 'HorizontalAlignment','center', 'FontSize',10)
    text(2,2,num2str(round(X.MIVals(ecnt,2),2)), 'HorizontalAlignment','center', 'FontSize',10)
    colormap jet;
    set(gca, 'xtick', [1,2], 'xticklabel', ''); %, 'xticklabel', {'mucosa', 'muscularis'});
    set(gca, 'ytick', [1,2], 'yticklabel', {'\rho', 'MI'}, 'FontSize', 10);
end

if(savefigs)
    set(gcf, 'PaperSize', [12,0.5], 'PaperPosition', [0,0,12,0.5]);
    print(gcf, [otag, '_ByMu_OverRndsHeatLegend.pdf'], '-dpdf', '-painters');
    print(gcf, [otag, '_ByMu_OverRndsHeatLegend.png'], '-dpng', '-painters');
end

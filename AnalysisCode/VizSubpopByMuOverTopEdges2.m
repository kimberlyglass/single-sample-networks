    function VizSubpopByMuOverTopEdges2(vizdatafile, otag, AgValFile);

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
		q1=quantile(Gmu1{mcnt,ecnt}',3);
		q2=quantile(Gmu2{mcnt,ecnt}',3);
		subplot(length(MethNames), NumEdges, pcnt),
		errorbar(2*[1:1:length(PerLabs)]-0.25, q1(2,:), q1(2,:)-q1(1,:), q1(3,:)-q1(2,:), 'b.', 'markersize', 10);
		hold on
		errorbar(2*[1:1:length(PerLabs)]+0.25, q2(2,:), q2(2,:)-q2(1,:), q2(3,:)-q2(2,:), 'r.', 'markersize', 10);
		if(mcnt<5)
            		plot([0,2*length(PerLabs)+2],[X.CorrVals(ecnt,1),X.CorrVals(ecnt,1)], 'b-', 'linewidth', 1);
		    	plot([0,2*length(PerLabs)+2],[X.CorrVals(ecnt,2),X.CorrVals(ecnt,2)], 'r-', 'linewidth', 1);
        	else
            		plot([0,2*length(PerLabs)+2],[X.MIVals(ecnt,1),X.MIVals(ecnt,1)], 'b-', 'linewidth', 1);
            		plot([0,2*length(PerLabs)+2],[X.MIVals(ecnt,2),X.MIVals(ecnt,2)], 'r-', 'linewidth', 1);
	        end
		hold off
		% imagesc([Gmu1{mcnt,ecnt}'; Gmu2{mcnt,ecnt}'], [-2,2])
		% colormap jet;
		set(gca, 'xtick', 2*[1:1:length(PerLabs)], 'xticklabel', PerLabs, 'FontSize', 10, 'xlim', [0,2*length(PerLabs)+2]);
		xtickangle(45)
	        ylabel('Edge Weight', 'FontSize', 10)
		% title([MethNames{mcnt}, ' - Edge#', num2str(ecnt)])
        	title([MethNames{mcnt}], 'FontSize', 12)
		pcnt=pcnt+1;
    end
end

if(savefigs)
	set(gcf, 'PaperSize', [12,12], 'PaperPosition', [0,0,12,12]);
	print(gcf, [otag, '_ByMu_OverRnds.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_ByMu_OverRnds.png'], '-dpng', '-painters');
end

figure(2),
for(ecnt=1:NumEdges)
	subplot(1,NumEdges,ecnt),
	plot([1,2*length(PerLabs)+1],[X.CorrVals(ecnt,1),X.CorrVals(ecnt,1)], 'b-', 'linewidth', 1);
	hold on
	plot([1,2*length(PerLabs)+1],[X.CorrVals(ecnt,2),X.CorrVals(ecnt,2)], 'r-', 'linewidth', 1);
	hold off
	lt1=['mucosa: \rho=', num2str(round(X.CorrVals(ecnt,1),2)), '; MI=', num2str(round(X.MIVals(ecnt,1),2))];
	lt2=['muscularis: \rho=', num2str(round(X.CorrVals(ecnt,2),2)), '; MI=', num2str(round(X.MIVals(ecnt,2),2))];
	legend({lt1,lt2}, 'location', 'northoutside')
	set(gca, 'ylim', [2,5])
	axis off
end


if(savefigs)
	set(gcf, 'PaperSize', [12,2], 'PaperPosition', [0,0,12,2]);
	print(gcf, [otag, '_ByMu_OverRndsLegend.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_ByMu_OverRndsLegend.png'], '-dpng', '-painters');
end

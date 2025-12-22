function VizSubpopByMuOverTopEdges(vizdatafile, otag, AgValFile);

if(strcmp(otag, ''))
	savefigs=0;
else
	savefigs=1;
end

MethNames={'LIONESS::PCC', 'SSN', 'SWEET', 'BONOBO', 'LIONESS::MI', 'CSN'};
load(vizdatafile); % PerLabs, Gmu1, Gmu2
[~,NumEdges]=size(Gvar1);
X=load(AgValFile);
AgVals=[X.CorrVals; X.MIVals];

% Visualization
pcnt=1;
figure(1),
for(mcnt=1:length(MethNames))
	for(ecnt=1:NumEdges)
		q1=quantile(Gmu1{mcnt,ecnt}',3);
		q2=quantile(Gmu2{mcnt,ecnt}',3);
		subplot(length(MethNames), NumEdges, pcnt),
		errorbar(2*[1:1:length(PerLabs)]-0.25, q1(2,:), q1(2,:)-q1(1,:), q1(3,:)-q1(2,:), 'b.', 'markersize', 10);
		hold on
		errorbar(2*[1:1:length(PerLabs)]+0.25, q2(2,:), q2(2,:)-q2(1,:), q2(3,:)-q2(2,:), 'r.', 'markersize', 10);
		plot([1,2*length(PerLabs)+1],[AgVals(ecnt,1),AgVals(ecnt,1)], 'b--', 'linewidth', 1);
		plot([1,2*length(PerLabs)+1],[AgVals(ecnt,2),AgVals(ecnt,2)], 'r--', 'linewidth', 1);
		hold off
		ylabel('Edge Weight')
		% imagesc([Gmu1{mcnt,ecnt}'; Gmu2{mcnt,ecnt}'], [-2,2])
		% colormap jet;
		set(gca, 'xtick', 2*[1:1:length(PerLabs)], 'xticklabel', PerLabs);
		xtickangle(45)
		title([MethNames{mcnt}, ' - Edge#', num2str(ecnt)])
		pcnt=pcnt+1;
	end
end

if(savefigs)
	set(gcf, 'PaperSize', [10,10], 'PaperPosition', [0,0,10,10]);
	print(gcf, [otag, '_ByMu_OverRnds.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_ByMu_OverRnds.png'], '-dpng', '-painters');
end

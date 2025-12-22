function CheckK2(datafile, otag, enum, K, kcols, showlegend);

if(strcmp(otag, ''))
	savefigs=0;
else
	savefigs=1;
end

% load the data
load(datafile);
sidx=1:2*numsampOrg;

% calculate SWEET networks as a function of K
SSNetsK=GenerateSWEETK(ExpData, K);

Sedge=zeros(length(K),length(sidx));
Ledge=zeros(length(K),length(sidx));
for(dcnt=1:length(K))
	Sedge(dcnt,:)=SSNetsK(enum,sidx,dcnt);
	Ledge(dcnt,:)=SSNets{1}(enum,sidx);
end

[~,ridx]=sort(Ledge(1,:));
plot(Sedge(1,ridx)', Ledge(1,ridx)', '.', 'MarkerSize', 10, 'Color', kcols(1,:))
hold on
for(dcnt=2:length(K))
	plot(Sedge(dcnt,ridx)', Ledge(dcnt,ridx)', '.', 'MarkerSize', 10, 'Color', kcols(dcnt,:));
end
hold off
set(gca, 'FontSize', 18)
xlabel('SWEET', 'FontSize', 24); % predicted edge-weight');
ylabel('LIONESS::PCC', 'FontSize', 24); % predicted edge-weight');

if(showlegend==1)
	lvals=cell(length(K),1);
	for(dcnt=1:length(K))
        	lvals{dcnt}=[num2str(K(dcnt))];
	end
	lvals{3}='0.1 (default)';

	legend(lvals, 'location', 'southeast', 'interpreter', 'latex', 'FontSize', 14)
end

if(savefigs)
        set(gcf, 'PaperSize', [6,5], 'PaperPosition', [0,0,6,5]);
        print(gcf, [otag, '_KRange.pdf'], '-dpdf', '-painters');
        print(gcf, [otag, '_KRange.png'], '-dpng', '-painters');
end

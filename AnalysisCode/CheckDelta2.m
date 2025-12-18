function CheckDelta2(datafile, otag, enum, delta, dcols, showlegend);

if(strcmp(otag, ''))
	savefigs=0;
else
	savefigs=1;
end

% load the data
load(datafile);
[~, numsamp]=size(ExpData);
sidx=1:2*numsampOrg;

% calculate BONOBO networks as a function of delta
SSNetsD=GenerateBONOBOdelta(ExpData, delta);

% pull out values of selected edge
Bedge=zeros(length(delta),length(sidx));
Ledge=zeros(length(delta),length(sidx));
for(dcnt=1:length(delta))
	Bedge(dcnt,:)=SSNetsD(enum,sidx,dcnt);
	Ledge(dcnt,:)=SSNets{1}(enum,sidx);
end
SSNetsD=GenerateBONOBOdelta(ExpData, 1/3);
SSNetsD=SSNetsD(:,sidx); 
BedgeT=SSNetsD(enum,:);
LedgeT=SSNets{1}(enum,:);

% Visualization
[~,ridx]=sort(Ledge(1,:));
plot(Bedge(1,ridx)', Ledge(1,ridx)', '.', 'MarkerSize', 10, 'Color', dcols(1,:));
hold on
for(dcnt=2:length(delta))
	plot(Bedge(dcnt,ridx)', Ledge(dcnt,ridx)', '.', 'MarkerSize', 10, 'Color', dcols(dcnt,:));
end
plot(BedgeT(ridx), LedgeT(ridx), '.', 'MarkerSize', 10, 'Color', [0.5,0.5,0.5]);
plot(SSNets{4}(enum,ridx), SSNets{1}(enum,ridx), 'k.', 'MarkerSize', 10);
hold off
set(gca, 'FontSize', 18)
xlabel('BONOBO', 'FontSize', 24); % predicted edge-weight');
ylabel('LIONESS::PCC', 'FontSize', 24); % predicted edge-weight');
% set(gca, 'xlim', [-1,1]); %, 'ylim', [-8,8]);

if(showlegend==1)
	lvals=cell(length(delta),1);
	for(dcnt=1:length(delta))
        	lvals{dcnt}=[num2str(delta(dcnt))];
	end
	lvals{length(delta)+1}='1/3 (theo. max)';
	lvals{length(delta)+2}='$\delta_q^{tuned}$ (default)';

	legend(lvals, 'location', 'southeast', 'interpreter', 'latex', 'FontSize', 14)
end

if(savefigs)
        set(gcf, 'PaperSize', [6,5], 'PaperPosition', [0,0,6,5]);
        print(gcf, [otag, '_DeltaRange.pdf'], '-dpdf', '-painters');
        print(gcf, [otag, '_DeltaRange.png'], '-dpng', '-painters');
end


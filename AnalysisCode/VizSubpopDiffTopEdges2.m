function VizSubpopDiffTopEdges2(vizdatafile, otag);

if(strcmp(otag, ''))
	savefigs=0;
else
	savefigs=1;
end

MethNames={'LIONESS::PCC', 'SSN', 'SWEET', 'BONOBO', 'LIONESS::MI', 'CSN'};
load(vizdatafile); % PerLabs, NodeIdx, Tstats, Pvals, NumTop1, NumTop2
[NumMeth, NumEdges]=size(Tstats);
[~,NumRnds]=size(Tstats{1,1});

% Visualization
methidx=repmat([1:1.5:1.5*length(MethNames)]',1,length(PerLabs), NumRnds);
peridx=repmat(1:length(PerLabs),length(MethNames),1,NumRnds);
elims={[-5,5], [-5,5], [-15,15], [-25,25]};
xtickvals=[1:1.5:1.5*length(MethNames)];
FalseRate=cell(NumEdges,1);
tcutoffs=[2.626, 3.390, 4.053, 4.654]; % t-statistic cutoffs, t-tailed test, df=100, p=[0.01, 0.001, 0.0001, 0.00001]

for(ecnt=1:NumEdges)
    FalseRate{ecnt}=zeros(NumMeth, length(PerLabs), length(tcutoffs));
    ratdata=zeros(length(MethNames), length(PerLabs), NumRnds);
    for(mcnt=1:length(MethNames))
        ratdata(mcnt,:,:)=Tstats{mcnt,ecnt};
        for(tcnt=1:length(tcutoffs))
            FalseRate{ecnt}(mcnt,:,tcnt)=sum(abs(Tstats{mcnt,ecnt})'>tcutoffs(tcnt));
            if(ecnt==3 && mcnt<5)
                FalseRate{ecnt}(mcnt,:,tcnt)=sum(Tstats{mcnt,ecnt}'<tcutoffs(tcnt)); % for edge-type 3 in the paper and for the PCC-based methods, anything not significantly higher in pop 1 vs 2 is a false negative
            end
        end
    end
    figure(1), subplot(NumEdges+1,1,ecnt),
    bc=boxchart(methidx(:), ratdata(:), 'groupbycolor', peridx(:));
    hold on
    plot([0,1.5*length(MethNames)+1],[0,0], 'k--');
    plot([0,1.5*length(MethNames)+1],[2.625, 2.625], 'k:') % this line assumes p<0.01 as sig cutoff
    plot([0,1.5*length(MethNames)+1],[-2.625, -2.625], 'k:') % this line assumes p<0.01 as sig cutoff
    hold off
    set(gca, 'linewidth', 1, 'FontSize', 14, 'xlim', [min(xtickvals)-0.75,max(xtickvals)+0.75]);
    set(gca, 'xtick', xtickvals, 'xticklabel', MethNames);
    set(gca, 'ylim', elims{ecnt});
    % legend(PerLabs, 'location', 'eastoutside', 'FontSize', 10)
    ylabel('T-statistic')

    col=turbo(length(PerLabs));
    for(pcnt=1:length(PerLabs))
        bc(pcnt).BoxFaceColor=col(pcnt,:);
        bc(pcnt).MarkerColor=col(pcnt,:);
    end

    for(mcnt=1:NumMeth)
        for(percnt=1:length(PerLabs))
            if(FalseRate{ecnt}(mcnt, percnt)>50)
                % text(xtickvals(mcnt)+1.35*(percnt/10-0.5), elims{ecnt}(2)*1.1, 'X', 'HorizontalAlignment','center');
                text(xtickvals(mcnt)+1.35*(percnt/10-0.5), elims{ecnt}(1), 'X', 'HorizontalAlignment','center');
            end
        end
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
set(gca, 'ylim', [22,42])
axis off

if(savefigs)
    set(gcf, 'PaperSize', [10,12], 'PaperPosition', [0,0,10,12]);
    print(gcf, [otag, '_Tstat_OverRnds.pdf'], '-dpdf', '-painters');
    print(gcf, [otag, '_Tstat_OverRnds.png'], '-dpng', '-painters');
end

figure(2),
pcnt=1;
for(ecnt=1:NumEdges)
    for(mcnt=1:NumMeth)
        dat=squeeze(FalseRate{ecnt}(mcnt,:,:))';
        subplot(NumEdges, NumMeth, pcnt), imagesc(dat,[0,100])
        for(percnt=1:length(PerLabs))
            for(tcnt=1:length(tcutoffs))
                if(dat(tcnt,percnt)>50)
                    text(percnt,tcnt,'X','HorizontalAlignment','center');
                end
            end
        end
        title(MethNames{mcnt})
        set(gca, 'xtick', 1:1:length(PerLabs), 'xTickLabel', PerLabs);
        % set(gca, 'ytick', 1:1:length(tcutoffs), 'yticklabel', {'0.01', '0.001', '0.0001', '0.00001'})
        set(gca, 'ytick', 1:1:length(tcutoffs), 'yticklabel', {'10^{-2}', '10^{-3}', '10^{-4}', '10^{-5}'})
        % colorbar
        pcnt=pcnt+1;
    end
end
figure(2), colormap("bone")

if(savefigs)
    set(gcf, 'PaperSize', [16,8], 'PaperPosition', [0,0,16,8]);
    print(gcf, [otag, '_FalseCalls_OverRnds.pdf'], '-dpdf', '-painters');
    print(gcf, [otag, '_FalseCalls_OverRnds.png'], '-dpng', '-painters');
end

figure(3), imagesc(dat,[0,100])
set(gca, 'xtick', 1:1:length(PerLabs), 'xTickLabel', PerLabs);
set(gca, 'ytick', 1:1:length(tcutoffs), 'yticklabel', {'10^{-2}', '10^{-3}', '10^{-4}', '10^{-5}'})
colorbar('FontSize', 24)
figure(3), colormap("bone")

if(savefigs)
    set(gcf, 'PaperSize', [6,6], 'PaperPosition', [0,0,6,6]);
    print(gcf, [otag, '_FalseCalls_OverRnds_colorbar.pdf'], '-dpdf', '-painters');
    print(gcf, [otag, '_FalseCalls_OverRnds_colorbar.png'], '-dpng', '-painters');
end

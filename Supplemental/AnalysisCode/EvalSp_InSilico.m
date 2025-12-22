function EvalSp_InSilico(datafile, PairIdx, otag);

if(strcmp(otag, ''))
        savefigs=0;
else
        savefigs=1;
end

% load the data
load(datafile);
[~, NumSamp]=size(ExpData);

% Calculate Sp for SWEET
x=0.01; % epsilon value in SWEET; note that 0.01 corresponds to SWEET's default value, updating this will only change the local calculation of Sq
PCCmat=corr(ExpData);
muPCC=(sum(PCCmat)-1)/(NumSamp-1); minPCC=min(muPCC); maxPCC=max(muPCC);
Sp=(muPCC-repmat(minPCC, 1, NumSamp)+x)./repmat(maxPCC-minPCC+x, 1, NumSamp); 

% indices for each 'sample-set'
ss1=1:numsampOrg;
ss2a=(numsampOrg+1):(1.5*numsampOrg);
ss2b=(1.5*numsampOrg):(2*numsampOrg);
ss3=(2*numsampOrg+1):NumSamp;
sidx={ss1,ss2a,ss2b,ss3};
scolor=[0.9,0.5,0.5;0.75,0.1,0.9;0.1,0.1,0.9;0.5,0.5,0.5];
slegend={'Sample Set 1', 'Sample Set 2a', 'Sample Set 2b', 'Sample Set 3'};
shortlegend={'Set 1', 'Set 2a', 'Set 2b', 'Set 3'};

% Visualization
figure(1),
imagesc(corr(ExpData),[-1,1]);
set(gca, 'xtick', [0,NumSamp+1], 'ytick', [0,NumSamp+1], 'FontSize', 18)
% title({'Pearson Correlation'; 'between Samples'});
xlabel('Samples \rightarrow', 'FontSize', 24);
ylabel('Samples \rightarrow', 'FontSize', 24);
cmap1=summer(32);
cmap2=spring(32);
cmap12=[cmap1;cmap2(end:-1:1,:)];
colormap(cmap12);
colorbar('FontSize', 16);

if(savefigs)
	set(gcf, 'PaperSize', [6,5], 'PaperPosition', [0,0,6,5]);
	print(gcf, [otag, '_SampleCorr.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_SampleCorr.png'], '-dpng', '-painters');
end


figure(2),
plot(sidx{1}, Sp(sidx{1}), '.', 'MarkerSize', 10, 'color', scolor(1,:));
hold on
for(scnt=2:length(sidx))
	plot(sidx{scnt}, Sp(sidx{scnt}), '.', 'MarkerSize', 10, 'color', scolor(scnt,:));
end
plot([numsampOrg, numsampOrg]+0.05, [0,1], 'k-');
plot(1.5*[numsampOrg, numsampOrg]+0.05, [0,1], '--', 'color', [0.5, 0.5, 0.5]);
plot(2*[numsampOrg, numsampOrg]+0.05, [0,1], 'k-');
hold off
% title('SWEET S_q parameter')
set(gca, 'xtick', [0,NumSamp+1], 'FontSize', 18);
xlabel('Samples \rightarrow', 'FontSize', 24);
% ylabel('S_q (\epsilon = 0.01)', 'FontSize', 24)
colormap(scolor(end:-1:1,:));
colorbar('ytick', 0.125:0.25:0.875, 'yticklabel', shortlegend(end:-1:1), 'FontSize', 16);

if(savefigs)
	set(gcf, 'PaperSize', [6,5], 'PaperPosition', [0,0,6,5]);
	print(gcf, [otag, '_SqValues.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_SqValues.png'], '-dpng', '-painters');
end


figure(3),
plot(SSNets{3}(PairIdx(1),sidx{1}), SSNets{1}(PairIdx(1),sidx{1}), '.', 'MarkerSize', 10, 'color', scolor(1,:));
hold on
for(scnt=2:length(sidx)-1)
	plot(SSNets{3}(PairIdx(1),sidx{scnt}), SSNets{1}(PairIdx(1),sidx{scnt}), '.', 'MarkerSize', 10, 'color', scolor(scnt,:));
end
hold off
set(gca, 'FontSize', 18)
xlabel('SWEET', 'FontSize', 24); % predicted edge-weight');
ylabel('LIONESS::PCC', 'FontSize', 24);% predicted edge-weight');
% legend(slegend(1:length(sidx)-1), 'location', 'southeast', 'FontSize', 12);

if(savefigs)
	set(gcf, 'PaperSize', [6,5], 'PaperPosition', [0,0,6,5]);
	print(gcf, [otag, '_LinearEdge.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_LinearEdge.png'], '-dpng', '-painters');
end


figure(4),
plot(SSNets{3}(PairIdx(2),sidx{1}), SSNets{1}(PairIdx(2),sidx{1}), '.', 'MarkerSize', 10, 'color', scolor(1,:));
hold on
for(scnt=2:length(sidx)-1)
	plot(SSNets{3}(PairIdx(2),sidx{scnt}), SSNets{1}(PairIdx(2),sidx{scnt}), '.', 'MarkerSize', 10, 'color', scolor(scnt,:));
end
hold off
set(gca, 'FontSize', 18)
xlabel('SWEET', 'FontSize', 24); % predicted edge-weight');
ylabel('LIONESS::PCC', 'FontSize', 24); % predicted edge-weight');
legend(slegend(1:length(sidx)-1), 'location', 'southeast', 'FontSize', 16);

if(savefigs)
	set(gcf, 'PaperSize', [6,5], 'PaperPosition', [0,0,6,5]);
	print(gcf, [otag, '_NonlinearEdge.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_NonlinearEdge.png'], '-dpng', '-painters');
end


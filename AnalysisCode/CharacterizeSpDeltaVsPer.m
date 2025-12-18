function CharacterizeSpDeltaVsPer(datafile, otag, selTissue, NumSamp, PerVals, PerLabs, NumRnds);

if(strcmp(otag, ''))
	savefigs=0;
else
	savefigs=1;
end


rng(5);
MixData=load(datafile);
% selTissue={'esophagus_mucosa', 'esophagus_muscularis'};

InSilicoExp1=MixData.ExpData(:,strcmp(MixData.Tissue, selTissue{1}));
InSilicoExp2=MixData.ExpData(:,strcmp(MixData.Tissue, selTissue{2}));
[NumGenes, Num1]=size(InSilicoExp1); Num2=size(InSilicoExp2,2);

ecmap=hsv2rgb([repmat([39/360, 0.9], 64,1), [0.9:-(0.9/63):0]']); % orange to black
ecmap=ecmap(end:-1:1,:);

cmap1=summer(32);
cmap2=spring(32);
cmap12=[cmap1;cmap2(end:-1:1,:)];

for(pcnt=1:length(PerVals))
	PopPer=PerVals(pcnt);
	n1=round(NumSamp*PopPer);
	n2=NumSamp-n1;

	ridx1=randperm(Num1); ridx1=ridx1(1:n1);
	ridx2=randperm(Num2); ridx2=ridx2(1:n2);
	ExpData=[InSilicoExp1(:,ridx1), InSilicoExp2(:,ridx2)];

	[~,eridx]=sort([mean(ExpData(:,1:n1),2)-mean(ExpData(:,(n1+1):end),2)].*[mean(ExpData(:,1:n1),2)+mean(ExpData(:,(n1+1):end),2)], 'descend');

	x=0.01;
	PCCmat=corr(ExpData);
	muPCC=(sum(PCCmat)-1)/(NumSamp-1); minPCC=min(muPCC); maxPCC=max(muPCC);
	Sp=(muPCC-repmat(minPCC, 1, NumSamp)+x)./repmat(maxPCC-minPCC+x, 1, NumSamp);

	delta=zeros(NumSamp,1);
	for(idx=1:NumSamp)
		idxvec=[1:(idx-1),(idx+1):NumSamp];
		cov_diag=var(ExpData(:,idxvec),[],2)';
		delta(idx)=1/(3+2*mean(sqrt(cov_diag))/var(cov_diag)); % BONOBO: delta
	end

	figure(pcnt),
	ax(1)=subplot(7,1,1:3); imagesc(ExpData(eridx,:),[0,15]);
	colormap(ax(1), ecmap);
	% set(gca, 'xtick', [n1/2,n2/2+n1], 'xticklabel', regexprep(selTissue, '_', ' '), 'yticklabel', '', 'ytick', [0,NumGenes+1]);
	set(gca, 'xticklabel','', 'yticklabel', '', 'ytick', [0,NumGenes+1], 'FontSize', 16);
	ylabel('Genes', 'FontSize', 20); xlabel('Samples', 'FontSize', 20);
	colorbar('FontSize', 16)
	% title(['% Samples in Group 1: ', PerLabs{pcnt}], 'FontSize', 20);

	ax(2)=subplot(7,1,4:5); imagesc(corr(ExpData), [0.6,1]);
	colormap(ax(2), cmap12);
	set(gca, 'xticklabel', '', 'yticklabel', '', 'FontSize', 16);
	xlabel('Samples', 'FontSize', 20); ylabel('Samples', 'FontSize', 20);
	colorbar('FontSize', 16);

	ax(3)=subplot(7,1,6); plot(1:NumSamp, Sp, 'k.');
	set(gca, 'ylim', [0,1], 'ytick', [0,0.5,1], 'xticklabel', '', 'FontSize', 16);
	xlabel('Samples', 'FontSize', 20); % ylabel('S_q (\epsilon=0.01)', 'FontSize', 20);
	colorbar('FontSize', 16)

	ax(4)=subplot(7,1,7), plot(1:NumSamp, delta, 'k.');
	set(gca, 'xticklabel', '', 'FontSize', 16);
	xlabel('Samples', 'FontSize', 20); % ylabel('\delta_q^{tuned}', 'FontSize', 20);

	% y1=ylim(gca);
	% set(gca, 'ylim', [0.23,0.33], 'ytick', [0.23, 0.28, 0.33], 'xticklabel', '');
	% text(0.5*NumSamp, y1(1)+0.15*(y1(2)-y1(1)), ['<\delta_q^{tuned}>=', num2str(mean(delta))], 'color', 'red', 'horizontalalignment', 'center');
	% title(['<\delta_q^{tuned}>=', num2str(mean(delta))]);
	colorbar('FontSize', 16)

	if(savefigs)
		set(gcf, 'PaperSize', [4,12], 'PaperPosition', [0,0,4,12]);
		print(gcf, [otag, '_PopPer', num2str(pcnt), '.pdf'], '-dpdf', '-painters');
		print(gcf, [otag, '_PopPer', num2str(pcnt), '.png'], '-dpng', '-painters');
	end
end


muSp=zeros(length(PerVals), NumRnds);
muD=zeros(length(PerVals), NumRnds);
muSp1=zeros(length(PerVals), NumRnds);
muD1=zeros(length(PerVals), NumRnds);
muSp2=zeros(length(PerVals), NumRnds);
muD2=zeros(length(PerVals), NumRnds);

for(pcnt=1:length(PerVals))
	PopPer=PerVals(pcnt);
	n1=round(NumSamp*PopPer);
	n2=NumSamp-n1;

	for(rcnt=1:NumRnds)
		ridx1=randperm(Num1); ridx1=ridx1(1:n1);
		ridx2=randperm(Num2); ridx2=ridx2(1:n2);
		ExpData=[InSilicoExp1(:,ridx1), InSilicoExp2(:,ridx2)];

		x=0.01;
		PCCmat=corr(ExpData);
		muPCC=(sum(PCCmat)-1)/(NumSamp-1); minPCC=min(muPCC); maxPCC=max(muPCC);
		Sp=(muPCC-repmat(minPCC, 1, NumSamp)+x)./repmat(maxPCC-minPCC+x, 1, NumSamp);

		delta=zeros(NumSamp,1);
		for(idx=1:NumSamp)
			idxvec=[1:(idx-1),(idx+1):NumSamp];
			cov_diag=var(ExpData(:,idxvec),[],2)';
			delta(idx)=1/(3+2*mean(sqrt(cov_diag))/var(cov_diag)); % BONOBO: delta
		end

		muSp(pcnt,rcnt)=mean(Sp);
		muD(pcnt,rcnt)=mean(delta);

		muSp1(pcnt,rcnt)=mean(Sp(1:n1));
		muSp2(pcnt,rcnt)=mean(Sp((n1+1):end));

		muD1(pcnt,rcnt)=mean(delta(1:n1));
		muD2(pcnt,rcnt)=mean(delta((n1+1):end));
	end
	disp(pcnt);
end

AvgSp=mean(muSp,2); StdSp=std(muSp,[],2);
AvgSp1=mean(muSp1,2); StdSp1=std(muSp1,[],2);
AvgSp2=mean(muSp2,2); StdSp2=std(muSp2,[],2);
figure(length(PerVals)+1), errorbar(PerVals, AvgSp, StdSp, StdSp, 'k.-', 'linewidth', 1.5, 'markersize', 15);
hold on
errorbar(PerVals, AvgSp1, StdSp1, StdSp1, 'b.-', 'linewidth', 1.5, 'markersize', 15);
errorbar(PerVals, AvgSp2, StdSp2, StdSp2, 'r.-', 'linewidth', 1.5, 'markersize', 15);
hold off
set(gca, 'xtick', PerVals, 'xticklabel', PerLabs, 'ylim', [0,1], 'FontSize', 16)
xlabel('% Samples in Group 1', 'FontSize', 20);
ylabel('<S_q>', 'FontSize', 20);
if(savefigs)
	set(gcf, 'PaperSize', [15,3], 'PaperPosition', [0,0,15,3]);
	print(gcf, [otag, '_AvgSpAcrossPopPer.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_AvgSpAcrossPopPer.png'], '-dpng', '-painters');
end

AvgD=mean(muD,2); StdD=std(muD,[],2);
AvgD1=mean(muD1,2); StdD1=std(muD1,[],2);
AvgD2=mean(muD2,2); StdD2=std(muD2,[],2);
figure(length(PerVals)+2), errorbar(PerVals, AvgD, StdD, StdD, 'k.-', 'linewidth', 1.5, 'markersize', 15);
hold on
errorbar(PerVals, AvgD1, StdD1, StdD1, 'b.-', 'linewidth', 1.5, 'markersize', 15);
errorbar(PerVals, AvgD2, StdD2, StdD2, 'r.-', 'linewidth', 1.5, 'markersize', 15);
hold off
set(gca, 'xtick', PerVals, 'xticklabel', PerLabs, 'ylim', [0.18,0.32], 'FontSize', 16)
ylabel('<\delta_q^{tuned}>', 'FontSize', 20);
xlabel('% Samples in Group 1', 'FontSize', 20);
legend(['All Samples', regexprep(selTissue, '_', ' ')], 'location', 'south');
if(savefigs)
	set(gcf, 'PaperSize', [15,3], 'PaperPosition', [0,0,15,3]);
	print(gcf, [otag, '_AvgDeltaAcrossPopPer.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_AvgDeltaAcrossPopPer.png'], '-dpng', '-painters');
end

function EvalSp_HumanData(datafile, NumRGenes, seltype, otag);

if(strcmp(otag, ''))
        savefigs=0;
else
        savefigs=1;
end

rng(1);
% read in the data
MixData=load(datafile);

% select genes from original data to use for the analysis
sigma=std(MixData.ExpData,[],2);
if(strcmp(seltype, 'rand'))
    gidx=randperm(length(MixData.GeneSymbols)); gidx=gidx(1:NumRGenes); % random genes
elseif(strcmp(seltype, 'high'))
    [~,gidx]=sort(sigma, 'descend'); gidx=gidx(1:NumRGenes); % most variable genes
elseif(strcmp(seltype, 'low'))
    [~,gidx]=sort(sigma, 'ascend'); gidx=gidx(1:NumRGenes); % least variable genes
else
    disp('sel variable not correctly defined, defaulting to selecting random genes');
    gidx=randperm(length(MixData.GeneSymbols)); gidx=gidx(1:NumRGenes); % random genes
end
ExpData=MixData.ExpData(gidx,:);

[~, NumSamp]=size(ExpData);

% Calculate Sp for SWEET
x=0.01; % epsilon value in SWEET; note that 0.01 corresponds to SWEET's default value, updating this will only change the local calculation of Sq
PCCmat=corr(ExpData);
muPCC=(sum(PCCmat)-1)/(NumSamp-1); minPCC=min(muPCC); maxPCC=max(muPCC);
Sp=(muPCC-repmat(minPCC, 1, NumSamp)+x)./repmat(maxPCC-minPCC+x, 1, NumSamp); 

% Calculate delta_tuned for BONOBO
delta=zeros(NumSamp,1);
for(idx=1:NumSamp)
    idxvec=[1:(idx-1),(idx+1):NumSamp];
    cov_diag=var(ExpData(:,idxvec),[],2)';
    delta(idx)=1/(3+2*mean(sqrt(cov_diag))/var(cov_diag)); % BONOBO: delta
end

disp(['Mean S_q: ', num2str(mean(Sp))]);
disp(['Mean \delta_q: ', num2str(mean(delta))]);

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
plot(1:NumSamp, Sp, 'k.', 'MarkerSize', 10);
% title('SWEET S_q parameter')
set(gca, 'xtick', [0,NumSamp+1], 'FontSize', 18);
xlabel('Samples \rightarrow', 'FontSize', 24);
% ylabel('S_q (\epsilon = 0.01)', 'FontSize', 24)

if(savefigs)
	set(gcf, 'PaperSize', [6,5], 'PaperPosition', [0,0,6,5]);
	print(gcf, [otag, '_SqValues.pdf'], '-dpdf', '-painters');
	print(gcf, [otag, '_SqValues.png'], '-dpng', '-painters');
end

figure(3),
plot(1:NumSamp, delta, 'k.', 'MarkerSize', 10);
% title('SWEET S_q parameter')
set(gca, 'xtick', [0,NumSamp+1], 'FontSize', 18);
xlabel('Samples \rightarrow', 'FontSize', 24);
% ylabel('delta_q^{(tuned)}', 'FontSize', 24)

if(savefigs)
    set(gcf, 'PaperSize', [6,5], 'PaperPosition', [0,0,6,5]);
    print(gcf, [otag, '_SqValues.pdf'], '-dpdf', '-painters');
    print(gcf, [otag, '_SqValues.png'], '-dpng', '-painters');
end

figure(4),
plot(Sp, delta, 'k.', 'MarkerSize', 10);
% title('SWEET S_q parameter')
ylabel('\delta_q^{(tuned)}', 'FontSize', 24)
xlabel('S_q', 'FontSize', 24)

if(savefigs)
    set(gcf, 'PaperSize', [6,5], 'PaperPosition', [0,0,6,5]);
    print(gcf, [otag, '_SqValues.pdf'], '-dpdf', '-painters');
    print(gcf, [otag, '_SqValues.png'], '-dpng', '-painters');
end

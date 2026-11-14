addpath('../method_code/')

numrnds=10;
MethNames={'LIONESS', 'SSN', 'SWEET', 'BONOBO'};
rng(1);

% vary number of genes, hold number of samples fixed
NumG=[50,100,200,400,800,1600]; NumSamples=100;
runT1=zeros(length(MethNames), length(NumG), numrnds);
for(gcnt=1:length(NumG))
    for(rndcnt=1:numrnds)
        ExpData=normrnd(0,1, NumG(gcnt), NumSamples);
        ExpData=ExpData-min(ExpData(:))+eps;
        tic 
        [~]=GenerateLIONESS(ExpData);
        runT1(1,gcnt,rndcnt)=toc;
        tic 
        [~]=GenerateSSN(ExpData);
        runT1(2,gcnt,rndcnt)=toc;
        tic 
        [~]=GenerateSWEET(ExpData);
        runT1(3,gcnt,rndcnt)=toc;
        tic 
        [~]=GenerateBONOBO(ExpData);
        runT1(4,gcnt,rndcnt)=toc;
    end
    disp(['Evaluating data with ', num2str(NumSamples), ' samples and ', num2str(NumG(gcnt)), ' genes.']);
end

mu=mean(runT1,3);
sigma=std(runT1,[],3);
figure(1), errorbar(NumG,mu,sigma,sigma,'-', 'linewidth',2)
set(gca, 'box', 'off', 'LineWidth', 2, 'FontSize', 14)
title(['Number of Samples = ', num2str(NumSamples)])
ylabel('Runtime (s)');
xlabel('Number of Genes')
legend(MethNames, 'FontSize', 12, 'location', 'northwest')
figure(1), set(gcf, 'theme', 'light')


% vary number of samples, hold number of genes fixed
NumS=[50,100,200,400,800,1600]; NumGenes=200;
runT2=zeros(length(MethNames), length(NumS), numrnds);
for(scnt=1:length(NumS))
    for(rndcnt=1:numrnds)
        ExpData=normrnd(0,1, NumGenes, NumS(scnt));
        ExpData=ExpData-min(ExpData(:))+eps;
        tic 
        [~]=GenerateLIONESS(ExpData);
        runT2(1,scnt,rndcnt)=toc;
        tic 
        [~]=GenerateSSN(ExpData);
        runT2(2,scnt,rndcnt)=toc;
        tic 
        [~]=GenerateSWEET(ExpData);
        runT2(3,scnt,rndcnt)=toc;
        tic 
        [~]=GenerateBONOBO(ExpData);
        runT2(4,scnt,rndcnt)=toc;
    end
    disp(['Evaluating data with ', num2str(NumS(scnt)), ' samples and ', num2str(NumGenes), ' genes.']);
end

mu=mean(runT2,3);
sigma=std(runT2,[],3);
figure(2), errorbar(NumS,mu,sigma,sigma,'-', 'linewidth',2)
set(gca, 'box', 'off', 'LineWidth', 2, 'FontSize', 14)
title(['N = ', num2str(NumGenes)])
ylabel('Runtime (s)');
xlabel('S')
legend(MethNames, 'FontSize', 12, 'location', 'northwest')
figure(2), theme('light')
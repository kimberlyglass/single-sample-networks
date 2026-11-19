addpath('AnalysisCode/')
addpath('NetworkCode/')

%%% IN SILICO DATA ANALYSIS %%%

% Create In Silico Expression Data and Corresponding Single-Sample Networks
numsampOrg=250; numsampRnd=100; numGenes=6;
datafile='./Data/InSilicoData.mat';
% MakeInSilicoData(numGenes,numsampOrg,numsampRnd,datafile);

% Visualize In Silico Data (Figure 2, Supplemental Figures 1-3)
otag='Figures/Figure1/VizPatterns'; % '' triggers the program to not save the figures
NodePairs={[1,4],[1,2]}; % two edges to visualize: 1-4 is nonlinear (4th edge calculated by functions), 1-2 is linear (1st edge calculated by functions)
VizSSNPatterns4(datafile, otag, NodePairs);
close all;
SingleSampleNetDists(datafile, NodePairs, otag);
close all;

% Evaluate results when an "add-one-in" approach is taken (Supplemental File)
netfile='./Data/InSilicoDataBG.mat';
CalcBGNets(datafile, netfile);
otag='./Figures/Figure1/VizPatternsBG'; % '' triggers the program to not save the figures
for(bgcnt=1:3)
	VizSSNPatternsBG(netfile, otag, NodePairs, bgcnt);
	close all;
end

% Evaluate BONOBO's delta parameter (Figure 3A-B)
delta=[0.01, 0.05, 0.1, 0.25, 0.5]; dcols=[228,26,28; 55,126,184; 77,175,74; 152,78,163; 255,127,0]/255;
enum=4; otag='Figures/Figure2/VizXPattern2'; showlegend=1;
CheckDelta2(datafile, otag, enum, delta, dcols, showlegend);
close all;
enum=1; otag='Figures/Figure2/VizLPattern2'; showlegend=0;
CheckDelta2(datafile, otag, enum, delta, dcols, showlegend);
close all;

% Evaluate SWEET's K parameter (Figure 3C-D)
K=[0.01, 0.05, 0.1, 0.25, 0.50]; kcols=[228,26,28; 55,126,184; 0,0,0; 152,78,163; 255,127,0]/255;
enum=4; otag='Figures/Figure2/VizXPattern2'; showlegend=1;
CheckK2(datafile, otag, enum, K, kcols, showlegend);
close all;
enum=1; otag='Figures/Figure2/VizLPattern2'; showlegend=0;
CheckK2(datafile, otag, enum, K, kcols, showlegend);
close all;

% Evalute SWEET's Sq parameter (Figure 3E-H)
otag='Figures/Figure2/VizSqPatterns';
PairIdx=[1,4];
EvalSp_InSilico(datafile, PairIdx, otag);
close all;



%%% GTEX / HUMAN DATA ANALYSIS %%%

% relevant GTEx data files used for the analysis published in [https://www.sciencedirect.com/science/article/pii/S2211124717314183] are copied to Data/GTExData/
MinSamp=250; minMedval=2; pcfilter=1;
datafile='./Data/GTEx_MixData3.mat';
% MakeGTEx_MixData3(MinSamp, minMedval, pcfilter, datafile);

% tissues and variables for imbalanced population analysis
selTissue={'esophagus_mucosa', 'esophagus_muscularis'}; % two tissues being investigated
NumSamp=100; % total number of samples selected to make an input dataset
NumRnds=100;

% Figure 4 & Supplemental Figure 4C
% show how delta and Sq vary based on percentage of samples from two populations / tissues
PerVals=[0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.95, 1];
PerLabs={'0%', '5%', '10%', '20%', '35%', '50%', '65%', '80%', '90%', '95%', '100%'};
otag='./Figures/Figure3/CharSpDelta';
CharacterizeSpDeltaVsPer(datafile, otag, selTissue, NumSamp, PerVals, PerLabs, NumRnds);
close all;

% Figure 5A & Supplemental Figure 4A-B
% for two tissues, compare tissue-specific correlation and select four 'edges' to systematically evaluate (Supplemental Figure 4)
otag='Figures/Figure3/PopulationVar'; % output tag for saving the figures
numRedges=100000; % number of edges to visualize (visualizing all the edges is too many, so we subset)
NodeIdx=SelectTopEdges2(datafile, otag, selTissue, numRedges);
close all;

% Figure 5B & Supplemental Figures 5-7
% for each selected edge, determine the mean, variance, and statistical difference in the predicted weight across samples in tissue 1 versus tissue 2, systematically varying the ratio of the tissues
NumRnds=100; % number of times to do this for each percentage value
PerVals=[0.05, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.95]; % percentage of samples selected from tissue 1
PerLabs={'5%', '10%', '20%', '35%', '50%', '65%', '80%', '90%', '95%'};
vizdatafile='./Data/GTExTopEdgeData.mat';
EvalSubpopTopEdges(datafile, selTissue, NumSamp, PerVals, PerLabs, NodeIdx, NumRnds, vizdatafile);

% visualize the mean (Supplemental Figure 5-6)
otag='Figures/Figure3/PopDiff'; % output tag for saving the figures
AgValFile='./Data/AgEdgeVals.mat';
AgValsOfTopEdges(datafile, selTissue, NodeIdx, AgValFile); % separately calculated the overall network values
VizSubpopByMuOverTopEdgesRnds(vizdatafile, otag, AgValFile);
VizSubpopByMuOverTopEdges2(vizdatafile, otag, AgValFile);
close all;

% Visualize the statistical difference (Figure 5B, Supplemental Figure 7)
otag='Figures/Figure3/PopDiff'; % output tag for saving the figures
VizSubpopDiffTopEdges2(vizdatafile, otag);
close all;

% Figure 6A-B
NumTissues=[2,4,8,15]; NumRnds=100; NumSel=100;
TissueLegend={'2 tissues', '4 tissues', '8 tissues', '15 tissues'};
RandIdxFile='./Data/GTEx_RandIdx.mat';
SubsetDataForMixtureAnalysis(datafile, NumTissues, NumRnds, NumSel, RandIdxFile);

NumRGenes=500;
dtag='./Data/MixCorrVals';
seltype={'low', 'rand', 'high'};
for(scnt=1:length(seltype))
 	EvalMethodsMixture2(datafile, NumRGenes, seltype{scnt}, RandIdxFile, dtag);
end

% visualize the results
otag='./Figures/Figure4/MixEvalBG'; 
VizDataMixture2(datafile, NumRGenes, seltype, otag); % Figure 6A, Supplemental Figure 8
close all;
VizMethodsMixture2(dtag, seltype, RandIdxFile, TissueLegend, otag); % Figure 6B
close all;

% evaluate Sp and dq values assuming all samples: GTEx data
for(scnt=1:length(seltype))
    disp([seltype{scnt}, ' variance genes:'])
    EvalSp_HumanData(datafile, NumRGenes, seltype{scnt}, '');
    close all;
end

% summarize results for Table 1
MethodsMixture2_Summary(dtag, seltype([1,3]), RandIdxFile);

% visualize this analysis pipeline for NumTissues=2 (1 randomization) and for seltype='low' (low variance genes)
% Supplemental Figure 9
VizIdxFile='./Data/GTEx_VizIdx.mat';
SubsetDataForMixtureAnalysis(datafile, 2, 1, NumSel, VizIdxFile);
otag=['./Figures/Figure4/VizBenchmark_', seltype{1}];
VizBenchmarkData(datafile, NumRGenes, seltype{1}, VizIdxFile, otag);
close all;

%{
% for non-linear methods / Supplemental Figure 10
% note this takes a *VERY* long time, so it is commented out
NumRGenes=50;
dtag='./Data/MixNLVals';
seltype={'low', 'rand', 'high'};
for(scnt=1:length(seltype))
	EvalMethodsMixtureNL(datafile, NumRGenes, seltype{scnt}, RandIdxFile, dtag);
end

% visualize the results
otag='./Figures/Figure4/MixEvalNL'; % otag='';
VizMethodsMixtureNL(dtag, seltype, RandIdxFile, TissueLegend, otag);
close all;
%}


%%% TCGA ANALYSIS %%%
% relevant TCGA data files used for the analysis published in [https://www.biorxiv.org/content/10.1101/2025.10.03.680219v1] are copied to Data/TGCAData/
MinSamp=80; minMedval=2;
datafile='./Data/TCGA_MixData.mat';
% MakeTCGA_MixData(MinSamp, minMedval, datafile);

% TCGA Evaluation
NumTissues=[4]; NumRnds=100; NumSel=100;
TissueLegend={'4 tissues'};
RandIdxFile='./Data/TGCA_RandIdx.mat';
SubsetDataForMixtureAnalysis(datafile, NumTissues, NumRnds, NumSel, RandIdxFile);

NumRGenes=500;
dtag='./Data/TCGAMixCorrVals';
seltype={'low', 'rand', 'high'};
for(scnt=1:length(seltype))
    EvalMethodsMixture2(datafile, NumRGenes, seltype{scnt}, RandIdxFile, dtag);
end

% visualize the results; Figure 6C, Supplemental Figure 11
otag='Figures/Figure4/TCGAMixEval';
VizMethodsMixtureSub3(dtag, seltype, RandIdxFile, TissueLegend, otag, datafile);
close all;

% evaluate Sp and dq values assuming all samples: TCGA data
for(scnt=1:length(seltype))
    disp([seltype{scnt}, ' variance genes:'])
    EvalSp_HumanData(datafile, NumRGenes, seltype{scnt}, '');
    close all;
end

% summarize results for Table 1
MethodsMixtureSub_Summary(dtag, seltype([2]), RandIdxFile, datafile);


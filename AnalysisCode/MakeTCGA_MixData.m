function MakeTCGA_MixData(MinSamp,minMedval,datafile);

gtex_data_path='./Data/TCGAData/'; % place where the original GTEx data files are stored.
rng(1);

disp('Reading in the data');
% read in sample info
fid=fopen([gtex_data_path, 'TCGA_BRCA_sample_subtypes.tsv'], 'r');
S=textscan(fid, '%s%s', 'delimiter', '\t', 'headerlines', 1);
fclose(fid);
SampleID=S{1};
SampleSubtype=S{2};
clear S;
[SubtypeNames,~,Sloc]=unique(SampleSubtype);
NumSubtypes=length(SubtypeNames);

% read in expression data
fid=fopen([gtex_data_path, 'TCGA_BRCA_normalized_expression.tsv'], 'r');
SampleIDs=fgetl(fid);
SampleIDs=strsplit(SampleIDs, '\t');
frewind(fid);
ExpData=textscan(fid, ['%s', repmat('%f', 1, length(SampleIDs))], 'delimiter', '\t', 'headerlines', 1);
fclose(fid);
GeneSymbols=ExpData{1};
ExpData=cat(2, ExpData{2:end});

% keep only the subtypes with at least the minimum number of associated samples
[uSubtype,~,Sloc]=unique(SampleSubtype);
nums=hist(Sloc,1:1:length(uSubtype));
snums=sort(nums, 'descend');
LocSubtypeNames=uSubtype(nums>=MinSamp); 
sfilter=ismember(SampleSubtype, LocSubtypeNames);
% filter ExpData
ExpData=ExpData(:,sfilter);
SampleSubtype=SampleSubtype(sfilter);

%{
% select a random MinSamp number of samples from the remaining tissues
[uSubtype,~,Sloc]=unique(SampleSubtype);
SubExp=cell(1, length(uSubtype));
SubIdx=cell(1, length(uSubtype));
gMeds=zeros(length(GeneSymbols), length(uSubtype));
for(scnt=1:length(uSubtype))
	sidx=find(scnt==Sloc);
	ridx=randperm(length(sidx));
	sidx=sidx(ridx(1:MinSamp));
	SubExp{scnt}=ExpData(:,sidx);
	SubIdx{scnt}=sidx';
	gMeds(:,scnt)=median(SubExp{scnt},2);
end
ExpData=cat(2, SubExp{:});
SampleSubtype=SampleSubtype(cat(2, SubIdx{:}));
%}

% Determine the median expression across each subtype
[uSubtype,~,Sloc]=unique(SampleSubtype);
gMeds=zeros(length(GeneSymbols), length(uSubtype));
for(scnt=1:length(uSubtype))
    gMeds(:,scnt)=median(ExpData(:,Sloc==scnt),2);
end

% remove genes that don't have some minimal expression across all tissues
minMed=min(gMeds,[],2);
gfilter=minMed>minMedval;
ExpData=ExpData(gfilter,:);
GeneSymbols=GeneSymbols(gfilter);
Tissue=SampleSubtype; % this makes it compatible with code written for GTEx

save(datafile, 'ExpData', 'GeneSymbols', 'Tissue');

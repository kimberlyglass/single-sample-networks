function MakeGTEx_MixData3(MinSamp,minMedval,pcfilter,datafile);

gtex_data_path='./Data/GTExData/'; % place where the original GTEx data files are stored.
rng(1);

disp('Reading in the data');
% read in sample info
fid=fopen([gtex_data_path, 'GTEx_SampleInfo.txt'], 'r');
headings=fgetl(fid); headings=strsplit(headings, '\t');
S=textscan(fid, repmat('%s', 1, 98), 'delimiter', '\t', 'headerlines', 1);
fclose(fid);
SampleID=S{2};
Tissue=S{98};
clear S;
[TissueNames,~,Tloc]=unique(Tissue);
NumTissues=length(TissueNames);

% read in expression data
fid=fopen([gtex_data_path, 'GTEx_NormExpressionData.txt'], 'r');
SampleIDs=fgetl(fid);
SampleIDs=strsplit(SampleIDs, '\t');
SampleIDs=SampleIDs(6:end);
frewind(fid);
ExpData=textscan(fid, ['%s%s%u%u%s', repmat('%f', 1, length(SampleIDs))], 'delimiter', '\t', 'headerlines', 1);
fclose(fid);
GeneNames=ExpData{1}; GeneSymbols=ExpData{2};
ExpData=cat(2, ExpData{6:end});

if(pcfilter)
	% remove 'dots' from ENSG
	delloc=strfind(GeneNames, '.');
	for(gcnt=1:length(GeneNames))
		GeneNames{gcnt}=GeneNames{gcnt}(1:delloc{gcnt}-1);
	end

	% keep only genes that are annotated as protein_coding
	[ENSG,~,BIOTYPE]=textread([gtex_data_path, 'GTEx_AnnInfo_AllGenes.txt'], '%s%s%s', 'headerlines', 1);
	gfilter=ismember(GeneNames,ENSG(strcmp(BIOTYPE,'protein_coding')));
	GeneNames=GeneNames(gfilter); GeneSymbols=GeneSymbols(gfilter);
	ExpData=ExpData(gfilter,:);
end

% keep only the tissues in GTEx with at least the minimum number of associated samples
[uTissue,~,loc]=unique(Tissue);
nums=hist(loc,1:1:length(uTissue));
snums=sort(nums, 'descend');
LocTissueNames=uTissue(nums>=MinSamp); 
tfilter=ismember(Tissue, LocTissueNames);
% filter ExpData
ExpData=ExpData(:,tfilter);
Tissue=Tissue(tfilter);

% select a random MinSamp number of samples from the remaining tissues
[uTissue,~,loc]=unique(Tissue);
SubExp=cell(1, length(uTissue));
SubIdx=cell(1, length(uTissue));
gMeds=zeros(length(GeneNames), length(uTissue));
for(tcnt=1:length(uTissue))
	tidx=find(tcnt==loc);
	ridx=randperm(length(tidx));
	tidx=tidx(ridx(1:MinSamp));
	SubExp{tcnt}=ExpData(:,tidx);
	SubIdx{tcnt}=tidx';
	gMeds(:,tcnt)=median(SubExp{tcnt},2);
end
ExpData=cat(2, SubExp{:});
Tissue=Tissue(cat(2, SubIdx{:}));

% remove genes that don't have some minimal expression across all tissues
minMed=min(gMeds,[],2);
gfilter=minMed>minMedval;
ExpData=ExpData(gfilter,:);
GeneSymbols=GeneSymbols(gfilter);
GeneNames=GeneNames(gfilter);

save(datafile, 'ExpData', 'GeneSymbols', 'Tissue');

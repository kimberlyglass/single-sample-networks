function EvalSubpopByVarOverTopEdges(datafile, selTissue, NumSamp, PerVals, PerLabs, NodeIdx, NumRnds, vizdatafile);

MethNames={'LIONESS::PCC', 'SSN', 'SWEET', 'BONOBO', 'LIONESS::MI', 'CSN'};

rng(1);
MixData=load(datafile);

InSilicoExp1=MixData.ExpData(:,strcmp(MixData.Tissue, selTissue{1}));
InSilicoExp2=MixData.ExpData(:,strcmp(MixData.Tissue, selTissue{2}));
[NumGenes, Num1]=size(InSilicoExp1); Num2=size(InSilicoExp2,2);

NumEdges=size(NodeIdx,1);

Gvar1=cell(length(MethNames),NumEdges);
Gvar2=cell(length(MethNames),NumEdges);
Gmu1=cell(length(MethNames),NumEdges);
Gmu2=cell(length(MethNames),NumEdges);
for(mcnt=1:length(MethNames))
	Gvar1{mcnt}=zeros(length(PerVals),NumRnds);
	Gvar2{mcnt}=zeros(length(PerVals),NumRnds);
	Gmu1{mcnt}=zeros(length(PerVals),NumRnds);
	Gmu2{mcnt}=zeros(length(PerVals),NumRnds);
end

for(pcnt=1:length(PerVals))
	disp([pcnt]);

	PopPer=PerVals(pcnt);
	n1=round(NumSamp*PopPer);
	n2=NumSamp-n1;

	for(rcnt=1:NumRnds)
		ridx1=randperm(Num1); ridx1=ridx1(1:n1);
		ridx2=randperm(Num2); ridx2=ridx2(1:n2);
		ExpData=[InSilicoExp1(:,ridx1), InSilicoExp2(:,ridx2)];

		for(ecnt=1:NumEdges)
			SSNets=GenerateSSCorrByEdge(ExpData,NodeIdx(ecnt,:)); % only generate SSCorr for this edge
			[SSNetsNL,~]=GenerateSSNL(ExpData(NodeIdx(ecnt,:),:));
			SSNets{5}=SSNetsNL{1};
			SSNets{6}=SSNetsNL{2};

			for(mcnt=1:length(MethNames))
				Gvar1{mcnt,ecnt}(pcnt,rcnt)=var(SSNets{mcnt}(1,1:n1),[],2); % variability of edge in top row (pop 1)
				Gvar2{mcnt,ecnt}(pcnt,rcnt)=var(SSNets{mcnt}(1,(n1+1):end),[],2); % variability of edge in top row (pop 2)
				Gmu1{mcnt,ecnt}(pcnt,rcnt)=mean(SSNets{mcnt}(1,1:n1),2); % mean of edge in top row (pop 1)
				Gmu2{mcnt,ecnt}(pcnt,rcnt)=mean(SSNets{mcnt}(1,(n1+1):end),2); % mean of edges in top row (pop 2)
			end
		end
	end
end

save(vizdatafile, 'Gvar1', 'Gvar2', 'Gmu1', 'Gmu2', 'PerLabs', 'NodeIdx');


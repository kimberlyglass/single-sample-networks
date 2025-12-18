function EvalSubpopByVarOverTopEdges_LWPCC(datafile, selTissue, NumSamp, PerVals, PerLabs, NodeIdx, NumRnds, vizdatafile);

MethNames={'LIONESS::WPCC'};

rng(1);
MixData=load(datafile);

InSilicoExp1=MixData.ExpData(:,strcmp(MixData.Tissue, selTissue{1}));
InSilicoExp2=MixData.ExpData(:,strcmp(MixData.Tissue, selTissue{2}));
[NumGenes, Num1]=size(InSilicoExp1); Num2=size(InSilicoExp2,2);

NumEdges=size(NodeIdx,1);

Gvar1=cell(1,NumEdges);
Gvar2=cell(1,NumEdges);
Gmu1=cell(1,NumEdges);
Gmu2=cell(1,NumEdges);
Gvar1{1}=zeros(length(PerVals),NumRnds);
Gvar2{1}=zeros(length(PerVals),NumRnds);
Gmu1{1}=zeros(length(PerVals),NumRnds);
Gmu2{1}=zeros(length(PerVals),NumRnds);

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
			w=[repmat(1/n1, 1, n1), repmat(1/n2, 1, n2)]/2;
			SSNets=GenerateLIONESS_WPCC(ExpData(NodeIdx(ecnt,:),:),w);

			Gvar1{1,ecnt}(pcnt,rcnt)=var(SSNets(1,1:n1),[],2); % variability of edge in top row (pop 1)
			Gvar2{1,ecnt}(pcnt,rcnt)=var(SSNets(1,(n1+1):end),[],2); % variability of edge in top row (pop 2)
			Gmu1{1,ecnt}(pcnt,rcnt)=mean(SSNets(1,1:n1),2); % mean of edge in top row (pop 1)
			Gmu2{1,ecnt}(pcnt,rcnt)=mean(SSNets(1,(n1+1):end),2); % mean of edges in top row (pop 2)
		end
	end
end

save(vizdatafile, 'Gvar1', 'Gvar2', 'Gmu1', 'Gmu2', 'PerLabs', 'NodeIdx');


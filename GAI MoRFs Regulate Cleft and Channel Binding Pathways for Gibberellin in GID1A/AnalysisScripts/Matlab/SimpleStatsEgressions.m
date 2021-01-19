
%load the data as a cell matrix. Each row is addative in counting raw
% # of times the residue of that column interacts with GAx
GAActivelist=load('./GA3minus/HbondTableEgressGA3.mat');
GAInactivelist=load('./GID1A/HbondTableEgressGID1A.mat');
name1="GA3";
name2="GID1A";
% SO this turns each row into raw samples
% This means each row now represents a single egression run instead of the
% addative count I use to plot histograms

GAInactiveResContactNumList=cell2mat(GAInactivelist.ResContactNumList);
GAActiveResContactNumList=cell2mat(GAActivelist.ResContactNumList);

DataList={GAActiveResContactNumList,GAInactiveResContactNumList};
newNumCount=[];
for k=1:length(DataList)
    RowLen=size(DataList{k},1);
    for i=2:size(DataList{k},1)
%         newNumCount=(DataList{k}(i:RowLen,:)-DataList{k}(i-1,:));
        newNumCount=bsxfun(@minus,DataList{k}(i:RowLen,:),DataList{k}(i-1,:));
        DataList{k}(i:RowLen,:)=newNumCount;
    end
end

GAInactiveresidues={GAInactivelist.ResContactNameList{end,:}};
GAActiveresidues={GAActivelist.ResContactNameList{end,:}};

% mutualdiff=setdiff({GA20list.ResContactNameList{end,:}},{GA4list.ResContactNameList{end,:}});
% returns the values in A that are not in B with no repetitions
GAInactivediff=setdiff(GAInactiveresidues,GAActiveresidues);
GAActivediff=setdiff(GAActiveresidues,GAInactiveresidues);

commonResidues=intersect(GAActiveresidues,GAInactiveresidues);

resultTable=cell(length(commonResidues),6);%table('VariableNames',['Name','KSp','Tp','Wp','AD1p','AD2p']);

for j=1:length(commonResidues)
    currentRes=commonResidues{j};
    for iA=1:length(GAActiveresidues)
        if isequal(currentRes, GAActiveresidues{iA})
            index1=iA;
        for iB=1:length(GAInactiveresidues)
            if isequal(currentRes, GAInactiveresidues{iB})
            index2=iB;
            [KSresults,KSpvalues]=kstest2(DataList{1}(:,index1),DataList{2}(:,index2));
            [Tresults, Tpvalues]=ttest2(DataList{1}(:,index1),DataList{2}(:,index2));
            [Wpvalues,Wresults]=ranksum(DataList{1}(:,index1),DataList{2}(:,index2));
            [AD1results, AD1pvalues]=adtest(DataList{1}(:,index1));
            [AD2results, AD2pvalues]=adtest(DataList{2}(:,index2));
            end
        end
        end
    end
    disp({currentRes,KSpvalues,Tpvalues,Wpvalues,AD1pvalues,AD2pvalues});
    resultTable{j,1}=currentRes;
    resultTable{j,2}=KSpvalues;
    resultTable{j,3}=Tpvalues;
    resultTable{j,4}=Wpvalues;
    resultTable{j,5}=AD1pvalues;
    resultTable{j,6}=AD2pvalues;
end
outputTable=cell2table(resultTable);
outputTable=sortrows(outputTable,3);
outputTable.Properties.VariableNames = {'Name','KSp','Tp','Wp','AD1p','AD2p'};
file="ResidueHbond"+name1+"vs"+name2+".txt";
writetable(outputTable,file, 'Delimiter','\t');

fid = fopen(file,'A');
fprintf(fid,"Common residues and significance of difference\n\n");
fprintf(fid, string(name1+" Unique Resiudes\n"));
fprintf(fid,'%s  ', GAActivediff{:});
fprintf(fid, '\n')
fprintf(fid, string(name2+" Unique Resiudes\n"));
fprintf(fid,'%s  ', GAInactivediff{:});
fprintf(fid, '\n')
fclose(fid);
disp("fin");

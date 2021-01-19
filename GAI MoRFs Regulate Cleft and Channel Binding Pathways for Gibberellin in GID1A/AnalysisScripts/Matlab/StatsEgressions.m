
%load the data as a cell matrix. Each row is addative in counting raw
% # of times the residue of that column interacts with GAx

%'MoRF1','MoRF2','MoRF3','lig','ligWater'
% prep field, change for different analysis
referenceTags={'lig'}; 
bond='H0bond';% 'Hbond'  'H0bond'
simulation='TableEgress'; % 'Table'   'TableUncap' 'TableEgress'
analName='GID1A-Actives_Hao';
finputStruct="GID1A-GA3_last.pdb";
%"GID1A-GA3_last.pdb";%ref.pdb Alt_GID1A_example  Hao_GID1A_example
% ----------------------------------------------------------------
%'/HbondTableEgress'   '/HbondTableUncap' '/HbondTable'
%'/H0bondTableEgress'   '/H0bondTableUncap' '/H0bondTable'
bondType=strcat('/',bond,simulation); 
% % % % "GA1216ox_Alt","GA1216oxminus"
ActiveNames={"GID1A-GA1216oxminus"};
ashortNames={"GID1A-GA1216ox_Hao"};

InactiveNames={"GID1A-GA3","GID1A-GA4"};
ishortNames={"GID1A-GA3_Hao","GID1A-GA4_Hao"};
% ActiveNames={"GA7minus","GA1minus","GA3minus","GA4minus","GA1216oxminus"};
% ashortNames={"GA7_Alt","GA1_Alt","GA3_Alt","GA4_Alt","GA1216ox_Alt"};
% 
% InactiveNames={"GA20minus","GA9minus","GA8minus","GA34minus","GA41617oxminus",...
%     "GA12minus"};
% ishortNames={"GA20_Alt","GA9_Alt","GA8_Alt","GA34_Alt","GA41617ox_Alt","GA12_Alt"};

% ActiveNames={"GID1A-GA7","GID1A-GA3","GID1A-GA4","GID1A-GA1216oxminus"};
% ashortNames={"GID1A-GA7_Alt","GID1A-GA3_Alt","GID1A-GA4_Alt","GID1A-GA1216ox_Alt"};
% 
% InactiveNames={"GID1A-GA20","GID1A-GA4MeO","GID1A-GA12"};
% ishortNames={"GID1A-GA20_Alt","GID1A-GA4MeO_Alt","GID1A-GA12_Alt"};

% ActiveNames={"GID1A-GA7","GID1A-GA3","GID1A-GA4","GID1A-GA1216oxminus"};
% ashortNames={"GID1A-GA7","GID1A-GA3","GID1A-GA4","GID1A-GA1216ox"};
% 
% InactiveNames={"GID1A-GA20","GID1A-GA4MeO","GID1A-GA12"};
% ishortNames={"GID1A-GA20","GID1A-GA4MeO","GID1A-GA12"};
% 
% ActiveNames={"GID1A-GA7","GID1A-GA3","GID1A-GA4","GID1A-GA1216oxminus"};
% ashortNames={"GID1A-GA7_Hao","GID1A-GA3_Hao","GID1A-GA4_Hao","GID1A-GA1216ox_Hao"};
% 
% InactiveNames={"GID1A-GA20","GID1A-GA4MeO","GID1A-GA12"};
% ishortNames={"GID1A-GA20_Hao","GID1A-GA4MeO_Hao","GID1A-GA12_Hao"};
% "GA1216oxminus","GA1216ox_Hao"
% ActiveNames={"GA7minus","GA1minus","GA3minus","GA4minus","GA1216oxminus"};
% ashortNames={"GA7_Hao","GA1_Hao","GA3_Hao","GA4_Hao","GA1216ox_Hao"};
% 
% InactiveNames={"GA20minus","GA9minus","GA8minus","GA34minus","GA41617oxminus",...
%     "GA12minus"};
% ishortNames={"GA20_Hao","GA9_Hao","GA8_Hao","GA34_Hao","GA41617ox_Hao","GA12_Hao"};
%                            Active Dir names
% "GA7minus","GA1minus","GA3minus","GA4minus","GA1216oxminus"
%                           Inactive Dir names
%  "GA20minus","GA9minus","GA8minus","GA34minus","GA41617oxminus"
%  "GA4MeO","GA12minus"
%                           Short tags for all 
% "GA7","GA1","GA3","GA4","GA20","GA9","GA8","GA34","GA41617ox","GA4MeO","GA12"
% "APO","GA1216ox"
% 
analysisTag={};
for i=1:length(referenceTags)
    analysisTag{i}=strcat('_',referenceTags{i},'.mat');
end
%This loads the text into cells, cells 2 contain egression time in frames
for ref=1:length(referenceTags)
uniqueTag=strcat(bond,simulation,"_",referenceTags{ref},"_",analName);
file=strcat(uniqueTag,".txt");


GAActiveresidues={};
GAInactiveresidues={};
% Some of these have different names inside the struct, idk why, im stupid
for i=1:length(ActiveNames)
    tempActiveList=load(strcat('./',ActiveNames{i},bondType,ashortNames{i},analysisTag{ref}));
    try
    GAActiveResContactNumList{i}=cell2mat(tempActiveList.ResContactNumList);
    GAActiveresidues{i}={tempActiveList.ResContactNameList{end,:}};
    catch
    GAActiveResContactNumList{i}=cell2mat(tempActiveList.ResContactList.values);
    GAActiveresidues{i}=tempActiveList.ResContactList.keys;
    end
end
for i=1:length(InactiveNames)
    tempInactiveList=load(strcat('./',InactiveNames{i},bondType,ishortNames{i},analysisTag{ref})); 
    try 
        GAInactiveResContactNumList{i}=cell2mat(tempInactiveList.ResContactNumList);
        GAInactiveresidues{i}={tempInactiveList.ResContactNameList{end,:}};
    catch
        GAInactiveResContactNumList{i}=cell2mat(tempInactiveList.ResContactList.values);
        GAInactiveresidues{i}=tempInactiveList.ResContactList.keys;
    end
end
% Each row now represents a single egression run 
% Data matrix of DGG systems, and egressions, mixed sizes
DataList={GAActiveResContactNumList,GAInactiveResContactNumList};
newNumCount=[];
for k=1:length(DataList)
    for m=1:length(DataList{k})
        RowLen=size(DataList{k}{m},1);
        for i=2:size(DataList{k}{m},1)
    %         newNumCount=(DataList{k}(i:RowLen,:)-DataList{k}(i-1,:));
            newNumCount=bsxfun(@minus,DataList{k}{m}(i:RowLen,:),DataList{k}{m}(i-1,:));
            DataList{k}{m}(i:RowLen,:)=newNumCount;
        end
    end
end
% Unique and common residues 
% So this is all or nothing, either every egression/system has it to be
% common or every subgroup has it to be unique
% Done per row, egression/simulation,
for i=1:length(GAActiveresidues)  
    if isequal(i,1)
        commonResidues=intersect(GAActiveresidues{1},GAInactiveresidues{1});
        % setdiff returns pos. 1 unique w/ resp to pos. 2
        GAInactivediff=setdiff(GAInactiveresidues{1},GAActiveresidues{1});
        GAActivediff=setdiff(GAActiveresidues{1},GAInactiveresidues{1});
    else
        commonResidues=intersect(intersect(GAActiveresidues{i},GAInactiveresidues{1}),commonResidues);
        GAInactivediff=intersect(setdiff(GAInactiveresidues{1},GAActiveresidues{i}),GAInactivediff);
        GAActivediff=intersect(setdiff(GAActiveresidues{i},GAInactiveresidues{1}), GAActivediff);
    end
    endNum=i;
end
for i=1:length(GAInactiveresidues)  
    commonResidues=intersect(intersect(GAActiveresidues{endNum},GAInactiveresidues{i}),commonResidues);
    GAInactivediff=intersect(setdiff(GAInactiveresidues{i},GAActiveresidues{endNum}),GAInactivediff);
    GAActivediff=intersect(setdiff(GAActiveresidues{endNum},GAInactiveresidues{i}), GAActivediff);
end

resultTable=cell(length(commonResidues),8);%table('VariableNames',['Name','KSp','Tp','Wp','AD1p','AD2p']);

for j=1:length(commonResidues)
    currentRes=commonResidues{j};
    activeCurrentResDist=[];
    inactiveCurrentResDist=[];
    for i=1:length(GAActiveresidues)
        disp(i);
        for iA=1:length(GAActiveresidues{i})
%             if isequal(currentRes, GAActiveresidues{i}{iA})
%                 tempList=[activeCurrentResDist;GAActiveResContactNumList{i}(:,iA)];
%                 activeCurrentResDist=cat(1,tempList);
%                 break;
%             end
             % For tRAMD runs
            if isequal(currentRes, GAActiveresidues{i}{end,iA})
                tempList=[activeCurrentResDist;GAActiveResContactNumList{i}(:,iA)];
                activeCurrentResDist=cat(1,tempList);
%                 disp(activeCurrentResDist);
                break;
            end
        end
        disp("For looped");
    end
    for i=1:length(GAInactiveresidues)
        disp(i);
        for iA=1:length(GAInactiveresidues{i})
%             if isequal(currentRes, GAInactiveresidues{i}{iA})
%                 tempList=[inactiveCurrentResDist;GAInactiveResContactNumList{i}(:,iA)];
%                 inactiveCurrentResDist=cat(1,tempList);
%                 continue
%             end
              if isequal(currentRes, GAInactiveresidues{i}{end,iA})
                tempList=[inactiveCurrentResDist;GAInactiveResContactNumList{i}(:,iA)];
                inactiveCurrentResDist=cat(1,tempList);
                break;
              end
        end
        disp("For looped");
    end
    [KSresults,KSpvalues]=kstest2(activeCurrentResDist,inactiveCurrentResDist);
    [Tresults, Tpvalues]=ttest2(activeCurrentResDist,inactiveCurrentResDist);
    [Wpvalues,Wresults]=ranksum(activeCurrentResDist,inactiveCurrentResDist);
    try 
    [AD1results, AD1pvalues]=adtest(activeCurrentResDist);
    [AD2results, AD2pvalues]=adtest(inactiveCurrentResDist);
    catch
        AD1pvalues=0;
        AD2pvalues=0;
    end
    disp({currentRes,KSpvalues,Tpvalues,Wpvalues,AD1pvalues,AD2pvalues});
    resultTable{j,1}=currentRes;
    resultTable{j,2}=KSpvalues;
    resultTable{j,3}=Tpvalues;
    resultTable{j,4}=Wpvalues;
    resultTable{j,5}=AD1pvalues;% 0;
    resultTable{j,6}=AD2pvalues;% 0;
    resultTable{j,7}=mean(activeCurrentResDist);
    resultTable{j,8}=mean(inactiveCurrentResDist);
end
outputTable=cell2table(resultTable);
outputTable=sortrows(outputTable,3);
outputTable.Properties.VariableNames = {'Name','KSp','Tp','Wp','AD1p','AD2p','Active','Inactive'};

writetable(outputTable,file, 'Delimiter','\t');

fid = fopen(file,'A');
fprintf(fid,"Common residues and significance of difference\n\n");
fprintf(fid, strcat(ashortNames{:}," Active GA Unique Resiudes\n"));
fprintf(fid,'%s  ', GAActivediff{:});
fprintf(fid, '\n')
fprintf(fid, strcat(ishortNames{:}," Inactive GA Unique Resiudes\n"));
fprintf(fid,'%s  ', GAInactivediff{:});
fprintf(fid, '\n')
fclose(fid);
disp("fin");

% Plot a B factor overlay to .pdb here. Give 1 and 0.5 for above
% 0.05 and below on common res only
PDBdata = pdb2mat(finputStruct);
% normContactFreq=arrayfun(@(x) x/max(max(A)),A(k,:));
% tempNames=char(ResContactNameList{k,:});
normResContactList = containers.Map(outputTable.Name,outputTable.Tp);
seqLength=PDBdata.resNum(length(PDBdata.resNum))-PDBdata.resNum(1);
tempResNum = PDBdata.resNum();
% if tempResNum(1) ~= 1
%     diffFromOne = tempResNum(1);
%     disp(diffFromOne);
%     tempResNum = tempResNum - diffFromOne;
% end
%disp(unique(tempResNum));
tempBFactor = PDBdata.betaFactor();
for i=2:8000
    tempResID=strcat(PDBdata.resName{i},num2str(tempResNum(i)+1));%For DGG
    %if( strcmp(tempResNum{i},tempResNum{i-1}))
    if isKey(normResContactList,tempResID)
%         disp(tempResID);
        if normResContactList(tempResID) < 0.05
        tempBFactor(i) =  1;
        else
        tempBFactor(i) =  0.5;
        end
    else
        tempBFactor(i)=0;
    end
end

%-------------------------------------------------------------------------
PDBdata.betaFactor = tempBFactor;
PDBdata.outfile = "ref_"+uniqueTag+".pdb";
mat2pdb(PDBdata);
end

function [PDBdata] = pdb2mat(readFile)
%%%  -- pdb2mat.m --
% This program is the most speedy way to read a PDB file that I could come
% up with. It's function is simple: give it a PDB file and out comes a
% matlab-friendly data structure. In cumbersomely large PDB's (such as those that 
% include solvent), this can shave off a good amount of time relative to
% many programs. Unfortunately there is no easy way to hasten the slowest
% step, which is turning strings into doubles.
%
% The output format is as given in online documentation 
% (as of July 2012 when writing this program)
% http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
%
% It outputs 14 pieces total of information about the PDB. 
% 
% -- mandatory information (11) --
% 
% outfile    (the name of the PDB, this is the only input on the command line)
% 
% recordName (the class or type of atom, such as ATOM, HETATM, SOL, etc)
% atomNum    (serial number of the atom)
% atomName   (elemental identification of the atom)
% altLoc     (alt. location indicator)
% resName    (name of the amino acid/residue)
% 
% chainID    (protein chain identifier)
% resNum     (index number of the amino acid)
% X          (X position of atom)
% Y          (Y position of atom)
% Z          (Z position of atom)
% 
% -- optional information (4) --
% These are extra data about the atoms. In PDBQT's they hold the partial
% charge, for CHARMM this is the chain name, and so on. 
% 
% occupancy
% betaFactor
% element
% charge
%
% 
% -- example usage: plot the atoms of 3IJU.pdb --
% 
% 
% PDBdata = pdb2mat('3IJU.pdb');               % read in data from PDB file
% plot3(PDBdata.X, PDBdata.Y, PDBdata.Z, '.'); % make a 3D plot of data
% 
% -- example usage: translate the atoms of 3IJU.pdb by 10 angstroms in x direction --
% 
% PDBdata = pdb2mat('3IJU.pdb');               % read in data from PDB file
% PDBdata.X = PDBdata.X + 10;                  % translate coordinates
% PDBdata.outfile = '3IJU_tran10angXdir.pdb';  % update file name
% mat2pdb(PDBdata);                            % output data in PDB format
% 
%%% --- HOW TO MAKE THIS CODE FASTER! --- >> COMMENT OUT WHAT YOU DON'T USE!!
%
% This program reads everything about the PDB by default. If you want a
% faster code for whatever reason, you can comment out the lines you don't
% need. Each numeric data removed (such at resNum, or betaFactor) speeds it 
% up by 7-8%. Each string data removed (such as resName or atomName) speeds
% it up by 1-2%.
%% -- OUTPUT --

%tic;

PDBdata.outfile = readFile;

% initialize file
FileID = fopen(readFile);
rawText = fread(FileID,inf,'*char');

% parse lines by end-of-lines
splitLines = strread(rawText, '%s', 'delimiter', '\n');

% initialize variables
numLines = length(splitLines);

recordName = cell(1,numLines);
atomNum    = cell(1,numLines);
atomName   = cell(1,numLines);
altLoc     = cell(1,numLines);
resName    = cell(1,numLines);

chainID    = cell(1,numLines);
resNum     = cell(1,numLines);
X          = cell(1,numLines);
Y          = cell(1,numLines);
Z          = cell(1,numLines);

comment    = cell(1,numLines);

% read each line
m = 1;
for n = 1:numLines
    
    thisLine = cell2mat(splitLines(n));
    
    if length(thisLine) > 53 && sum(isstrprop(thisLine(23:53), 'alpha')) == 0
        
        recordName(m) = {thisLine(1:6)};
        atomNum(m)    = {thisLine(7:11)};
        atomName(m)   = {thisLine(13:16)};
        altLoc(m)     = {thisLine(17)};
        resName(m)    = {thisLine(18:20)};
        
        chainID(m)    = {thisLine(22)};
        resNum(m)     = {thisLine(23:26)};
        X(m)          = {thisLine(31:38)};
        Y(m)          = {thisLine(39:46)};
        Z(m)          = {thisLine(47:54)};
        
        comment(m)            = {thisLine(55:end)};
        
        m = m + 1;
    end
    
end

% trim exess
keepData = logical(strcmp(recordName,'ATOM  ') + strcmp(recordName,'HETATM'));

recordName = recordName(keepData);
atomNum    = atomNum(keepData);
atomName   = atomName(keepData);
altLoc     = altLoc(keepData);
resName    = resName(keepData);

chainID    = chainID(keepData);
resNum     = resNum(keepData);
X          = X(keepData);
Y          = Y(keepData);
Z          = Z(keepData);

comment    = comment(keepData);

% parse out "comment" section
occupancy  = cell(1, length(recordName));
betaFactor = cell(1, length(recordName));
element    = cell(1, length(recordName));
charge     = cell(1, length(recordName));

% fix spacing
for n = 1:length(recordName)
    thisLine = sprintf('%-26s',cell2mat(comment(n)));
    occupancy(n)  = {thisLine(1:6)};
    betaFactor(n) = {thisLine(7:12)};
    element(n)    = {thisLine(13:24)};
    charge(n)     = {thisLine(25:26)};
end

% reformat data for convenience
PDBdata.recordName = strtrim(recordName);
PDBdata.atomNum    = str2double(atomNum);
PDBdata.atomName   = strtrim(atomName);
PDBdata.altLoc     = altLoc;
PDBdata.resName    = strtrim(resName);

PDBdata.chainID    = chainID;
PDBdata.resNum     = str2double(resNum);
PDBdata.X          = str2double(X);
PDBdata.Y          = str2double(Y);
PDBdata.Z          = str2double(Z);

PDBdata.occupancy  = str2double(occupancy);
PDBdata.betaFactor = str2double(betaFactor);
PDBdata.element    = strtrim(element);
PDBdata.charge     = strtrim(charge);

% I commented these lines out, since they cause more problems than they
% solve. They do clean up the output for certain situations.

% if isnan(PDBdata.occupancy(1))
%     PDBdata.occupancy = strtrim(PDBdata.occupancy);
% end
% if isnan(PDBdata.betaFactor(1))
%     PDBdata.occupancy = strtrim(PDBdata.betaFactor);
% end



% close file
fclose(FileID);

%toc;

end
function mat2pdb(input)
%% review XYZ coordinate data 

% coordinate data is required! Checking XYZ input
if ~isfield(input, 'X') || ~isfield(input, 'Y') || ~isfield(input, 'Z')
    fprintf('we need xyz coordinate data to make a PDB!!\n\texiting...\n');
    return;
end
X = input.X;
Y = input.Y;
Z = input.Z;
if length(X) ~= length(Y) || length(X) ~= length(Z)
    fprintf('xyz coordinate data is not of equal lengths!\n\texiting...\n');
    return;
end

%% review optional data inputs

% in case optional data data not given, fill in blanks
if ~isfield(input, 'outfile')
    input.outfile = 'mat2PDB.pdb';
end
if ~isfield(input, 'recordName')
    input.recordName = cell(1,length(X));
    input.recordName(1:end) = {'ATOM'};
end
if ~isfield(input, 'atomNum')
    input.atomNum = 1:length(X);
end
if ~isfield(input, 'atomName')
    input.atomName = cell(1,length(X));
    input.atomName(1:end) = {'OW'};
end
if ~isfield(input, 'altLoc')
    input.altLoc = cell(1,length(X));
    input.altLoc(1:end) = {' '};
end
if ~isfield(input, 'resName')
    input.resName = cell(1,length(X));
    input.resName(1:end) = {'SOL'};
end
if ~isfield(input, 'chainID')
    input.chainID = cell(1,length(X));
    input.chainID(1:end) = {'A'};
end
if ~isfield(input, 'resNum')
    input.resNum = 1:length(X);
end
if ~isfield(input, 'occupancy')
    input.occupancy = ones(1,length(X));
end
if ~isfield(input, 'betaFactor')
    input.betaFactor = zeros(1, length(X));
end
if ~isfield(input, 'element')
    input.element = cell(1,length(X));
    input.element(1:end) = {'O'};
end
if ~isfield(input, 'charge')
    input.charge = cell(1,length(X));
    input.charge(1:end) = {' '};
end

outfile    = input.outfile;
recordName = input.recordName;
atomNum    = input.atomNum;
atomName   = input.atomName;
altLoc     = input.altLoc;
resName    = input.resName;
chainID    = input.chainID;
resNum     = input.resNum;
occupancy  = input.occupancy;
betaFactor = input.betaFactor;
element    = input.element;
charge     = input.charge;

%% remove faulty inputs

if length(recordName) ~= length(X)
    fprintf('recordName input is not the correct length!\n\tignoring user input\n');
    recordName = cell(1,length(X));
    recordName(1:end) = {'ATOM'};
end
if length(atomNum) ~= length(X)
    fprintf('atom serial number input is not the correct length!\n\tignoring user input\n');
    atomNum = 1:length(X);
end
if length(atomName) ~= length(X)
    fprintf('atom name input is not the correct length!\n\tignoring user input\n');
    atomName = cell(1,length(X));
    atomName(1:end) = {'OW'};
end
if length(altLoc) ~= length(X)
    fprintf('alternate location input is not the correct length!\n\tignoring user input\n');
    altLoc = cell(1,length(X));
    altLoc(1:end) = {' '};
end
if length(resName) ~= length(X)
    fprintf('residue name input is not the correct length!\n\tignoring user input\n');
    resName = cell(1,length(X));
    resName(1:end) = {'SOL'};
end
if length(chainID) ~= length(X)
    fprintf('chain ID input is not the correct length!\n\tignoring user input\n');
    chainID = cell(1,length(X));
    chainID(1:end) = {'A'};
end
if length(resNum) ~= length(X)
    fprintf('residue number input is not the correct length!\n\tignoring user input\n');
    resNum = 1:length(X);
end
if length(occupancy) ~= length(X)
    fprintf('occupancy input is not the correct length!\n\tignoring user input\n');
    occupancy = ones(1,length(X));
end
if length(betaFactor) ~= length(X)
    fprintf('beta factor input is not the correct length!\n\tignoring user input\n');
    betaFactor = zeros(1, length(X));
end
if length(element) ~= length(X)
    fprintf('element symbol input is not the correct length!\n\tignoring user input\n');
    element = cell(1,length(X));
    element(1:end) = {'O'};
end
if length(charge) ~= length(X)
    fprintf('charge input is not the correct length!\n\tignoring user input\n');
    charge = cell(1,length(X));
    charge(1:end) = {' '};
end

% fix atomName spacing
for n = 1:length(atomName)
    atomName(n) = {sprintf('%-3s',cell2mat(atomName(n)))};
end


%% create PDB

% open file
% fprintf('outputting PDB in file %s\n', outfile);
FILE = fopen(outfile, 'w');

% output data
for n = 1:length(atomNum)
    
    % standard PDB output line
    fprintf( FILE, '%-6s%5u%5s%1.1s%3s %1.1s%4u%12.3f%8.3f%8.3f%6.2f%6.2f%12s%2s\n', ...
        cell2mat(recordName(n)), atomNum(n), cell2mat(atomName(n)), ...
        cell2mat(altLoc(n)), cell2mat(resName(n)), cell2mat(chainID(n)), ...
        resNum(n), X(n), Y(n), Z(n), occupancy(n), betaFactor(n), ...
        cell2mat(element(n)), cell2mat(charge(n)));
    
%     % output progress in terminal
%     if ~mod(n,400)
%         fprintf('   %6.2f%%', 100*n / length(atomNum));
%         if ~mod(n, 4000)
%             fprintf('\n');
%         end
%     end
    
end
fprintf( FILE, 'END\n');

% close file
% fprintf('   %6.2f%%\n    done! closing file...\n', 100);

fclose(FILE);

end
%%%
%Author: JP
%Adaption of first4_4 and analysis work up for Hbond list produced in vmd
%The file format expected:Frame #Hb #H20b following pattern[3let ResNum]
%ResNum comes from the file put into NAMD, so numbering specific
%%%
set(0,'DefaultFigureColor','white')
fig.InvertHardcopy = 'off';
figure('Color','white','Visible','off');

width = 8.37;     % Width in inches - adjust as necessary
height = 3.84;    % Height in inches - adjust as necessary
alw = 1.5;%0.75;    % AxesLineWidth 
fsz = 14;           % Fontsize 
lw = 1.5;           % LineWidth dir (fullfile(
msz = 8;            % MarkerSize 

set(0,'defaultAxesFontSize',fsz); 
set(0,'defaultLineLineWidth',lw);   
set(0,'defaultLineMarkerSize',msz); 
set(0,'defaultAxesLineWidth',alw); 

% finputData=input("Match name of data files\n",'s');
% finputStruct=input("Name of pdb file\n",'s');
% uniqueTag=input("Name of unique tag\n",'s');

sysName=input("Name of system\n",'s');


function hbondEgressionMapper(sysName)
finputData=["MoRF1_Hbond_3.5.dat","MoRF2_Hbond_3.5.dat",...
    "MoRF3_Hbond_3.5.dat","lig_Hbond_3.5.dat",...
    "ligWater_Hbond_3.5.dat"];
finputStruct="ref.pdb";
uniqueTag=[strcat(sysName,"_MoRF1"),strcat(sysName,"_MoRF2"),...
    strcat(sysName,"_MoRF3"),strcat(sysName,"_lig"),strcat(sysName,"_ligWater")];

for a=1: legnth(finputData)
hold on;
pdbFrame = pdb2mat(finputStruct);
FileList = dir(strcat('./VMD_output/egression*',finputData{a}));
colourList=hsv(length(FileList));
TotalResList={};
TotalWaterList={};
for i=1 : length(FileList)
legnendList(i)="Egression "+i;
if( isfile(strcat('./VMD_output/',FileList(i).name)) )
   fid = fopen(strcat('./VMD_output/',FileList(i).name));
   textLine = fgets(fid); % Read first line.
   lineCounter = 1;

 while ischar(textLine)
  fprintf('\nLine #%d of text = %s\n', lineCounter, textLine);
  % get into numbers array.
  numbers = sscanf(textLine, '%f ');
  for k = 1 : length(numbers)
    ca{lineCounter, k} = numbers(k);
  end
  TotalResList{i,lineCounter} = ca{lineCounter,2};
  TotalWaterList{i,lineCounter}=ca{lineCounter,3};
% ALternate way where the whole array is in one cell.
  %ca2{lineCounter} = numbers;
  tempLine = split(textLine)';
  Names{i,lineCounter} = {tempLine{4:length(tempLine)-1}};
% Read the next line.
  textLine = fgets(fid);
  lineCounter = lineCounter + 1;
  end
   fclose(fid);
else
   error(['data file not found!  Looked for: ',finputData{a}]);
end
%TotalResList{i,:}=[ca{:,2}];
%TotalWaterList{i,:}=[ca{:,3}];

%Plot int of Hbonds over frames
figure(1);

plot([ca{:,1}],[ca{:,2}],':o','LineWidth',length(FileList)+1-i,...
                     'Color',colourList(i,:),...
                     'MarkerEdgeColor','g',...
                     'MarkerFaceColor','g',...
                     'DisplayName',legnendList(i),...
                     'MarkerSize',3)

plot([ca{:,1}],[ca{:,3}],':o','LineWidth',length(FileList)+1-i,...
                     'Color',colourList(i,:),...
                     'MarkerEdgeColor','b',...
                     'MarkerFaceColor','b',...
                     'HandleVisibility','off',...
                     'MarkerSize',3)
end
hold off;
legend();
savefig(figure(1),"egressHbondFreq"+uniqueTag{a});  

maxLengthCell=size(TotalResList,2); %finding the longest vector in the cell array
for i=1:length(FileList)
    for j=sum(cellfun('size',TotalResList(i,:),2)):maxLengthCell
         TotalResList{i,j}=0;   %zeropad the elements in each cell array with a length shorter than the maxlength
    end
    for j=sum(cellfun('size',TotalWaterList(i,:),2)):maxLengthCell
        TotalWaterList{i,j}=0;
    end
end
AvgResFreq=mean(cell2mat(TotalResList));
StdResFreq=std(cell2mat(TotalResList));
AvgWaterFreq=mean(cell2mat(TotalWaterList));
StdWaterFreq=std(cell2mat(TotalWaterList));
figure(2);
hold on;
x=1:1:length(AvgResFreq);
%plot(x,AvgResFreq,':o','LineWidth',5,...
%                      'Color','k',...
%                      'MarkerEdgeColor','k',...
%                      'MarkerFaceColor','g',...
%                      'MarkerSize',3)
errorbar(x,AvgResFreq, StdResFreq)
%plot(x,AvgWaterFreq,':o','LineWidth',5,...
%                      'Color','k',...
%                      'MarkerEdgeColor','b',...
%                      'MarkerFaceColor','b',...
%                      'HandleVisibility','off',...
%                      'MarkerSize',3)
errorbar(x,AvgWaterFreq,StdWaterFreq)
legend(["residue","water"]);
hold off;
savefig(figure(2),"egressAvgHbondFreq"+uniqueTag{a});
%%
% Make the propensity scale onto pdb structure
%plot Histogram of H-bond contact residueID
ResContactNameList={};%containers.Map; cant use, fucking sorts the keys
ResContactNumList={};
ResidueCountCorrection = 110; %DELLA length
%The NAMD file takes all residiues for counting; 
NumOfRes = 1;
figure(3);
hold on;
for k=1 : length(FileList)
    if k>=2
        for z=1:size(ResContactNameList,2)
            ResContactNameList{k,z}=ResContactNameList{k-1,z};
            ResContactNumList{k,z}=ResContactNumList{k-1,z};
        end
        %ResContactNameList{k,:}=ResContactNameList{k-1,:};
    end
    for i=1 : length(Names)
        fold = length(Names{k,i})/2;
        for j=1 : length(Names{k,i})/2 %half array is 3 letter, other res number
            ResidueID2 = strcat(Names{k,i}{j},Names{k,i}{j+fold});
            %ResNum=str2double(Names{i}{j+fold})-ResidueCountCorrection;
            %ResidueID = strcat(Names{i}{j},num2str(ResNum));
            if isempty(ResContactNameList)
                ResContactNameList{1,NumOfRes}= ResidueID2;
                ResContactNumList{k,NumOfRes}= 1;
            end
            for z=1:size(ResContactNameList,2)
                if isequal(ResContactNameList{k,z},ResidueID2)
                    logicPresent=1;%disp ('YES')
                    idx=z;
                    break;
                else
                    logicPresent=0;%disp ('NO')
                end
            end
            if ~logicPresent%any(strcmp(ResContactNameList{k,:},ResidueID2))%ismember(ResidueID2,ResContactNameList{k,:})
                NumOfRes = NumOfRes+1;
                ResContactNameList{k,NumOfRes}= ResidueID2;
                ResContactNumList{k,NumOfRes}= 1;
            else
                %idx = find(contains(ResContactNameList{k,:}, ResidueID2));
                temp = ResContactNumList{k,idx};
                temp = temp + 1;
                ResContactNumList{k,idx}= temp;
            end
        end
    end
    
    disp({ResContactNameList{k,:}});
end
maxLengthCell=size(ResContactNameList,2); %finding the longest vector in the cell array
for i=1:length(FileList)
    for j=sum(cellfun('size',ResContactNumList(i,:),2)):maxLengthCell
         ResContactNumList{i,j}=0;   %zeropad the elements in each cell array with a length shorter than the maxlength
    end
end
A=log(double(cell2mat(ResContactNumList))); %A is your matrix
% basic histogram of h-bond residues
bar( A','stacked');%,'FaceColor',colourList
set(gca,'XTick',[1:length(ResContactNameList)]);
xtickangle(90);
set(gca,'xticklabel', char(ResContactNameList{k,:}) );
hold off;

savefig(figure(3),"HbondHistEgress"+uniqueTag{a});
save("HbondTableEgress"+uniqueTag{a}, 'ResContactNameList', 'ResContactNumList');


    % % ------------------- Plot colmn mean in FlxFlx onto PDB file
PDBdata = pdb2mat(finputStruct);
normContactFreq=arrayfun(@(x) x/max(max(A)),A(k,:));
tempNames=char(ResContactNameList{k,:});
normResContactList = containers.Map(string(tempNames),normContactFreq');
seqLength=PDBdata.resNum(length(PDBdata.resNum))-PDBdata.resNum(1);
tempResNum = PDBdata.resNum();
if tempResNum(1) ~= 1
    diffFromOne = tempResNum(1);
    disp(diffFromOne);
    tempResNum = tempResNum - diffFromOne;
end
%disp(unique(tempResNum));
tempBFactor = PDBdata.betaFactor();
for i=2:length((PDBdata.resNum))
    tempResID=strcat(PDBdata.resName{i},num2str(tempResNum(i)));
    %if( strcmp(tempResNum{i},tempResNum{i-1}))
    if isKey(normResContactList,tempResID)
        tempBFactor(i) =  normResContactList(tempResID);
    else
        tempBFactor(i)=0;
    end
end

%-------------------------------------------------------------------------
PDBdata.betaFactor = tempBFactor;
%disp(tempBFactor);
% PDBdata.betaFactor = tempFlxFlxIndx;
PDBdata.outfile = "ref_Hcontacts_egress"+uniqueTag{a}+".pdb";
mat2pdb(PDBdata);
%--------------------------------------------------------------------------
%%
% plot of h bond residues over time (attendence over time?)
figure(4);
ylim([0, length(tempNames)]);
logicIndex = string(tempNames);
yticklabels(logicIndex);
yticks([1:1:length(logicIndex)]);
hold on;
grid on;
colList=hsv(length(FileList));
%add third demention here for FT per egression...
logicMatrix=zeros(length(tempNames),length(Names));
for k=1 : length(FileList)
    for i=1:1:length(Names)
        fold = length(Names{k,i})/2;
        present=[];
        logicArray=zeros(length(tempNames),1);
        %present = zeros(1,length(logicIndex));
        for j=1 : length(Names{k,i})/2 %half array is 3 letter, other res number
            ResidueID2 = strcat(Names{k,i}{j},Names{k,i}{j+fold});
            for z=1:size(ResContactNameList,2)
                if isequal(ResContactNameList{k,z},ResidueID2)
                    logicPresent=1;%disp ('YES')
                    idx=z;
                    logicArray(z)=1;
                    break;
                else
                    logicPresent=0;%disp ('NO')
                end
            end
            if logicPresent
               logicPresentNum = find(strcmp(logicIndex, ResidueID2));
               present = [present; logicPresentNum];
               %uses actual numbers for position on y axis
            else
               error(['residue not found!  Looked for: ',ResidueID2]);
           end
        end
        x = i * ones(1, length(present));
        plot(x, present', 'o','Color',colList(k,:),'MarkerSize', 12-k, 'LineWidth', 3);
        logicMatrix(:,i)= logical(logicMatrix(:,i))+logical(logicArray);
    end
end
savefig(figure(4),"HbondAttendenceEgress"+uniqueTag{a});
hold off;
%Plot Fourier transformation of Attendence
%Note that currently this is a all sum attendence. Not egression by
%egression like it probably should be
figure(5);
ax1 = subplot(1,2,1);
plot(abs(fft(logicMatrix(5:10,:))));
lgh = legend(string(tempNames(5:10,:))');
ax2 = subplot(1,2,2);
set(lgh,'position',get(ax2,'position'));
axis(ax2,'off');
title('legend')
%%
pause(20);
close all
clear
end
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
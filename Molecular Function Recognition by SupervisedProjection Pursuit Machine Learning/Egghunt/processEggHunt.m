function processEggHunt(Ntype, egg, opv, dof)
% stripped down version of layEggs
% The program creates eggs, the concealment environment, places the eggs,
% searches for the eggs using PCA, and SPLOC, and writes the results to a
% single output file. 
%
% INPUT
% User defined parameters
% 
% PROCESS
% defines type of concealment
% create egg
% run PCA for finding the egg
% prepare and run SPLOC
% 
% OUTPUT
% single file data dump
% df  sDd  sDu  sPd  sPu  sEpdf  jS pDd  pDu  pPd  pPu  pEpdf  T1  T2  T3 
% df = degree of freedom
% sDd = sploc discriminant subspace dimension
% sDu = sploc undetermined subspace dimension
% sPd = sploc percent of egg reconstructed within discriminant subspace
% sPu = sploc percent of egg reconstructed within undetermined subspace
% sEpdf = sploc efficacy per degree of freedom
% jS = PCA result selection = {1, 2, 3} => F, N, B 
% pDd = pca discriminant subspace dimension
% pDu = pca undetermined subspace dimension
% pPd = pca percent of egg reconstructed within discriminant subspace
% pPu = pca percent of egg reconstructed within undetermined subspace
% pEpdf = pca efficacy per degree of freedom
% T1 = cputime to generate concealing-environment, place egg, setup traits
% T2 = cputime to get initial conditions for pca results and run sploc
% T3 = wall time to get initial conditions for pca results and run sploc
% ------------------------------------------------------------------------
%%                                                parameterize calculation
scaleFactor = 4; 
GNtype = Ntype;  % 'SGN','CGN' => structureless, correleted Guassian noise
placeEgg = egg;                                         % (Big,Small,None)
%k1 = 1;                               % large egg: k1=1   small egg: k1=5
%k2 = 2;                               % large egg: k2=2   small egg: k2=6
%opv = 100;        % 4 20  100                   observations per variable
d = dof;                                              % dimension of space
% ------------------------------------------------------------------------
%%                                                           sanity checks
   if( d < 10 )
   error('minimum dimension = 10');
   end
   if( scaleFactor < 0.01 )
   error('minimum scaleFactor is set at 0.01');
   end
   if( placeEgg=="Big" )
   eggDIM = 2;
   k1 = 1;                                    % Don's egg:       k1 = 1
   k2 = 2;                                    % Don's egg:       k2 = 2
   strEgg = 'B';
   elseif( placeEgg=="Small" )
   eggDIM = 2;
   k1 = 5;                                    % Don's egg:       k1 = 5
   k2 = 6;                                    % Don's egg:       k2 = 6
   strEgg = 'S';
   elseif( placeEgg=="None" )
   eggDIM = 0;
   k1 = 1;                                    % no egg, control
   k2 = 2;       
   strEgg = 'N';
   else
   error("Egg Type Not Recognized. Choose 'Big' 'Small' or 'None'");
   end
%%                                            build concealing-environment
t0 = cputime();                   % for tracking time to build data stream
tic
% -------------------------------------- build concealment characteristics
   switch GNtype
       case 'SGN'
       ave = zeros(1,d);
       sig0 = ones(1,d);
       sig = diag(sig0);
       strEnv = 'sgn';
       case 'CGN'
       ave = zeros(1,d);
       sig0 = ( sqrt(d)./(1:d) );
       sig = diag(sig0);
          for i=1:d
             for j=i+1:d
             g = sig0(j);
             sig(i,j) = g/( 1 + abs(i-j) ).^0.5;     % adjustable exponent
             sig(j,i) = sig(i,j);         % ^^^--> nice at 1/2  
             end
          end
       strEnv = 'cgn';
       otherwise
       error('unknown concealment');
   end
% --------------------- generate nonfunctional and functional data streams
n = d*opv;                                             % number of samples
nsplit = floor(n/3);              % split samples into 3 separate sections
%%                                     construct nonfunctional data matrix
xN = mvnrnd(ave,sig,n);
dataRefNameN = ['emptyness',GNtype,num2str(d,'%04i'),'OPV', ...
                 num2str(opv,'%02i')];
matrixNname = ['XNeggAbsent',GNtype,num2str(d,'%04i'),'OPV', ...
                num2str(opv,'%02i'),'.abc'];
dataMatrixN = struct;
dataMatrixN.dataRefName = dataRefNameN;                           % string
dataMatrixN.dataMatrixName = cell(1,1);            % cell array of strings
dataMatrixN.dataMatrixName{1} = matrixNname;       % cell array of strings
dataMatrixN.nVariables = d;                           % numerical variable
dataMatrixN.pType = 'notFormated';           % string defines packing type
dataMatrixN.dim = 1;                     % # of components in local vector
dataMatrixN.n = 1;                                    % number of matrices
dataMatrixN.A = cell(1,1);                             % define cell array
dataMatrixN.A{1} = xN';                                 % numerical matrix
dataMatrixN.nSamples = n;                                % numerical array
% -------------------------------------------- build trait data structures
traitNa = getMultivariateStats4sploc(dataMatrixN,n,0,'cov');         % all
traitNb = getBoostedMVStats4sploc(dataMatrixN,1,0,'cov');   % bootstrapped
traitNs = getMultivariateStats4sploc(dataMatrixN,nsplit,0,'cov');  % split
clear dataMatrixN
%%                                        construct functional data matrix
xF = mvnrnd(ave,sig,n);
%%                                                               place egg
iMAX = floor( 0.8*d + 0.01 );
iMIN = iMAX - 5;                     % consider a six dimensional subspace
Qsub = sig(iMIN:iMAX,iMIN:iMAX);
[V,D] = eig(Qsub);
[lambda,indx] = sort( diag(D), 'descend' );
lambda = sqrt(lambda);       % convert from variance to standard deviation
V = V(:,indx);
% -------------------------------------------------------- skip if control
   if( eggDIM > 0 )
   xsub = xF(:,iMIN:iMAX);
   ysub = xsub*V;
% ------------------------------------------------------------ quick check
   zsub = ysub*V'; % =>  zsub = xsub*(V*V') = xsub* I = xsub
   maxDiff = max( max( abs(xsub - zsub) ) );
      if( maxDiff > 1.0e-12 )
      error('inverse of sub-covariance matrix is incorrect');
      end
% --------------------------------------------------------- transform data
   xs = 0; %2*a;
   ys = 0; %2*b;
   ysub(:,k1) = xs + scaleFactor*lambda(k1)*randn(n,1); 
   ysub(:,k2) = ys + scaleFactor*lambda(k2)*randn(n,1);
   new_xsub = ysub*V';
   xF(:,iMIN:iMAX) = new_xsub;
   end
% -------------------------------------------------- construct dataMatrixF
   if( eggDIM > 0 )
   dataRefNameF = ['eggPlaced',GNtype,num2str(d,'%04i'),'OPV', ...
                   num2str(opv,'%02i')];
   matrixFname = ['XFeggPlaced',GNtype,num2str(d,'%04i'),'OPV', ...
                   num2str(opv,'%02i'),'.abc'];
   else
   dataRefNameF = ['eggAbsent',GNtype, ...
                   num2str(d,'%04i'),'OPV',num2str(opv,'%02i')];
   matrixFname = ['XFeggAbsent',GNtype,num2str(d,'%04i'),'OPV', ...
                  num2str(opv,'%02i'),'.abc'];
   end
dataMatrixF = struct;
dataMatrixF.dataRefName = dataRefNameF;                           % string
dataMatrixF.dataMatrixName = cell(1,1);            % cell array of strings
dataMatrixF.dataMatrixName{1} = matrixFname;       % cell array of strings
dataMatrixF.nVariables = d;                           % numerical variable
dataMatrixF.pType = 'notFormated';           % string defines packing type
dataMatrixF.dim = 1;                     % # of components in local vector
dataMatrixF.n = 1;                                    % number of matrices
dataMatrixF.A = cell(1,1);                             % define cell array
dataMatrixF.A{1} = xF';                                 % numerical matrix
dataMatrixF.nSamples = n;                                % numerical array
% -------------------------------------------- build trait data structures
traitFa = getMultivariateStats4sploc(dataMatrixF,n,0,'cov');         % all
traitFb = getBoostedMVStats4sploc(dataMatrixF,1,0,'cov');   % bootstrapped
traitFs = getMultivariateStats4sploc(dataMatrixF,nsplit,0,'cov');  % split
clear dataMatrixF
%%                                                     find eggs using PCA
muF = mean(xF);
AF = (xF - muF);
QF = (AF'*AF)/(n-1);
% --------------------------
muN = mean(xN);
AN = (xN - muN);
QN = (AN'*AN)/(n-1);
% --------------------------
xB = [xF; xN];
muB = mean(xB);
AB = (xB - muB);
QB = (AB'*AB)/(2*n-1);
% --------------------------
% % xD = xF - xN;
% % muD = mean(xD);
% % AD = (xD - muD);
% % QD = (AD'*AD)/(n-1);
% ---------------------------------------------------- compare three cases
clear AF AN AB xF xN muF  muN  muB
[sbvF,~] = eig(QF);
[sbvN,~] = eig(QN);
[sbvB,~] = eig(QB);
% % [sbvD,~] = eig(QD);
% ---------------------------------------- stage 1: generate data & do PCA
cpuGenerate = cputime() - t0;
wallTime1 = toc;
% ----------------------------------------- stage2: other sploc operations
t0 = cputime();
tic 
% ----------------------------------------------- get SPLOC basis spectrum
[splocResultF,~,~] = getBasisVecSpectrum(sbvF,'pcF',traitFs,traitNs);
[splocResultN,~,~] = getBasisVecSpectrum(sbvN,'pcN',traitFs,traitNs);
[splocResultB,~,~] = getBasisVecSpectrum(sbvB,'pcB',traitFs,traitNs);
%[splocResultD,~,~] =getBasisVecSpectrum(sbvD,'pcD',traitFs,traitNs);
clear QF QN QB
% ------------------------------------------------- pick the best PCA case
eeeF = splocResultF.efficacy;
eeeN = splocResultN.efficacy;
eeeB = splocResultB.efficacy;
%eeeD = splocResultD.efficacy;
eeeD = -dof;
[~,jsort] = sort( [eeeF,eeeN,eeeB,eeeD] );
breakDegeneracy = zeros(1,4);
breakDegeneracy( jsort(1) ) = 0.3;
breakDegeneracy( jsort(2) ) = 0.2;
breakDegeneracy( jsort(3) ) = 0.1;
%breakDegeneracy( jsort(4) )= 0.0;
eeeF = splocResultF.Dd + breakDegeneracy(1);
eeeN = splocResultN.Dd + breakDegeneracy(2);
eeeB = splocResultB.Dd + breakDegeneracy(3);
%eeeD = splocResultD.Dd + breakDegeneracy(4);
eeeD = -dof;
[~,jsort] = max( [eeeF,eeeN,eeeB,eeeD] ); 
jSelect = jsort(1);
      switch jSelect
          case 1                                % functional PCA basis set
          splocResultsPCA = splocResultF;
          case 2                             % nonfunctional PCA basis set
          splocResultsPCA = splocResultN;
          case 3                            % pooled average PCA basis set
          splocResultsPCA = splocResultB;
          case 4                                % functional PCA basis set
          splocResultsPCA = splocResultF;
      end
%%                                                   output of PCA results
Dd = splocResultsPCA.Dd; 
% ----------------------------------- calculate reconstruction percentages
reconstruct = 0;
   if( Dd > 0 )
   Ud = getDiscriminantSBV(splocResultsPCA);     % discriminant modes only
      for k=1:Dd
      reconstruct = reconstruct + sum( Ud(iMIN:iMAX,k).*V(:,k1) )^2;
      reconstruct = reconstruct + sum( Ud(iMIN:iMAX,k).*V(:,k2) )^2;
      % note that V = 0 outside of range:  iMIN:iMAX
      end
   end  
pcaPd = round( reconstruct*500 )/10;            % good to the tenths place
Du = splocResultsPCA.Du; 
reconstruct = 0;
   if( Du > 0 )
   Uu = getUndeterminedSBV(splocResultsPCA);     % discriminant modes only
      for k=1:Du
      reconstruct = reconstruct + sum( Uu(iMIN:iMAX,k).*V(:,k1) )^2;
      reconstruct = reconstruct + sum( Uu(iMIN:iMAX,k).*V(:,k2) )^2;
      % note that V = 0 outside of range:  iMIN:iMAX
      end
   end  
pcaPu = round( reconstruct*500 )/10;            % good to the tenths place
pcaDd = Dd;
pcaDu = Du;
pcaEpdf = splocResultsPCA.efficacy/dof;       % => Efficacy per dof = Epdf
%%                                                   find eggs using SPLOC
sbv = splocResultsPCA.SBV;                      % same as internal process
splocResult = sploc(2,sbv,'eggHunt',traitFa,traitNa,0);         % SPLOC it
%                     ^^^---------> same as the default initial conditions
sbv = splocResult.SBV;
splocResult = sploc(1,sbv,'eggHunt',traitFb,traitNb,0);
sbv = splocResult.SBV;
splocResult = sploc(0,sbv,'eggHunt',traitFs,traitNs,0);
%%                                                    look at SPLOC output
Dd = splocResult.Dd; 
% ----------------------------------- calculate reconstruction percentages
Du = splocResult.Du;
reconstruct = 0;
   if( Dd > 0 )
   Ud = getDiscriminantSBV(splocResult);         % discriminant modes only
      for k=1:Dd
      reconstruct = reconstruct + sum( Ud(iMIN:iMAX,k).*V(:,k1) )^2;
      reconstruct = reconstruct + sum( Ud(iMIN:iMAX,k).*V(:,k2) )^2;
      % note that V = 0 outside of range:  iMIN:iMAX
      end
   end  
splocPd = round( reconstruct*500 )/10;          % good to the tenths place
reconstruct = 0;
   if( Du > 0 )
   Uu = getUndeterminedSBV(splocResult);         % undetermined modes only
      for k=1:Du
      reconstruct = reconstruct + sum( Uu(iMIN:iMAX,k).*V(:,k1) )^2;
      reconstruct = reconstruct + sum( Uu(iMIN:iMAX,k).*V(:,k2) )^2;
      % note that V = 0 outside of range:  iMIN:iMAX
      end
   end  
splocPu = round( reconstruct*500 )/10;          % good to the tenths place
splocDd = Dd;
splocDu = Du;
splocEpdf = splocResult.efficacy/dof;         % => Efficacy per dof = Epdf
wallTime2 = toc;
cpuSPLOC = cputime() - t0;
%%                                              append information to file
% df  sDd  sDu  sPd  sPu  sEpDOF  j pDd  pDu  pPd  pPu  pEpDOF  T1  T2  T3 
% REMARK: see OUTPUT description above in header of program
outFileName = [strEnv,strEgg,num2str(opv,'%03i'),'.txt'];
fileID = fopen(outFileName,'a');
outputFormat = ['%5d',' %5d %5d %4.1f %4.1f %6.1f','   %1d', ...
                      ' %5d %5d %4.1f %4.1f %6.1f', ...
                      ' %11.1f %11.1f %11.1f %11.1f\n'];
writeLine = [dof,splocDd,splocDu,splocPd,splocPu,splocEpdf,jSelect, ...
                 pcaDd,  pcaDu,  pcaPd  ,pcaPu  ,pcaEpdf, ...
                 cpuGenerate,wallTime1,cpuSPLOC,wallTime2];
fprintf(fileID,outputFormat,writeLine); 
fclose(fileID);
% disp(writeLine);
% disp('   ');     
end
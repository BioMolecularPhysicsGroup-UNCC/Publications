% Multivariate Gaussian Distribution (MGD) plus Egg or no Egg
% The program creates a MGD to hide an egg and performs an egg hunt.
% This is a test benchmark program that serves as a stepping stone.
% Create an egg, and for 1 case at a time, try to find it with SPLOC
%
% INPUT
% User defined parameters
% 
% PROCESS
% defines type of concealment
% create egg
% prepare and run SPLOC
% 
% OUTPUT
% sploc - spectrum outputs
% d-mode projections
% feature space classification
% 
% ------------------------------------------------------------------------
close all
clear 
%%                                                parameterize calculation
GNtype = 'CGN';  % 'SGN','CGN' => structureless, correleted Guassian noise
k1 = 1;                                % large egg: k1=1   small egg: k1=5
k2 = 2;                                % large egg: k2=2   small egg: k2=6
scaleFactor = 4; 
placeEgg = true;               % (true,false) => (egg present, egg absent)
opv = 20;        % 4 20  100
d = 1000;                                             % dimension of space
% ------------------------------------------------------------------------
vb = 2;                                            % verbosity for sploc()
%%                                                           sanity checks
   if( d < 10 )
   error('minimum dimension = 10');
   end
   if( scaleFactor < 0.01 )
   error('minimum scaleFactor is set at 0.01');
   end
   if( placeEgg )
   eggDIM = 2;
   else
   eggDIM = 0;
   end
%%                                         set default graphics parameters
% REMARKS: These set functions are to make the plots ~publication quality.
%          Using the default set functions for axes, font size and figure
%          color are for convenience only.
%fig.InvertHardcopy = 'off';         
set(0,'DefaultFigureColor','white')
set(0,'defaultAxesFontSize',14); 
set(0,'defaultLineLineWidth',1.8);   
set(0,'defaultLineMarkerSize',7); 
set(0,'defaultAxesLineWidth',1.5);
%set(0,'defaultGeographicAxes','on');
%%                                                        initialize sploc
initializeSPLOC(1,'fName','FriedEgg','gType','png');
setDataMatrixFormat('notVectored',2);               % by-pass screen input
dof = d;
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
xD = xF - xN;
muD = mean(xD);
AD = (xD - muD);
QD = (AD'*AD)/(n-1);
% ----------------------------------------------------- compare four cases
[sbvF,lambdaF] = eig(QF);
[sbvN,lambdaN] = eig(QN);
[sbvB,lambdaB] = eig(QB);
[sbvD,lambdaD] = eig(QD);
lambdaF = diag(lambdaF);
lambdaN = diag(lambdaN);
lambdaB = diag(lambdaB);
lambdaD = diag(lambdaD);
lambdaF = sort(lambdaF,'descend');
lambdaN = sort(lambdaN,'descend');
lambdaB = sort(lambdaB,'descend');
lambdaD = sort(lambdaD,'descend');
% ---------------------------------------- stage 1: generate data & do PCA
cpuGenerate = cputime() - t0;
wallTime1 = toc;
% ----------------------------------------- stage2: other sploc operations
t0 = cputime();
tic 
% ----------------------------------------------- get SPLOC basis spectrum
[splocResultF,~,mS2U_F] = getBasisVecSpectrum(sbvF,'pcF',traitFs,traitNs);
[splocResultN,~,mS2U_N] = getBasisVecSpectrum(sbvN,'pcN',traitFs,traitNs);
[splocResultB,~,mS2U_B] = getBasisVecSpectrum(sbvB,'pcB',traitFs,traitNs);
[splocResultD,~,mS2U_D] = getBasisVecSpectrum(sbvD,'pcD',traitFs,traitNs);
DF = splocResultF.Dd;
DN = splocResultN.Dd;
DB = splocResultB.Dd;
DD = splocResultD.Dd;
ddd4 = [DF,DN,DB,DD];
DFu = splocResultF.Du;
DNu = splocResultN.Du;
DBu = splocResultB.Du;
DDu = splocResultD.Du;
uuu4 = [DFu,DNu,DBu,DDu];
DFi = splocResultF.Di;
DNi = splocResultN.Di;
DBi = splocResultB.Di;
DDi = splocResultD.Di;
iii4 = [DFi,DNi,DBi,DDi];
eeeF = splocResultF.efficacy;
eeeN = splocResultN.efficacy;
eeeB = splocResultB.efficacy;
eeeD = splocResultD.efficacy;
eee4 = [eeeF,eeeN,eeeB,eeeD]/dof;             % => Efficacy per dof = Epdf
disp('  ');
disp(dividerLine( 'comparison of PCA cases: F,N,B,D') );
disp( [ddd4; uuu4; iii4; eee4] );
disp('  ');
% ------------------------------------------------- pick the best PCA case
ddd4(4) = -dof;                      % this eliminates option 4 altogether
[~,jsort] = sort(eee4,'descend');
breakDegeneracy = zeros(1,4);
breakDegeneracy( jsort(1) ) = 0.3;
breakDegeneracy( jsort(2) ) = 0.2;
breakDegeneracy( jsort(3) ) = 0.1;
%breakDegeneracy( jsort(4) )= 0.0;
DFe = ddd4(1) + breakDegeneracy(1);
DNe = ddd4(2) + breakDegeneracy(2);
DBe = ddd4(3) + breakDegeneracy(3);
DDe = ddd4(4) + breakDegeneracy(4);
[~,jsort] = max( [DFe,DNe,DBe,DDe] ); 
jSelect = jsort(1);
      switch jSelect
          case 1                                % functional PCA basis set
          lambda = lambdaF;
          splocResultsPCA = splocResultF;
          mapS2U = mS2U_F;
          case 2                             % nonfunctional PCA basis set
          mapS2U = mS2U_N;
          splocResultsPCA = splocResultN;
          lambda = lambdaN;
          case 3                            % pooled average PCA basis set
          lambda = lambdaB;
          splocResultsPCA = splocResultB;
          mapS2U = mS2U_B;
          case 4                                % difference PCA basis set
          lambda = lambdaD;
          splocResultsPCA = splocResultD;
          mapS2U = mS2U_D;
      end
% ---------------------------------- determine initial condition for sploc
% ttt = eee4 + ddd4 - uuu4;
% [~,iSelect] = max(ttt);
%       switch iSelect
%           case 1                              % functional PCA basis set
%           splocResult0 = splocResultF;
%           case 2                           % nonfunctional PCA basis set
%           splocResult0 = splocResultN;
%           case 3                          % pooled average PCA basis set
%           splocResult0 = splocResultB;
%           case 4                              % difference PCA basis set
%           splocResult0 = splocResultD;
%       end
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
% ----------------------------------------------------- determine loadings
pcaLoadings = zeros(d,1);
   if( Dd > 0 )
   Ud = getDiscriminantSBV(splocResultsPCA);     % discriminant modes only
      for k=1:Dd
      pcaLoadings = pcaLoadings + Ud(:,k).^2;
      end
   end
%%                                                   find eggs using SPLOC
sbv = splocResultsPCA.SBV;                      % same as internal process
splocResult = sploc(2,sbv,'eggHunt',traitFa,traitNa,vb);        % SPLOC it
%                     ^^^---------> same as the default initial conditions
sbv = splocResult.SBV;
splocResult = sploc(1,sbv,'eggHunt',traitFb,traitNb,vb);
sbv = splocResult.SBV;
splocResult = sploc(0,sbv,'eggHunt',traitFs,traitNs,vb);
%{
% ------------------------------------------------------------------ try 1
sbv0 = splocResult0.SBV;
splocResult1 = sploc(2,sbv0,'eggHunt',traitFa,traitNa,0);       % SPLOC it
%                    ^---> use default initial conditions for benchmarking
sbv = splocResult1.SBV;
splocResult1 = sploc(0,sbv,'eggHunt',traitFs,traitNs,0);
efficacy1 = splocResult1.efficacy;
% ------------------------------------------------------------------ try 2
splocResult2 = sploc(0,sbv0,'eggHunt',traitFs,traitNs,0);       % SPLOC it
efficacy2 = splocResult2.efficacy;
   if( efficacy1 > efficacy2 )
   splocResult = splocResult1;
   bias = 1;
   else
   splocResult = splocResult2;
   bias = 0;
   end
%}
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
% ----------------------------------------------------- determine loadings
splocLoadings = zeros(d,1);
   if( Dd > 0 )
   Ud = getDiscriminantSBV(splocResult);         % discriminant modes only
      for k=1:Dd
      splocLoadings = splocLoadings + Ud(:,k).^2;
      end
   end
wallTime2 = toc;
cpuSPLOC = cputime() - t0;
%%                                              append information to file
disp( dividerLine('egg hunt setup') );
disp(['  environment = ',GNtype]);
   if( placeEgg )
   disp('  egg present = true');
   disp(['           k1 = ',num2str(k1)]);
   disp(['           k2 = ',num2str(k2)]);
   else
   disp('  egg present = false');
   end
disp(['           df = ',num2str(dof)]);
disp(['          opv = ',num2str(opv)]);
disp( dividerLine('SPLOC results') );
%disp(['SPLOC-jSelect = ',num2str(iSelect)]);
%disp(['         bias = ',num2str(bias)]);
disp(['           Dd = ',num2str(splocDd)]);
disp(['           Du = ',num2str(splocDu)]);
disp(['           Di = ',num2str(dof - splocDd - splocDu)]);
disp(['        %eggD = ',num2str(splocPd)]);
disp(['        %eggU = ',num2str(splocPu)]);
disp(['        %eggI = ',num2str(100 - splocPd - splocPu)]);
disp(['         Epdf = ',num2str(splocEpdf)]);
disp( dividerLine('PCA results') );
disp(['  PCA-jSelect = ',num2str(jSelect)]);
disp(['           Dd = ',num2str(pcaDd)]);
disp(['           Du = ',num2str(pcaDu)]);
disp(['           Di = ',num2str(dof - pcaDd - pcaDu)]);
disp(['        %eggD = ',num2str(pcaPd)]);
disp(['        %eggU = ',num2str(pcaPu)]);
disp(['        %eggI = ',num2str(100 - pcaPd - pcaPu)]);
disp(['         Epdf = ',num2str(pcaEpdf)]);
disp( dividerLine('timing') );
disp(['setup cputime = ',num2str(cpuGenerate)]);
disp(['    walltime1 = ',num2str(wallTime1)]);
disp(['SPLOC cputime = ',num2str(cpuSPLOC)]);
disp(['    walltime2 = ',num2str(wallTime2)]);
disp(['total cputime = ',num2str(cpuGenerate + cpuSPLOC)]);
disp([' net walltime = ',num2str(wallTime1 + wallTime2)]);
% ----------------------------------------------------------- plot figures
figure(1);
clf;
plot(1:d,log10(lambdaF),'b');
hold on;
plot(1:d,log10(lambdaN),'r');
plot(1:d,log10(lambdaB),'k');
plot(1:d,log10(lambdaD),'g','linewidth',2.2);
xlabel('PCA mode rank');
ylabel('PCA mode variance');
legend('Functional','Nonfunctional','Pooled','Difference', ...
       'location','northeast');
plotCongruencySpectrum(splocResult,3);
plotCongruencySpectrum(splocResultsPCA,3);
figure(6)
clf;
hold on;
plot(1:d,pcaLoadings,'b');
plot(1:d,splocLoadings,'r');
xlabel('variables');
ylabel('loadings');
legend('PCA','SPLOC','location','northwest');
%%                                                   plot further analysis
%{
figure(10)
clf;
hold on;
iii = 1:d;
scatter(iii,pMode,60,'k','o');                             % missed target
scatter(iii(iMIN:iMAX),pMode(iMIN:iMAX),80,'m','o','filled'); % hit target                                      
xlabel('variable');
ylabel('egg likelihood');
title('PCA modes');
% ----------------------------------------------- show congruency spectrum
plotCongruencySpectrum(splocResultsPCA,3);
   if( Dd > 0 )
   Ud = getDiscriminantSBV(splocResultsPCA);     % discriminant modes only
   qF = xF*Ud;
   qN = xN*Ud;
   figure(11);
   clf;
   hold on;
      if( Dd == 1 )
      stdN = std(qN);
      stdF = std(qF);
      ts = 1:n;
         if( scaleFactor > 1 )
         plot(ts,qF,'bo');
         plot(ts,qN,'ro');
         else
         plot(ts,qN,'ro');
         plot(ts,qF,'bo');
         end
      xlabel('samples');
      ylabel('PCA projection 1');
      title(['stdDevF = ',num2str(stdF),'   stdDevN = ',num2str(stdN)]);
      else
         if( scaleFactor > 1 )
         plot(qF(:,1),qF(:,2),'bo');
         plot(qN(:,1),qN(:,2),'ro');
         else
         plot(qN(:,1),qN(:,2),'ro');
         plot(qF(:,1),qF(:,2),'bo');
         end
      xlabel('PCA projection 1');
      ylabel('PCA projection 2');
      end
% ------------------------------------------ show feature space clustering
   featureF = getFeatureVectors(traitFs,Ud);
   featureN = getFeatureVectors(traitNs,Ud);
% featureX. <-- data structure
% dataRefName = reference name for data with similar traits for sploc
% mMatrixName = cell array for file names that store the mMatrix data
% cMatrixName = cell array for file names that store the cMatrix data
%   nXsystems = number of systems being projected into feature space
%               X represents any system (0,1,unknown) that is ranked
%                 1-system => from training set labeled as functional
%                 0-system => from training set labeled as nonfunctional
%      nModes = # of discriminant modes contained in U.
%   nFeatures = 2xnModes = number of distinct features
%     Fmatrix = nFeatures x nXsystems 
   FmatrixF = featureF.Fmatrix;
   FmatrixN = featureN.Fmatrix;
   jStDv = 0;
         for k=1:featureF.nModes
         figure(20 + k);
         clf;
         hold on;
         jMean = jStDv + 1;
         jStDv = jMean + 1;
         % -----------------------------------
         x1 = FmatrixF(jMean,:);
         x2 = FmatrixF(jStDv,:);
         scatter(x1,x2,75,'filled','b');
         % -----------------------------------
         x1 = FmatrixN(jMean,:);
         x2 = FmatrixN(jStDv,:);
         scatter(x1,x2,75,'filled','r');
         xlabel('mean');
         ylabel('standard deviation');
         title(['Feature space for PCA mode: ',num2str(k)]);
         %pause
         end
   end
% -------------------------------------------------------- plot scree plot
figure(12)
clf;
yy = max( -12 , log10(lambda) );
xx = 1:d;
plot(xx,yy,'k','linewidth',2.0);
   if( Dd > 0 )
   hold on;
   pcaModeList = mapS2U(1:Dd);
   plot(pcaModeList, yy(pcaModeList),'ro');
   end
xlabel('PCA mode index');
ylabel('log10(lambda)');
%%                                                   find eggs using SPLOC
nTrials = 1;
splocName = [dataRefNameF,'sploc',num2str(nTrials,'%02i')];
% ------------------------------------------------------------------ logic
% start with initial conditions that are invariant w.r.t. rotation
% drive the initial basis set to get as many d-modes as possible.
% use normal scoring function to finalize the results.
% These 2 steps is essentially the same procedure previously hard coded
% Note: sploc() performs the above PCA calcualtion internally. Since it is
% done outside of sploc() to get the results as an intermediate step for
% PCA approach, and because it gives the same results, it is best to 
% start the initial condition as that result, fed into sploc() explicitly.
% --------------------------------------------------------------- SPLOC it
sbv = splocResultsPCA.SBV;
splocResult = sploc(2,sbv,splocName,traitFa,traitNa,vb);
sbv = splocResult.SBV;
splocResult = sploc(0,sbv,splocName,traitFs,traitNs,vb);
%%                                                    look at SPLOC output
Dd = splocResult.Dd; 
% ---------------------- calc. prob to be within the targeted 6 dimensions
pMode = zeros(d,1);
reconstruct2 = 0;
   if( Dd > 0 )
   Ud = getDiscriminantSBV(splocResult);         % discriminant modes only
      for k=1:Dd
      pMode = pMode + Ud(:,k).^2;
      reconstruct2 = reconstruct2 + sum( Ud(iMIN:iMAX,k).*V(:,k1) )^2;
      reconstruct2 = reconstruct2 + sum( Ud(iMIN:iMAX,k).*V(:,k2) )^2;
      % V = 0 outside of range:  iMIN:iMAX
      end
   end
reconstruct2 = round( reconstruct2*500 )/10;
disp(['PCA  results: %-reconstruct= ',num2str(reconstruct), ...
      '  nWrong= ',num2str(mWrong), ...
      '  overCount= ',num2str(mE1),'  underCount= ',num2str(mE2)]);
nCorrect = min( sum( pMode(iMIN:iMAX) ) , eggDIM);
nE2 = splocResult.Du + eggDIM - nCorrect;                 % under counting
nE1 = max( splocResult.Di - d + eggDIM, 0) ...             % over counting
    + Dd - nCorrect;  
nWrong = eggDIM - nCorrect;
disp(['SPLOC results: %-reconstruct= ',num2str(reconstruct2), ...
      '  nWrong = ',num2str(nWrong), ...
      '  overCount = ',num2str(nE1),'  underCount = ',num2str(nE2)]);
figure(30)
clf;
hold on;
iii = 1:d;
scatter(iii,pMode,60,'k','o');                                % off target
scatter(iii(iMIN:iMAX),pMode(iMIN:iMAX),80,'m','o','filled');  % on target
xlabel('variable');
ylabel('egg likelihood');
title('SPLOC modes');
% ----------------------------------------------- show congruency spectrum
plotCongruencySpectrum(splocResult,3);
   if( Dd > 0 )
   Ud = getDiscriminantSBV(splocResult);         % discriminant modes only
   qF = xF*Ud;
   qN = xN*Ud;
   figure(31);
   clf;
   hold on;
      if( Dd == 1 )
      stdN = std(qN);
      stdF = std(qF);
      ts = 1:n;
         if( scaleFactor > 1 )     
         plot(ts,qF,'bo');
         plot(ts,qN,'ro');
         else
         plot(ts,qN,'ro');
         plot(ts,qF,'bo');
         end
      xlabel('samples');
      ylabel('SPLOC projection 1');
      title(['stdDevF = ',num2str(stdF),'   stdDevN = ',num2str(stdN)]);
      else
         if( scaleFactor > 1 )
         plot(qF(:,1),qF(:,2),'bo');
         plot(qN(:,1),qN(:,2),'ro');
         else
         plot(qN(:,1),qN(:,2),'ro');
         plot(qF(:,1),qF(:,2),'bo');
         end
      xlabel('SPLOC projection 1');
      ylabel('SPLOC projection 2');
      end
% ------------------------------------------ show feature space clustering
   featureF = getFeatureVectors(traitFs,Ud);
   featureN = getFeatureVectors(traitNs,Ud);
% featureX. <-- data structure
% dataRefName = reference name for data with similar traits for sploc
% mMatrixName = cell array for file names that store the mMatrix data
% cMatrixName = cell array for file names that store the cMatrix data
%   nXsystems = number of systems being projected into feature space
%               X represents any system (0,1,unknown) that is ranked
%                 1-system => from training set labeled as functional
%                 0-system => from training set labeled as nonfunctional
%      nModes = # of discriminant modes contained in U.
%   nFeatures = 2xnModes = number of distinct features
%     Fmatrix = nFeatures x nXsystems 
   FmatrixF = featureF.Fmatrix;
   FmatrixN = featureN.Fmatrix;
   jStDv = 0;
         for k=1:featureF.nModes
         figure(40 + k);
         clf;
         hold on;
         jMean = jStDv + 1;
         jStDv = jMean + 1;
         % -----------------------------------
         x1 = FmatrixF(jMean,:);
         x2 = FmatrixF(jStDv,:);
         scatter(x1,x2,75,'filled','b');
         % -----------------------------------
         x1 = FmatrixN(jMean,:);
         x2 = FmatrixN(jStDv,:);
         scatter(x1,x2,75,'filled','r');
         xlabel('mean');
         ylabel('standard deviation');
         title(['Feature space for SPLOC mode: ',num2str(k)]);
         %pause
         end
   end
%}
%%          play
%{
Usploc = splocResult0id.SBV;
Upca = splocResultsPCA.SBV;
qPCAmode = zeros(d,1);
qSPLOCmode = zeros(d,1);
for k=90:100
qPCAmode   = qPCAmode   + Upca(:,k).^2;
qSPLOCmode = qSPLOCmode + Usploc(:,k).^2;
end
figure(200);
clf;
hold on;
plot(1:d,qSPLOCmode,'or');
plot(1:d,qPCAmode,'*b');
totalSPLOC = sum( qSPLOCmode(iMIN:iMAX) );
totalPCA = sum( qPCAmode(iMIN:iMAX) );
disp( [totalPCA, totalSPLOC] );
% ------------------------------------------------
pPCAmode = zeros(d,1);
pSPLOCmode = zeros(d,1);
for k=1:1
pPCAmode   = pPCAmode   + Upca(:,k).^2;
pSPLOCmode = pSPLOCmode + Usploc(:,k).^2;
end
figure(201);
clf;
hold on;
plot(1:d,pSPLOCmode,'or');
plot(1:d,pPCAmode,'*b');
totalSPLOC = sum( pSPLOCmode(iMIN:iMAX) );
totalPCA = sum( pPCAmode(iMIN:iMAX) );
disp( [totalPCA, totalSPLOC] );
% ------------------------------------------------
pPCAmode = zeros(d,1);
pSPLOCmode = zeros(d,1);
for k=1:2
pPCAmode   = pPCAmode   + Upca(:,k).^2;
pSPLOCmode = pSPLOCmode + Usploc(:,k).^2;
end
figure(202);
clf;
hold on;
plot(1:d,pSPLOCmode,'or');
plot(1:d,pPCAmode,'*b');
totalSPLOC = sum( pSPLOCmode(iMIN:iMAX) );
totalPCA = sum( pPCAmode(iMIN:iMAX) );
disp( [totalPCA, totalSPLOC] );
%}

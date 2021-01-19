% automate the Egg Hunt
% User specifies a job # from 1 to 18.
% 18 different batches can be run in parallel without any clashing.
% This will be a depth first filling of the data generation process.
% 
% INPUT:
% jobNumber     1 through  18       => runs the specified job
%              -1 through -18       => test run that bypasses the egg hunt
%
% PROCESS
% processEggHunt(Ntype, egg, dof, opv)
%
% OUTPUT
% print to single file
% ------------------------------------------------------- clear everything
clear all
%%                                                             specify job
job = 1;      % a number from 1 to 18    OR    -1 to -18 for testing only
% --------------------------------------------- specify system size series
ddAll = round(logspace(1,4,31));               % set dof (entire egg hunt)
%dd = ddAll;                                   % default is to do all sizes
%dd=ddAll(1:4);                                             % select range
%dd=ddAll(1:4:17);                                % pick desired selection
%dd=ddAll(31);                                    % pick desired selection
dd=[1000,2500];                                                % custom :)
nT=3;                                                   % number of trials
%%                                                    decompose job number
   if( job < 0 )
   testRun = true;
   job = -job;
   else
   testRun = false;
   end
% ------------------------------------------------------------ check scope
   if( (job < 1) || (job > 18) )
   error('job number is out of range!');
   end
% ----------------------------------------------- decompose the job number
% job = 9*B6 + 3*T3 + 1*T1;
%                       T1 = 1,2,3 for opv = 4, 20, 100 respectively
%                T3 = 0,1,2 for ee = Big, Small, None
%         B9 = 0,1 for nn = SGN, CGN
%  job    B9    T3    T1
%    1     0     0     1
%    2     0     0     2
%    3     0     0     3
%    4     0     1     1
%    5     0     1     2
%    6     0     1     3
%    7     0     2     1
%    8     0     2     2
%    9     0     2     3
% -----------------------
%   10     1     0     1
%   11     1     0     2
%   12     1     0     3
%   13     1     1     1
%   14     1     1     2
%   15     1     1     3
%   16     1     2     1
%   17     1     2     2
%   18     1     2     3
% ----------------------------------------- build job decomposition arrays
getB9 = zeros(1,18);
getT3 = zeros(1,18);
getT1 = zeros(1,18);
   for B9=0:1
     for T3=0:2
       for T1=1:3
       j = 9*B9 + 3*T3 + T1;
       getB9(j) = B9;
       getT3(j) = T3;
       getT1(j) = T1;
       end
     end
   end
%%                                                                set case
B9 = getB9(job);
T3 = getT3(job);
T1 = getT1(job);
% ------------------------------------------------- concealing environment
   switch B9
       case 0
       nn = "SGN";
       case 1
       nn = "CGN";
   end
% --------------------------------------------------------------- egg type
   switch T3
       case 0
       ee = "Big";
       case 1
       ee = "Small";
       case 2
       ee = "None";
   end
% ---------------------------------------------- observations per variable
   switch T1
       case 1
       oo = 4;
       case 2
       oo = 20;
       case 3
       oo = 100;
   end
tt = 1:nT;                                                        % trials
%%                                                      run egg hunt cases
% ------------------------------------------------------- initialize sploc
initializeSPLOC(1,'fName','FriedEgg','gType','png');
setDataMatrixFormat('notVectored',2);               % by-pass screen input
mType = 'cov';
% ------------------------------------------------------- sweep over cases
for dof = dd
   
   for Ntype = nn
       
      for egg = ee
       
         for opv = oo
              
            for trial = tt 
                
                if( testRun )                        % this is just a test
                str = ['dof= ',num2str(dof),'  ', ...
                       'job= ',num2str(job), ...
                       ' B9= ',num2str(B9), ...
                       ' T3= ',num2str(T3), ...
                       ' T1= ',num2str(T1)];
                disp(str);
                else                                   %=> do it for real!
                processEggHunt(Ntype, egg, opv, dof)
                end
                
            end
              
         end
          
      end
      
   end
   
end
% ----------------------------------------------- verify job decomposition
  if( testRun )
  disp('  ');
  disp('---------------------------------- decomposition of job numbers');
  j=1:18;
  disp( [j', getB9', getT3', getT1'] );
  disp('  ');
  end

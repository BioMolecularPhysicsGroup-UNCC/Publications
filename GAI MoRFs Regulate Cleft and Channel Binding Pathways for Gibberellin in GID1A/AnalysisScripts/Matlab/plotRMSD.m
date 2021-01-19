% FileList = dir(strcat('./VMD_output/egression*',finputData));
ActiveNames={"GA1minus", "GA3minus","GA4minus","APO","GA7minus","GA1216oxminus"};
InactiveNames={"GA9minus","GA34minus","GA41617oxminus","GA12minus","GA20minus","GA4MeO","GA8minus"};
%"GA1minus", "GA3minus","GA4minus","APO","GA7minus","GA1216oxminus"
% "GA9minus","GA34minus","GA41617oxminus","GA12minus","GA20minus","GA4MeO","GA8minus"
% "GID1A-GA3","GID1A-GA4","GID1A-GA1216oxminus","GID1A-GA7"
% "GID1A-GA4MeO","GID1A-GA12","GID1A-GA20"
skipInt=300;
avgTable=containers.Map();
hold on;
for i=1:length(ActiveNames)
if isequal( ActiveNames{i},  "GA3minus")
   RMSDdata=readtable(strcat(ActiveNames{i},"/VMD_output/out-MoRF3_rmsd.dat"));
   avgTable(ActiveNames{i})=median(RMSDdata.Var2);
   x=[];
   y=[];
   z=[];
   for i=skipInt:skipInt:length(RMSDdata.Var2)
   x(i/skipInt)=mean(RMSDdata.Var2(i-(skipInt-1):i));
   y(i/skipInt)=i/skipInt;
   z(i/skipInt)=std(RMSDdata.Var2(i-(skipInt-1):i));
   end
   errorbar(y,x,z,"color",'k');
   %plot(RMSDdata.Var1,RMSDdata.Var2,"color",'k');
elseif isequal( ActiveNames{i},  "")
   RMSDdata=readtable(strcat(ActiveNames{i},"/VMD_output/out-MoRF3_rmsd.dat"));
   avgTable(ActiveNames{i})=median(RMSDdata.Var2);
   x=[];
   y=[];
   z=[];
   for i=skipInt:skipInt:length(RMSDdata.Var2)
   x(i/skipInt)=mean(RMSDdata.Var2(i-(skipInt-1):i));
   y(i/skipInt)=i/skipInt;
   z(i/skipInt)=std(RMSDdata.Var2(i-(skipInt-1):i));
   end
   errorbar(y,x,z,"color",'b');
else
   RMSDdata=readtable(strcat(ActiveNames{i},"/VMD_output/out-MoRF3_rmsd.dat"));
   avgTable(ActiveNames{i})=median(RMSDdata.Var2);
   x=[];
   y=[];
   z=[];   
   for i=skipInt:skipInt:length(RMSDdata.Var2)
   x(i/skipInt)=mean(RMSDdata.Var2(i-(skipInt-1):i));
   y(i/skipInt)=i/skipInt;
   z(i/skipInt)=std(RMSDdata.Var2(i-(skipInt-1):i));
   end
   errorbar(y,x,z,"color",'g');
%plot(RMSDdata.Var1,RMSDdata.Var2,"color",'g');
end
end
for i=1:length(InactiveNames)
   if isequal( InactiveNames{i},  "GA12minus")
   RMSDdata=readtable(strcat(InactiveNames{i},"/VMD_output/out-MoRF3_rmsd.dat"));
   avgTable(InactiveNames{i})=median(RMSDdata.Var2);
   x=[];
   y=[];
   z=[];
   for i=skipInt:skipInt:length(RMSDdata.Var2)
   x(i/skipInt)=mean(RMSDdata.Var2(i-(skipInt-1):i));
   y(i/skipInt)=i/skipInt;
   z(i/skipInt)=std(RMSDdata.Var2(i-(skipInt-1):i));
   end
   errorbar(y,x,z,"color",'m');
   %plot(RMSDdata.Var1,RMSDdata.Var2,"color",'k');
   elseif isequal( InactiveNames{i},  "")
   RMSDdata=readtable(strcat(InactiveNames{i},"/VMD_output/out-MoRF3_rmsd.dat"));
   avgTable(InactiveNames{i})=median(RMSDdata.Var2);
   x=[];
   y=[];
   z=[];
   for i=skipInt:skipInt:length(RMSDdata.Var2)
   x(i/skipInt)=mean(RMSDdata.Var2(i-(skipInt-1):i));
   y(i/skipInt)=i/skipInt;
   z(i/skipInt)=std(RMSDdata.Var2(i-(skipInt-1):i));
   end
   errorbar(y,x,z,"color",'b');
   else
   RMSDdata=readtable(strcat(InactiveNames{i},"/VMD_output/out-MoRF3_rmsd.dat"));
   avgTable(InactiveNames{i})=median(RMSDdata.Var2);
   x=[];
   y=[];
   z=[];   
   for i=skipInt:skipInt:length(RMSDdata.Var2)
   x(i/skipInt)=mean(RMSDdata.Var2(i-(skipInt-1):i));
   y(i/skipInt)=i/skipInt;
   z(i/skipInt)=std(RMSDdata.Var2(i-(skipInt-1):i));
   end
   errorbar(y,x,z,"color",'r');
%plot(RMSDdata.Var1,RMSDdata.Var2,"color",'g');
   end
end
title("RMSD Prod. Run MoRF 3");
xlabel(strcat("frames/",num2str(skipInt)));
ylabel("RMSD");
legend([ActiveNames,InactiveNames],'Location','northeastoutside');
hold off;
tempTable=table(avgTable.keys',avgTable.values');

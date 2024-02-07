% Playing around with MOR density results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add paths
baseloc="C:\Users\dk818\OneDrive - Imperial College London\NCORE\miFCdFC";

% For outputs
str='\Outputs';
outputpath=strcat(baseloc,str);
addpath(outputpath);

% For Inputs
str='\AtlasAndOtherInputs';
atlaspath=strcat(baseloc,str);
addpath(atlaspath);

% For miFC
str='\Functions\mi';
addpath(strcat(baseloc,str));

% For ICC, violinplot, and colors
str='\Functions';
addpath(strcat(baseloc,str));

%% 
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of sig res in VMN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Split results into edges w/ higher FC in MD vs HC pars, and
vice versa
GroupResTbl should be the go-to - cols 1 and 2 are the roi nums, col 3 =
pval, col 4 = zstat. Pos = MD>HC, neg = HC>MD. 
%}
% Ok, can see many more ROIs w/ higher miFC in HC as compared to MD
MD=GroupResTbl(GroupResTbl(:,4)>0,:);
HC=GroupResTbl(GroupResTbl(:,4)<0,:);
VMNROIs=[51,52,80,81,82,157,158,159,160,187,188,189,201,203,206,207,209,210,213,214];

MDLongCol=[MD(:,1);MD(:,2)];
MDNumVMNROIs=ismember(MDLongCol,VMNROIs);
disp(strcat('MD>HC: There are ',num2str(sum(MDNumVMNROIs)),'VMN ROIs out of ',num2str(size(MD,1)),'edges'))
disp(strcat('Put another way, ',num2str((sum(MDNumVMNROIs)/size(MD,1)*100)) ,'% of edges higher in MD pars contain a VMN ROI'))

HCLongCol=[HC(:,1);HC(:,2)];
HCNumVMNROIs=ismember(HCLongCol,VMNROIs);
disp(strcat('HC>MD: There are ',num2str(sum(HCNumVMNROIs)),'VMN ROIs out of ',num2str(size(HC,1)),'edges'))
disp(strcat('Put another way, ',num2str((sum(HCNumVMNROIs)/size(HC,1)*100)) ,'% of edges higher in HC pars contain a VMN ROI'))
% Interesting - despite there being a greater number of sig higher edges in
% HC as compared to MD pars, those that are higher in MD pars contain a
% greater number of VMN ROIs


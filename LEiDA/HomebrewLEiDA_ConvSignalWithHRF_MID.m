%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What I need to do - convolve each person's task timeseries with the HRF.
% I can read in each peron's EV text file, make a boxcar function with it,
% convolve it with the HRF, then look at brain states during periods of
% task. 
addpath('C:\Users\dk818\OneDrive - Imperial College London\NCORE\miFCdFC\AtlasAndOtherInputs\EVs');
addpath('C:\Users\dk818\OneDrive - Imperial College London\LilThingsAndMATLABPrograms\spm12')
baseloc="C:\Users\dk818\OneDrive - Imperial College London\NCORE\miFCdFC";

% Add paths
currloc=string(pwd);
homeDir=("C:\Users\dk818\OneDrive - Imperial College London\NCORE\miFCdFC");

% For inputs
str='\AtlasAndOtherInputs';
addpath(strcat(homeDir,str));

str='\Outputs';
outputpath=strcat(homeDir,str);
addpath(outputpath);

% For LEiDA
str='\LEiDA';
addpath(strcat(homeDir,str));

% For ICC, violinplot, colors, Dunn score, and more
str='\Functions';
addpath(strcat(homeDir,str));

%%
load df_forLEiDA.mat
load templateEigAndMat.mat
load HomeBrewLEiDA_MID.mat 
num_states=size(templateEig,2);

df=df_mid;
n_Subjects=22;
[N_areas, Tmax]=size(df{1,1});
num_condi=n_Task; % number of experimental conditions
TR=1.5;
NumVols=size(df{1,1},2);
TMax=NumVols*TR;  
timeser=zeros(TMax,1);
canonical_hrf=spm_hrf(TR);
n_sess=2;
DF={};

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify RA trials per subj 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MD_PLC
str='\AtlasAndOtherInputs\EVs\MD_PLC_MID';
addpath(strcat(baseloc,str));
pth=strcat(baseloc,str);
cd (pth)
files=dir('*txt');

for subj=1:size(files,1)
% read in EV 
tmpdf=readtable(files(subj).name);

% Now I need to use the EV to find the onset of each win anticipation
% period, and for the length specified, add 1 to those values
for ii=1:size(tmpdf,1)
    start=table2array(tmpdf(ii,1));
    fin=start+table2array(tmpdf(ii,2));
    timeser(start:fin,1)=timeser(start:fin)+1;
end

hrf_convolved_with_stim_time_series(subj,:) = conv(timeser,canonical_hrf);

% Sweeeeeet! And for our ease of use, let's now binarise this for all
% seconds greater than 0
bin_hrf(subj,:)=hrf_convolved_with_stim_time_series(subj,:)>0;

% And put this in TR format
bin_hrf_inTR{1}(subj,:)=bin_hrf(subj,1:1.5:end);

% Nice!!! So what I could do is pull out all the brain states where reward
% anticipation is, and see if there's a different probability of occurrence
% during that period
% In order to do that in a way that reliably matches the df used in LEiDA,
end

% HC_PLC
str='\AtlasAndOtherInputs\EVs\HC_PLC_MID';
addpath(strcat(baseloc,str));
pth=strcat(baseloc,str);
cd (pth)
files=dir('*txt');

for subj=1:size(files,1)
tmpdf=readtable(files(subj).name);
for ii=1:size(tmpdf,1)
    start=table2array(tmpdf(ii,1));
    fin=start+table2array(tmpdf(ii,2));
    timeser(start:fin,1)=timeser(start:fin)+1;
end
hrf_convolved_with_stim_time_series(subj,:) = conv(timeser,canonical_hrf);
bin_hrf(subj,:)=hrf_convolved_with_stim_time_series(subj,:)>0;
bin_hrf_inTR{2}(subj,:)=bin_hrf(subj,1:1.5:end);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract RA parts of state timeseries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s=1:n_Subjects
for task=1:num_condi
Ctime=StateTimeDF{s,task}(:,1);
tmpdf=bin_hrf_inTR{task}(s,:)';
% Then remove edges
fin=size(tmpdf,1);
tmpdf([1:5,fin-4:fin],:)=[];
% Make sure tmpdf is the same time as the Ctime series
cutoff=size(Ctime,1);
if size(tmpdf,1) > cutoff
    tmpdf(cutoff:end,:)=[];
end
% Now only keep vols where RA occurred
Ctime=Ctime(tmpdf,1);
RAStateTimeDF{s,task}(:,1)=Ctime;
end
end

% For ease, switch back naming
StateTimeDF=RAStateTimeDF;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 - Analyse the clustering results between states
% This assess the lifetime (in seconds) and probability each state will occur 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=zeros(n_Task,n_Subjects,num_states); 
LT=zeros(n_Task,n_Subjects,num_states); 

for state=1:num_states
for task=1:n_Task   % for each condition
for s=1:n_Subjects % for each subject
     for c=1:num_states
        Ctime_bin=StateTimeDF{s,task}(:,1)==c;
        P(task,s,c)=mean(Ctime_bin);   
        LT(task,s,c)=sum(Ctime_bin)*TR;
    end                
end
end
end
% save('Checkpoint2_LEiDA_MID.mat','-v7.3')
% save('Checkpoint2_LEiDA_CUE.mat','-v7.3')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyse the effect of condition on states 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%e%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:num_states
% For LT
LTStatePerTask=squeeze(LT(:,:,ii));
[LT_P(ii,1),~,stats]=signrank(LTStatePerTask(1,:),LTStatePerTask(2,:));
if isfield(stats,'zval')==1
LT_ETA(ii,1)=stats.zval;
else
LT_ETA(ii,1)=0;
end

% For prob
PStatePerTask=squeeze(LT(:,:,ii));
[Prob_P(ii,1),~,stats]=signrank(PStatePerTask(1,:),PStatePerTask(2,:));
if isfield(stats,'zval')==1
Prob_ETA(ii,1)=stats.zval;
else
Prob_ETA(ii,1)=0;
end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDR correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,LT_P_FDR(:,1)] = fdr(LT_P(:,1));
[~,~,Prob_P_FDR(:,1)] = fdr(Prob_P(:,1));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Organise Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We want p, then eta, for Prob and Lt, then all complexity 
StandardRes=[LT_P_FDR,LT_ETA,Prob_P_FDR,Prob_ETA];

%%
% cd(outputpath)
% save('HomeBrewLEiDA_MID.mat','-v7.3')
% save('HomeBrewLEiDA_CUE.mat','-v7.3')
save('HomeBrewLEiDA_MID_Eig.mat','-v7.3')
% Just a note - state order, 1-9, is the same for Eigs and Mats. is Vis,
% Somat, DorsAttn, SalAttn, Lim, Cont, DMN, Tmp, VMN
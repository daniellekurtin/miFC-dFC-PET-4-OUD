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

% For data
str='\ExtractedTimeseries';
addpath(strcat(baseloc,str));

%%
%{
Organise data from how it's extracted from fslmeants to a friendlier format
This will load each condition, combine cortex and subcortex files, and put
it all in a BigStruct
%}

load('df_forLEIDA.mat')
df=df_mid;
% df=df_cue;
[n_Subjects, n_Task]=size(df);
n_Subjects=22;
n_subj=n_Subjects;
TR=1.5;
NumVols=size(df_mid{1,1},2);
TMax=NumVols*TR;  
timeser=zeros(TMax,1);
canonical_hrf=spm_hrf(TR);
n_sess=2;
DF={};

%%
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
end

%%
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
% ADDED SECTION FOR NORMALISATION OF TIMESERIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zscore everything to have centre mean 0 and 1 SD
for sess=1:n_sess
for subj=1:n_subj
for rr=1:num_rois
    BigStructNorm{subj,sess}(:,rr)=zscore(df{subj,sess}(rr,:)');
end
end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute mi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for sess=1:n_sess
for subj=1:n_subj
    
% This is this subj and session's miFC for each region (col) across time
% (rows)
% We only want rows where RA occurrs
tmp=BigStructNorm{subj,sess};

% Now I need to get the chunks of reward anticipation - first
% pull out this subj's data
tmpdf=bin_hrf_inTR{sess}(subj,:)';

% Make sure tmpdf is the same size
cutoff=size(tmp,1);
if size(tmpdf,1) > cutoff
    tmpdf(cutoff+1:end,:)=[];
end

% Now only keep vols where RA occurred
tmp=tmp(tmpdf,:);

for rr1=1:num_rois
    for rr2=1:num_rois
        MI{subj,sess}(rr1,rr2) = mutualinfo(tmp(:,rr1),tmp(:,rr2));
    end
end  

end
end

% Since this is an expensive computational bit, let's make a checkpoint
cd (outputpath)
save("Checkpoint1_ComputeMI_RA.mat",'-v7.3')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assess reliability of each region's MI across subject per session
% Helpful link: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4913118/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ComboNames=readtable('ComboNames.xlsx');

for sesh=1:n_sess
for rr=1:num_rois

% Each col is a subj, and each row is the miFC of a region to region rr 
for subj=1:n_subj
ROIMI(:,subj)=MI{subj,sesh}(rr,:)';
end

% This is checking the reliability of interregional miFC across subjs for
% each condition
[BetSubjRelR(sesh,rr), ~, ~, ~, ~, ~, BetSubjRelP(sesh,rr)] = ICC(ROIMI,'A-1');
end
end 

BetSubjRelR=BetSubjRelR';

% Now we want to get average ICC
for ii=1:n_sess
meanICC(ii,1)=mean(BetSubjRelR(:,ii));
meanICC(ii,2)=std(BetSubjRelR(:,ii));
% Adding median
meanICC(ii,3)=median(BetSubjRelR(:,ii));
end

% Compute number of outliers -  what % are outside 3x SD
for nn=1:n_sess
    benchmark=meanICC(nn,1)-(3*meanICC(nn,2));
    outliercount=0;
    h=0;
    for ii=1:length(BetSubjRelR)
    if BetSubjRelR(ii,nn)<benchmark
        h=h+1;
        outliercount=outliercount+1;
        nameRec(h,nn)=SchaeferNames(ii,2);
    end
    end
    OutlierCount(1,nn)=outliercount;
end

C=linspecer(2);
figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
violins = violinplot(BetSubjRelR);
% names={'MD_PLC','MD_APR','HC_PLC','HC_AP'};
% xticklabels=names;
for ii=1:n_sess
violins(1,ii).ShowData = 1;
violins(1,ii).ShowMean = 'yes';
violins(1,ii).ViolinColor = C(ii,:);
end
hold off
saveas(gcf,'RA_MI_ICCScoreViolins.png');
close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the effect of condition on miFC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For each pair of ROIs . . .
for rr1=1:num_rois
for rr2=1:rr1
%  . . . get a table of the pairwise miFC between for each subj (row) for each session (col)    
for subj=1:n_subj
for sesh=1:n_sess
ROIMI(subj,sesh)=MI{subj,sesh}(rr1,rr2);
end
end
[p,~,stats]=signrank(ROIMI(:,1),ROIMI(:,2));
if p < 0.05 
GroupSigP(rr1,rr2)=p;
if isfield(stats,'zval')==1
GroupSigETA(rr1,rr2)=stats.zval;
else
GroupSigETA(rr1,rr2)=0;
end
end
clear ROIMI
end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the FDR-corrected pvals and save Eta vals 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GroupSigP_FDR=zeros(num_rois,num_rois);
maxw=size(GroupSigP,2);
for ii=1:length(GroupSigP)
if ii<maxw
[~,~,tmpP]=fdr(GroupSigP(ii,1:ii));
tmpETA=GroupSigETA(ii,1:ii);
pad=zeros(num_rois-size(tmpP,2),1)';
GroupSigP_FDR(ii,:)=[tmpP,pad];
GroupSigETASurvive(ii,:)=[tmpETA,pad];
else
[~,~,tmpP]=fdr(GroupSigP(ii,1:end));
pad=zeros(num_rois-size(tmpP,2),1)';
GroupSigP_FDR(ii,:)=[tmpP,pad];    
tmpETA=GroupSigETA(ii,1:end);
GroupSigETASurvive(ii,:)=[tmpETA,pad];
end
end

GroupSigETASurvivepadded=GroupSigETASurvive;
clear GroupSigETASurvive


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the FDR-corrected pvals and save Eta vals 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brilliant! Now to see how many sig effects we have.
GroupSigCountSurviveC4MC=0;
for ii=1:num_rois
    for jj=1:ii
        if GroupSigP_FDR(ii,jj)<0.05 && GroupSigP_FDR(ii,jj)>0
            GroupSigCountSurviveC4MC=GroupSigCountSurviveC4MC+1;
            GroupSigETASurvive(ii,jj)=GroupSigETA(ii,jj);
        end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN EFFECT OF DRUG 
% Compute post hocs for effects suriving corr for mult comp
% This will show the direction and size of effect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GroupPostHocP={};
GroupPostHocPSig={};
GroupPostHocStats={};
GroupPostHocStatsSig={};
Groupposthocsigcount=zeros(4,4);
h=0;

% For each ROI pair. . . 
for rr1=1:num_rois
for rr2=1:rr1
    if GroupSigP_FDR(rr1,rr2)<0.05 && GroupSigP_FDR(rr1,rr2)>0
    
        % Make the main results table
        % Schaefer ROI label 1, Schaefer ROI label 2, 
        % main p val, main effect size, 
        % and for each of the comparisons, a pval and 
        % effect size.
        h=h+1;
        GroupResTblName(h,1)=ComboNames(rr1,1);
        GroupResTblName(h,2)=ComboNames(rr2,1);
        GroupResTbl(h,1)=rr1;
        GroupResTbl(h,2)=rr2;
        GroupResTbl(h,3)=GroupSigP_FDR(rr1,rr2);
        GroupResTbl(h,4)=GroupSigETASurvive(rr1,rr2);
    end   
end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ok, to plot stuff - I would like to have a circle graph with diff
% hemispheres on the top and bottom of the plot, and sections of the plot
% colored according to network

% % To make circgraph with the proper alignment, I've got to reorganize the
% % data and the labels. 
% SchaeferNames=readtable('Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv');
% SchaeferLabels=SchaeferNames(:,2);
% SchaeferLabels2(1:(num_rois/2),1)=SchaeferLabels(1:(num_rois/2),1);
% Flipped=flip(SchaeferLabels((num_rois/2)+1:num_rois,1));
% SchaeferLabels2((num_rois/2)+1:num_rois,1)=Flipped;
% Flipped=[];
% 
% % Cool, now just gotta flip results the same way
% PartEta2Survive(end:num_rois,end:num_rois)=0;
% PartEta2SurviveFlip=PartEta2Survive;
% PartEta2SurviveFlip((num_rois/2)+1:num_rois,1:(num_rois/2))=flip(PartEta2SurviveFlip((num_rois/2)+1:num_rois,1:(num_rois/2)));
% PartEta2SurviveFlip(1:(num_rois/2),(num_rois/2)+1:num_rois)=flip(PartEta2SurviveFlip(1:(num_rois/2),(num_rois/2)+1:num_rois),2);
% PartEta2SurviveFlip((num_rois/2)+1:num_rois,(num_rois/2)+1:num_rois) = flip(flip(PartEta2SurviveFlip((num_rois/2)+1:num_rois,(num_rois/2)+1:num_rois),2));

% % Now I have to make a color map where the first 50 vals map
% C=linspecer(17);
% load Clist
% for ii=1:height(Clist)
%    colr(1,:)= C(Clist(ii,1),:);  
%    Clist2(ii,:)=colr ;
% end

% % Plot some sheeeet
% figure()
% circularGraph(GroupSigETASurvive)
% 
% figure()
% circularGraph(DrugSigETASurvive)
% 
% figure()
% circularGraph(GrpDrgSigETASurvive)
% 

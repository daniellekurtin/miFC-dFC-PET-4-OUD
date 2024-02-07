% Helpful link = http://web.mit.edu/hst.583/www/LABS/lab5b_manual.pdf

% What I need to do - convolve each person's task timeseries with the HRF.
% I can read in each peron's EV text file, make a boxcar function with it,
% convolve i
% t with the HRF, then look at brain states during periods of
% task. 
addpath('C:\Users\dk818\OneDrive - Imperial College London\NCORE\miFCdFC\AtlasAndOtherInputs\EVs');
addpath('C:\Users\dk818\OneDrive - Imperial College London\LilThingsAndMATLABPrograms\spm12')
baseloc="C:\Users\dk818\OneDrive - Imperial College London\NCORE\miFCdFC";

load 'Checkpoint3_LEiDA_MID.mat'
% load 'df_forLEida.mat';
% df=df_mid;

% We know the task is 292 frames long, and each frame is 1.5 seconds, so
% total seconds = 292 * 1.5
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
% In order to do that in a way that reliably matches the df used in LEiDA,
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
% Time to see the number of states during each period, and whether that's
% affected by condition!

RA_P=zeros(n_Task,n_Subjects,maxk-1,maxk); 
RA_LT=zeros(n_Task,n_Subjects,maxk-1,maxk); 

for k=mink-1:maxk-1   % for each cluster
    for task=1:n_Task   % for each condition
        for s=1:n_Subjects % for each subject
            
            % Select the time points representing this subject and task               
            T=(Time_all(1,:)==s&Time_all(2,:)==task); 
            Ctime=Kmeans_results{k}.IDX(T); % The cluster each timecourse component belongs to

            % Now I need to get the chunks of reward anticipation - first
            % pull out this subj's data
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
            
             for c=1:k+1  
                RA_P(task,s,k,c)=mean(Ctime==c);   
                Ctime_bin=Ctime==c;
                a=find(diff(Ctime_bin)==1); 
                b=find(diff(Ctime_bin)==-1);
                
                % We discard the cases where state starts or ends ON
                if length(b)>length(a)
                    b(1)=[];
                elseif length(a)>length(b)
                    a(end)=[];
                elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
                    b(1)=[];
                    a(end)=[];
                end

                if ~isempty(a) && ~isempty(b)
                    C_Durations=b-a;
                else
                    C_Durations=0;
                end
                RA_LT(task,s,k,c)=mean(C_Durations)*TR;
            end                
        end
    end
end   

%%
% Compute the effect of group, drug, and group and drug on LT and P during
% RA periods
clear ind_sort

for KK=mink:maxk

Best_Clusters=Kmeans_results{KK-1};
ProbC=zeros(1,KK);

for c=1:KK
    ProbC(c)=mean(Best_Clusters.IDX==c);
end

[~, ind_sort{KK}(1,:)]=sort(ProbC,'descend'); 

end

% Reminder- LT(num_coni, num_subj, k, num_clus_in_K
for KK=1:maxk-1
    for ii=1:KK+1
        
        % For LT 
        LTStatePerTask=squeeze(RA_LT(:,:,KK,ind_sort{KK+1}(1,ii)))';
        [RA_LT_P{KK}(ii,1),~,stats]=signrank(LTStatePerTask(:,1),LTStatePerTask(:,2));
        if isfield(stats,'zval')==1
        RA_LT_ETA{KK}(ii,1)=stats.zval;
        else
        RA_LT_ETA{KK}(ii,1)=0;
        end

        % % For prob 
        PStatePerTask=squeeze(RA_P(:,:,KK,ind_sort{KK+1}(1,ii)))';
        [RA_Prob_P{KK}(ii,1),~,stats]=signrank(PStatePerTask(:,1),PStatePerTask(:,2));
        if isfield(stats,'zval')==1
        RA_Prob_ETA{KK}(ii,1)=stats.zval;
        else
        RA_Prob_ETA{KK}(ii,1)=0;
        end

    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post hocs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % P(task,subj,k,c)
% for KK=mink-1:maxk-1
% for cc=1:KK+1
% PStatePerTask=squeeze(RA_P(:,:,KK,ind_sort{KK+1}(1,cc)))';
% LTStatePerTask=squeeze(RA_LT(:,:,KK,ind_sort{KK+1}(1,cc)))';
% 
% for ii=1:n_Task
% for jj=1:ii
% 
%         % Non-parametric - Wilcoxon Signed Rank 
%         xp=PStatePerTask(:,ii);
%         yp=PStatePerTask(:,jj);
%         RA_Prob_Posthoc{KK,cc}(ii,jj)=signrank(xp,yp);
% 
%         if RA_Prob_Posthoc{KK,cc}(ii,jj)<0.05
%             disp(strcat('K=',num2str(KK+1),'s=',num2str(cc),'stimCond',num2str(ii),'vs stimCond',num2str(jj)))
%         end
% 
%         xlt=LTStatePerTask(:,ii);
%         ylt=LTStatePerTask(:,jj);
%         RA_LT_Posthoc{KK,cc}(ii,jj)=signrank(xlt,ylt);  
% 
%         if RA_LT_Posthoc{KK,cc}(ii,jj)<0.05
%             disp(strcat('K=',num2str(KK+1),'s=',num2str(cc),'stimCond',num2str(ii),'vs stimCond',num2str(jj)))
%         end
% end
% end
% 
% end
% end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDR correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:maxk-1
    [~,~,RA_LT_P_FDR{ii}(:,1)] = fdr(RA_LT_P{ii}(:,1));
    [~,~,RA_Prob_P_FDR{ii}(:,1)] = fdr(RA_Prob_P{ii}(:,1));
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Organising outputs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Put in nice table
RA_Prob_Fancy={};
RA_LT_Fancy={};

for k=mink-1:maxk-1
for c=1:k+1

tmp=strcat(num2str(RA_Prob_P_FDR{k}(c,1)),', ',num2str(round(RA_Prob_ETA{k}(c,1),2,'decimal')));
RA_Prob_Fancy{1}(c,k)=convertCharsToStrings(tmp);

tmp=strcat(num2str(RA_LT_P_FDR{k}(c,1)),', ',num2str(round(RA_LT_ETA{k}(c,1),2,'decimal')));
RA_LT_Fancy{1}(c,k)=convertCharsToStrings(tmp);
end
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot pvals for LT and Prob
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=linspecer(14);

% Probability
figure()
hold on
for kk=mink-1:maxk-1
    clear xx
    for jj=1:kk+1
    	xx(jj,1)=jj;
    end
    scatter(xx,log10(RA_Prob_P_FDR{1,kk}(:,1)),30,C(kk,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
%     swarmchart(xx,log10(Prob_P_FDR{1,kk}(:,1)),30,C(kk,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
end
yline(log10(0.05))
alpha(0.7)
% name='RA_SwarmChatPvals_Prob_Result.png';
% saveas(gcf,name);
% close all

% LT
figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
for kk=mink-1:maxk-1
    clear xx
    for jj=1:kk+1
    	xx(jj,1)=jj;
    end
     scatter(xx,log10(RA_LT_P_FDR{1,kk}(:,1)),30,C(kk,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
%         swarmchart(xx,log10(LT_P_FDR{1,kk}(:,jj)),30,C(kk,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
end
yline(log10(0.05))
alpha(0.7)
% name='RA_SwarmChatPvals_LT_Result.png';
% saveas(gcf,name);
% close all



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMENTED THIS OLD VERSION OUT 
% NEXT SECTION SHOULD WORK
% Make violins for LT and P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % For LT
% % MD_PLC=squeeze(RA_LT(1,:,k-1,1:k));
% % MD_APR=squeeze(RA_LT(2,:,k-1,1:k));
% % HC_PLC=squeeze(RA_LT(3,:,k-1,1:k));
% % HC_APR=squeeze(RA_LT(4,:,k-1,1:k));
% 
% % For Prob
% MD_PLC=squeeze(RA_P(1,:,k-1,1:k));
% HC_PLC=squeeze(RA_P(2,:,k-1,1:k));
% 
% ind_sort2=[];
% for KK=mink:maxk
% Best_Clusters=Kmeans_results{rangeK==KK};
% ProbC=zeros(1,K);
% for c=1:K
%     ProbC(c)=mean(Best_Clusters.IDX==c);
% end
% [~, ind_sort2(KK,:)]=sort(ProbC,'descend'); 
% end
% 
% State={};
% for ii=1:K
%    State{ii}(:,1)=MD_PLC(:,ind_sort2(K,ii));
%    State{ii}(:,3)=HC_PLC(:,ind_sort2(K,ii));
% end
% 
% C=linspecer(4);
% 
% for cc=1:K
% figure()
% hold on
% set(gca,'FontSize',20)
% set(gca,'FontName','Arial')
% tmp=State{cc};   % get the LT or the Prob
% violins = violinplot(tmp);
% 
% for ii=1:num_condi
% violins(1,ii).ShowData = 1;
% violins(1,ii).ShowMean = 'yes';
% violins(1,ii).ViolinColor = C(ii,:);
% end
% hold off
% 
% name=strcat('RA_ProbViolins_State',num2str(cc),'.png');
% % name=strcat('RA_LTViolins_State',num2str(cc),'.png');
% saveas(gcf,name);
% close all
% end

%%
% For LT
% MD_PLC=squeeze(LT(1,:,k-1,1:k));
% HC_PLC=squeeze(LT(2,:,k-1,1:k));

% For Prob
MD_PLC=squeeze(P(1,:,k-1,1:k));
HC_PLC=squeeze(P(2,:,k-1,1:k));

ind_sort2=[];
% for KK=mink:maxk
% Best_Clusters=Kmeans_results{rangeK==KK};
% ProbC=zeros(1,K);
% for c=1:K
%     ProbC(c)=mean(Best_Clusters.IDX==c);
% end
% [~, ind_sort2(KK,:)]=sort(ProbC,'descend'); 
% end

State={};
% for ii=1:K
%    State{ii}(:,1)=MD_PLC(:,ind_sort2(K,ii));
%    State{ii}(:,2)=HC_PLC(:,ind_sort2(K,ii));
% end
for ii=1:K
   State{ii}(:,1)=MD_PLC(:,ii);
   State{ii}(:,2)=HC_PLC(:,ii);
end

C=linspecer(2);

for cc=1:K
figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
tmp=State{cc};   % get the LT or the Prob
violins = violinplot(tmp);

for ii=1:num_condi
violins(1,ii).ShowData = 1;
violins(1,ii).ShowMean = 'yes';
violins(1,ii).ViolinColor = C(ii,:);
end
hold off
% 
% % name=strcat('CUE_ProbViolins_State',num2str(cc),'.png');
% name=strcat('MID_ProbViolins_State',num2str(cc),'.png');
% % name=strcat('CUE_LTViolins_State',num2str(cc),'.png');
% % name=strcat('MID_LTViolins_State',num2str(cc),'.png');
% saveas(gcf,name);
% close all


name=strcat('RA_ProbViolins_State',num2str(cc),'.png');
% name=strcat('RA_LTViolins_State',num2str(cc),'.png');
saveas(gcf,name);
close all
end

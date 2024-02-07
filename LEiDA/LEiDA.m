
%{
LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)
This function processes, clusters and analyses BOLD data using LEiDA in the following order:
1. Read the BOLD data from the folders and computes the BOLD phase 
- Calculate the instantaneous BOLD synchronization matrix
- Compute the Leading Eigenvector at each frame from all fMRI scans
2. Cluster the Leading Eigenvectors for k=2:12. This can be adapted via the
'mink' and 'maxk' vars
3. Compute the probability and lifetimes each cluster in each session
- Calculate signigifance between task
- Saves the Eigenvectors, Clusters and statistics into LEiDA_results.mat
4. Plots FC states and errorbars for each clustering solution
- Adds an asterisk when results are significantly different between blocks

Created by Joana Cabral, Oct 2017, joana.cabral@psych.ox.ac.uk
First use in: Cabral, et al. 2017 Scientific reports 7, no. 1 (2017): 5135.

Adapted by Danielle Kurtin, Oct 2021, d.kurtin@surrey.ac.uk 

%}

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add paths
currloc=string(pwd);
homeDir=("C:\Users\dk818\OneDrive - Imperial College London\NCORE\miFCdFC");

% For outputs
str='\Outputs';
% str='\Outputs\EdgeRemoved_CUE';
outputpath=strcat(homeDir,str);
addpath(outputpath);

% For LEiDA
str='\LEiDA';
addpath(strcat(homeDir,str));

% For ICC, violinplot, colors, Dunn score, and more
str='\Functions';
addpath(strcat(homeDir,str));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEiDA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load df_forLEiDA.mat

%%%%%%%%%%%%%
% NOTE TO SEELLLLFFFFFFFF - MUST CHANGE THIS LINE DEPENDING WHETHER YOU
% WANT TO RUN CUE OR NOT!!!!!!!!!
%%%%%%%%%%%%
% df=df_mid;
df=df_cue;

[n_Subjects, n_Task]=size(df);
n_Subjects=22;
[N_areas, Tmax]=size(df{1,1});
num_condi=n_Task; % number of experimental conditions
TR=1.5;

%%%%%%%%%%%%%%%%% DOING THIS TO REMOVE EDGE EFFECTS
% Leading_Eig=zeros(Tmax*n_Subjects,N_areas); % All leading eigenvectors- creates (29 frames*18 runs=)529 rows, and 166 columns.
% Time_all=zeros(2, n_Subjects*Tmax); % vector with subject nr and task at each t
Leading_Eig=zeros((Tmax-10)*n_Subjects,N_areas); % All leading eigenvectors- creates (29 frames*18 runs=)529 rows, and 166 columns.
Time_all=zeros(2, n_Subjects*(Tmax-10)); % vector with subject nr and task at each t

t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)
fnq=1/(2*TR);                 % Nyquist frequency.
flp = .02;                    % lowpass frequency of filter (Hz). This allows anything below 50 Hz. 
fhi = 0.1;                    % highpass- we've already applied one, though it is less inclusive (allowing anything over 0.01 to pass). I will keep this for now. 
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency. 
k=2;                          % 2nd order butterworth filter. This determines the steepness of the gain function (and 2 is pretty smooth). 
[bfilt,afilt]=butter(k,Wn);   % "Butter" is a MATLAB function than constructs the butterworth filter using defined cutoff frequencies.

for s=1:n_Subjects %for all subjects
    for task=1:n_Task %for all Blocks
        
        % Get the BOLD signals from this subject in this task
        BOLD = df{s,task};
        Phase_BOLD=zeros(N_areas,Tmax); 

        % Get the BOLD phase using the Hilbert transform
        for seed=1:N_areas
            BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:)); %for this region, demean the timecourse
            signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));                    
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
        end

        %%%%%%%%%%%%%%%%% DOING THIS TO REMOVE EDGE EFFECTS
        Phase_BOLD(:,[1:5,Tmax-4:Tmax])=[];
        
        %%%%%%%%%%%%%%%%% DOING THIS TO REMOVE EDGE EFFECTS
        for t=1:Tmax-10 %for each time point      
            iFC=zeros(N_areas); 
            for n=1:N_areas 
                for p=1:N_areas
                    iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
                end
            end
            
            [V1,~]=eigs(iFC,1);
 
            if mean(V1>0)>.5
                V1=-V1;
            elseif mean(V1>0)==.5 && sum(V1(V1>0))>-sum(V1(V1<0))
                V1=-V1;
            end

            % Save V1 from all frames in all fMRI sessions in Leading eig
            t_all=t_all+1; % Update time
            Leading_Eig(t_all,:)=V1;
            Time_all(:,t_all)=[s task]; % Information that at t_all, V1 corresponds to subject s in a given task
        end
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Squeeze leading Eigs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_sub=size(Leading_Eig,1);
X=[];
for s=1:N_sub    
    X=cat(1,X,squeeze(Leading_Eig(s,:,:))); 
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kmeans clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mink=2;
maxk=15;
rangeK=mink:maxk;
opt= statset('UseParallel',1); %,'UseSubstreams',1);    % The options may vary according to the Matlab version
Kmeans_results=cell(size(rangeK));

for k=mink:maxk 
    [IDX, C, SUMD, D]=kmeans(X,k,'Distance','cityblock','Replicates',20,'Display','off'); %,'Options',opt);   
    Kmeans_results{k-(mink-1)}.IDX=IDX;
    Kmeans_results{k-(mink-1)}.C=C; 
    Kmeans_results{k-(mink-1)}.SUMD=SUMD; 
    Kmeans_results{k-(mink-1)}.D=D; 
    Kmeans_results{k-(mink-1)}.CalHar = evalclusters(X,IDX,'CalinskiHarabasz');
    Kmeans_results{k-(mink-1)}.DavBoul = evalclusters(X,IDX,'DaviesBouldin');
    Kmeans_results{k-(mink-1)}.Sil = evalclusters(X,IDX,'silhouette');
end

cd (outputpath)
% save('Checkpoint1_LEiDA_MID.mat','-v7.3')
save('Checkpoint1_LEiDA_CUE.mat','-v7.3')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate Clustering performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distM_fcd=squareform(pdist(X,'cityblock'));
dunn_score=zeros(maxk,1);
for j=mink-1:maxk-1
   dunn_score(j)=dunns(j,distM_fcd,Kmeans_results{j}.IDX);
   disp(['Performance for ' num2str(j) ' clusters'])
end

figure()
hold on
subplot(2,2,1)
plot(dunn_score)
title('Higher Dunn score = optimal solution')
xticklabels=[mink:1:maxk];
xlabel('Number of clusters')
ylabel('Dunn score')

ch=zeros(maxk,1);
db=zeros(maxk,1);
sil=zeros(maxk,1);
for j=mink-1:maxk-1
    ch(j)=Kmeans_results{j}.CalHar.CriterionValues;
    db(j)=Kmeans_results{j}.DavBoul.CriterionValues;
    sil(j)=Kmeans_results{j}.Sil.CriterionValues;
end

[~,ind_max]=max(ch);
disp(['Best clustering solution: ' num2str(ind_max) ' clusters'])

subplot(2,2,2)
xticklabels=[mink:1:maxk];
plot(ch(1:maxk-1,1))
title('Higher CH score = optimal solution')
xlabel('Number of clusters')
ylabel('Calinski Harabasz score')

subplot(2,2,3)
xticklabels=[mink:1:maxk];
plot(db(1:maxk-1,1))
title('Smallest DB score = optimal solution')
xlabel('Number of clusters')
ylabel('Davies Bouldin score')

subplot(2,2,4)
xticklabels=[mink:1:maxk];
plot(sil(1:maxk-1,1))
title('Silhoutte score closest to 1 = optimal solution')
xlabel('Number of clusters')
ylabel('Silhouette score')

% saveas(gcf,'MID_ClusteringPerformance.png');
saveas(gcf,'ClusteringPerformance_CUE.png');
close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 - Analyse the clustering results between states
% This assess the lifetime (in seconds) and probability each state will occur 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=zeros(n_Task,n_Subjects,maxk-1,maxk); 
LT=zeros(n_Task,n_Subjects,maxk-1,maxk); 

for k=mink-1:maxk-1   % for each cluster
    for task=1:n_Task   % for each condition
        for s=1:n_Subjects % for each subject
            
            % Select the time points representing this subject and task               
            T=(Time_all(1,:)==s&Time_all(2,:)==task); 
   %         Ctime=Kmeans_results{k}.IDX(T(1:end-2)); % The cluster each timecourse component belongs to
            Ctime=Kmeans_results{k}.IDX(T); % The cluster each timecourse component belongs to

             for c=1:k+1  
                P(task,s,k,c)=mean(Ctime==c);   
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
                LT(task,s,k,c)=mean(C_Durations)*TR;
            end                
        end
    end
end   

% save('Checkpoint2_LEiDA_MID.mat','-v7.3')
save('Checkpoint2_LEiDA_CUE.mat','-v7.3')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyse the effect of condition on states 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%e%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Make sure data is properly sorted
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
        LTStatePerTask=squeeze(LT(:,:,KK,ind_sort{KK+1}(1,ii)))';
        % LTStatePerTask=array2table(LTStatePerTask);

        % % Order of cols = 'MD_PLC','MD_APR','HC_PLC','HC_AP'
        % % I renamed them C1-C4 for ease of running hte ANOVA
        % LTStatePerTask.Properties.VariableNames = {'C1','C2'};
        % 
        % % Create the effects design
        % withinDesign = table([1 1 2 2]',[1 2 1 2]','VariableNames',{'Group','Drug'});
        % withinDesign.Group = categorical(withinDesign.Group);
        % withinDesign.Drug = categorical(withinDesign.Drug);
        % 
        % % Create and run the repeated measures ANOVA model 
        % rm = fitrm(LTStatePerTask, 'C1-C4 ~ 1', 'WithinDesign', withinDesign);
        % AT = ranova(rm, 'WithinModel', 'Group*Drug');


        [LT_P{KK}(ii,1),~,stats]=signrank(LTStatePerTask(:,1),LTStatePerTask(:,2));

        if isfield(stats,'zval')==1
        LT_ETA{KK}(ii,1)=stats.zval;
        else
        LT_ETA{KK}(ii,1)=0;
        end


        % % Effect of drug
        % LT_P{KK}(ii,2)=table2array(AT(5,8));
        % Eta=table2array(AT(5,1))/(table2array(AT(5,1))+table2array(AT(6,1)));
        % LT_ETA{KK}(ii,2)=Eta;
        % 
        % % Effect of group and drug
        % LT_P{KK}(ii,3)=table2array(AT(7,8));
        % Eta=table2array(AT(7,1))/(table2array(AT(7,1))+table2array(AT(8,1)));
        % LT_ETA{KK}(ii,3)=Eta;
           

        % For prob 
        PStatePerTask=squeeze(P(:,:,KK,ind_sort{KK+1}(1,ii)))';
        [Prob_P{KK}(ii,1),~,stats]=signrank(PStatePerTask(:,1),PStatePerTask(:,2));
        
        if isfield(stats,'zval')==1
        Prob_ETA{KK}(ii,1)=stats.zval;
        else
        Prob_ETA{KK}(ii,1)=0;
        end
        % 
        % PStatePerTask=array2table(PStatePerTask);
        % PStatePerTask.Properties.VariableNames = {'C1','C2','C3','C4'};

        % % Create the effects design
        % withinDesign = table([1 1 2 2]',[1 2 1 2]','VariableNames',{'Group','Drug'});
        % withinDesign.Group = categorical(withinDesign.Group);
        % withinDesign.Drug = categorical(withinDesign.Drug);
        % 
        % % Create and run the repeated measures ANOVA model 
        % rm = fitrm(PStatePerTask, 'C1-C4 ~ 1', 'WithinDesign', withinDesign);
        % AT = ranova(rm, 'WithinModel', 'Group*Drug');

        % Effect of group
        % Prob_P{KK}(ii,1)=table2array(AT(3,8));
        % Eta=table2array(AT(3,1))/(table2array(AT(3,1))+table2array(AT(4,1)));
        % Prob_ETA{KK}(ii,1)=Eta;
        
        % % Effect of drug
        % Prob_P{KK}(ii,2)=table2array(AT(5,8));
        % Eta=table2array(AT(5,1))/(table2array(AT(5,1))+table2array(AT(6,1)));
        % Prob_ETA{KK}(ii,2)=Eta;
        % 
        % % Effect of group and drug
        % Prob_P{KK}(ii,3)=table2array(AT(7,8));
        % Eta=table2array(AT(7,1))/(table2array(AT(7,1))+table2array(AT(8,1)));
        % Prob_ETA{KK}(ii,3)=Eta;
           

    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post hocs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P(task,subj,k,c)
for KK=mink-1:maxk-1
for cc=1:KK+1
PStatePerTask=squeeze(P(:,:,KK,ind_sort{KK+1}(1,cc)))';
LTStatePerTask=squeeze(LT(:,:,KK,ind_sort{KK+1}(1,cc)))';

for ii=1:n_Task
for jj=1:ii
        
        % Non-parametric - Wilcoxon Signed Rank 
        xp=PStatePerTask(:,ii);
        yp=PStatePerTask(:,jj);
        Prob_Posthoc{KK,cc}(ii,jj)=signrank(xp,yp);
        
        if Prob_Posthoc{KK,cc}(ii,jj)<0.05
            disp(strcat('K=',num2str(KK+1),'s=',num2str(cc),'stimCond',num2str(ii),'vs stimCond',num2str(jj)))
        end
        
        xlt=LTStatePerTask(:,ii);
        ylt=LTStatePerTask(:,jj);
        LT_Posthoc{KK,cc}(ii,jj)=signrank(xlt,ylt);  
        
end
end

end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDR correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:maxk-1
    [~,~,LT_P_FDR{ii}(:,1)] = fdr(LT_P{ii}(:,1));
    [~,~,Prob_P_FDR{ii}(:,1)] = fdr(Prob_P{ii}(:,1));
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Organising outputs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Put in nice table
Prob_Fancy={};
LT_Fancy={};

for k=mink-1:maxk-1
for c=1:k+1
tmp=strcat(num2str(Prob_P_FDR{k}(c,1)),', ',num2str(round(Prob_ETA{k}(c,1),2,'decimal')));
Prob_Fancy{1}(c,k)=convertCharsToStrings(tmp);
tmp=strcat(num2str(LT_P_FDR{k}(c,1)),', ',num2str(round(LT_ETA{k}(c,1),2,'decimal')));
LT_Fancy{1}(c,k)=convertCharsToStrings(tmp);
end
end

% save('Checkpoint3_LEiDA_MID.mat','-v7.3');
save('Checkpoint3_LEiDA_CUE.mat','-v7.3');
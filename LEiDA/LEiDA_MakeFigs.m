%{
This script creates the following figures
1. Updated connectivity plots
2. Correlation matrix between regions per state
3. Radarplots for probability and lifetime
4. Matrices and Digraphs showing transition probabilities
5. Switches per condition
6. The proportion each state occurs per subject
7. Cluster time series per subject and task timeseries

Credits
LEiDA was first created by Joana Cabral, Oct 2017, joana.cabral@psych.ox.ac.uk. First use in: Cabral, et al. 2017 Scientific reports 7, no. 1 (2017): 5135.
Many of these plots were made by Henry Hebron, Oct 2021, h.hebron@surrey.ac.uk
Adapted by Danielle Kurtin, Jan 2022, d.kurtin@surrey.ac.uk 
%}

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add paths
currloc=string(pwd);
homeDir=('C:\Users\dk818\OneDrive - Imperial College London\NCORE\miFCdFC');

% For outputs
str='\Outputs';
outputpath=strcat(homeDir,str);
addpath(outputpath);

% For atlas and other data
str='\AtlasAndOtherInputs';
addpath(strcat(homeDir,str));

% For LEiDA
str='\LEiDA';
addpath(strcat(homeDir,str));

% For ICC, violinplot, colors, Dunn score, and more
str='\Functions';
addpath(strcat(homeDir,str));

% For spm
% Because it's such a large package I'm just linking it directly to my local version, instead
% of moving it to a diff folder
addpath('C:\Users\dk818\OneDrive - Imperial College London\LilThingsAndMATLABPrograms\spm12');

% For colorbrewer - is only used to visualise the state prob per subj, may
% not need it in general, so I've linked it directly in hopes of not using
% it eventually
addpath('C:\Users\dk00549\OneDrive - University of Surrey\Documents\Surrey\MATLAB\PluginsProgramsToolboxes\ColorBrewer')

% KEY PARAMETER - SET VALUE FOR REPRESENTATIVE VALUE OF K
K = 5;
k = K;
ind_max = K;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOK DOWN

% LOOK DOWN 

% LOOK DOWN DON'T IGNOREEEEEE!!!!!!!!!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WHCIH TASK YOU WANNA DO?!
% THIS WILL CHANGE WHAT ALL THE SAVE NAMES ARE - BE SURE TO UNCOMMENT
% APPROPRIATELY!!!!!!!!!!!!!!!

load('Checkpoint1_LEiDA_MID.mat')
% load('Checkpoint1_LEiDA_CUE.mat')

% Loading and general params
colors = linspecer(num_condi);
timeser = size(Phase_BOLD,2); % this calculates how long the extracted timeseries are
h=hist(Best_Clusters.IDX,ind_max);
[y, ind]=sort(h,'descend');
SchaeferNames=readtable('ComboNames.xlsx');
% SchaeferLabels=SchaeferNames(:,2);

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
    scatter(xx,log10(Prob_P_FDR{1,kk}(:,1)),30,C(kk,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
%     swarmchart(xx,log10(Prob_P_FDR{1,kk}(:,1)),30,C(kk,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
end
yline(log10(0.05))
alpha(0.7)
% name=strcat('CUE_SwarmChatPvals_Prob_Result',num2str(hh),'.png');
name='MID_SwarmChatPvals_Prob.png';
saveas(gcf,name);
close all

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
     scatter(xx,log10(LT_P_FDR{1,kk}(:,1)),30,C(kk,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
%         swarmchart(xx,log10(LT_P_FDR{1,kk}(:,jj)),30,C(kk,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
end
yline(log10(0.05))
alpha(0.7)
% name=strcat('CUE_SwarmChatPvals_LT_Result',num2str(hh),'.png');
% name=strcat('MID_SwarmChatPvals_LT_Result',num2str(hh),'.png');
name='MID_SwarmChatPvals_LT.png';
saveas(gcf,name);
close all


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make violins for LT and P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For LT
% MD_PLC=squeeze(LT(1,:,k-1,1:k));
% MD_APR=squeeze(LT(2,:,k-1,1:k));
% HC_PLC=squeeze(LT(3,:,k-1,1:k));
% HC_APR=squeeze(LT(4,:,k-1,1:k));

% For Prob
MD_PLC=squeeze(P(1,:,k-1,1:k));
MD_APR=squeeze(P(2,:,k-1,1:k));
HC_PLC=squeeze(P(3,:,k-1,1:k));
HC_APR=squeeze(P(4,:,k-1,1:k));

ind_sort2=[];
for KK=mink:maxk
Best_Clusters=Kmeans_results{rangeK==KK};
ProbC=zeros(1,K);
for c=1:K
    ProbC(c)=mean(Best_Clusters.IDX==c);
end
[~, ind_sort2(KK,:)]=sort(ProbC,'descend'); 
end

State={};
for ii=1:K
   State{ii}(:,1)=MD_PLC(:,ind_sort2(K,ii));
   State{ii}(:,2)=MD_APR(:,ind_sort2(K,ii));
   State{ii}(:,3)=HC_PLC(:,ind_sort2(K,ii));
   State{ii}(:,4)=HC_APR(:,ind_sort2(K,ii));
end

C=linspecer(4);

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
% name=strcat('CUE_ProbViolins_State',num2str(cc),'.png');
% name=strcat('MID_ProbViolins_State',num2str(cc),'.png');
% name=strcat('LTViolins_State',num2str(cc),'.png');
% saveas(gcf,name);
% close all
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot pvals for Complexity - effect of group
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:15
   ComplexSig(ii,3)=ii; 
end

for jj=1:3
C=linspecer(7);
figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
% LZC
scatter(ComplexSig(2:15,3),log10(ComplexSig(2:15,1,jj)),50,C(1,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
% BDMC
scatter(ComplexSig(2:15,3),log10(ComplexSig(2:15,2,jj)),50,C(2,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
% 0TE
scatter(ComplexSig(2:15,3),log10(TransSig(2:15,1,jj)),50,C(3,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
% 1TE
scatter(ComplexSig(2:15,3),log10(TransSig(2:15,2,jj)),50,C(4,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
% 2TE
scatter(ComplexSig(2:15,3),log10(TransSig(2:15,3,jj)),50,C(5,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
% 3TE
scatter(ComplexSig(2:15,3),log10(TransSig(2:15,4,jj)),50,C(6,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
% 4TE
scatter(ComplexSig(2:15,3),log10(TransSig(2:15,5,jj)),50,C(7,:),'filled','jitter', 'on', 'jitterAmount', 0.2,'MarkerEdgeColor',[0 0 0])
yline(log10(0.05))
alpha(0.7)

% name=strcat('MID_SwarmChart_Complexity_Result',num2str(jj),'.png');
% name=strcat('CUE_SwarmChart_Complexity_Result',num2str(jj),'.png');
% saveas(gcf,name);
% close all

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Violins for each complexity metric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=linspecer(7);

% First, I need a table w/ everyone's results: Results + TransH tbl
for CC=K
for tt=1:num_condi

VizTbl(:,1,tt)=Results{CC}(:,1,tt);
VizTbl(:,2,tt)=Results{CC}(:,2,tt);
VizTbl(:,3,tt)=TransitionsH{CC}(1,:,tt)';
VizTbl(:,4,tt)=TransitionsH{CC}(2,:,tt)';
VizTbl(:,5,tt)=TransitionsH{CC}(3,:,tt)';
VizTbl(:,6,tt)=TransitionsH{CC}(4,:,tt)';
VizTbl(:,7,tt)=TransitionsH{CC}(5,:,tt)';

for ii=1:7
   VizTblNorm(:,ii,tt)=InvRankTrans(VizTbl(:,ii,tt));
end

end
end

% Now our table is subj:ComplexMetrix:task, adn we want the same metric,
% but across tasks
for ii=1:7
figure
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
tmp=squeeze(VizTbl(:,ii,:));
violins = violinplot(tmp);
%xticklabels({'LZC','BDMC','0TE','1TE','2TE','3TE','4TE'})
% xticklabels({'MD PLC','MD APR','HC PLC','HC APR'})
% ylabel({'Normalized complexity vals'})

for i = 1:num_condi
    violins(i).ShowMean = 'yes';
    violins(i).ViolinColor = C(i,:);
end

% for sub = 1:n_Subjects
%     pl_ = plot(1:num_condi,tmp(sub,:),'Color',[.5 .5 .5],'LineStyle','-.');
%     pl_.Color(4) = .4;
% end


% name=strcat('MID_ComplexityViolin_Metric',num2str(ii),'_Result',num2str(jj),'.png');
% name=strcat('CUE_ComplexityViolin_Metric',num2str(ii),'_Result',num2str(jj),'.png');
% saveas(gcf,name);
% close all

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sexy brain plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for K=4 %mink:maxk
disp(strcat('Now running for K=',num2str(K)))
Best_Clusters=Kmeans_results{K-1};
[N_Cl, N_ba]=size(Best_Clusters.C);
h=hist(Best_Clusters.IDX,N_Cl);
[y, ind]=sort(h,'descend');
V=Best_Clusters.C(ind,:);

if K==mink
    for c=1:K
        plot_nodes_in_cortex_SchaeferShen_216ROI(V(c,:))
        title(['#' num2str(c)])
        name=strcat('SexyBrainPlot_K',num2str(K),'_State',num2str(c),".png");
        saveas(gcf,name)
        close all
    end
else
   for c=2:K
    plot_nodes_in_cortex_SchaeferShen_216ROI(V(c,:))
    title(['#' num2str(c)])
    name=strcat('SexyBrainPlot_K',num2str(K),'_State',num2str(c),".png");
    saveas(gcf,name)
    close all
   end
end
end


%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % All in one - sexy brain plots
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for K=mink:maxk
% disp(strcat('Now running for K=',num2str(K)))
% Best_Clusters=Kmeans_results{K-1};
% [N_Cl, N_ba]=size(Best_Clusters.C);
% h=hist(Best_Clusters.IDX,N_Cl);
% [y, ind]=sort(h,'descend');
% V=Best_Clusters.C(ind,:);
% 
% figure
% hold on
% % 
% % if K==mink
%     for c=1:K
% %         subplot(1,K,c)
%         subaxis(1,K,c, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
%         plot_nodes_in_cortex_SchaeferShen_216ROI(V(c,:))
% %         title(['#' num2str(c)])
% %         name=strcat('SexyBrainPlot_K',num2str(K),'_State',num2str(c),".png");
% %         saveas(gcf,name)
% %         close all
% %     end
% % else
% %    for c=2:K
% %     subplot(1,K,c)
% %     plot_nodes_in_cortex_SchaeferShen_216ROI(V(c,:))
% %     title(['#' num2str(c)])
% %     name=strcat('SexyBrainPlot_K',num2str(K),'_State',num2str(c),".png");
% %     saveas(gcf,name)
% %     close all
%    end
% % end
% 
% axis tight
% axis off
% hold off
% % name=strcat('SexyBrainPlot_K',num2str(K),'_State',num2str(c),".png");
% name=strcat('SexyBrainPlot_K',num2str(K),".png");
% saveas(gcf,name)
% close all
% end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOLD projection plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Best_Clusters=Kmeans_results{K-1};
[N_Cl, N_ba]=size(Best_Clusters.C);
h=hist(Best_Clusters.IDX,N_Cl);
[y, ind]=sort(h,'descend');
V=Best_Clusters.C(ind,:);
Order=[1:2:N_ba N_ba:-2:2];
SchaeferNamesRed=table2array(SchaeferNames);
% SchaeferNamesRed(1:10:end,:)=[];

figure
for c=1:K
    subplot(1,K,c)
    Vc=V(c,:);
    hold on
    barh(Vc.*(Vc<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh(Vc.*(Vc>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 N_areas])
    xlim([-.15 .15])
    set(gca,'YTick',1:10:N_areas,'Fontsize',14)
    if c==1
%         set(gca,'YTickLabel',SchaeferNames(end:-1:1,:),'Fontsize',8)
        set(gca,'YTickLabel',SchaeferNamesRed(1:10:N_areas),'Fontsize',8)
    else
        set(gca,'YTickLabel',[])
    end
    ylim([0 N_areas])
    xlim([-.15 .15])
    set(gca,'Ydir','reverse')
    title(['State ' num2str(c)])
    grid on
end
% cd (outputpath)
% saveas(gcf,'MID_BOLDProjectionPlot.png');
% saveas(gcf,'CUE_BOLDProjectionPlot.png');

%%
% Make a list of ROIs that are in 2ndary community per state
Best_Clusters=Kmeans_results{K-1};
[N_Cl, N_ba]=size(Best_Clusters.C);
h=hist(Best_Clusters.IDX,N_Cl);
[y, ind]=sort(h,'descend');
V=Best_Clusters.C(ind,:);
clear nameList
h=0;
for c=2:K
    Vc=V(c,:);
    for ii=1:length(Vc)
        if Vc(1,ii)>0
            h=h+1;
            nameList(h,c-1)=SchaeferNamesRed(ii,1);
        end
    end
    h=0;
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proportion each state occurs per subj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Best_Clusters=Kmeans_results{ind_max};

% Compute the cluster time series per person
scan_length = length(Best_Clusters.IDX)/n_Subjects;
new_index = [];
I = 1;

for sub = 1:n_Subjects
    new_index = [new_index; Best_Clusters.IDX(I:I+scan_length-1); NaN];
    subject_states(sub,:) = Best_Clusters.IDX(I:I+scan_length-1);
    I = I+scan_length;       
end

% Compute the proportion each state has per each subject's cluster time
% series
for state = 1:ind_max
    prob_per_sub(:,state) = 100.*(sum(subject_states==state,2)./length(subject_states(1,:)));  
end

% this sorts all rows by col 1
newProbPerSub=sortrows(prob_per_sub,1,'descend');

% check to see which 
newProbPerSub2=[];
for ii=[3,4,2,1]
    w=size(newProbPerSub2,2)+1;
    newProbPerSub2(:,w)=newProbPerSub(:,ii);
end

figure()
hold on
imagesc(newProbPerSub2)
% [rows,cols] = size(prob_per_sub);
axis tight
% colormap(brewermap([],"YlGnBu"));
colorbar
ylabel('Subject')
xlabel('State')
title('State probability per subject')
% saveas(gcf,'MID_StateProbPerSubj.png');
% saveas(gcf,'CUE_StateProbPerSubj.png');
% close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ICC for state occurrence per subj 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk=mink-1:maxk-1

Best_Clusters=Kmeans_results{kk};

% Compute the cluster time series per person
scan_length = length(Best_Clusters.IDX)/n_Subjects;
new_index = [];
I = 1;

for sub = 1:n_Subjects
    new_index = [new_index; Best_Clusters.IDX(I:I+scan_length-1); NaN];
    subject_states(sub,:) = Best_Clusters.IDX(I:I+scan_length-1);
    I = I+scan_length;       
end

% Compute the proportion each state has per each subject's cluster time
% seres
for state = 1:kk
    prob_per_sub(:,state) = 100.*(sum(subject_states==state,2)./length(subject_states(1,:)));  
end

% Compute the ICC among subjs - gotta flip the matrix to get it
flip_prob_per_sub=prob_per_sub';
[ICCr(kk,1), LB, UB, F, df1, df2, ICCp(kk,1)] = ICC(flip_prob_per_sub, 'A-1');

end

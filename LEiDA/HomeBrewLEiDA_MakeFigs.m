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


%%
% Need to make LT, P, and Results for both Cue and MID task, so I can plot
% everything nicely
% LT & P = cond x subj x network
% Results =  cond x subj x complexity measure

load('HomeBrewLEiDA_MID_Matrix.mat')
% For Prob
Prob_MD_PLC=squeeze(P(1,:,:));
Prob_HC_PLC=squeeze(P(2,:,:));
Prob_State={};
for ii=1:num_states
   Prob_State{ii}(:,1)=Prob_MD_PLC(:,ii);
   Prob_State{ii}(:,2)=Prob_HC_PLC(:,ii);
end

% For LT
LT_MD_PLC=squeeze(LT(1,:,:));
LT_HC_PLC=squeeze(LT(2,:,:));
LT_State={};
for ii=1:num_states
   LT_State{ii}(:,1)=LT_MD_PLC(:,ii);
   LT_State{ii}(:,2)=LT_HC_PLC(:,ii);
end

% For Complexity
Cm_MD_PLC=squeeze(Results(1,:,:));
Cm_HC_PLC=squeeze(Results(2,:,:));
Cm_State={};
for ii=1:7 % num of complexity res
   Cm_State{ii}(:,1)=Cm_MD_PLC(:,ii);
   Cm_State{ii}(:,2)=Cm_HC_PLC(:,ii);
end

clear LT P Results

load('HomeBrewLEiDA_CUE_Matrix.mat')
% For Prob
Prob_MD_PLC=squeeze(P(1,:,:));
Prob_HC_PLC=squeeze(P(2,:,:));
for ii=1:num_states
   Prob_State{ii}(:,3)=Prob_MD_PLC(:,ii);
   Prob_State{ii}(:,4)=Prob_HC_PLC(:,ii);
end

% For LT
LT_MD_PLC=squeeze(LT(1,:,:));
LT_HC_PLC=squeeze(LT(2,:,:));
for ii=1:num_states
   LT_State{ii}(:,3)=LT_MD_PLC(:,ii);
   LT_State{ii}(:,4)=LT_HC_PLC(:,ii);
end

% For Complexity
Cm_MD_PLC=squeeze(Results(1,:,:));
Cm_HC_PLC=squeeze(Results(2,:,:));
for ii=1:7 % num of complexity res
   Cm_State{ii}(:,3)=Cm_MD_PLC(:,ii);
   Cm_State{ii}(:,4)=Cm_HC_PLC(:,ii);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make violins for LT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=linspecer(2);

for cc=1:num_states
figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
tmp=LT_State{cc};   % get the LT or the Prob
violins = violinplot(tmp);
for i = 1:4
    if i == 1 || i ==3
    violins(i).ShowMean = 'yes';
    violins(i).ViolinColor = C(1,:);
    else
    violins(i).ShowMean = 'yes';
    violins(i).ViolinColor = C(2,:);
    end 
end
hold off
name=strcat('HomeBrew_ViolinsLT_State',num2str(cc),'.png');
saveas(gcf,name);
close all
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make violins for Prob
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=linspecer(2);

for cc=1:num_states
figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
tmp=Prob_State{cc};   % get the LT or the Prob
violins = violinplot(tmp);
for i = 1:4
    if i == 1 || i ==3
    violins(i).ShowMean = 'yes';
    violins(i).ViolinColor = C(1,:);
    else
    violins(i).ShowMean = 'yes';
    violins(i).ViolinColor = C(2,:);
    end 
end
hold off
name=strcat('HomeBrew_ViolinsProb_State',num2str(cc),'.png');
saveas(gcf,name);
close all
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make violins for Complexity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=linspecer(2);

for cc=1:7
figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
tmp=Cm_State{cc};   % get the LT or the Prob
violins = violinplot(tmp);
for i = 1:4
    if i == 1 || i ==3
    violins(i).ShowMean = 'yes';
    violins(i).ViolinColor = C(1,:);
    else
    violins(i).ShowMean = 'yes';
    violins(i).ViolinColor = C(2,:);
    end 
end
hold off
name=strcat('HomeBrew_ViolinsComplexity_Result',num2str(cc),'.png');
saveas(gcf,name);
close all
end


%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Sexy brain plots
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for K=5 %mink:maxk
% disp(strcat('Now running for K=',num2str(K)))
% Best_Clusters=Kmeans_results{K-1};
% [N_Cl, N_ba]=size(Best_Clusters.C);
% h=hist(Best_Clusters.IDX,N_Cl);
% [y, ind]=sort(h,'descend');
% V=Best_Clusters.C(ind,:);
% 
% if K==mink
%     for c=1:K
%         plot_nodes_in_cortex_SchaeferShen_216ROI(V(c,:))
%         title(['#' num2str(c)])
%         name=strcat('CUE_SexyBrainPlot_K',num2str(K),'_State',num2str(c),".png");
%         % name=strcat('MID_SexyBrainPlot_K',num2str(K),'_State',num2str(c),".png");
%         saveas(gcf,name)
%         close all
%     end
% else
%    for c=2:K
%     plot_nodes_in_cortex_SchaeferShen_216ROI(V(c,:))
%     title(['#' num2str(c)])
%     name=strcat('CUE_SexyBrainPlot_K',num2str(K),'_State',num2str(c),".png");
%     % name=strcat('MID_SexyBrainPlot_K',num2str(K),'_State',num2str(c),".png");
%     saveas(gcf,name)
%     close all
%    end
% end
% end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proportion each state occurs per subj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Best_Clusters=Kmeans_results{ind_max};
% 
% % Compute the cluster time series per person
% scan_length = length(Best_Clusters.IDX)/n_Subjects;
% new_index = [];
% I = 1;
% 
% for sub = 1:n_Subjects
%     new_index = [new_index; Best_Clusters.IDX(I:I+scan_length-1); NaN];
%     subject_states(sub,:) = Best_Clusters.IDX(I:I+scan_length-1);
%     I = I+scan_length;       
% end
% 
% % Compute the proportion each state has per each subject's cluster time
% % series
% for state = 1:ind_max
%     prob_per_sub(:,state) = 100.*(sum(subject_states==state,2)./length(subject_states(1,:)));  
% end
% 
% % this sorts all rows by col 1
% newProbPerSub=sortrows(prob_per_sub,1,'ascend');
% 
% % check to see which 
% newProbPerSub2=[];
% for ii=[1,2,4,3,5] 
%     w=size(newProbPerSub2,2)+1;
%     newProbPerSub2(:,w)=newProbPerSub(:,ii);
% end
% 
% figure()
% hold on
% imagesc(newProbPerSub2)
% % [rows,cols] = size(prob_per_sub);
% axis tight
% % colormap(brewermap([],"YlGnBu"));
% colorbar
% ylabel('Subject')
% xlabel('State')
% title('State probability per subject')
% set(gca,'FontSize',20)
% set(gca,'FontName','Arial')
% % saveas(gcf,'MID_StateProbPerSubj.png');
% saveas(gcf,'CUE_StateProbPerSubj.png');
% close all
% 
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ICC for state occurrence per subj 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for kk=mink-1:maxk-1
% 
% Best_Clusters=Kmeans_results{kk};
% 
% % Compute the cluster time series per person
% scan_length = length(Best_Clusters.IDX)/n_Subjects;
% new_index = [];
% I = 1;
% 
% for sub = 1:n_Subjects
%     new_index = [new_index; Best_Clusters.IDX(I:I+scan_length-1); NaN];
%     subject_states(sub,:) = Best_Clusters.IDX(I:I+scan_length-1);
%     I = I+scan_length;       
% end
% 
% % Compute the proportion each state has per each subject's cluster time
% % seres
% for state = 1:kk
%     prob_per_sub(:,state) = 100.*(sum(subject_states==state,2)./length(subject_states(1,:)));  
% end
% 
% % Compute the ICC among subjs - gotta flip the matrix to get it
% flip_prob_per_sub=prob_per_sub';
% [ICCr(kk,1), LB, UB, F, df1, df2, ICCp(kk,1)] = ICC(flip_prob_per_sub, 'A-1');
% 
% end

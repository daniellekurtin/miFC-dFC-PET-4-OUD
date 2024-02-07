%{ 
Code for complexity, from Greg Scott
(gregory.scott99@imperial.ac.uk)

Adapted for use w/ cluster time series computed using Leading Eigenvector
Dynamic Analysis by Danielle Kurtin

Clear LZ ref: https://information-dynamics.github.io/complexity/information/2019/06/26/lempel-ziv.html
Greg & Rich's preprint: https://www.biorxiv.org/content/biorxiv/early/2020/10/15/2020.06.17.156570.full.pdf

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

% For complexity
str='\Functions\complexity';
addpath(strcat(homeDir,str));

% For ICC, violinplot, colors, Dunn score, and more
str='\Functions';
addpath(strcat(homeDir,str));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading, general params, and organizing data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('Checkpoint1_LEiDA_v1.mat')
num_subj=n_Subjects;
timeser = size(Phase_BOLD,2); % this calculates how long the extracted timeseries are

for CC=mink:maxk

Best_Clusters=Kmeans_results{CC-1};
K=CC;
disp(strcat('Now running for ',num2str(K),' clusters'));

ProbC=zeros(1,K);
for c=1:K
    ProbC(c)=mean(Best_Clusters.IDX==c);
end
[~, ind_sort]=sort(ProbC,'descend'); 


for clust = 1:K
    cluster_time_series(Best_Clusters.IDX==ind_sort(clust),:) = clust;
end

ClustTimeSer = reshape(cluster_time_series,timeser,num_condi,n_Subjects);
clear transition_matrix
for sub = 1:n_Subjects
    for condi = 1:num_condi
        transition_matrix.cond{condi}(:,:,sub) = transitions_v2(ClustTimeSer(:,condi,sub),K).*100;
    end
end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute complexity
% Note - future versions should delete complexity metrics that aren't used
% May be as simple as only using the first two cols in the res table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results=[];
Transitions = {};
TransitionsH = [];

for CC=mink:maxk

    Best_Clusters=Kmeans_results{CC-1};
    K=CC;
    disp(strcat('Now running for ',num2str(K),' clusters'));
    ProbC=zeros(1,K);

    for c=1:K
        ProbC(c)=mean(Best_Clusters.IDX==c);
    end
    [~, ind_sort]=sort(ProbC,'descend'); 

    for clust = 1:K
        cluster_time_series(Best_Clusters.IDX==ind_sort(clust),:) = clust;
    end
    ClustTimeSer = reshape(cluster_time_series,timeser,num_condi,n_Subjects);

    for subj=1:num_subj
        fprintf('.');

        for i=1:num_condi
            S = squeeze(ClustTimeSer(:,i,subj));

            % Convert the string into a binary string with n-bit encoding,
            % where n=the number of bits to compute K
            % For a group of n bits, it is possible to represent 2^n values
            S0=S-1;  % This gives us state 0 (instead of starting at 1)
            if K==2 
                nbits=2;
            else
                
                nbits=ceil(sqrt(max(S0)))+1;
                %nbits=ceil(log102(max(S0)));
                % The original is below
                % nbits=ceil(log102(max(K-1)));  
            end
          
            Sb = ctx_dec2bi1d(S0, nbits);
            % Convert to state switches - this now is a vector starting from
            % the first switch. I'll add a zero as the first val, so it's the
            % length of everything else, and can be added to the same Results
            % matric
            tmp = diff(S)~=0;
            Ss=zeros(length(tmp)+1,1);
            Ss(2:length(tmp)+1,1)=tmp;
            Ss=Ss';

            % Convert the state timseries sequence into a 2D form - 1-hop encoding
            Sc = ctx_dec2col(S);

            % Convert the switch timseries sequence into a 2D form - 1-hop
            % encoding. Need to have Ss not start from 0
            SS=Ss+1;
            SsSc = ctx_dec2col(SS);

            % Calculate LZC, normalise by LZC on a single shuffled version of the data
            % (could divide by mean of N shuffled versions) - Greg isn't sure
            % whether there's a need to shuffle. 
            Results{CC}(subj,1, i) = calc_lz_complexity(Sb, 'exhaustive', false) ./ ...
                calc_lz_complexity(Sb(randperm(length(Sb))), 'exhaustive', false);

            % Convert the binary string into a Matlab string format for BDM1d
            Sstr = ctx_bi2str(Sb);

            % Calculate BDM1d
            Results{CC}(subj,2,i) = ctx_stringbdm1d(Sstr);

            % Calculate BDM 2d version (needs a 2D non-string input)
             Results{CC}(subj,3,i) = ctx_bdm2d(Sc);

            % Use switching timeseries, not state timeseries
            % Calculate LZC on switching timeseries (not state timeseries) 
             Results{CC}(subj,4, i) = calc_lz_complexity(Ss, 'exhaustive', false) ./ ...
                calc_lz_complexity(Ss(randperm(length(Ss))), 'exhaustive', false);

            % Convert the binary string into a Matlab string format for BDM1d
            SsSstr = ctx_bi2str(Ss);

            % Calculate BDM1d using state switching string
            Results{CC}(subj,5,i) = ctx_stringbdm1d(SsSstr);

            % Calculate BDM 2d version (needs a 2D non-string input)
            Results{CC}(subj,6,i) = ctx_bdm2d(SsSc);

            % Calculate transition entropy
            for trorder = 0:4
                Transitions{subj,i} = markovmodel(ClustTimeSer(:,i,subj), trorder);
                n = Transitions{subj,i}(:,end);
                p = n./sum(n);
                p = p(p>0);
                TransitionsH{CC}(trorder + 1, subj,i) = -sum(p .*log2(p));
    
            end
        end
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the effect of session on complexity metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TableC=[];

% For LZC and BDMC
for CC=mink:maxk
disp(strcat('Computing the effect of condition on complexity metric when K=', num2str(CC)))
clear Table TMP
TMP=Results{CC};

for cidx=1:2
    res = squeeze(TMP(:,cidx,:));   %this gives a participant x condition table        
    [ComplexSig(CC,cidx),~,stats]=signrank(res(:,1),res(:,2));
    if isfield(stats,'zval')==1
    TableC(CC,cidx)=stats.zval;
    else
    TableC(CC,cidx)=0;
    end
end
end

%%
% For trans ent
for CC=mink:maxk

disp(strcat('Now running for K=',num2str(CC)))
% Below = 0-4th order transition matrices
% Transitions{K}(transorder,subj,task)
Trans0=squeeze(TransitionsH{CC}(1,:,:));
Trans1=squeeze(TransitionsH{CC}(2,:,:));
Trans2=squeeze(TransitionsH{CC}(3,:,:));
Trans3=squeeze(TransitionsH{CC}(4,:,:));
Trans4=squeeze(TransitionsH{CC}(5,:,:));

% Trans0
[TransSig(CC,1),~,stats]=signrank(Trans0(:,1),Trans0(:,2));
if isfield(stats,'zval')==1
TableTE(CC,1)=stats.zval;
else
TableTE(CC,1)=0;
end

% Trans1
[TransSig(CC,2),~,stats]=signrank(Trans1(:,1),Trans1(:,2));
if isfield(stats,'zval')==1
TableTE(CC,2)=stats.zval;
else
TableTE(CC,2)=0;
end

% Trans2
[TransSig(CC,3),~,stats]=signrank(Trans2(:,1),Trans2(:,2));
if isfield(stats,'zval')==1
TableTE(CC,3)=stats.zval;
else
TableTE(CC,3)=0;
end

% Trans3
[TransSig(CC,4),~,stats]=signrank(Trans3(:,1),Trans3(:,2));
if isfield(stats,'zval')==1
TableTE(CC,4)=stats.zval;
else
TableTE(CC,4)=0;
end

% Trans4
[TransSig(CC,5),~,stats]=signrank(Trans4(:,1),Trans4(:,2));
if isfield(stats,'zval')==1
TableTE(CC,5)=stats.zval;
else
TableTE(CC,5)=0;
end

end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Organising pvals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put pvals in a nice format with ChiSq val
% ComplexSig(CC,cidx)
% TableC(CC,cidx)

% Put in nice table
% FancyFriedmanC={};

for CC=mink:maxk
for cidx=1:2
%tmp=strcat(num2str(round(FriedmanNice_P_padj(c,k),2,'decimal')),', ',num2str(round(FriedmanNiceForm_P_ChiSq(c,k),2,'decimal')));
tmp=strcat(num2str(ComplexSig(CC,cidx)),', ',num2str(round(TableC(CC,cidx),2,'decimal')));
FancyFriedmanC(CC,cidx)=convertCharsToStrings(tmp);
end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post hocs for LZC and BDMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TMP=Results{5};
res = squeeze(TMP(:,1,:));   %this gives a participant x condition table
for ii=1:num_condi
   for jj=1:num_condi
       ComplexPostHocLZC(ii,jj)=signrank(res(:,ii),res(:,jj));
   end
end    

res = squeeze(TMP(:,2,:));   %this gives a participant x condition table
for ii=1:num_condi
   for jj=1:num_condi
       ComplexPostHocBDMC(ii,jj)=signrank(res(:,ii),res(:,jj));
   end
end    


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post hocs for Transition Entropy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for CC=mink:maxk
for oo=1:5 % order of trans ent
for ii=1:num_condi
for jj=1:num_condi
        TransSigPostHocs{CC,oo}(ii,jj)=signrank(squeeze(TransitionsH{CC}(oo,:,ii))',squeeze(TransitionsH{CC}(oo,:,jj))');
end
end
end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Organise Pvals 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TransSig(CC,torder)
% TableTE(CC,torder)
% Put in nice table
% FancyFriedmanC={};


for CC=mink:maxk
for torder=1:5
tmp=strcat(num2str(TransSig(CC,torder)),', ',num2str(round(TableTE(CC,torder),2,'decimal')));
FancyFriedmanC(CC,2+torder)=convertCharsToStrings(tmp);
end
end

%%
% save('Checkpoint4_LEiDA_MID.mat','-v7.3')
save('Checkpoint4_LEiDA_CUE.mat','-v7.3')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEiDA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load df_forLEiDA.mat
load templateEigAndMat.mat
num_states=size(templateEig,2);
%%%%%%%%%%%%%
% NOTE TO SEELLLLFFFFFFFF - MUST CHANGE THIS LINE DEPENDING WHETHER YOU
% WANT TO RUN CUE OR NOT!!!!!!!!!
%%%%%%%%%%%%
df=df_cue;
% df=df_mid;
[~, n_Task]=size(df);
n_Subjects=22;
[N_areas, Tmax]=size(df{1,1});
num_condi=n_Task; % number of experimental conditions
TR=1.5;
Leading_Eig=zeros((Tmax-10)*n_Subjects,N_areas); % All leading eigenvectors, accounting for removal of edges
Time_all=zeros(2, n_Subjects*(Tmax-10)); % vector with subject nr and task at each t, except the edges
t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)
fnq=1/(2*TR);                 % Nyquist frequency.
flp = .02;                    % lowpass frequency of filter (Hz). This allows anything below 50 Hz. 
fhi = 0.1;                    % highpass- we've already applied one, though it is less inclusive (allowing anything over 0.01 to pass). I will keep this for now. 
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency. 
k=2;                          % 2nd order butterworth filter. This determines the steepness of the gain function (and 2 is pretty smooth). 
[bfilt,afilt]=butter(k,Wn);   % "Butter" is a MATLAB function than constructs the butterworth filter using defined cutoff frequencies.

for s=1:n_Subjects %for all subjects
    for task=1:n_Task %for all Blocks
        disp(strcat('Running New LEiDA subj',num2str(s),'task',num2str(task)))
        
        % Get the BOLD signals from this subject in this task
        BOLD = df{s,task};
        Phase_BOLD=zeros(N_areas,Tmax); 

        % Get the BOLD phase using the Hilbert transform
        for seed=1:N_areas
            BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:)); %for this region, demean the timecourse
            signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));                    
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
        end

        % Remove edge effects
        Phase_BOLD(:,[1:5,Tmax-4:Tmax])=[];
        for t=1:Tmax-10 %for each time point      
            iFC=zeros(N_areas); 
            for n=1:N_areas 
                for p=1:N_areas
                    iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % THIS THAT NEW SHIT
            % Here I'll compare the Eig to a set of vectors that represent
            % each of the 8 subcortical networks, and make a state
            % timeseries where the network w/ the strongest relationship
            % to the eig
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            % To assess strength of relationship between the iFC and the
            % template network matrices, I can pull the iFC values for the
            % intranetwork edges, and then rather than correlate them, find
            % the mean iFC among all regions in the network - that way it
            % also accounts for the number of regions. I'm also tempted to
            % go mean/std, though - would be a better measure of signal to
            % noise. However, std does not work with B, so we're going to
            % move on for now.
            
            for ii=1:num_states
            tmpMat=squeeze(templateMat(:,:,ii));
            B = iFC .* tmpMat;  % All values except those in the network template are set to zero
            B = triu(B,1); % Returns upper trangle without the diag
            B(B == 0) = NaN;  % makes 0's nan so we can ignore them and it doesn't drag down the overall mean
            storeMatStr(ii,1)=mean(B,[1 2],"omitmissing"); % this gives the mean over all elements in the matrix
            end
            [~,StateTimeDFMat{s,task}(t,1)]=max(abs(storeMatStr));

            % Ok, back to Eig extraction now
            [V1,~]=eigs(iFC,1);
            
            % Since we now don't care about communities, we can take the
            % abs value - this now represent each region's strength of
            % contribution to the FC state 
            V1=abs(V1);

            storeCor=zeros(num_states,1);
            for ii=1:size(templateEig,2)
            [storeCor(ii,1),~]=corr(V1,templateEig(:,ii),'type','Spearman');
            end

            [~,StateTimeDFEig{s,task}(t,1)]=max(abs(storeCor));

            % Save V1 from all frames in all fMRI sessions in Leading eig
            t_all=t_all+1; % Update time
            Time_all(:,t_all)=[s task]; % Information that at t_all, V1 corresponds to subject s in a given task
        end
    end
end

%%
% Decide here whether you want to do it with Matrix or Eig - I think I'll
% do matrix for now 
StateTimeDF=StateTimeDFMat;


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Complexity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results=[];
Transitions = {};
TransitionsH = [];
nbits=4; % rather than have a variable number of states that needed diff bit encoding, 4-bit encoding can represent our 9 numbers

for subj=1:n_Subjects
    fprintf('.');
    for i=1:num_condi
        S = StateTimeDF{subj,i};

        % Convert the string into a binary string with n-bit encoding,
        % where n=the number of bits to compute K
        % For a group of n bits, it is possible to represent 2^n values
        S0=S-1;  % This gives us state 0 (instead of starting at 1)

        Sb = ctx_dec2bi1d(S0, nbits);

        % Calculate LZC, normalise by LZC on a single shuffled version of the data
        % (could divide by mean of N shuffled versions) - Greg isn't sure
        % whether there's a need to shuffle. 
        Results(i,subj,1) = calc_lz_complexity(Sb, 'exhaustive', false) ./ ...
            calc_lz_complexity(Sb(randperm(length(Sb))), 'exhaustive', false);

        % Convert the binary string into a Matlab string format for BDM1d
        Sstr = ctx_bi2str(Sb);
        % Calculate BDM1d
        
        Results(i,subj,2) = ctx_stringbdm1d(Sstr);
        % Calculate transition entropy
        for trorder = 0:4
            Transitions{subj,i} = markovmodel(S, trorder);
            n = Transitions{subj,i}(:,end);
            p = n./sum(n);
            p = p(p>0);
            Results(i,subj,trorder+3)= -sum(p .*log2(p));
            % TransitionsH{CC}(trorder + 1, subj,i) = -sum(p .*log2(p));

        end
    end
end
%%
% Assess whether there's a sig effect of group on Information theroetic
% metrics - 1 = LZC, 2 = BDMC, 3-7 = 0-4th order TE
for ii=1:7
tmp=squeeze(Results(:,:,ii)');
[ComplexRes(ii,1),~,stats]=signrank(tmp(:,1),tmp(:,2));
if isfield(stats,'zval')==1
ComplexRes(ii,2)=stats.zval;
else
ComplexRes(ii,2)=0;
end
end


%%
% cd(outputpath)
% save('HomeBrewLEiDA_MID_Matrix.mat','-v7.3')
% save('HomeBrewLEiDA_CUE_Matrix.mat','-v7.3')
% save('HomeBrewLEiDA_MID_Eig.mat','-v7.3')
save('HomeBrewLEiDA_CUE_Eig.mat','-v7.3')
% Just a note - state order, 1-9, is the same for Eigs and Mats. is Vis,
% Somat, DorsAttn, SalAttn, Lim, Cont, DMN, Tmp, VMN
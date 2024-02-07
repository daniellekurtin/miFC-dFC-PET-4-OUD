%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are vars that should probs be input when running this as a function
% md_subj=26;
% hc_subj=22;
n_sess=2;
num_rois=214;

% Add paths
baseloc="C:\Users\dk818\OneDrive - Imperial College London\NCORE\miFCdFC";

% For outputs
str='\Outputs';
outputpath=strcat(baseloc,str);
addpath(outputpath);

% For inputs
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CUE
Load and organise data
Organise data from how it's extracted from fslmeants to a friendlier format
This will load each condition, combine cortex and subcortex files, and put
it all in a BigStruct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

tmpdf={};
cueFiles_Cor=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*Cortex*'));
n_subj=size(cueFiles_Cor,1);

cueFiles_l_acc=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*l_acc*'));
cueFiles_l_amy=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*l_amy*'));
cueFiles_l_cau=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*l_cau*'));
cueFiles_l_hip=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*l_hip*'));
cueFiles_l_pal=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*l_pal*'));
cueFiles_l_put=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*l_put*'));
cueFiles_l_thl=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*l_thl*'));
cueFiles_r_nac=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*r_nac*'));
cueFiles_r_amy=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*r_amy*'));
cueFiles_r_cau=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*r_cau*'));
cueFiles_r_hip=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*r_hip*'));
cueFiles_r_pal=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*r_pal*'));
cueFiles_r_put=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*r_put*'));
cueFiles_r_thl=dir(strcat(baseloc,'\ExtractedTimeseries\*cue*r_thl*'));

for subj=1:n_subj
% Load and transpose cortex 
Cor = load(cueFiles_Cor(subj).name);
Cor = rows2vars(array2table(Cor)) ;
Cor(:,1) = [];
    
% Now the l_acc
l_acc=load(cueFiles_l_acc(subj).name);
l_acc = rows2vars(array2table(l_acc)) ;
l_acc(:,1) = [];

% Now the l_amy
l_amy=load(cueFiles_l_amy(subj).name);
l_amy = rows2vars(array2table(l_amy)) ;
l_amy(:,1) = [];

% Now the l_cau
l_cau=load(cueFiles_l_cau(subj).name);
l_cau = rows2vars(array2table(l_cau)) ;
l_cau(:,1) = [];

% Now the l_hip
l_hip=load(cueFiles_l_hip(subj).name);
l_hip = rows2vars(array2table(l_hip)) ;
l_hip(:,1) = [];

% Now the l_pal
l_pal=load(cueFiles_l_pal(subj).name);
l_pal = rows2vars(array2table(l_pal)) ;
l_pal(:,1) = [];

% Now the l_put
l_put=load(cueFiles_l_put(subj).name);
l_put = rows2vars(array2table(l_put)) ;
l_put(:,1) = [];

% Now the l_thl
l_thl=load(cueFiles_l_thl(subj).name);
l_thl = rows2vars(array2table(l_thl)) ;
l_thl(:,1) = [];

% Now the r_nac
r_nac=load(cueFiles_r_nac(subj).name);
r_nac = rows2vars(array2table(r_nac)) ;
r_nac(:,1) = [];

% Now the r_amy
r_amy=load(cueFiles_r_amy(subj).name);
r_amy = rows2vars(array2table(r_amy)) ;
r_amy(:,1) = [];

% Now the r_cau
r_cau=load(cueFiles_r_cau(subj).name);
r_cau = rows2vars(array2table(r_cau)) ;
r_cau(:,1) = [];

% Now the r_hip
r_hip=load(cueFiles_r_hip(subj).name);
r_hip = rows2vars(array2table(r_hip)) ;
r_hip(:,1) = [];

% Now the r_pal
r_pal=load(cueFiles_r_pal(subj).name);
r_pal = rows2vars(array2table(r_pal)) ;
r_pal(:,1) = [];

% Now the r_put
r_put=load(cueFiles_r_put(subj).name);
r_put = rows2vars(array2table(r_put)) ;
r_put(:,1) = [];

% Now the r_thl
r_thl=load(cueFiles_r_thl(subj).name);
r_thl = rows2vars(array2table(r_thl)) ;
r_thl(:,1) = [];

% % Stack 'em
stack = [table2array(Cor);table2array(l_acc);table2array(l_amy);table2array(l_cau);table2array(l_hip);table2array(l_pal);table2array(l_put);table2array(l_thl);table2array(r_nac);table2array(r_amy);table2array(r_cau);table2array(r_hip);table2array(r_pal);table2array(r_put);table2array(r_thl)] ;
tmpdf{subj,1}=stack;

end

% Now divide df into MD and HC pars - first 25 = MD pars
df_cue={};
h1=0;
h2=0;
for ii=1:n_subj
if ii < 26
    h1=h1+1;
    df_cue{h1,1}=tmpdf{ii,1};
else
    h2=h2+1;
    df_cue{h2,2}=tmpdf{ii,1};
end
end

%%

tmpdf={};
midFiles_Cor=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*Cortex*'));
n_subj=size(midFiles_Cor,1);

midFiles_l_acc=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*l_acc*'));
midFiles_l_amy=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*l_amy*'));
midFiles_l_cau=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*l_cau*'));
midFiles_l_hip=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*l_hip*'));
midFiles_l_pal=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*l_pal*'));
midFiles_l_put=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*l_put*'));
midFiles_l_thl=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*l_thl*'));
midFiles_r_nac=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*r_nac*'));
midFiles_r_amy=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*r_amy*'));
midFiles_r_cau=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*r_cau*'));
midFiles_r_hip=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*r_hip*'));
midFiles_r_pal=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*r_pal*'));
midFiles_r_put=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*r_put*'));
midFiles_r_thl=dir(strcat(baseloc,'\ExtractedTimeseries\*mid*r_thl*'));

for subj=1:n_subj
% Load and transpose cortex 
Cor = load(midFiles_Cor(subj).name);
Cor = rows2vars(array2table(Cor)) ;
Cor(:,1) = [];
    
% Now the l_acc
l_acc=load(midFiles_l_acc(subj).name);
l_acc = rows2vars(array2table(l_acc)) ;
l_acc(:,1) = [];

% Now the l_amy
l_amy=load(midFiles_l_amy(subj).name);
l_amy = rows2vars(array2table(l_amy)) ;
l_amy(:,1) = [];

% Now the l_cau
l_cau=load(midFiles_l_cau(subj).name);
l_cau = rows2vars(array2table(l_cau)) ;
l_cau(:,1) = [];

% Now the l_hip
l_hip=load(midFiles_l_hip(subj).name);
l_hip = rows2vars(array2table(l_hip)) ;
l_hip(:,1) = [];

% Now the l_pal
l_pal=load(midFiles_l_pal(subj).name);
l_pal = rows2vars(array2table(l_pal)) ;
l_pal(:,1) = [];

% Now the l_put
l_put=load(midFiles_l_put(subj).name);
l_put = rows2vars(array2table(l_put)) ;
l_put(:,1) = [];

% Now the l_thl
l_thl=load(midFiles_l_thl(subj).name);
l_thl = rows2vars(array2table(l_thl)) ;
l_thl(:,1) = [];

% Now the r_nac
r_nac=load(midFiles_r_nac(subj).name);
r_nac = rows2vars(array2table(r_nac)) ;
r_nac(:,1) = [];

% Now the r_amy
r_amy=load(midFiles_r_amy(subj).name);
r_amy = rows2vars(array2table(r_amy)) ;
r_amy(:,1) = [];

% Now the r_cau
r_cau=load(midFiles_r_cau(subj).name);
r_cau = rows2vars(array2table(r_cau)) ;
r_cau(:,1) = [];

% Now the r_hip
r_hip=load(midFiles_r_hip(subj).name);
r_hip = rows2vars(array2table(r_hip)) ;
r_hip(:,1) = [];

% Now the r_pal
r_pal=load(midFiles_r_pal(subj).name);
r_pal = rows2vars(array2table(r_pal)) ;
r_pal(:,1) = [];

% Now the r_put
r_put=load(midFiles_r_put(subj).name);
r_put = rows2vars(array2table(r_put)) ;
r_put(:,1) = [];

% Now the r_thl
r_thl=load(midFiles_r_thl(subj).name);
r_thl = rows2vars(array2table(r_thl)) ;
r_thl(:,1) = [];

% % Stack 'em
stack = [table2array(Cor);table2array(l_acc);table2array(l_amy);table2array(l_cau);table2array(l_hip);table2array(l_pal);table2array(l_put);table2array(l_thl);table2array(r_nac);table2array(r_amy);table2array(r_cau);table2array(r_hip);table2array(r_pal);table2array(r_put);table2array(r_thl)] ;
tmpdf{subj,1}=stack;

end

% Now divide df into MD and HC pars - first 25 = MD pars
df_mid={};
h1=0;
h2=0;
for ii=1:n_subj
if ii ~= 14
    if ii < 27
        h1=h1+1;
        df_mid{h1,1}=tmpdf{ii,1};
    else
        h2=h2+1;
        df_mid{h2,2}=tmpdf{ii,1};
    end
end
end

%%
% need to even out the playing field for both HC and MD pars
for ii=23:26
df_mid{ii,1}=[];
df_cue{ii,1}=[];
end

%%
% % Save all this sheeeeeeet
cd(outputpath)
clear cue*
clear mid*
clear l_*
clear r_*
clear stack str tmpdf h1 h2 ii n_subj subj
save('df_forLEIDA.mat','-v7.3');

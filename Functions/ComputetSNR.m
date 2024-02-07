%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup to compute tSNR
% tSNR was computed using a shell script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currloc="C:\Users\dk818\OneDrive - Imperial College London\NCORE\miFCdFC";

% For ICC, violinplot, and colors
str='\Functions';
addpath(strcat(currloc,str));

% For outputs
str='\Outputs';
outputpath=strcat(currloc,str);
addpath(outputpath);

% For atlas and other data
str='\AtlasAndOtherInputs';
addpath(strcat(currloc,str));

% For tSNR
str='\tSNR';
addpath(strcat(currloc,str));

% Read in the table
tmp=readtable('tSNROutput_Final.csv');
mid=tmp(1:2:end,:);
cue=tmp(2:2:end,:);

% make df w/ MID md, hc, then Cue md, hc
midPT=table2array(mid([1:11,23:33],1));
midHC=table2array(mid([12:22,34:44],1));
cuePT=table2array(cue([1:11,23:33],1));
cueHC=table2array(cue([12:22,34:44],1));

df=[midPT(:,1),midHC(:,1),cuePT(:,1),cueHC(:,1)];

[p,h,stats] = signrank(df(:,1),df(:,2))

[p,h,stats] = signrank(df(:,3),df(:,4))

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=linspecer(2);

figure()
hold on
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
violins = violinplot(df);
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
sigstar({[1,2]})
% 
% saveas(gcf,'tSNR.png')
% close all

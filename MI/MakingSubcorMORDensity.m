%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MakingSubcorMORDensity
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copy the FSLmeants output into an array
Output=[0.3902335925,1.534402468,1.470731698,1.957733244,2.637560248,0.6445972312,1.438231096,1.602469022,0.4547576417,1.526716828,1.649634442,1.992754721,2.462172308,0.6413373909,1.474345219,1.404098817];
% Now need to avg the thalamus ones - 
RHThl=(Output(1,3)+Output(1,4))/2;
LHThl=(Output(1,11)+Output(1,12))/2;

% L And R Thl = 0
dictionary=[13,10,16,9,14,15,0,5,2,8,1,6,7,0];

% Now let's make an array with our numbers 
for ii=1:14
if ii ~= 7 && ii~=14
SubCorMOR(ii,1)=Output(1,dictionary(1,ii));
elseif ii==7
SubCorMOR(ii,1)=LHThl;
else 
SubCorMOR(ii,1)=RHThl;
end
end    

% Now to scale it - 
% The mean and std for the other ROIs = 21 and 8, w/ the max = 38, just
% over 2 std from the mean. 
% This is a bit hacky, but I can rescale all subcor ROIs to be within 1 or
% 2 stds of the mean; i.e., between 29-37
ReScaleSubCorMOR=rescale(SubCorMOR,29,37);

% Now let's load the MOR table 
MORTbl=readtable('MU_carfentanil_hc204_kantonen.csv');
MORTbl=table2array(MORTbl);
MORTbl=MORTbl';

% And add the new cols
MORTbl(201:214,1)=ReScaleSubCorMOR;

% Zscore the whole thing
MORTbl=zscore(MORTbl);

% Rescale
MORTbl=rescale(MORTbl);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Output=[1.075796774,2.388936295,2.150033539,2.950152966,9.870969236,6.057843045,10.31921246,9.732880904,1.011703029,2.237421249,2.279978015,3.007171287,8.562547102,4.069452433,9.907230421,8.45111662];
% Now need to avg the thalamus ones - 
RHThl=(Output(1,3)+Output(1,4))/2;
LHThl=(Output(1,11)+Output(1,12))/2;

% L And R Thl = 0
dictionary=[13,10,16,9,14,15,0,5,2,8,1,6,7,0];

% Now let's make an array with our numbers 
for ii=1:14
if ii ~= 7 && ii~=14
SubCorMOR(ii,1)=Output(1,dictionary(1,ii));
elseif ii==7
SubCorMOR(ii,1)=LHThl;
else 
SubCorMOR(ii,1)=RHThl;
end
end    

% Now to scale it - 
% The mean and std for the other ROIs = 21 and 8, w/ the max = 38, just
% over 2 std from the mean. 
% This is a bit hacky, but I can rescale all subcor ROIs to be within 1 or
% 2 stds of the mean; i.e., between 29-37
ReScaleSubCorMOR=rescale(SubCorMOR,0.3,2.58);

D2Tbl1=readtable('D2_fallypride_hc49_jaworska.csv');
D2Tbl1=table2array(D2Tbl1);
D2Tbl1=D2Tbl1';
D2Tbl2=readtable('D2_flb457_hc37_smith.csv');
D2Tbl2=table2array(D2Tbl2);
D2Tbl2=D2Tbl2';
D2Tbl3=readtable('D2_flb457_hc55_sandiego.csv');
D2Tbl3=table2array(D2Tbl3);
D2Tbl3=D2Tbl3';

D2Tbl=(D2Tbl1+D2Tbl2+D2Tbl3)/3;

% And add the new cols
D2Tbl(201:214,1)=ReScaleSubCorMOR;

% Zscore the whole thing
D2Tbl=zscore(D2Tbl);

% Rescale
D2Tbl=rescale(D2Tbl);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:num_rois
    for jj=1:ii
        HeatMapSum(ii,jj)=MORTbl(ii,1)+MORTbl(jj,1);
    end
end

for ii=1:num_rois
    for jj=1:ii
        HeatMapRatio1(ii,jj)=MORTbl(ii,1)/MORTbl(jj,1);
    end
end

for ii=1:num_rois
    for jj=1:ii
        HeatMapRatio2(ii,jj)=MORTbl(jj,1)/MORTbl(ii,1);
    end
end

SigmiFC=zeros(num_rois,num_rois);
for ii=1:size(GroupResTbl,1)
SigmiFC(GroupResTbl(ii,1),GroupResTbl(ii,2))=1;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:num_rois
    for jj=1:ii
        D2HeatMapSum(ii,jj)=D2Tbl(ii,1)+D2Tbl(jj,1);
    end
end

for ii=1:num_rois
    for jj=1:ii
        D2HeatMapRatio1(ii,jj)=D2Tbl(ii,1)/D2Tbl(jj,1);
    end
end

for ii=1:num_rois
    for jj=1:ii
        D2HeatMapRatio2(ii,jj)=D2Tbl(jj,1)/D2Tbl(ii,1);
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOR SANITY CHECK 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Because I ran the below w/ HeatMatSum and D2HeatMapSum, I want to make
% sure my results replicate for the ratio ones. Ratio1 results are same as
% Ratio 2.
% HeatMapSum=HeatMapRatio1;
% D2HeatMapSum=HeatMapRatio1;

% Also, this needs to be added a few sections above, bc if the MOR or D2 rescaled = 0, 
% then the ratio becomes infinity
% for the sanity check - make anything 0 a 0.01
MORTbl(7,1)=0.01;
D2Tbl(204,1)=0.01;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make fig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now to color this - 
addpath('C:\Users\dk818\OneDrive - Imperial College London\LilThingsAndMATLABPrograms\ColorBrewer')

c=brewermap([],"RdPu");
h=size(c,1);
C=[1,1,1];
C(2:h+1,:)=c;

% Cols 1, 2, and 3 = miFC, MOR, D2
HighMD=[];
LowMD=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure()
% subplot(1,2,1)
% imagesc(D2HeatMapSum,'AlphaData',0.5) 
% colormap(C) 
% hold on
for ii=1:num_rois
    for jj=1:ii
        % Red means MD miFC > HC miFC, blue means HC mifc > MD mifc
        if GroupSigETASurvivepadded(ii,jj) > 0 
            % scatter(jj,ii,30,[0.69 0 0.16],'filled','square')
            h=size(HighMD,1)+1;

            %%%%%%%%%%%%%%%%%%%
            % NOTE!!!! 
            % I AM CHANGING ETA TO BE A PROPER ETA! 
            % Before, I just had the z-stat; now, I'm dividing it by the
            % square root of the total sample size 
            HighMD(h,1)=(GroupSigETASurvivepadded(ii,jj))/sqrt(n_subj*2);
            HighMD(h,2)=HeatMapSum(ii,jj);
            HighMD(h,3)=D2HeatMapSum(ii,jj);

        elseif GroupSigETASurvivepadded(ii,jj) < 0
            % scatter(jj,ii,30,[0.02 0.12 1],'filled','square')
            h=size(LowMD,1)+1;
            LowMD(h,1)=(GroupSigETASurvivepadded(ii,jj))/sqrt(n_subj*2);
            LowMD(h,2)=HeatMapSum(ii,jj);
            LowMD(h,3)=D2HeatMapSum(ii,jj);

        end
    end
end
% % Now to add lines - see end of this script to delineate network boundaries
% for ii=[13,29,40,51,57,75,99,101,113,131,142,157,165,184,197,201]
%     YLine=ones(ii,1)*ii;
%     XLine=[1:1:ii];
%     plot(XLine,YLine,'Color',[62/256 36/256 100/256],'LineWidth',1)
% 
%     XLine2=ones(num_rois-ii,1)*ii;
%     YLine2=[1+ii:1:num_rois];
%     plot(XLine2,YLine2,'Color',[62/256 36/256 100/256],'LineWidth',1)
% end
% axis tight
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% subplot(1,2,2)
% imagesc(D2HeatMapSum,'AlphaData',0.5) 
% colormap(C) 
% hold on
% for ii=1:num_rois
%     for jj=1:ii
%         % Red means MD miFC > HC miFC, blue means HC mifc > MD mifc
%         if GroupSigETASurvivepadded(ii,jj) > 0 
%             scatter(jj,ii,30,[0.69 0 0.16],'filled','square')
% 
%         elseif GroupSigETASurvivepadded(ii,jj) < 0
%             scatter(jj,ii,30,[0.02 0.12 1],'filled','square')
% 
%         end
%     end
% end
% % Now to add lines - see end of this script to delineate network boundaries
% for ii=[13,29,40,51,57,75,99,101,113,131,142,157,165,184,197,201]
%     YLine=ones(ii,1)*ii;
%     XLine=[1:1:ii];
%     plot(XLine,YLine,'Color',[62/256 36/256 100/256],'LineWidth',1)
% 
%     XLine2=ones(num_rois-ii,1)*ii;
%     YLine2=[1+ii:1:num_rois];
%     plot(XLine2,YLine2,'Color',[62/256 36/256 100/256],'LineWidth',1)
% end
% axis tight
% hold off


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  relationship between sig effects and MOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now make plots for the
% FOR CUE CB - For some reason the eta val for HighMD, first res, is nuts - gonna just
% make it a little higher than the max
HighMD(1,1)=0.5;

% MOR
figure()
hold on
scatter(HighMD(:,1),HighMD(:,2),[],'filled','MarkerFaceColor',[0.69 0 0.16], 'MarkerFaceAlpha',0.5);
scatter(abs(LowMD(:,1)),LowMD(:,2),[],'filled','MarkerFaceColor',[0.02 0.12 1],'MarkerFaceAlpha',0.5)
h=lsline();
h(2).Color=[0.69 0 0.16];
h(1).Color=[0.02 0.12 1];
h(1).LineWidth=1;
h(2).LineWidth=1;

% % Stats time! 
disp('Strength of rel btwn MOR Receptors and miFC when MD>HC')
[r,p]=corr(HighMD(:,1),HighMD(:,2),'type','Spearman')

% When you rerun this for sanity checks, load and rerun for MID, RA, Cue,
% and CB 
% SanityTable(1,1)=p;
% SanityTable(1,2)=r;

disp('Strength of rel btwn MOR Receptors and miFC when HC>MD')
[r,p]=corr(LowMD(:,1),LowMD(:,2),'type','Spearman')

% SanityTable(2,1)=p;
% SanityTable(2,2)=r;

% saveas(gcf,'CB_RelBtwnMORAndMIFC.png')
% saveas(gcf,'CUE_RelBtwnMORAndMIFC.png')
% saveas(gcf,'MID_RelBtwnMORAndMIFC.png')
% saveas(gcf,'RA_RelBtwnMORAndMIFC.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  relationship between sig effects and D2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% D2
figure()
hold on
scatter(HighMD(:,1),HighMD(:,3),[],'filled','MarkerFaceColor',[0.69 0 0.16], 'MarkerFaceAlpha',0.5);
scatter(abs(LowMD(:,1)),LowMD(:,3),[],'filled','MarkerFaceColor',[0.02 0.12 1],'MarkerFaceAlpha',0.5)
h=lsline();
h(2).Color=[0.69 0 0.16];
h(1).Color=[0.02 0.12 1];
h(1).LineWidth=1;
h(2).LineWidth=1;

% % Stats time! 
disp('Strength of rel btwn D2 Receptors and miFC when MD>HC')
[r,p]=corr(HighMD(:,1),HighMD(:,3),'type','Spearman')
% SanityTable(1,3)=p;
% SanityTable(1,4)=r;


disp('Strength of rel btwn D2 Receptors and miFC when HC>MD')
[r,p]=corr(LowMD(:,1),LowMD(:,3),'type','Spearman')
% SanityTable(2,3)=p;
% SanityTable(2,4)=r;

% saveas(gcf,'CB_RelBtwnD2AndMIFC.png')
% saveas(gcf,'CUE_RelBtwnD2AndMIFC.png')
% saveas(gcf,'MID_RelBtwnD2AndMIFC.png')
% saveas(gcf,'RA_RelBtwnD2AndMIFC.png')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Old version of making many subplots, and the legend 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath('C:\Users\dk818\OneDrive - Imperial College London\LilThingsAndMATLABPrograms\ColorBrewer')
% 
% c=brewermap([],"YlGn");
% h=size(c,1);
% C=[1,1,1];
% C(2:h+1,:)=c;
% 
% figure()
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % MOR
% % subplot(2,3,1)
% subplot(1,2,1)
% imagesc(HeatMapSum,'AlphaData',0.5) 
% colormap(C) 
% % colorbar()
% hold on
% for ii=1:num_rois
%     for jj=1:ii
%         % Red means MD miFC > HC miFC, blue means HC mifc > MD mifc
%         if GroupSigETASurvivepadded(ii,jj) > 0 
%             scatter(jj,ii,30,[0.69 0 0.16],'filled','square')
%         elseif GroupSigETASurvivepadded(ii,jj) < 0
%             scatter(jj,ii,30,[0.02 0.12 1],'filled','square')
% 
%         end
%     end
% end
% % Now to add lines - see end of this script to delineate network boundaries
% for ii=[13,29,40,51,57,75,99,101,113,131,142,157,165,184,197,201]
%     YLine=ones(ii,1)*ii;
%     XLine=[1:1:ii];
%     plot(XLine,YLine,'Color',[62/256 36/256 100/256],'LineWidth',1)
% 
%     XLine2=ones(num_rois-ii,1)*ii;
%     YLine2=[1+ii:1:num_rois];
%     plot(XLine2,YLine2,'Color',[62/256 36/256 100/256],'LineWidth',1)
% end
% % colorbar()
% axis tight
% hold off
% 
% subplot(2,3,2)
% imagesc(HeatMapRatio1,'AlphaData',0.5) 
% colormap(C) 
% hold on
% for ii=1:num_rois
%     for jj=1:ii
%         if SigmiFC(ii,jj) > 0
%             scatter(jj,ii,30,[205/256 92/256 92/256],'filled','square')
%         end
%     end
% end
% for ii=[13,29,40,51,57,75,99,101,113,131,142,157,165,184,197,201]
%     YLine=ones(ii,1)*ii;
%     XLine=[1:1:ii];
%     plot(XLine,YLine,'Color',[62/256 36/256 100/256],'LineWidth',1)
% 
%     XLine2=ones(num_rois-ii,1)*ii;
%     YLine2=[1+ii:1:num_rois];
%     plot(XLine2,YLine2,'Color',[62/256 36/256 100/256],'LineWidth',1)
% end
% % colorbar()
% axis tight
% hold off
% 
% subplot(2,3,3)
% imagesc(HeatMapRatio2,'AlphaData',0.5) 
% colormap(C) 
% hold on
% for ii=1:num_rois
%     for jj=1:ii
%         if SigmiFC(ii,jj) > 0
%             scatter(jj,ii,30,[205/256 92/256 92/256],'filled','square')
%         end
%     end
% end
% for ii=[13,29,40,51,57,75,99,101,113,131,142,157,165,184,197,201]
%     YLine=ones(ii,1)*ii;
%     XLine=[1:1:ii];
%     plot(XLine,YLine,'Color',[62/256 36/256 100/256],'LineWidth',1)
% 
%     XLine2=ones(num_rois-ii,1)*ii;
%     YLine2=[1+ii:1:num_rois];
%     plot(XLine2,YLine2,'Color',[62/256 36/256 100/256],'LineWidth',1)
% end
% % colorbar()
% axis tight
% hold off
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % D2
% % subplot(2,3,4)
% subplot(1,2,2)
% imagesc(D2HeatMapSum,'AlphaData',0.5) 
% colormap(C) 
% hold on
% for ii=1:num_rois
%     for jj=1:ii
%         if SigmiFC(ii,jj) > 0
%             scatter(jj,ii,30,[205/256 92/256 92/256],'filled','square')
%         end
%     end
% end
% % Now to add lines - see end of this script to delineate network boundaries
% for ii=[13,29,40,51,57,75,99,101,113,131,142,157,165,184,197,201]
%     YLine=ones(ii,1)*ii;
%     XLine=[1:1:ii];
%     plot(XLine,YLine,'Color',[62/256 36/256 100/256],'LineWidth',1)
% 
%     XLine2=ones(num_rois-ii,1)*ii;
%     YLine2=[1+ii:1:num_rois];
%     plot(XLine2,YLine2,'Color',[62/256 36/256 100/256],'LineWidth',1)
% end
% % colorbar()
% axis tight
% hold off
% 
% subplot(2,3,5)
% imagesc(D2HeatMapRatio1,'AlphaData',0.5) 
% colormap(C) 
% hold on
% for ii=1:num_rois
%     for jj=1:ii
%         if SigmiFC(ii,jj) > 0
%             scatter(jj,ii,30,[205/256 92/256 92/256],'filled','square')
%         end
%     end
% end
% for ii=[13,29,40,51,57,75,99,101,113,131,142,157,165,184,197,201]
%     YLine=ones(ii,1)*ii;
%     XLine=[1:1:ii];
%     plot(XLine,YLine,'Color',[62/256 36/256 100/256],'LineWidth',1)
% 
%     XLine2=ones(num_rois-ii,1)*ii;
%     YLine2=[1+ii:1:num_rois];
%     plot(XLine2,YLine2,'Color',[62/256 36/256 100/256],'LineWidth',1)
% end
% % colorbar()
% axis tight
% hold off
% 
% subplot(2,3,6)
% imagesc(D2HeatMapRatio2,'AlphaData',0.5) 
% colormap(C) 
% hold on
% for ii=1:num_rois
%     for jj=1:ii
%         if SigmiFC(ii,jj) > 0
%             scatter(jj,ii,30,[205/256 92/256 92/256],'filled','square')
%         end
%     end
% end
% for ii=[13,29,40,51,57,75,99,101,113,131,142,157,165,184,197,201]
%     YLine=ones(ii,1)*ii;
%     XLine=[1:1:ii];
%     plot(XLine,YLine,'Color',[62/256 36/256 100/256],'LineWidth',1)
% 
%     XLine2=ones(num_rois-ii,1)*ii;
%     YLine2=[1+ii:1:num_rois];
%     plot(XLine2,YLine2,'Color',[62/256 36/256 100/256],'LineWidth',1)
% end
% % colorbar()
% axis tight
% hold off
% 
% % name=('CUE_MORD2.png');
% % saveas(gcf, name)
% % close all


%{
ROINums per network to show edges
LH
Vis 1:12
SoMat 13:28
DorsAttn 29:39
SalVent 40:50
Limbic 51:56
Cont 57:74
Default 75:98
TmpPar 99:100

RH
Vis 101:112
SoMat 113:130
DorsAttn 131:141
SalVent 142:156
Limb 157:164
Cont 165:183
DMN 184:196
TmpPar 197:200
Subcortical 201:214
%}

%%
%{
Shen ROI label link = https://github.com/yetianmed/subcortex/blob/master/Group-Parcellation/3T/Subcortex-Only/Tian_Subcortex_S1_3T_label.txt#L1

1 HIP-rh 
2 AMY-rh
3 pTHA-rh
4 aTHA-rh
5 NAc-rh
6 GP-rh
7 PUT-rh
8 CAU-rh
9 HIP-lh
10 AMY-lh
11 pTHA-lh
12 aTHA-lh
13 NAc-lh
14 GP-lh
15 PUT-lh
16 CAU-lh

My order = Shen place num
1 LH Nac = 13
2 LH Amy = 10 
3 LH Cau = 16
4 LH Hip = 9 
5 LH Pal = 14
6 LH Put = 15
7 LH Thl = 3, LHThl, made below
8 RH Nac = 5
9 RH Amy = 2 
10 RH Cau = 8
11 RH Hip = 1
12 RH Pal = 6
13 RH Put = 7
14 RH Thl = 11, RHThl

dictionary = [13,10,16,9,14,15,LHThl,5,2,8,1,6,7,RHthl]
%}
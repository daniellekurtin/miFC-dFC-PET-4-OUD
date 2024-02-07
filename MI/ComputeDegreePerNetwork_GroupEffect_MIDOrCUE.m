%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This gives you ROI-Num network key, to identify which networks ROIs are
% in, the empty arrays to count number of results. It counts all subcor
% regions as SubCor
load UseToComputeDegPerNet.mat
load templateEigAndMat.mat

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making a table to store results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All I actually need is ROI1, ROI2
% and then what's computed first is 
% NetworkID 17 Networks ROI, NetworkID 17 Networks ROI, 
% and then the same but for Network ID 7 Networks next
WorkingTable(:,1:2)=GroupResTbl(:,1:2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Counting ROI1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:length(WorkingTable)
   
    if ismember(WorkingTable(ii,1),VisCent) == 1
       VisCentCount=VisCentCount+1 ;
%        VisCentCxd(height(VisCentCxd)+1,1) =  WorkingTable(ii,4);
        % Let's add a col with the Network ID for 17 networks and Reduced
        % Netowkrs
       WorkingTable(ii,3)=1;
       WorkingTable(ii,5)=1;
       
   elseif ismember(WorkingTable(ii,1),VisPeri) == 1
       VisPeriCount=VisPeriCount+1 ;
%        VisPeriCxd(height(VisPeriCxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=2;
       WorkingTable(ii,5)=1;
       
   elseif ismember(WorkingTable(ii,1),SoMatA) == 1
       SoMatACount=SoMatACount+1 ;
%        SoMatACxd(height(SoMatACxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=3;
       WorkingTable(ii,5)=2;
        
   elseif ismember(WorkingTable(ii,1),SoMatB) == 1
       SoMatBCount=SoMatBCount+1 ;
%        SoMatBCxd(height(SoMatBCxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=4;
       WorkingTable(ii,5)=2;
        
   elseif ismember(WorkingTable(ii,1),DorsAttnA) == 1
       DorsAttnACount=DorsAttnACount+1  ;
%        DorsAttnACxd(height(DorsAttnACxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=5;
       WorkingTable(ii,5)=3;
        
   elseif ismember(WorkingTable(ii,1),DorsAttnB) == 1
       DorsAttnBCount=DorsAttnBCount+1 ;
%        DorsAttnBCxd(height(DorsAttnBCxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=6;
       WorkingTable(ii,5)=3;
        
   elseif ismember(WorkingTable(ii,1),SalVentAttnA) == 1
       SalVentAttnACount=SalVentAttnACount+1 ;
%        SalVentAttnACxd(height(SalVentAttnACxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=7;
       WorkingTable(ii,5)=4;
        
   elseif ismember(WorkingTable(ii,1),SalVentAttnB) == 1
       SalVentAttnBCount=SalVentAttnBCount+1  ;
%        SalVentAttnBCxd(height(SalVentAttnBCxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=8;
       WorkingTable(ii,5)=4;
        
   elseif ismember(WorkingTable(ii,1),LimB) == 1
       LimBCount=LimBCount+1 ;
%        LimBCxd(height(LimBCxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=9;
       WorkingTable(ii,5)=5;
        
   elseif ismember(WorkingTable(ii,1),LimA) == 1
       LimACount=LimACount+1 ;
%        LimACxd(height(LimACxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=10;
       WorkingTable(ii,5)=5;
        
   elseif ismember(WorkingTable(ii,1),ContA) == 1
       ContACount=ContACount+1 ;    
%        ContACxd(height(ContACxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=11;
       WorkingTable(ii,5)=6;
        
   elseif ismember(WorkingTable(ii,1),ContB) == 1
       ContBCount=ContBCount+1 ;
%        ContBCxd(height(ContBCxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=12;
       WorkingTable(ii,5)=6;
        
   elseif ismember(WorkingTable(ii,1),ContC) == 1
       ContCCount=ContCCount+1 ;
%        ContCCxd(height(ContCCxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=13;
       WorkingTable(ii,5)=6;
        
   elseif ismember(WorkingTable(ii,1),DMNA) == 1
       DMNACount=DMNACount+1  ;
%        DMNACxd(height(DMNACxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=14;
       WorkingTable(ii,5)=7;
        
   elseif ismember(WorkingTable(ii,1),DMNB) == 1
       DMNBCount=DMNBCount+1 ;
%        DMNBCxd(height(DMNBCxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=15;
       WorkingTable(ii,5)=7;
        
   elseif ismember(WorkingTable(ii,1),DMNC) == 1
       DMNCCount=DMNCCount+1 ;
%        DMNCCxd(height(DMNCCxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=16;
       WorkingTable(ii,5)=7;
        
   elseif ismember(WorkingTable(ii,1),TmpPar) == 1
       TmpParCount=TmpParCount+1 ;    
%        TmpParCxd(height(TmpParCxd)+1,1) =  WorkingTable(ii,4);
       WorkingTable(ii,3)=17;
       WorkingTable(ii,5)=8;
       
    elseif ismember(WorkingTable(ii,1),SubCor) == 1
       SubCorCount=SubCorCount+1 ; 
       WorkingTable(ii,3)=18;
       WorkingTable(ii,5)=9;       

    else 
        disp(strcat('no match for ',num2str(ii)))
              
   end
       
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Counting ROI2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:length(WorkingTable)
   
    if ismember(WorkingTable(ii,2),VisCent) == 1
       VisCentCount=VisCentCount+1 ;
       WorkingTable(ii,4)=1;
       WorkingTable(ii,6)=1;
       
   elseif ismember(WorkingTable(ii,2),VisPeri) == 1
       VisPeriCount=VisPeriCount+1 ;
       WorkingTable(ii,4)=2;
       WorkingTable(ii,6)=1;
       
   elseif ismember(WorkingTable(ii,2),SoMatA) == 1
       SoMatACount=SoMatACount+1 ;
       WorkingTable(ii,4)=3;
       WorkingTable(ii,6)=2;
       
   elseif ismember(WorkingTable(ii,2),SoMatB) == 1
       SoMatBCount=SoMatBCount+1 ;
       WorkingTable(ii,4)=4;
       WorkingTable(ii,6)=2;
       
   elseif ismember(WorkingTable(ii,2),DorsAttnA) == 1
       DorsAttnACount=DorsAttnACount+1  ;
       WorkingTable(ii,4)=5;
       WorkingTable(ii,6)=3;
       
   elseif ismember(WorkingTable(ii,2),DorsAttnB) == 1
       DorsAttnBCount=DorsAttnBCount+1 ;
       WorkingTable(ii,4)=6;
       WorkingTable(ii,6)=3;
       
   elseif ismember(WorkingTable(ii,2),SalVentAttnA) == 1
       SalVentAttnACount=SalVentAttnACount+1 ;
       WorkingTable(ii,4)=7;
       WorkingTable(ii,6)=4;
       
   elseif ismember(WorkingTable(ii,2),SalVentAttnB) == 1
       SalVentAttnBCount=SalVentAttnBCount+1  ;
       WorkingTable(ii,4)=8;
       WorkingTable(ii,6)=4;
       
   elseif ismember(WorkingTable(ii,2),LimB) == 1
       LimBCount=LimBCount+1 ;
       WorkingTable(ii,4)=9;
       WorkingTable(ii,6)=5;
       
   elseif ismember(WorkingTable(ii,2),LimA) == 1
       LimACount=LimACount+1 ;
       WorkingTable(ii,4)=10;
       WorkingTable(ii,6)=5;
       
   elseif ismember(WorkingTable(ii,2),ContA) == 1
       ContACount=ContACount+1 ;    
       WorkingTable(ii,4)=11;
       WorkingTable(ii,6)=6;
       
   elseif ismember(WorkingTable(ii,2),ContB) == 1
       ContBCount=ContBCount+1 ;
       WorkingTable(ii,4)=12;
       WorkingTable(ii,6)=6;
       
   elseif ismember(WorkingTable(ii,2),ContC) == 1
       ContCCount=ContCCount+1 ;
       WorkingTable(ii,4)=13;
       WorkingTable(ii,6)=6;
       
   elseif ismember(WorkingTable(ii,2),DMNA) == 1
       DMNACount=DMNACount+1  ;
       WorkingTable(ii,4)=14;
       WorkingTable(ii,6)=7;
       
   elseif ismember(WorkingTable(ii,2),DMNB) == 1
       DMNBCount=DMNBCount+1 ;
       WorkingTable(ii,4)=15;
       WorkingTable(ii,6)=7;
       
   elseif ismember(WorkingTable(ii,2),DMNC) == 1
       DMNCCount=DMNCCount+1 ;
       WorkingTable(ii,4)=16;
       WorkingTable(ii,6)=7;
       
   elseif ismember(WorkingTable(ii,2),TmpPar) == 1
       TmpParCount=TmpParCount+1 ;    
       WorkingTable(ii,4)=17;
       WorkingTable(ii,6)=8;
       
    elseif ismember(WorkingTable(ii,2),SubCor) == 1
       SubCorCount=SubCorCount+1 ; 
       WorkingTable(ii,4)=18;
       WorkingTable(ii,6)=9;         
       
    else 
        disp(strcat('no match for ',num2str(ii)))
              
   end
       
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute degree per 17 network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we have our counts per network. Let's name them. 
NetList=[VisCentCount,VisPeriCount,SoMatACount,SoMatBCount,DorsAttnACount,DorsAttnBCount,SalVentAttnACount,SalVentAttnBCount,LimBCount,LimACount,ContACount,ContBCount,ContCCount,DMNACount,DMNBCount,DMNCCount,TmpParCount,SubCorCount]';
NetNames=["Visual Central","Visual Peripheral","Somatomotor A",'Somatomotor B','Dorsal Attention A','Dorsal Attention B','Salience/Ventral Attention A','Salience/Ventral Attention B','Limbic B','Limbic A','Control A','Control B','Control C','DMN A','DMN B','DMN C','Temporal Parietal','Subcortical']';

% Order the networks by degree - SortNetList
[Order,SortNetList]=sort(NetList);

% Use the order of the network degree to sort the names of the networks
% SortNetName=[];
% for ii=1:width(NetList)
%     SortNetName(1,ii)=NetNames(1,Order(1,ii));
% end

mm=mean(NetList)
md=median(NetList)
sd=std(NetList)

PropNetList=NetList/sum(NetList);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute degree per 7 network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RedNetList=[VisCentCount+VisPeriCount,SoMatACount+SoMatBCount,DorsAttnACount+DorsAttnBCount,SalVentAttnACount+SalVentAttnBCount,LimBCount+LimACount,ContACount+ContBCount+ContCCount,DMNACount+DMNBCount+DMNCCount,TmpParCount,SubCorCount]';
RedNetNames=["Visual","Somatomotor","Dorsal Attention","Salience/Ventral Attention","Limbic","Control","DMN","Temporal Parietal","Subcortical"]';

Redmm=mean(RedNetList)
Redmd=median(RedNetList)
Redsd=std(RedNetList)

PropRedNetList=RedNetList/sum(RedNetList);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make plots for 17 network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hmmm.... I'd like to make a matrix with the number of connections between
% networks listed. . . 
% WorkingTable cols  5 adn 6 are coded with the following NetNames
A=zeros(18);

for ii=1:length(WorkingTable)
for aa=1:18
for bb=1:18
    
    T=isequal([WorkingTable(ii,3),WorkingTable(ii,4)],[aa,bb]);
    
    if  T == logical(1);
        A(WorkingTable(ii,3),WorkingTable(ii,4))= A(WorkingTable(ii,3),WorkingTable(ii,4)) + 1; 
    end
    
end
end
end

% Because the above script treats a transition of 1-->2 and 2-->1 as
% unique, we get an asymmetric matrix. To combat this we can add the 
% upper and lower diags (after flipping and rotating them). 
AA=A;
Aup=triu(A);
Alow=tril(A,-1);
AlowFlip=rot90(flip(Alow),-1);

A=[];
A=AlowFlip+Aup;

C=linspecer(max(NetList));
C(1,:)=[1,1,1];
figure()
imagesc(A)
colormap(C)
colorbar()
xticks([1:1:18]);
yticks=([1:1:18]');
% xticklabels(NetNames)
% yticklabels(NetNames)
alpha(0.5)
% cd (outputpath)
% saveas(gcf,'Matrix_DrgEffect_17Networks.png');
% close all



% Note - removing within-network connections 
% C=linspecer(17);
% G = digraph(A);  % Can remove self-connections via G = digraph(A,'omitselfloops');
% pl_ = plot(G,'Layout','layered','EdgeAlpha',0.5);
% G.Edges.LWidths = (G.Edges.Weight/sum(G.Edges.Weight))*100;
% pl_.LineWidth = G.Edges.LWidths;
% pl_.EdgeColor=[0 0 0];
% pl_.ShowArrows=0;
% pl_.NodeColor=C;
% pl_.NodeLabel=["","","","","","","","","","","","","","","","",""];
% % pl_.MarkerSize=(PropNetList*100)+15;

% Now I have to make a color map where the first 50 vals map
C=linspecer(18);
NetNames2={"Visual Central","Visual Peripheral","Somatomotor A",'Somatomotor B','Dorsal Atten. A','Dorsal Atten. B','Sal/Vent. Atten. A','Sal/Vent. Atten. B','Limbic B','Limbic A','Control A','Control B','Control C','DMN A','DMN B','DMN C','Temp. Parietal','Subcortical'}';

% Plot some sheeeet
figure()
circularGraph(A,'Label',NetNames2,'Colormap',C)
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
% saveas(gcf,'MID_CirclePlot_GrpEffect_17Networks.png');
% saveas(gcf,'CUE_CirclePlot_GrpEffect_17Networks.png');
% saveas(gcf,'CB_CirclePlot_GrpEffect_17Networks.png');
% saveas(gcf,'RA_CirclePlot_GrpEffect_17Networks.png');
close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make plots for 7 network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=zeros(9);

for ii=1:length(WorkingTable)
for aa=1:9
for bb=1:9
    
    T=isequal([WorkingTable(ii,5),WorkingTable(ii,6)],[aa,bb]);
    
    if  T == logical(1);
        A(WorkingTable(ii,5),WorkingTable(ii,6))= A(WorkingTable(ii,5),WorkingTable(ii,6)) + 1; 
    end
    
end
end
end

AA=A;
Aup=triu(A);
Alow=tril(A,-1);
AlowFlip=rot90(flip(Alow),-1);

A=[];
A=AlowFlip+Aup;

% C=linspecer(max(RedNetList));
% C(1,:)=[1,1,1];
% figure()
% imagesc(A)
% colormap(C)
% colorbar()
% xticks([1:1:9]);
% yticks=([1:1:9]');
% xticklabels(RedNetNames)
% yticklabels(RedNetNames)
% alpha(0.5)
% cd (outputpath)
% saveas(gcf,'Matrix_DrgEffect_7Networks.png');
% close all

% Used this: https://stackoverflow.com/questions/17308565/draw-network-or-graph-from-matrix-in-matlab
% And this: https://uk.mathworks.com/help/matlab/ref/gplot.html
% To build the below
% N=9;
% coords = [cos(2*pi*(1:N)/N); sin(2*pi*(1:N)/N)]';
% [x,y]=gplot(A,coords)
% figure()
% a=plot(x,y)

% Note - removing within-network connections (20% of connections) 
% C=linspecer(8);
% G = digraph(A,'omitselfloops');  % Can remove self-connections via G = digraph(A,'omitselfloops');
% pl_ = plot(G,'Layout','layered','EdgeAlpha',0.5);
% G.Edges.LWidths = (G.Edges.Weight/sum(G.Edges.Weight))*100;
% pl_.LineWidth = G.Edges.LWidths;
% pl_.EdgeColor=[0 0 0];
% pl_.ShowArrows=0;
% pl_.NodeColor=C;
% pl_.NodeLabel=["","","","","","","",""];
% pl_.MarkerSize=(PropRedNetList*100)+15;

% C=linspecer(9);
C=zeros([9,3]);

RedNetNames2={"Vis","SoMat","DorsAttn","SalVent","Limbic","Control","DMN","TempPar","Subcor"}';
% Plot some sheeeet
figure()
circularGraph(A,'Label',RedNetNames2,'Colormap',C)
set(gca,'FontSize',20)
set(gca,'FontName','Arial')
% saveas(gcf,'MID_CirclePlot_GrpEffect_7Networks.png');
% saveas(gcf,'CUE_CirclePlot_GrpEffect_7Networks.png');
% saveas(gcf,'RA_CirclePlot_GrpEffect_7Networks.png');
% saveas(gcf,'CB_CirclePlot_GrpEffect_7Networks.png');
% close all

%%
% save("Checkpoint2_ComputeMI_MID.mat",'-v7.3')
% save("Checkpoint2_ComputeMI_CUE.mat",'-v7.3')
% save("Checkpoint2_ComputeMI_RA.mat",'-v7.3')
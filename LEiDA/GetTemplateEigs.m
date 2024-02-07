%%
% Other sandbox - making eigs for each of the networks. 
load UseToComputeDegPerNet.mat

% Get rid of what we don't need
clear ContACount ContACxd ContBCount ContBCxd ContCCount ContCCxd DMNACount DMNACxd DMNBCxd DMNBCount DMNCCount DMNCCxd LimACxd LimACount LimBCxd LimBCount SalVentAttnBCxd SalVentAttnBCount SalVentAttnACxd SalVentAttnACount SoMatBCxd SoMatBCount SoMatACxd SoMatACount TmpParCxd TmpParCount VisPeriCxd VisPeriCount VisCentCxd VisCentCount

num_rois=214;

Vis=[VisCent;VisPeri];
SoMat=[SoMatA;SoMatB];
DorsAttn=[DorsAttnA;DorsAttnB];
SalVent=[SalVentAttnA;SalVentAttnB];
Limb=[LimB;LimA];
Cont=[ContA;ContB;ContC];
DMN=[DMNA;DMNB;DMNC];
TmpPar;
VMN=[51,52,80,81,82,157,158,159,160,187,188,189,201,203,206,207,209,210,213,214]';

VisCentEig=zeros(num_rois,1);
for ii=1:size(Vis)
VisCentEig(Vis(ii,1),1)=1;
end

SoMatEig=zeros(num_rois,1);
for ii=1:size(SoMat)
SoMatEig(SoMat(ii,1),1)=1;
end

DorsAttnEig=zeros(num_rois,1);
for ii=1:size(DorsAttn)
DorsAttnEig(DorsAttn(ii,1),1)=1;
end

SalVentEig=zeros(num_rois,1);
for ii=1:size(SalVent)
SalVentEig(SalVent(ii,1),1)=1;
end

LimbEig=zeros(num_rois,1);
for ii=1:size(Limb)
LimbEig(Limb(ii,1),1)=1;
end

ContEig=zeros(num_rois,1);
for ii=1:size(Cont)
ContEig(Cont(ii,1),1)=1;
end

DMNEig=zeros(num_rois,1);
for ii=1:size(DMN)
DMNEig(DMN(ii,1),1)=1;
end

TmpParEig=zeros(num_rois,1);
for ii=1:size(TmpPar)
TmpParEig(TmpPar(ii,1),1)=1;
end

VMNEig=zeros(num_rois,1);
for ii=1:size(VMN)
VMNEig(VMN(ii,1),1)=1;
end

% now make a nice template array 
templateEig=[VisCentEig,SoMatEig,DorsAttnEig,SalVentEig,LimbEig,ContEig,DMNEig,TmpParEig,VMNEig];

%%
% Make template matrices 
VisMat=zeros(num_rois,num_rois);
for ii=1:size(Vis)
    for jj=1:size(Vis)
        VisMat(Vis(ii,1),Vis(jj,1))=1;
    end 
end

SoMatMat=zeros(num_rois,num_rois);
for ii=1:size(SoMat)
    for jj=1:size(SoMat)
        SoMatMat(SoMat(ii,1),SoMat(jj,1))=1;
    end 
end

DorsAttnMat=zeros(num_rois,num_rois);
for ii=1:size(DorsAttn)
    for jj=1:size(DorsAttn)
        DorsAttnMat(DorsAttn(ii,1),DorsAttn(jj,1))=1;
    end 
end

SalVentMat=zeros(num_rois,num_rois);
for ii=1:size(SalVent)
    for jj=1:size(SalVent)
        SalVentMat(SalVent(ii,1),SalVent(jj,1))=1;
    end 
end

LimbMat=zeros(num_rois,num_rois);
for ii=1:size(Limb)
    for jj=1:size(Limb)
        LimbMat(Limb(ii,1),Limb(jj,1))=1;
    end 
end

ContMat=zeros(num_rois,num_rois);
for ii=1:size(Cont)
    for jj=1:size(Cont)
        ContMat(Cont(ii,1),Cont(jj,1))=1;
    end 
end

DMNMat=zeros(num_rois,num_rois);
for ii=1:size(DMN)
    for jj=1:size(DMN)
        DMNMat(DMN(ii,1),DMN(jj,1))=1;
    end 
end

TmpParMat=zeros(num_rois,num_rois);
for ii=1:size(TmpPar)
    for jj=1:size(TmpPar)
        TmpParMat(TmpPar(ii,1),TmpPar(jj,1))=1;
    end 
end

VMNMat=zeros(num_rois,num_rois);
for ii=1:size(VMN)
    for jj=1:size(VMN)
        VMNMat(VMN(ii,1),VMN(jj,1))=1;
    end 
end


templateMat(:,:,1)=VisMat;
templateMat(:,:,2)=SoMatMat;
templateMat(:,:,3)=DorsAttnMat;
templateMat(:,:,4)=SalVentMat;
templateMat(:,:,5)=LimbMat;
templateMat(:,:,6)=ContMat;
templateMat(:,:,7)=DMNMat;
templateMat(:,:,8)=TmpParMat;
templateMat(:,:,9)=VMNMat;

%%
save('templateEigAndMat.mat')

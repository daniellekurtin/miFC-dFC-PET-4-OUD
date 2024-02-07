
%{
MID WHOLE TASK
Higher in HC - .46
	- 0.04 = VMN ROIs, 0.42 = no VMN ROI
Higher in MD - .54
	- 0.30 = VMN ROIs, 0.24 = no VMN ROI

REWARD ANTICIPATION TRIALS
Higher in HC - 0.54
	- 0.05 VMN ROIs, 0.49 no VMN ROI
Higher in MD - .46
	- 0.25 VMN ROIs, 0.21 no VMN

CUE REACTIVITY WHOLE TASK
Higher in HC - 0.77
	- 0.07 VMN ROIs, 0.70 non VMN
Higher in MD - 0.23
	- 0.06 VMN, 0.17 non VMN

CUE BLOCKS
Higher in HC - 0.84
	- 0.08 VMN, 0.76 nonVMN
Higher in MD - 0.16
	- 0.06 VMN, 0.10 non VMN
%}

close all

% MID whole task
C=linspecer(2);
C2=C;
C2(1:2,:)=[C(1,:);C(1,:)];
C2(3:4,:)=[C(2,:);C(2,:)];

% MID whole task
figure()
h=pie([0.04 0.42 0.30 0.24]); %,{'HC>MD, VMN','HC>MD, no VMN','MD>HC, VMN','MD>HC, no VMN'});
patchHand = findobj(h, 'Type', 'Patch'); 
set(patchHand, {'FaceColor'}, mat2cell(C2, ones(size(C2,1),1), 3));
set(patchHand, {'FaceAlpha'},[{0.5};{1};{0.5};{1}]);


% MID reward anticipation
figure()
h=pie([0.05 0.49 0.25 0.21]); %,{'HC>MD, VMN','HC>MD, no VMN','MD>HC, VMN','MD>HC, no VMN'});
patchHand = findobj(h, 'Type', 'Patch'); 
set(patchHand, {'FaceColor'}, mat2cell(C2, ones(size(C2,1),1), 3));
set(patchHand, {'FaceAlpha'},[{0.5};{1};{0.5};{1}]);

% Cue Reactivity task
figure()
h=pie([0.07 0.70 0.06 0.17]); %,{'HC>MD, VMN','HC>MD, no VMN','MD>HC, VMN','MD>HC, no VMN'});
patchHand = findobj(h, 'Type', 'Patch'); 
set(patchHand, {'FaceColor'}, mat2cell(C2, ones(size(C2,1),1), 3));
set(patchHand, {'FaceAlpha'},[{0.5};{1};{0.5};{1}]);

% Cue blocks
figure()
h=pie([0.08 0.76 0.06 0.10]); %,{'HC>MD, VMN','HC>MD, no VMN','MD>HC, VMN','MD>HC, no VMN'});
patchHand = findobj(h, 'Type', 'Patch'); 
set(patchHand, {'FaceColor'}, mat2cell(C2, ones(size(C2,1),1), 3));
set(patchHand, {'FaceAlpha'},[{0.5};{1};{0.5};{1}]);


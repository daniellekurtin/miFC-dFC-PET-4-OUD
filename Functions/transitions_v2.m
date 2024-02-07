
function transition_matrix = transitions(time_series,num_substates)

% Function to compute transition matrix for a time-series

% example 

% 4 sub-states, 22 participants, 3 conditions

%{

num_substates = 4;

for cond = 1:3

for sub = 1:22
    
    time_series = randi(num_substates,1,230);
    transition_matrix.cond{cond}(:,:,sub) = transitions(time_series,num_substates);

end

end

%plotting

%}


%%



%difference = time_series(2:end)-time_series(1:end-1);
clear switches
switches(1,:) = time_series(1:end-1);
switches(2,:) = time_series(2:end);

transition_matrix = zeros(num_substates,num_substates);

for frame = 1:length(switches(1,:))
   
    transition_matrix(switches(1,frame),switches(2,frame)) = transition_matrix(switches(1,frame),switches(2,frame))+1;%Counts transition types
    
end

trans_mat=[];

for ii=1:size(transition_matrix,1)
    for jj=1:size(transition_matrix,2)
        trans_mat(ii,jj)=transition_matrix(ii,jj)/ sum(transition_matrix(ii,:))     ;   
        
    end
end


clear transition_matrix
transition_matrix = trans_mat;
% for substate = 1:num_substates
% 
%    % BELOW IS HENRY'S ORIGINAL SCRIPT
%    % It computes the proportion of directed transitions as a proportion of
%    % time that state occurred - not out of all transitions that state could
%    % make!
%    % transition_matrix(substate,:) = transition_matrix(substate,:)./sum(switches(1,:)==substate);%Express transitions from each state as a proportion of total occurences of that state
%     
% end

end



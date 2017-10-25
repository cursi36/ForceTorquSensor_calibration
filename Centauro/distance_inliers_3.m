%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inliers Computation
%
% - n_in = [2 x n_comp] # inliers on all data and last data set for each
% component;
% - Inliers = Vector of inliers for each component;
% - Pop = [ 1 x n_comp] Population of the inliers;
% - Spread = [ 1 x n_comp] Spread of the inliers.
% - F_sample_ref = [m x 6] Reference wrench;
% - S = Acquired data;
% - r = [1 x 2]  Dimension of latest acquired data set = 15;
% - S_sample = [n_comp x 6] rows of the shape matrix;
% - limits = sensor's range.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n_in,Inliers,Pop,Spread] = distance_inliers_3(F_sample_ref,S_sample,S,limits,res,r)

% Threshold for inliers is: error in V smaller than 2*average of errors

[sr,~] = size(S_sample); % number of components to analyse


V_ref = S(:,2:end).'; % refrence voltage

V = S_sample*F_sample_ref;  % fitted voltage (computed one)

[sS,~] = size(S);
Points = [1:sS]';

n_in = zeros(2,sr);  %[tot inliers; inliers on new set]

if sr > 1 % More than one component analysed
    
    Inliers = cell(1,sr);
    Pop = zeros(1,sr);
    Spread = zeros(1,sr);
    
    pop = zeros(1,6);
    spread = zeros(1,6);
    
    for i = 1:sr
        
        err = (V_ref(i,:)-V(i,:));
        
        
        thresh = 1/sqrt(length(err))*norm(err); % Inliers threshold
        w_lim = 2*thresh;
        
        f = find( abs(err) <= w_lim); % inliers
        n_in(1,i) = length(f);
        Inliers{1,i} = Points(f,1);
        
        if isempty(f) == 0 %There are inliers
            
            for ii = 1:6
                cnts = hist(F_sample_ref(ii,Points(f,1)),-limits(ii):res(ii):limits(ii));
                f_cnts = find(cnts ~= 0);
                pop(ii) = length(f_cnts)/length(cnts);
                
                spread(ii) = (max(f_cnts)-min(f_cnts))/(length(cnts)-1);
            end
            Pop(i) = 1/6*sum(pop);
            Spread(i) = 1/6*sum(spread);
            
        end
        
        % Inliers on new set
        f = find( abs(err(:,r(end-1):r(end))) <= w_lim); % inliers
        n_in(2,i) = length(f);
        
    end
    
else % only one voltage component analysed
    
    pop = zeros(1,6);
    spread = zeros(1,6);
    
    Pop = 0;
    Spread = 0;
    
    
    err = (V_ref-V); % [0 .... 0 err]
    
    thresh = 1/sqrt(length(err))*norm(err);
    w_lim = 2*thresh;
    
    f = find( abs(err) <= w_lim); % inliers
    n_in(1) = length(f);
    Inliers = Points(f,1);
    
    if isempty(f) == 0 %There are inliers
        
        for ii = 1:6
            cnts = hist(F_sample_ref(ii,Points(f,1)),-limits(ii):res(ii):limits(ii));
            f_cnts = find(cnts ~= 0);
            pop(ii) = length(f_cnts)/length(cnts);
            
            spread(ii) = (max(f_cnts)-min(f_cnts))/(length(cnts)-1);
        end
        Pop = 1/6*sum(pop);
        Spread = 1/6*sum(spread);
        
    end
    
    % Inliers on new set
    f = find( err(:,r(end-1):r(end)) <= w_lim); % inliers
    n_in(2) = length(f);
    
    
    
end





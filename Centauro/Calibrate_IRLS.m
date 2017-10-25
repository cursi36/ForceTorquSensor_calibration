%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute shape matrix with Iteratively Constrained Iteratively Reweighted Least
% Squares and inliers evaluation
%
% - Values = [4 x 6] cell array of inliers for each voltage component;
% - S_sample = [6 x6] Shape matrix;
% - N_in = [ 2 x n] Numer of inliers for each component;
% - S = Acquired data;
% - r = [1 x 2]  Dimension of latest acquired data set = 15;
% - F_sample_ref = [m x 6] Reference wrench;
% - options = solver options;
% - A_eq,b_eq,A_ineq,b_ineq,lb,ub = equality, inequality constraints and
% bounds;
% - limits = [6 x 1] Range of each wreanch component;
% res = [6 x 1] Resolution for each force component.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Values,S_sample,N_in] = Calibrate_IRLS(S,T,Values,r,F_sample_ref,S_sample,options,A_eq,b_eq,A_ineq,b_ineq,lb,ub,limits,res)


% Tuning Parameters
w_pop = 0.4;
w_spread = 0.4;
w_in = 0.2;

P_lim = 0.9; % Limit percentage of inliers

S_mat = zeros(6,6,2); % Inizialization of matrix on new data set

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Old model evaluation
% Compute inliers of old model on all dta set and on new set, and compute
% function f_old.
% - N_in = [# inliers on all data; # on new set ] [2 x n]

[N_in,Inliers,Pop,Spread] = distance_inliers_3(F_sample_ref,S_sample,S(:,1:7),limits,res,r);

F_choice = w_pop*Pop+w_spread*Spread+w_in*N_in(1,:)/length(S); % f_old

p = N_in(2,:)/(r(end)-r(end-1)+1); % Percentage of inliers on new set

fp_p = find(p >= P_lim); % components with enough inliers
fp_n = find(p < P_lim ); % components with few inliers

% If there are enough inliers, model for component fp_p(ii) is good ----> refine on total inliers

if isempty(fp_p) == 0
    
    
    for ii = 1:length(fp_p)
        
        if rank(F_sample_ref(:,Inliers{1,fp_p(ii),1})) == 6
            
            [S_sample(fp_p(ii),:)] = optimal_lsq(A_eq,b_eq,A_ineq,b_ineq,S_sample(fp_p(ii),:),...
                fp_p(ii),S(Inliers{1,fp_p(ii),1},fp_p(ii)+1),F_sample_ref(:,Inliers{1,fp_p(ii),1}).');
            
        end
    end
    %
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Not enough inliers for component fp_n on new set ------> Either old model wrong or wrong data
% Compute new models from new data set and from all data set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(fp_n) == 0
    
    % Model on new set
    
    for ii = 1:length(fp_n)
        
        if rank(F_sample_ref(:,r(end-1):r(end))) == 6
            
            [S_mat(fp_n(ii),:,1)] = optimal_lsq(A_eq,b_eq,A_ineq,b_ineq,zeros(1,6),...
                fp_n(ii),S(r(end-1):r(end),fp_n(ii)+1),F_sample_ref(:,r(end-1):r(end)).');
            
        end
        
        % Model on entire data set
        
        if rank(F_sample_ref) == 6
            
            [S_mat(fp_n(ii),:,2)] = optimal_lsq(A_eq,b_eq,A_ineq,b_ineq,zeros(1,6),...
                fp_n(ii),S(:,fp_n(ii)+1),F_sample_ref.');
            
        end
        
        % Compute inliers for both models on all data set
        
        [n_in,inliers{1,1},pop,spread] = distance_inliers_3(F_sample_ref,S_mat(fp_n(ii),:,1),[S(:,1) S(:,fp_n(ii)+1)],limits,res,r);
        
        f_choice = w_pop*pop+w_spread*spread+w_in*n_in(1,:)/length(S);
        
        
        [n_in_tot,inliers{1,2},pop,spread] = distance_inliers_3(F_sample_ref,S_mat(fp_n(ii),:,2),[S(:,1) S(:,fp_n(ii)+1)],limits,res,r);
        
        f_choice(2,:) = w_pop*pop+w_spread*spread+w_in*n_in_tot(1,:)/length(S);
        
        CHOICE = ([f_choice; F_choice(fp_n(ii))]); % 3 x 1
        
        f = find(CHOICE == max(CHOICE));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Choose best model for each component and refine on inliers.
        
        if f(1) == 3  % Old model
            
            if rank([F_sample_ref(:,Inliers{1,(fp_n(ii))}) ]) == 6
                
                [S_sample(fp_n(ii),:)] = optimal_lsq(A_eq,b_eq,A_ineq,b_ineq,S_sample(fp_n(ii),:),...
                    fp_n(ii),S(Inliers{1,(fp_n(ii))},fp_n(ii)+1),F_sample_ref(:,Inliers{1,(fp_n(ii))}).');
                
            end
            
        else
            
            if rank([F_sample_ref(:,inliers{1,f(1)}) ]) == 6
                
                [S_sample(fp_n(ii),:)] = optimal_lsq(A_eq,b_eq,A_ineq,b_ineq,S_mat(fp_n(ii),:,f(1)),...
                    fp_n(ii),S(inliers{1,f(1)},fp_n(ii)+1),F_sample_ref(:,inliers{1,f(1)}).');
                
            else
                S_sample(fp_n(ii),:) = S_mat(fp_n(ii),:,f(1));
                
            end
            
            
            
            
        end
        
        
    end
    
    
end


[N_in,Inliers,~,~] = distance_inliers_3(F_sample_ref,S_sample,S(:,1:7),limits,res,r);


for i = 1:6
    
    Values{1,i} = S(Inliers{1,i},2:7);
    Values{2,i} = S(Inliers{1,i},8:13);
    Values{3,i} = T(i,:)*Values{2,i}.';
    Values{4,i} = S(Inliers{1,i},1);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model computation
%
% - s_i = [6 x1] i-th row of shape matrix;
% - f = component whose model has to be computed;
% - s0 = [1 x 6];
% - V = [m x 1] voltage component;

function [s_i] = optimal_lsq(A_eq,b_eq,A_ineq,b_ineq,s0,f,V,F_ref)

% x = [6x1];
% y = [nx1] = V;
% A = [nx6] = F_ref;
% z = [x;u]; u = Lagrange multipliers;

A = F_ref;
y = V;

constr = any(A_eq(:,:,f) ~= 0, 2);
k = find( constr ~= 0); % number of constraints;
D = A_eq(k,:,f);
h = b_eq(k,1);

% A_big *z = b_big;

A_big = zeros(6+length(k));
A_big(1:6,1:6) = 2*(A)'*A;    A_big(1:6,7:end) = D.';
A_big(7:end,1:6) = D;

b_big = [2*(A)'*y;h];

z = A_big\b_big;
x(:,1) = z(1:6);

i = 0;
err = Inf;
eps = 1e-03;
n_iter = 3000;


while  err > eps  || i > n_iter
    
    
    d =  (y-A*x(:,1)); % [m x 1] vector of residuals
    
    thresh_d = 2*1/sqrt(length(d))*norm(d); % Inliers threshold
    
    w = (abs(d) < thresh_d).*((1-d.^2/(thresh_d)^2)).^2; % [m x 1] vector of weights
    
    A_big(1:6,1:6) = 2*A'*((w.^2).*A );
    b_big(1:6) = 2*A'*((w.^2).*y);
    
    
    z = A_big\b_big;
    
    x(:,2) = z(1:6);
    
    
    err = norm(x(:,1)-x(:,2))/norm(x(:,1));
    
    x(:,1) = x(:,2);
    i = i+1;
    
    
end

s_i = x(:,2);

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute shape matrix with Constrained Least Squares and no inliers
% evaluation. All data set is considered.
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

function [Values,S_sample,N_in] = Calibrate_opt_tot(S,T,Values,r,F_sample_ref,S_sample,options,A_eq,b_eq,A_ineq,b_ineq,lb,ub)
%
for ii = 1:6
    
    [S_sample(ii,:)] = optimal_lsq(A_eq,b_eq,A_ineq,b_ineq,zeros(1,6),...
        ii,S(:,ii+1),F_sample_ref.',options,lb,ub);
    
    
end

N_in = [];

end

function [c_i] = optimal_lsq(A_eq,b_eq,A_ineq,b_ineq,x0,f,V,F_ref,options,lb,ub)


y = V;
x = F_ref;

[s,~,flag] = fmincon(@(s)cost(s,y,x),x0.',[],[],[],[],lb(:,f),ub(:,f),[],options);



if flag == -1 || flag == -2
    
    c_i = s0;
    
else

    c_i = s';
end

end

function err = cost(A,y,x)


d =  (y-[x ]*A);
err = d.'*d;


end






function [n_in,Inliers,Pop,Spread] = distance_inliers_normal(F_sample_ref,S_sample,S,sigmaV,limits,res)

% Distance in V smaller than 3*sigmaV
% sigmaV = n_components x 1
% S_sample = ncomp x 6;
% sigmaV = nx1;
% F_sample_ref = 6 x n

[sr,~] = size(S_sample);

n = [S_sample -ones(sr,1)];  % normals to planes (not unitary)

V_ref = S(:,2:end).'; %6 x sS

[sS,~] = size(S);
Points = [1:sS]';

w_lim = (3*sigmaV);

n_in = zeros(1,sr);
Inliers = cell(1,sr);
Pop = zeros(1,sr);
Spread = zeros(1,sr);

pop = zeros(1,6);
spread = zeros(1,6);


V = S_sample*F_sample_ref;  % sr x sS



for i = 1:sr
    
    err = (V_ref(i,:)-V(i,:)); % [0 .... 0 err]
    d = err;
    %                d = abs(err*n(i,7)/norm(n(i,:),2));
    
    f = find( d <= w_lim(i)); % inliers
    n_in(i) = length(f);
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
    
end

end





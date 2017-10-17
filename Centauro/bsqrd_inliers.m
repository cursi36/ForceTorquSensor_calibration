function [n_in,Inliers,pop,spread] = bsqrd_inliers(F_sample_ref,F_sample,tune,S,w_lim,limits,res)

X = S(:,2:7);
[sS,~] = size(S);
Points = [1:sS]';

resid = (F_sample_ref-F_sample).';

[~,sr] = size(resid);

h = leverage(X);

MAD = mad(resid,1,1);

s = MAD/0.6745;

if sr > 1
    
    s = repmat(s,length(h),1);
    
    h = repmat(h,1,sr);
    
end

r = resid./(tune*s.*sqrt(1-h));

w = (abs(r)<1) .* (1 - r.^2).^2;

n_in = zeros(1,sr);

Inliers = cell(1,sr);

pop = zeros(1,sr);
spread = zeros(1,sr);

for i = 1:sr
    
    f = find(w(:,i) >= w_lim);
    
    n_in(1,i) = length(f);
    
    Inliers{1,i} = Points(f,1);
    
    if isempty(f) == 0 %There are inliers
        cnts = hist(F_sample_ref(i,Points(f,1)),-limits(i):res(i):limits(i));
        f_cnts = find(cnts ~= 0);
        pop(i) = length(f_cnts)/length(cnts);
        %     spread(i) = var(f_cnts);
        spread(i) = (max(f_cnts)-min(f_cnts))/(length(cnts)-1);
        
        
    end
    
end

end

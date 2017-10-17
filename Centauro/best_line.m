function [A] = best_line(F_sample_ref,F_calib)

[r,c] = size(F_calib) ;

A = zeros(2,r); % [ang.coeff ; intercept] for each model/component

for i = 1:r
    
    F = [F_calib(i,:).' ];
    A(1,i) = pinv(F)*F_sample_ref(i,:).';
    
end

end
function [R_2]=R_sqr2(F_ref,F_calib)

[s,~]=size(F_ref);

sigma = s*var(F_ref); % Column vector

sqr_err = sum((F_calib-F_ref.').^2,2);  % Row vectors

R_2=1-(sqr_err./sigma.');


end
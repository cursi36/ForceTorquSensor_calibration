%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mean Squared Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mse] = immse(A,B)

err2 = (A-B).^2;
[r,c] = size(A); 

if r > 1
mse = sqrt(sum(err2,2))/c;

elseif r == 1 %(row vectors)
    
    mse = sqrt(sum(err2))/c;
end
end
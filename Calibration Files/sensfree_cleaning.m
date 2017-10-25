%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get info from unloaded sensor
%
% Min, Max = minimum and maximum values recorded;
% offset = offset of voltages and forces;
% Sigma = standard deviation of voltages and forces.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Min,Max,offset,Sigma]=sensfree_cleaning

sens_free = load('sens_free.txt');

[s,~] = size(sens_free);

offset = mean(sens_free(:,2:13));
offset = repmat(offset,s,1);

sens_free(:,2:13) = sens_free(:,2:13)-offset;

Min = min(sens_free(1:end,2:13));
Max = max(sens_free(1:end,2:13));

Sigma = sqrt(var(sens_free(:,2:13)));

end

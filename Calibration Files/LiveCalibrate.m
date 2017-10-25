%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Live Calibration
%
% Obtain Shape matrix with different methods:
% - with constrained optimization;
% - unconstrained robustfit;
% - constrained optimization without inliers.
%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  LiveCalibrate(f,hs,hati,ax_pop,ax_R2,limits,res,h_in,T,h_mat)

%% Initialization 

options = optimoptions('lsqlin','Display','off','algorithm','active-set');

options_fmincon = optimoptions('fmincon','Display','off','algorithm','sqp');

warning('off','optimlib:lsqlin:WillBeRemoved');

load 'C_ref'

m = max(abs(C_ref),[],2);
C_ref = C_ref./repmat(m,1,6);

S_ref = inv(C_ref);
m = max(abs(S_ref),[],2);
S_ref = S_ref./repmat(m,1,6);

% Inequality Constraints

A_constr = zeros(6,6,6);
b_constr = zeros(6,1);

% Equality Constraints for each voltage component
Aeq_constr = zeros(6,6,6);
beq_constr = 1e-06*zeros(6,1);

Aeq_constr(1,1,1) = 1;  Aeq_constr(2,2,1) = 1;  Aeq_constr(3,6,1) = 1; % V1
Aeq_constr(1,3,2) = 1;  Aeq_constr(2,4,2) = 1;  Aeq_constr(3,5,2) = 1; % V2
Aeq_constr(1,1,3) = 1;  Aeq_constr(2,2,3) = 1;  Aeq_constr(3,4,3) = 1;   Aeq_constr(4,6,3) = 1; % V3
Aeq_constr(1,1,4) = 1;  Aeq_constr(2,3,4) = 1;  Aeq_constr(3,4,4) = 1;   Aeq_constr(4,5,4) = 1; % V4
Aeq_constr(1,1,5) = 1;  Aeq_constr(2,2,5) = 1;  Aeq_constr(3,6,5) = 1; % V5
Aeq_constr(1,3,6) = 1;  Aeq_constr(2,4,6) = 1;  Aeq_constr(3,5,6) = 1; % V6

% Lowe and Upper bounds for each voltage
eps = -1e-06; % Set to Inf if unconstrained solution wanted
Eps = Inf;
%
lb(:,1) = -[eps;eps;Eps;Eps;Eps;eps];
lb(:,2) = -[Eps;Eps;eps;eps;eps;Eps];
lb(:,3) = -[eps;eps;Eps;eps;Eps;eps];
lb(:,4) = -[eps;Eps;eps;eps;eps;Eps];
lb(:,5) = -[eps;eps;Eps;Eps;Eps;eps];
lb(:,6) = -[Eps;Eps;eps;eps;eps;Eps];

ub(:,1) = -lb(:,1);
ub(:,2) = -lb(:,2);
ub(:,3) = -lb(:,3);
ub(:,4) = -lb(:,4);
ub(:,5) = -lb(:,5);
ub(:,6) = -lb(:,6);



i = 1;
r = ones(2,1);
s = zeros(2,1);

dir_sens=['Fx';'Fy';'Fz';'Mx';'My';'Mz'];

% Get noise values and offest from unloaded sensor
[Min,Max,off_free,~] = sensfree_cleaning;


packs = 15; % number of points to analize


Values = cell(4,6);
Values_opt = cell(4,6);
Values_tot = cell(4,6);
Values_2 = cell(4,6);
Values_line = cell(4,6);

C_sample = zeros(6);

S_sample = Inf*ones(6);
S_sample_tot = Inf*ones(6);

C_tot = zeros(6);
C_pinv = zeros(6);

S_sample_opt = S_sample;
S_sample_2 = S_sample;
offset = zeros(6,1);
N_in = zeros(1,6);
N_in_opt = zeros(1,6);
N_in_tot = zeros(1,6);
N_in_2 = zeros(1,6);
N_in_line = zeros(1,6);

j = 1;
R2 = zeros(6,2);
err = 100*ones(6,1);
R2_dot = 100*zeros(6,1);

k_stop = 0;

data_thresh = 0; % Data Threshold; if 0 ----> Take all data

%% Online Process 

while i >= 1
    
        try

    Sensors = load ('sens.txt'); % File being saved from sensor reading
    [s(2),~] = size(Sensors);
    
    Sensors(:,1) = Sensors(:,1)/(1e+09); % Time stamp
    
    Sensors(s(1)+1:s(2),2:13) = Sensors(s(1)+1:s(2),2:13)-repmat(off_free(1,:),s(2)-s(1),1); % Remove offset
    
    f_Thresh = find(any(Sensors(:,2:13) >= repmat(data_thresh*Max,s(2),1) | Sensors(:,2:13) <= repmat(data_thresh*Min,s(2),1),2));
    
    if isempty(f_Thresh) == 0
        
        S = Sensors(f_Thresh,:);
        [r(2),~] = size(S);
        
        F_sample_ref = T*S(:,8:13).'; 
        
        Population(F_sample_ref,limits,res,f,ax_pop); % Plot population of applied wrench

        if r(2)-r(1) >= packs  
            
        % Constrained Robust method 
            [Values,S_sample,N_in] = Calibrate_opt_dist(S,T,Values,r,F_sample_ref,...
                S_sample,options_fmincon,Aeq_constr,beq_constr,A_constr,b_constr,lb,ub,limits,res);

        % Robustfit (Unconstrained)
%              [Values,S_sample,N_in] = Calibrate_rob_fit(S,T,Values,r,F_sample_ref,...
%                 S_sample,options_fmincon,Aeq_constr,beq_constr,A_constr,b_constr,lb,ub,limits,res);

        % Constrained Least Squares for each row
      %     [Values,S_sample,N_in] = Calibrate_opt_tot(S,T,Values,r,F_sample_ref,...
     %          S_sample,options_fmincon,Aeq_constr,beq_constr,A_constr,b_constr,lb,ub);
            
 
            C_sample = inv(S_sample);
            
            F_calib = C_sample*S(:,2:7).';
            
            S_sample_tot = S(:,2:7).'*pinv(F_sample_ref); % OLS on whole Shape matrix
            
            C_tot = inv(S_sample_tot);
            
            if any(C_sample(:) == Inf) == 0
                
                R2(:,2) = R_sqr2(F_sample_ref.',F_calib);
                
            end
            
            j = j+1;
            r(1) = r(2);
            
        end
   
    end
    
    
    f_min = find(R2(:,2) == min(R2(:,2))); % Worst R^2
    lab(2,:) = dir_sens(f_min(1),:);
    
    if R2(f_min(1),2) >= 0.9
        
        k_stop = k_stop+1;
        
        if k_stop == 10

            figure(f(3))
            plot(ax_R2(1),0,0,'.g','MarkerSize',50)

        end
        
    else
        
        k_stop = 0;
        figure(f(3))
        plot(ax_R2(1),0,0,'.r','MarkerSize',50)
        
    end
    
            C_sample_opt_tot = inv(S_sample_opt);
    
     save('C_sample.mat','C_sample')
%             save('C_sample_opt_tot.mat','C_sample_opt_tot')
%     save('C_tot.mat','C_tot')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F_calib = C_sample*Sensors(:,2:7).';
    
    F_sample_ref = T*Sensors(:,8:13).';
    
    plotFigures_shape(f,hs,hati,h_in,F_calib,F_sample_ref,Sensors,ax_R2,min(R2(:,2)),min(R2_dot),max(err),lab,Values,h_mat,S_ref,S_sample)
    
    s(1) = s(2);
    R2(:,1) = R2(:,2);
    pause(0.1);
    
         catch

        end
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot applied wrench's population, mapping the range into interval [0 1]
%
% - F_sample_ref = reference wrench [m x 6];
% - limits = range of each wrench components [6 x 1];
% res = resolution for each wrench component  [6 x 1];
% f = figure handle;

function  Population(F_sample_ref,limits,res,f,ax_pop)

figure(f(4))

for i = 1:6
    
    x = -limits(i):res(i):limits(i);
    
    cnts = hist(F_sample_ref(i,:),x);
    f_cnts = find(cnts ~= 0);
    
    y = (f_cnts-1)/(length(cnts)-1); % mapping to [0 1]
    
    hold on
    plot(i*ones(length(f_cnts),1),y,'.b')
    hold off
    
end


end


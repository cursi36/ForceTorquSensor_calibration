%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Live Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  Simulate(f,hs,hati,ax_pop,ax_R2,limits,res,h_in,T,h_mat)

options = optimoptions('lsqlin','Display','off','algorithm','active-set');

options_fmincon = optimoptions('fmincon','Display','off','algorithm','sqp');

warning('off','optimlib:lsqlin:WillBeRemoved');

load 'C_ref';
% C_sim = load('calib_test6.txt'); % Matrix for simulation


for i = 1:6
    
    pop_lim(i) = length(-limits(i):res(i):limits(i));
    
    
end

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



% 
% sets = ['sens_1.txt';'sens_2.txt';'sens_3.txt';'sens_4.txt';...
%     'sens_5.txt';'sens_6.txt';'sens_7.txt';'sens_8.txt';'sens_9.txt';...
%     'sens10.txt';'sens11.txt';'sens12.txt';'sens13.txt';'sens14.txt'];


sets = ['sens21.txt';'sens22.txt';'sens23.txt';'sens24.txt';'sens25.txt';'sens26.txt';'sens27.txt';'sens28.txt';'sens29.txt';'sens30.txt'];



[s_sets,~] = size(sets);

sS = 0;


for II = 1:s_sets
    
 
    Sens = load(sets(II,:));
    
    [sS,~] = size(Sens);

    Sensors = [];

    Sens(:,1) = [1:sS].';

    i = 1;
    r = ones(2,1);
    s = zeros(2,1);
    
    dir_sens=['Fx';'Fy';'Fz';'Mx';'My';'Mz'];
    
    %% Take only ATI values
    [Min,Max,off_free,Sigma] = sensfree_cleaning;
     
    Sens(:,2:13) = Sens(:,2:13)-repmat(off_free(1,:),sS,1);
 
    
    Values = cell(4,6);
    Values_opt = cell(4,6);
    Values_tot = cell(4,6);
    Values_2 = cell(4,6);
    Values_opt_tot= cell(4,6);
    
    C_sample = zeros(6);
    
    S_sample = Inf*ones(6);
    S_sample_tot = Inf*ones(6);
    S_sample_opt_tot = Inf*ones(6);
    
    
    
    S_sample_opt = S_sample;
    S_sample_2 = S_sample;
    offset = zeros(6,1);
    N_in = zeros(1,6);
    N_in_opt = zeros(1,6);
    N_in_tot = zeros(1,6);
    N_in_2 = zeros(1,6);
    N_in_opt_tot = zeros(1,6);
    
    A_line = zeros(2,6);
    
    j = 1;
    R2 = zeros(6,2);
    err = 100*ones(6,1);
    R2_dot = 100*zeros(6,1);
    
    
        
    %  Add fake values
    
    Sens(:,2:7) = Sens(:,8:13)*T.'*(inv(C_ref)).';
    f_out = find(any(Sens(:,2:13) >= repmat(0*Max,sS,1) | Sens(:,2:13) <= repmat(0*Min,sS,1),2)); % Values above Max
    
%      Sens(f_out(300:349),2:7) = Sens(f_out(300:349),2:7)+10*randn(50,6); Fake values addded in interval
%     

        % Fake values spread
%     rng('default')
%     n_fake = round(length(f_out)*0.2);
%     fake = round(length(f_out)*rand(n_fake,1));
%     fake(fake == 0) = 1;
%     Sens(f_out(fake),2:7) = Sens(f_out(fake),2:7)+100*randn(n_fake,6);
    
    
    
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    packs = 0:15:length(Sens);
    
    time = zeros(5,1);
    
    I = 1;
    
    
    w_lim = Sigma(1:6).';
    
    for i = 2:length(packs)
        
        
        Sensors(packs(i-1)+1:packs(i),:) = Sens(packs(i-1)+1:packs(i),:);
        
        
        [s(2),~] = size(Sensors);

        
        f_Thresh = find(any(Sensors(:,2:13) >= repmat(0*Max,s(2),1) | Sensors(:,2:13) <= repmat(0*Min,s(2),1),2)); % Values above Max
        
        if isempty(f_Thresh) == 0
            
            I = I+1;

            S = Sensors(f_Thresh,:);
            
            r(2) = length(f_Thresh);
            
            F_sample_ref = T*S(:,8:13).';
            
            Population(F_sample_ref,limits,res,f,ax_pop); % Plot population of applied wrench
            
            if r(2)-r(1) >= 15% 

% Constrained Robust method 
         %   [Values,S_sample,N_in] = Calibrate_opt_dist(S,T,Values,r,F_sample_ref,...
          %      S_sample,options_fmincon,Aeq_constr,beq_constr,A_constr,b_constr,lb,ub,limits,res);

        % Robustfit (Unconstrained)
%              [Values,S_sample,N_in] = Calibrate_rob_fit(S,T,Values,r,F_sample_ref,...
%                 S_sample,options_fmincon,Aeq_constr,beq_constr,A_constr,b_constr,lb,ub,limits,res);

        %Constrained Least Squares for each row

                [Values_opt_tot,S_sample_opt_tot,N_in_opt_tot] = Calibrate_opt_tot(S,T,Values_opt_tot,r,F_sample_ref,...
                    S_sample_opt_tot,N_in_opt_tot,options,Aeq_constr,beq_constr,A_constr,b_constr,lb,ub,w_lim,limits,res);
                %%%

                
                C_sample = inv(S_sample_opt);
                
                

                      F_calib = C_sample*S(:,2:7).';
                      
                      S_sample_tot = S(:,2:7).'*pinv(F_sample_ref);
                      
                      C_tot = inv(S_sample_tot);
                      C_pinv = F_sample_ref*pinv(S(:,2:7).');
     
     if any(C_sample(:) == Inf) == 0
    
   R2(:,2) = R_sqr2(F_sample_ref.',F_calib);
   
   end

                      j = j+1;
                      r(1) = r(2);
                 
                   end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     
    
    
     end
     

                      f_min = find(R2(:,2) == min(R2(:,2)));
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
          
                      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F_calib = C_sample*Sensors(:,2:7).';
    
    F_sample_ref = T*Sensors(:,8:13).';
    
    plotFigures_shape(f,hs,hati,h_in,F_calib,F_sample_ref,Sensors,ax_R2,min(R2(:,2)),min(R2_dot),max(err),lab,Values,h_mat,S_ref,S_sample)
     
        s(1) = s(2);
        R2(:,1) = R2(:,2);
        pause(0.05);
        
        
        
    end
    
    C_dist(:,:,II) = C_sample;
    C_opt_tot(:,:,II) = C_sample_tot;
    C_Pinv(:,:,II) = C_pinv;
    S_opt_tot(:,:,II) = S_sample_opt_tot;
    
    Time(:,:,II) = time;
    %      save('C_opt.mat','C_opt')
    %      save('C_opt_tot.mat','C_opt_tot')
    %      save('C.mat','C')
    %      save('C_2.mat','C_2')
    %      save('C_line.mat','C_line')
    %      save('Time.mat','Time')
end

end

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



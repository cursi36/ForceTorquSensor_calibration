%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Live Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  Simulate(f,hs,hati,ax_pop,ax_R2,limits,res,h_in,T,h_mat)

options = optimoptions('lsqlin','Display','off','MaxIter',3000);
load 'C_ref'
% C_sim = load('calib_test6.txt');

options_fmincon = optimoptions('fmincon','Display','off','algorithm','sqp');

for i = 1:6
    
    pop_lim(i) = length(-limits(i):res(i):limits(i));
    
    
end


%load 'C_ref'
%C_sim = load('calib_test6.txt');
%load 'Sigma'

A_constr = zeros(6,6,6);
b_constr = zeros(6,1);

Aeq_constr = zeros(6,6,6);
beq_constr = 1e-03*zeros(6,1);

Aeq_constr(1,1,1) = 1;  Aeq_constr(2,3,1) = 1;  Aeq_constr(3,4,1) = 1;  Aeq_constr(4,5,1) = 1;
% Aeq_constr(6,2,1) = 1; Aeq_constr(6,6,1) = 1;
Aeq_constr(1,1,2) = 1;  Aeq_constr(2,3,2) = 1;  Aeq_constr(3,5,2) = 1;
% Aeq_constr(4,2,2) = 1;  Aeq_constr(4,6,2) = -1;   Aeq_constr(5,2,2) = 2;  Aeq_constr(5,4,2) = 1;
Aeq_constr(1,2,3) = 1;  Aeq_constr(2,4,3) = 1;  Aeq_constr(3,6,3) = 1;
% Aeq_constr(4,1,3) = 1;  Aeq_constr(4,3,3) = -1; Aeq_constr(5,1,3) = 1;  Aeq_constr(5,5,3) = -1;
Aeq_constr(1,2,4) = 1;  Aeq_constr(2,3,4) = 1;  Aeq_constr(3,4,4) = 1;  Aeq_constr(4,6,4) = 1;
% Aeq_constr(5,1,4) = 1;  Aeq_constr(5,5,4) = 1;
Aeq_constr(1,2,5) = 1;  Aeq_constr(2,4,5) = 1;  Aeq_constr(3,6,5) = 1;
% Aeq_constr(4,1,5) = 1;  Aeq_constr(4,5,5) = 1;  Aeq_constr(5,1,5) = 2;  Aeq_constr(5,3,5) = -1;
Aeq_constr(1,1,6) = 1;  Aeq_constr(2,3,6) = 1;  Aeq_constr(3,5,6) = 1;
% Aeq_constr(4,2,6) = 1;  Aeq_constr(4,4,6) = -1;  Aeq_constr(5,2,6) = 1;  Aeq_constr(5,6,6) = -1;


eps = 1e-6;
%
lb(:,1) = -[eps;eps;Inf;Inf;Inf;eps];
lb(:,2) = -[Inf;Inf;eps;eps;eps;Inf];
lb(:,3) = -[eps;eps;Inf;eps;Inf;eps];
lb(:,4) = -[eps;Inf;eps;eps;eps;Inf];
lb(:,5) = -[eps;eps;Inf;Inf;Inf;eps];
lb(:,6) = -[Inf;Inf;eps;eps;eps;Inf];

ub(:,1) = -lb(:,1);
ub(:,2) = -lb(:,2);
ub(:,3) = -lb(:,3);
ub(:,4) = -lb(:,4);
ub(:,5) = -lb(:,5);
ub(:,6) = -lb(:,6);
%
% lb = -Inf*eye(6);
% ub = Inf*eye(6);




% Aeq_constr = zeros(6,6,6);
% beq_constr = zeros(6,1);


% 
% sets = ['sens_1.txt';'sens_2.txt';'sens_3.txt';'sens_4.txt';...
%     'sens_5.txt';'sens_6.txt';'sens_7.txt';'sens_8.txt';'sens_9.txt';...
%     'sens10.txt';'sens11.txt';'sens12.txt';'sens13.txt';'sens14.txt'];


sets = ['sens21.txt';'sens22.txt';'sens23.txt';'sens24.txt';'sens25.txt';'sens26.txt';'sens27.txt';'sens28.txt';'sens29.txt';'sens30.txt'];



[s_sets,~] = size(sets);

sS = 0;


for II = 1:s_sets
    
    %Sens = load('sens_5.txt');
    
    Sens = load(sets(II,:));
    
    [sS,~] = size(Sens);
    
   
    
    Sensors = [];
    
    
    
    Sens(:,1) = [1:sS].';
    
    % Verify that at least 6 values are always added to have pinv working
    i = 1;
    r = ones(2,1);
    s = zeros(2,1);
    
    plane = nchoosek(1:6,2); % components combinations
    dir_sens=['Fx';'Fy';'Fz';'Mx';'My';'Mz'];
    
    %% Take only ATI values
    [Min,Max,off_free,Sigma] = sensfree_cleaning;
     
    Sens(:,2:13) = Sens(:,2:13)-repmat(off_free(1,:),sS,1);
    
    %% Transformation Matrix
    height = 0.16;
    T = eye(6,6);
    T(4,2) = -height;  T(5,1) = height;
    
    windowsize = 5;
    b = 1/windowsize*ones(1,windowsize);
    a = 1;
    
    
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
    
%      Sens(f_out(300:349),2:7) = Sens(f_out(300:349),2:7)+10*randn(50,6);
%     
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
            
            %          S(r(i):r(i)+length(f_Thresh)-1,:) = Sens(f_Thresh,:);
            S = Sensors(f_Thresh,:);
            
            r(2) = length(f_Thresh);
            
            F_sample_ref = T*S(:,8:13).';
            
            %                 [pop_ref] = Population(F_sample_ref,limits,res);
            
            if r(2)-r(1) >= 15% Either take packages of points, or don't care, since also all data will be considered
                
%                 tic
%                 
%                 [Values,S_sample,N_in] = Calibrate_opt_dist(j,S,T,Values,r,F_sample_ref,...
%                     S_sample,N_in,options,Aeq_constr,beq_constr,A_constr,b_constr,w_lim,lb,ub,limits,res,pop_lim);
%                 
%                 time(1) = time(1)+toc;
%                 
%                 tic
% [Values_opt,S_sample_opt,N_in_opt] = Calibrate_opt_dist_3(j,S,T,Values_opt,r,F_sample_ref,...
%                      S_sample_opt,N_in_opt,options_fmincon,Aeq_constr,beq_constr,A_constr,b_constr,w_lim,lb,ub,limits,res,pop_lim);
% time(2) = time(2)+toc;
% 
% tic
% [Values_2,S_sample_2,N_in_2] = Calibrate_rob_fit(j,S,T,Values_2,r,F_sample_ref,...
%                      S_sample_2,N_in_2,options_fmincon,Aeq_constr,beq_constr,A_constr,b_constr,w_lim,lb,ub,limits,res,pop_lim);
% 
% time(3) = time(3)+toc;
                
%                
               tic
                [Values_opt_tot,S_sample_opt_tot,N_in_opt_tot] = Calibrate_opt_tot(S,T,Values_opt_tot,r,F_sample_ref,...
                    S_sample_opt_tot,N_in_opt_tot,options,Aeq_constr,beq_constr,A_constr,b_constr,lb,ub,w_lim,limits,res);
                %%%
                
                time(4) = time(4)+toc;
                
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
                      
                      %% Good Calibration reached
                      
                      figure(f(3))
                      plot(ax_R2(1),0,0,'.g','MarkerSize',50)
                      
%                      save('C_sample.txt','C_sample')
%              save('C_tot.txt','C_tot')
                      
                      
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
                           
                     
        
        
        
        plotFigures_shape(f,hs,hati,h_in,F_calib,F_sample_ref,Sensors,ax_R2,min(R2(:,2)),min(R2_dot),max(err),lab,Values_opt,h_mat,C_ref,C_sample)
  
        s(1) = s(2);
        R2(:,1) = R2(:,2);
        pause(0.005);
        
        
        
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

function [pop] = Population(F_sample_ref,limits,res)

figure(10)

for i = 1:6
    subplot(2,3,i)
    hist(F_sample_ref(i,:),-limits(i):res(i):limits(i));
    
    cnts = hist(F_sample_ref(i,:),-limits(i):res(i):limits(i));
    f_cnts = find(cnts ~= 0);
    pop(i) = length(f_cnts)/length(cnts);
    
end

end

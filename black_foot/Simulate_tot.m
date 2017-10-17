%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Live Calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  Simulate_tot(f,hs,hati,limits,res,h_in)

options = optimoptions('lsqlin','Display','off','MaxIter',3000,'Algorithm','active-set'); 
load 'C_ref'
C_sim = load('calib_test6.txt');

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


eps = 1e-10;
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
%lb = -Inf*eye(6);
%ub = Inf*eye(6);


% Aeq_constr = zeros(6,6,6);
% beq_constr = zeros(6,1);


%
sets = ['sens_0.txt';'sens_1.txt';'sens_2.txt';'sens_3.txt';'sens_4.txt';...
    'sens_5.txt';'sens_6.txt';'sens_7.txt';'sens_8.txt';'sens_9.txt';...
    'sens10.txt';'sens11.txt';'sens12.txt';'sens13.txt';'sens14.txt'];

[s_sets,~] = size(sets);

sS = 0;

for II = 1:3
    
    [ss,~] = size(load(sets(II,:)));
    Sens(sS+1:sS+ss,:) = load(sets(II,:));
    [sS,~] = size(Sens);
end

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
    
    
    %% Transformation Matrix
    height = 0.16;
    T = eye(6,6);
    T(4,2) = -height;  T(5,1) = height;
    
    windowsize = 5;
    b = 1/windowsize*ones(1,windowsize);
    a = 1;
    
    
    Values = cell(4,6);
    Values_opt = cell(4,6);
    Values_opt_tot = cell(4,6);
    Values_2 = cell(4,6);
    Values_line = cell(4,6);
    
    C_sample = zeros(6);
    
        S_sample = 1e6*ones(6);
    
             

    
    C_sample_opt = zeros(6,6);
    C_sample_opt_tot = zeros(6,6);   
    C_sample_2 = zeros(6,6);
    C_sample_line = zeros(6,6);
    offset = zeros(6,1);
    N_in = zeros(1,6);
    N_in_opt = zeros(1,6);
    N_in_opt_tot = zeros(1,6);
    N_in_2 = zeros(1,6);
    N_in_line = zeros(1,6);

    A_line = zeros(2,6);
    
    j = 1;
    R2 = zeros(6,2);
    err = 100*ones(6,1);
    R2_dot = 100*zeros(6,1);
    
 
%  Add fake cvalib matrix
% 
%    Sens(:,2:7) = Sens(:,8:13)*T.'*(inv(C_ref)).';
% Sens(300:399,2:7) = Sens(300:399,8:13)*T.'*inv(C_sim);%F_ati is the simulated force with simulated Calib matrix
% Sens(300:399,8:13) = 100*rand(100,6);%F_ati is the simulated force with simulated Calib matrix


% %     Sens(1:100,2:7) = Sens(300:399,8:13)*T.'*inv(C_sim);%F_ati is the simulated force with simulated Calib matrix

    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
      Sens(:,2:13) = filter(b,a,Sens(:,2:13));
      
    packs = 0:15:length(Sens);
    
    time = zeros(5,1);
    
    I = 1;
    

w_lim = Sigma(1:6).';
    
for i = 2:length(packs)
     
        
        Sensors(packs(i-1)+1:packs(i),:) = Sens(packs(i-1)+1:packs(i),:);
       

        [s(2),~] = size(Sensors);
        
%         Sensors(s(1)+1:s(2),2:13) = Sensors(s(1)+1:s(2),2:13)-repmat(off_free(1,:),s(2)-s(1),1);
        
       


        f_Thresh = find(any(Sensors(:,2:13) > repmat(10*Max,s(2),1) | Sensors(:,2:13) < repmat(10*Min,s(2),1),2)); % Values above Max
        
        if isempty(f_Thresh) == 0
            
            I = I+1;
            
            %          S(r(i):r(i)+length(f_Thresh)-1,:) = Sens(f_Thresh,:);
            S = Sensors(f_Thresh,:);
            
            r(2) = length(f_Thresh);
                        
            F_sample_ref = T*S(:,8:13).';
            
%                 [pop_ref] = Population(F_sample_ref,limits,res);  
            
            if r(2)-r(1) >= 30 % Either take packages of points, or don't care, since also all data will be considered
              
tic

%if j == 1
%
%S_sample = S(:,2:7).'*pinv(F_sample_ref);
%end

%        
%                   [Values,S_sample,N_in] = Calibrate_opt_dist(j,S,T,Values,r,F_sample_ref,...
%                    S_sample,N_in,options,Aeq_constr,beq_constr,A_constr,b_constr,w_lim,lb,ub,limits,res,pop_lim);
                   
%                   [Values,C_sample,N_in] = Calibrate_opt(S,T,Values,r,F_sample_ref,...
%                   C_sample,N_in,options,Aeq_constr,beq_constr,A_constr,b_constr,lb,ub,limits,res);
%                   
%                   [A_line,Values,C_sample,N_in] = Calibrate_best_line(A_line,S,T,Values,r,F_sample_ref,...
%                   C_sample,N_in,options,Aeq_constr,beq_constr,A_constr,b_constr,w_lim,lb,ub,limits,res);
                   
%                   [Values,C_sample,N_in] = Calibrate_opt_dist2(S,T,Values,r,F_sample_ref,...
%                   C_sample,N_in,options,Aeq_constr,beq_constr,A_constr,b_constr,w_lim,lb,ub,limits,res);
                  
                  [Values,S_sample,N_in] = Calibrate_opt_tot(S,T,Values,r,F_sample_ref,...
                  S_sample,N_in,options,Aeq_constr,beq_constr,A_constr,b_constr,lb,ub,w_lim,limits,res);
%%%%                   
                   C_sample = inv(S_sample);
                   
                 time(3) = time(3)+toc;

                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                C_sample_tot = F_sample_ref*pinv(S(:,2:7).');


                                
                F_calib = C_sample*S(:,2:7).';

                R2(:,2) = R_sqr2(F_sample_ref.',F_calib);
                


                 j = j+1;
                 r(1) = r(2);
                 
                
            end
%             
            
            
        end
        
            
            F_calib = C_sample*Sensors(:,2:7).';

%        f_min = find(err == max(err));
%        lab(1,:) = dir_sens(f_min(1),:);
%        f_min = find(R2(:,2) == min(R2(:,2)));
%        lab(2,:) = dir_sens(f_min(1),:);
%        f_min = find(abs(R2_dot) == max(abs(R2_dot)));
%        lab(3,:) = dir_sens(f_min(1),:);

lab = dir_sens(1:3,:);
        
        
        F_sample_ref = T*Sensors(:,8:13).';
        
        
        
        plotFigures(f,hs,hati,h_in,F_calib,F_sample_ref,Sensors,min(R2(:,2)),min(R2_dot),max(err),lab,cell(4,6))

        s(1) = s(2);
        R2(:,1) = R2(:,2);
        pause(0.005);
   
        
        
    end
    C_opt(:,:,II) = C_sample_opt;
    C_opt_tot(:,:,II) = C_sample_opt_tot;
    C(:,:,II) = C_sample;
    C_2(:,:,II) = C_sample_2;
    C_line(:,:,II) = C_sample_line;
    Line(:,:,II) = A_line;
    Time(:,:,II) = time;
%      save('C_opt.mat','C_opt')
%      save('C_opt_tot.mat','C_opt_tot')
%      save('C.mat','C')
%      save('C_2.mat','C_2')
%      save('C_line.mat','C_line')
%      save('Time.mat','Time')
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

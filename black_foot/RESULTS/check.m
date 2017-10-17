


clear all
close all

load 'C_ref'


limits = [300;300;300;25; 25; 25];
  
  res = ones(6,1);


[Min,Max,off_free,Sigma] = sensfree_cleaning;

sets = ['sens21.txt';'sens22.txt';'sens23.txt';'sens24.txt';'sens25.txt';'sens26.txt';'sens27.txt';'sens28.txt';'sens29.txt';'sens30.txt'];
[s_sets,~] = size(sets);

 load ('C_sample_21.mat');
C(:,:,1) = C_sample;

 load ('C_sample_22.mat');
 C(:,:,2) = C_sample;

 load ('C_sample_23.mat');
C(:,:,3) = C_sample;

 load ('C_sample_24.mat');
C(:,:,4) = C_sample;

 load ('C_sample_25.mat');
C(:,:,5) = C_sample;

load ('C_sample_26.mat');
C(:,:,6) = C_sample;

load ('C_sample_27.mat');
C(:,:,7) = C_sample;

load ('C_sample_28.mat');
C(:,:,8) = C_sample;


load ('C_sample_29.mat');
C(:,:,9) = C_sample;

load ('C_sample_210.mat');
C(:,:,10) = C_sample;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 load ('C_sample_rob_fit_21.mat');
C_rob_fit(:,:,1) = C_sample;

 load ('C_sample_rob_fit_22.mat');
 C_rob_fit(:,:,2) = C_sample;

 load ('C_sample_rob_fit_23.mat');
C_rob_fit(:,:,3) = C_sample;

 load ('C_sample_rob_fit_24.mat');
C_rob_fit(:,:,4) = C_sample;

 load ('C_sample_rob_fit_25.mat');
C_rob_fit(:,:,5) = C_sample;

load ('C_sample_rob_fit_26.mat');
C_rob_fit(:,:,6) = C_sample;

load ('C_sample_rob_fit_27.mat');
C_rob_fit(:,:,7) = C_sample;

load ('C_sample_rob_fit_28.mat');
C_rob_fit(:,:,8) = C_sample;


load ('C_sample_rob_fit_29.mat');
C_rob_fit(:,:,9) = C_sample;

load ('C_sample_rob_fit_30.mat');
C_rob_fit(:,:,10) = C_sample;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 load ('C_tot_21.mat');
C_Tot(:,:,1) = C_tot;

 load ('C_tot_22.mat');
C_Tot(:,:,2) = C_tot;

 load ('C_tot_23.mat');
C_Tot(:,:,3) = C_tot;

 load ('C_tot_24.mat');
C_Tot(:,:,4) = C_tot;

 load ('C_tot_25.mat');
C_Tot(:,:,5) = C_tot;

 load ('C_tot_26.mat');
C_Tot(:,:,6) = C_tot;

 load ('C_tot_27.mat');
C_Tot(:,:,7) = C_tot;

 load ('C_tot_28.mat');
C_Tot(:,:,8) = C_tot;

 load ('C_tot_29.mat');
C_Tot(:,:,9) = C_tot;

 load ('C_tot_210.mat');
C_Tot(:,:,10) = C_tot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 load ('C_sample_opt_tot_21.mat');
C_opt_Tot(:,:,1) = C_sample_opt_tot;

 load ('C_sample_opt_tot_22.mat');
 C_opt_Tot(:,:,2) = C_sample_opt_tot;

 load ('C_sample_opt_tot_23.mat');
C_opt_Tot(:,:,3) = C_sample_opt_tot;

 load ('C_sample_opt_tot_24.mat');
C_opt_Tot(:,:,4) = C_sample_opt_tot;

 load ('C_sample_opt_tot_25.mat');
C_opt_Tot(:,:,5) = C_sample_opt_tot;

load ('C_sample_opt_tot_26.mat');
C_opt_Tot(:,:,6) = C_sample_opt_tot;

load ('C_sample_opt_tot_27.mat');
C_opt_Tot(:,:,7) = C_sample_opt_tot;

load ('C_sample_opt_tot_28.mat');
C_opt_Tot(:,:,8) = C_sample_opt_tot;


load ('C_sample_opt_tot_29.mat');
C_opt_Tot(:,:,9) = C_sample_opt_tot;

load ('C_sample_opt_tot_30.mat');
C_opt_Tot(:,:,10) = C_sample_opt_tot;

height=0.16;
T=eye(6,6);
T(4,2)=-height;  T(5,1)=height;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validation set


R = zeros(6,s_sets,s_sets); % R(:,ii,i) on set i by matrix ii

% for i = 1:s_sets
%     
%     S = load(sets(i,:));
%     F_ref = T*S(:,8:13).';
%     C_Pinv(:,:,i) = F_ref*pinv(S(:,2:7).');
% end

L = 22.5/1000;

for i = 1:s_sets
    
    S_mat(:,:,i) = inv(C(:,:,i));
    S_mat(:,4:6,i) = S_mat(:,4:6,i)*L;
    m = max(abs(S_mat(:,:,i)),[],2);
    
    S_scal(:,:,i) = S_mat(:,:,i)./repmat(m,1,6);
    
    figure(1)
    subplot(2,5,i)
    imagesc(abs(S_scal(:,:,i)))
    colormap(jet)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        S_mat_Tot(:,:,i) = inv(C_Tot(:,:,i));
    S_mat_Tot(:,4:6,i) = S_mat_Tot(:,4:6,i)*L;
    
     m = max(abs(S_mat_Tot(:,:,i)),[],2);
    
    S_scal_Tot(:,:,i) = S_mat_Tot(:,:,i)./repmat(m,1,6);
    
    figure(2)
    subplot(2,5,i)
    imagesc(abs(S_scal_Tot(:,:,i)))
    colormap(jet)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
     S_mat_opt_Tot(:,:,i) = inv(C_opt_Tot(:,:,i));
    S_mat_opt_Tot(:,4:6,i) = S_mat_opt_Tot(:,4:6,i)*L;
    
     m = max(abs(S_mat_opt_Tot(:,:,i)),[],2);
    
    S_scal_opt_Tot(:,:,i) = S_mat_opt_Tot(:,:,i)./repmat(m,1,6);
    
    figure(3)
    subplot(2,5,i)
    imagesc(abs(S_scal_opt_Tot(:,:,i)))
    colormap(jet)
    Sensor = load(sets(i,:));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    S_mat_rob_fit(:,:,i) = inv(C_rob_fit(:,:,i));
    S_mat_rob_fit(:,4:6,i) = S_mat_rob_fit(:,4:6,i)*L;
    m = max(abs(S_mat_rob_fit(:,:,i)),[],2);
    
    S_scal_rob_fit(:,:,i) = S_mat_rob_fit(:,:,i)./repmat(m,1,6);
    
    figure(4)
    subplot(2,5,i)
    imagesc(abs(S_scal_rob_fit(:,:,i)))
    colormap(jet)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
       Sensor(:,2:13) = Sensor(:,2:13)-repmat(off_free(1,:),length(Sensor),1);
      
      f_thresh = find(any(Sensor(:,2:13) >= repmat(5*Max,length(Sensor),1) | Sensor(:,2:13) <= repmat(5*Min,length(Sensor),1),2)); % Values above Max
        
      S = Sensor(f_thresh,:); 
    F_ref = T*S(:,8:13).';
    
    F_ref_col = reshape(F_ref',6*length(S),1);
    
    Det(i) = det(F_ref*F_ref.');
    
    
    for iii = 1:6
            cnts = hist(F_ref(iii,:),-limits(iii):res(iii):limits(iii));
            f_cnts = find(cnts ~= 0);
            pop(iii) = length(f_cnts)/length(cnts);
            
            spread(iii) = (max(f_cnts)-min(f_cnts))/(length(cnts)-1);
            end
            Pop(i) = 1/6*sum(pop);
            Spread(i) = 1/6*sum(spread);
            
    for ii = 1:s_sets
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Optimal Solution
        F_calib = C(:,:,ii)*S(:,2:7).';
        R(:,ii,i) = R_sqr2(F_ref.',F_calib);
        
        R_V(:,ii,i) = R_sqr2(S(:,2:7),inv(C(:,:,ii))*F_ref);
         R_adj(:,ii,i) = 1-(1-R(:,ii,i))*(length(F_calib)-1)/(length(F_calib)-5);

        
        F = blkdiag([F_calib(1,:).' ], [F_calib(2,:).' ], [F_calib(3,:).' ] ...
            ,[F_calib(4,:).' ], [F_calib(5,:).' ], [F_calib(6,:).' ]);
        A(:,ii,i) = pinv(F)*F_ref_col;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Optimal Solution
        F_calib = C_rob_fit(:,:,ii)*S(:,2:7).';
        R_rob_fit(:,ii,i) = R_sqr2(F_ref.',F_calib);
        
       
        
        F = blkdiag([F_calib(1,:).' ], [F_calib(2,:).' ], [F_calib(3,:).' ] ...
            ,[F_calib(4,:).' ], [F_calib(5,:).' ], [F_calib(6,:).' ]);
        A_rob_fit(:,ii,i) = pinv(F)*F_ref_col;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %% OPt Tot
        
        F_calib = C_opt_Tot(:,:,ii)*S(:,2:7).';
        R_opt_Tot(:,ii,i) = R_sqr2(F_ref.',F_calib);
        
        R_V_opt_Tot(:,ii,i) = R_sqr2(S(:,2:7),inv(C_opt_Tot(:,:,ii))*F_ref);
         R_adj_opt_Tot(:,ii,i) = 1-(1-R(:,ii,i))*(length(F_calib)-1)/(length(F_calib)-5);

        
        F = blkdiag([F_calib(1,:).' ], [F_calib(2,:).' ], [F_calib(3,:).' ] ...
            ,[F_calib(4,:).' ], [F_calib(5,:).' ], [F_calib(6,:).' ]);
        A_opt_Tot(:,ii,i) = pinv(F)*F_ref_col;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Tot
        F_calib = C_Tot(:,:,ii)*S(:,2:7).';
        R_Tot(:,ii,i) = R_sqr2(F_ref.',F_calib);
         R_V_Tot(:,ii,i) = R_sqr2(S(:,2:7),inv(C_Tot(:,:,ii))*F_ref);
         R_adj_Tot(:,ii,i) = 1-(1-R_Tot(:,ii,i))*(length(F_calib)-1)/(length(F_calib)-5);
         
%          Norm_Tot(:,ii) = norm(C_Tot(:,:,ii),2);
%          S_mat = inv(C_Tot(:,:,ii));
%          for II = 1:6
%         MSE_Tot(II,ii,i) = sqrt(immse(F_calib(II,:),F_ref(II,:)));
%         MSE_V_Tot(II,ii,i) = sqrt(immse(S(:,II+1).',S_mat(II,:)*F_ref));
%          
%          end
           
                F = blkdiag([F_calib(1,:).' ], [F_calib(2,:).' ], [F_calib(3,:).' ] ...
            ,[F_calib(4,:).' ], [F_calib(5,:).' ], [F_calib(6,:).' ]);
        A_Tot(:,ii,i) = pinv(F)*F_ref_col;
        
    end
    
end
S_ref = inv(C_ref);
m = max(abs(S_ref),[],2);
S_ref_scal = S_ref./repmat(m,1,6);
figure()
imagesc(abs(S_ref_scal))
colormap(jet)
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variances of matrices

C_mean = sum(C,3)/s_sets;
C_opt_Tot_mean = sum(C_opt_Tot,3)/s_sets;
C_Tot_mean = sum(C_Tot,3)/s_sets;

% Standard Deviations
var_C = sqrt(sum((C-repmat(C_mean,1,1,s_sets)).^2,3)/s_sets);
var_C_opt_Tot = sqrt(sum((C_opt_Tot-repmat(C_opt_Tot_mean,1,1,s_sets)).^2,3)/s_sets);
var_C_Tot = sqrt(sum((C_Tot-repmat(C_Tot_mean,1,1,s_sets)).^2,3)/s_sets);



  S_mat(:,4:6,:) = S_mat(:,4:6,:)/L;
    S_mat_opt_Tot(:,4:6,:) = S_mat_opt_Tot(:,4:6,:)/L;
     S_mat_Tot(:,4:6,:) = S_mat_Tot(:,4:6,:)/L; 
  
S_mean = sum(S_mat,3)/s_sets;
S_opt_Tot_mean = sum(S_mat_opt_Tot,3)/s_sets;
S_Tot_mean = sum(S_mat_Tot,3)/s_sets;

% Standard Deviations
var_S = sqrt(sum((S_mat-repmat(S_mean,1,1,s_sets)).^2,3)/s_sets);
var_S_opt_Tot = sqrt(sum((S_mat_opt_Tot-repmat(S_opt_Tot_mean,1,1,s_sets)).^2,3)/s_sets);
var_S_Tot = sqrt(sum((S_mat_Tot-repmat(S_Tot_mean,1,1,s_sets)).^2,3)/s_sets);

perc_var_opt = var_S./S_mean*100;
perc_var_opt_Tot = var_S_opt_Tot./S_opt_Tot_mean*100;
perc_var_Tot = var_S_Tot./S_Tot_mean*100;


m = max(abs(C_ref),[],2);
C_ref_scal = C_ref./repmat(m,1,6);

figure()
for i = 1:6
    R_comp = reshape(R(i,:,:),100,1);
    R_Tot_comp = reshape(R_Tot(i,:,:),100,1);
    R_opt_Tot_comp = reshape(R_opt_Tot(i,:,:),100,1);
    subplot(2,3,i)
    hold on
    plot(0:length(R_comp)-1,R_comp,'.b','MarkerSize',15)
     plot(0:length(R_comp)-1,R_Tot_comp,'.r','MarkerSize',15)
      plot(0:length(R_comp)-1,R_opt_Tot_comp,'.k','MarkerSize',15)
     set(gca,'xtick',[0:10:length(R_comp)-1])
     hold off
end

% for i = 1:s_sets
% m = max(abs(C(:,:,i)),[],2);
% C_scal(:,:,i) = C(:,:,i)./repmat(m,1,6);
% Err(1,i) = sqrt(immse(abs(C_ref_scal),abs(C_scal(:,:,i))));
% m = max(abs(C_Pinv(:,:,i)),[],2);        M(:,ii,i) = mean(F_calib.')-mean(F_ref.');
% C_Pinv_scal(:,:,i) = C_Pinv(:,:,i)./repmat(m,1,6);
% Err(2,i) = sqrt(immse(abs(C_ref_scal),abs(C_Pinv_scal(:,:,i))));
% m = max(abs(C_Tot(:,:,i)),[],2);
% C_Tot_scal(:,:,i) = C_Tot(:,:,i)./repmat(m,1,6);
% Err(3,i) = sqrt(immse(abs(C_ref_scal),abs(C_Tot_scal(:,:,i))));
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C_Pinv(:,:,6) =%%%%%%%%%%%%%%
% %% Best Matrix
% R_best = sum(sum(R,3),1);
% R_best = R_best/(6*s_sets);
% 
% R_min = min(R);
% R_min = sum(R_min,3)/s_sets;
% 
% col = find(R_best == max(R_best));
% 
% C_final = C(:,:,col);
% R_final = R(:,col,:);
% 
% S_sample = inv(C_final);
% 
% limits = [300;300;300;25; 25; 25];
% res = ones(6,1);
% 
% 
% sigmaV = Sigma(1:6).';
% 
% for i = 1:s_sets
% 
%     S = load(sets(i,:));
%  
%      S(:,2:13) = S(:,2:13)-repmat(off_free(1,:),length(S),1);
%     F_ref = T*S(:,8:13).';
%     F_calib = C_final*S(:,2:7).';
% % %     [n_in(:,:,i),Inliers,Pop,Spread] = distance_inliers_normal(F_ref,S_sample,S,sigmaV,limits,res);
% % figure()
% % for ii = 1:6
% % subplot(3,2,ii)
% % plot(F_ref(ii,:),F_calib(ii,:),'.')
% % hold on
% % plot(F_ref(ii,:),F_ref(ii,:),'.r')
% % axis equal
% % end
%    
% end
% 
% figure()
% imagesc(abs(C_scal(:,:,col)))
% colormap('jet')
% 
%     Sensor = load ('sens_Fz_2Kg.txt');
% % Sens = load ('sens_free.txt');
%     Sensor(:,2:13) = Sensor(:,2:13)-repmat(off_free(1,:),length(Sensor),1);
%     
%      f_thresh = find(any(Sensor(:,2:13) >= repmat(5*Max,length(Sensor),1) | Sensor(:,2:13) <= repmat(5*Min,length(Sensor),1),2)); % Values above Max
%         
%       Sens = Sensor(f_thresh,:); 
%     
% %     windowsize = 50;
% %     b = 1/windowsize*ones(1,windowsize);
% %     a = 1;
% %     Sens(:,2:13) = filter(b,a,Sens(:,2:13));
%     
% F_ref = T*Sens(:,8:13).';
% F_calib = C_mean*Sens(:,2:7).';
% 
% mean(F_calib.')
% % 
% % R_sqr2(F_ref.',F_calib)
% % R_sqr2(Sens(:,2:7),inv(C_final)*F_ref);
% 
% % sqrt(var(F_calib(:,100:end).'))
% 
% 
% 
% figure()
% for i = 1:6
%     subplot(3,2,i)
%     hold on
%     plot(Sens(:,1),F_ref(i,:),'r')
%     plot(Sens(:,1),F_calib(i,:),'b')
% end
% 
% F_ref = T*Sens(:,8:13).';
% F_calib = C_Tot_mean*Sens(:,2:7).';
% 
% mean(F_calib.')
% 
% % R_sqr2(F_ref.',F_calib)
% % R_sqr2(Sens(:,2:7),inv(C_Tot(:,:,col))*F_ref);
% % sqrt(var(F_calib(:,100:end).'))
% 
% 
% figure()
% for i = 1:6
%     subplot(3,2,i)
%     hold on
%     plot(Sens(:,1),F_ref(i,:),'r')
%     plot(Sens(:,1),F_calib(i,:),'b')
% end
% 
% 
% V = mean(Sens(100:end,2:7));
% 
% 
% for i = 1:s_sets
%     F_2kg(:,i) = C(:,:,i)*V.';
%     F_tot_2kg(:,i) = C_Tot(:,:,i)*V.';
% end
%     
% 
% 

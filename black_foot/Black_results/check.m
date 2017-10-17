


clear all
close all


sets = ['sens_1.txt';'sens_2.txt';'sens_3.txt'];
[s_sets,~] = size(sets);

C(:,:,1)= load ('C_sample_1.txt');

C(:,:,2) = load ('C_sample_2.txt');

C(:,:,3) = load ('C_sample_3.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_Tot(:,:,1) = load ('C_tot_1.txt');

C_Tot(:,:,2) = load ('C_tot_2.txt');

C_Tot(:,:,3) = load ('C_tot_3.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_Pinv(:,:,1) = load ('C_pinv_1.txt');

C_Pinv(:,:,2) = load ('C_pinv_2.txt');

C_Pinv(:,:,3) = load ('C_pinv_3.txt');

height=0.16;
T=eye(6,6);
T(4,2)=-height;  T(5,1)=height;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variances of matrices

C_mean = sum(C,3)/s_sets;
C_Pinv_mean = sum(C_Pinv,3)/s_sets;
C_Tot_mean = sum(C_Tot,3)/s_sets;

% Standard Deviations
var_C = sqrt(sum((C-repmat(C_mean,1,1,s_sets)).^2,3)/s_sets);
var_C_Pinv = sqrt(sum((C_Pinv-repmat(C_Pinv_mean,1,1,s_sets)).^2,3)/s_sets);
var_C_Tot = sqrt(sum((C_Tot-repmat(C_Tot_mean,1,1,s_sets)).^2,3)/s_sets);

perc_var_opt = var_C./C_mean*100;
perc_var_Pinv = var_C_Pinv./C_Pinv_mean*100;
perc_var_Tot = var_C_Tot./C_Tot_mean*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Validation set


R = zeros(6,s_sets,s_sets); % R(:,ii,i) on set i by matrix ii

for i = 1:s_sets
    
    S = load(sets(i,:));
    F_ref = T*S(:,8:13).';
    
    F_ref_col = reshape(F_ref',6*length(S),1);
    
    for ii = 1:s_sets
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Optimal Solution
        F_calib = C(:,:,ii)*S(:,2:7).';
        R(:,ii) = R_sqr2(F_ref.',F_calib);
        
        F = blkdiag([F_calib(1,:).' ], [F_calib(2,:).' ], [F_calib(3,:).' ] ...
            ,[F_calib(4,:).' ], [F_calib(5,:).' ], [F_calib(6,:).' ]);
        A(:,ii,i) = pinv(F)*F_ref_col;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Pinv
        F_calib = C_Pinv(:,:,ii)*S(:,2:7).';
        R_Pinv(:,ii,i) = R_sqr2(F_ref.',F_calib);
        
        F = blkdiag([F_calib(1,:).' ], [F_calib(2,:).' ], [F_calib(3,:).' ] ...
            ,[F_calib(4,:).' ], [F_calib(5,:).' ], [F_calib(6,:).' ]);
        A_Pinv(:,ii,i) = pinv(F)*F_ref_col;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Tot
        F_calib = C_Tot(:,:,ii)*S(:,2:7).';
        R_Tot(:,ii,i) = R_sqr2(F_ref.',F_calib);
        
        F = blkdiag([F_calib(1,:).' ], [F_calib(2,:).' ], [F_calib(3,:).' ] ...
            ,[F_calib(4,:).' ], [F_calib(5,:).' ], [F_calib(6,:).' ]);
        A_Tot(:,ii,i) = pinv(F)*F_ref_col;
        
    end
    
end


% Check Results
function Check

clear all
close all

sets = ['sens_0.txt';'sens_1.txt';'sens_2.txt';'sens_3.txt';'sens_4.txt';...
    'sens_5.txt';'sens_6.txt';'sens_7.txt';'sens_8.txt';'sens_9.txt';...
    'sens10.txt';'sens11.txt';'sens12.txt';'sens13.txt';'sens14.txt'];

[s_sets,~] = size(sets);
ss = zeros(1,s_sets);

    height = 0.16;
    T = eye(6,6);
    T(4,2) = -height;  T(5,1) = height;
    
windowsize = 5;
b = 1/windowsize*ones(1,windowsize);
a = 1;

load 'C_Pinv'
load 'C_ref'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unconstr
load 'C_dist_unconstr'
load 'C_opt_tot_unconstr'

[~,~,sc] = size(C_dist);

for i = 1:s_sets
    
    
    Sens = load(sets(i,:));
    Sens(:,2:13) = filter(b,a,Sens(:,2:13));
    
    F_sample_ref = T*Sens(:,8:13).';
    F_ref = reshape(F_sample_ref,6*length(Sens),1);
    
    for ii = 1:sc
        
        
        F_calib = C_dist(:,:,ii)*Sens(:,2:7).';
        R2_dist(:,i,ii) = R_sqr2(F_sample_ref.',F_calib);
        
        F = blkdiag([F_calib(1,:).' ones(length(F_calib(1,:)),1)], [F_calib(2,:).' ones(length(F_calib(1,:)),1)], [F_calib(3,:).' ones(length(F_calib(1,:)),1)] ...
            ,[F_calib(4,:).' ones(length(F_calib(1,:)),1)], [F_calib(5,:).' ones(length(F_calib(1,:)),1)], [F_calib(6,:).' ones(length(F_calib(1,:)),1)]);
        A_dist(:,i,ii) = pinv(F)*F_ref;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F_calib = C_opt_tot(:,:,ii)*Sens(:,2:7).';
        R2_tot(:,i,ii) = R_sqr2(F_sample_ref.',F_calib);
        
        F = blkdiag([F_calib(1,:).' ones(length(F_calib(1,:)),1)], [F_calib(2,:).' ones(length(F_calib(1,:)),1)], [F_calib(3,:).' ones(length(F_calib(1,:)),1)] ...
            ,[F_calib(4,:).' ones(length(F_calib(1,:)),1)], [F_calib(5,:).' ones(length(F_calib(1,:)),1)], [F_calib(6,:).' ones(length(F_calib(1,:)),1)]);
        A_tot(:,i,ii) = pinv(F)*F_ref;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F_calib = C_Pinv(:,:,ii)*Sens(:,2:7).';
        R2_pinv(:,i,ii) = R_sqr2(F_sample_ref.',F_calib);
        
        F = blkdiag([F_calib(1,:).' ones(length(F_calib(1,:)),1)], [F_calib(2,:).' ones(length(F_calib(1,:)),1)], [F_calib(3,:).' ones(length(F_calib(1,:)),1)] ...
            ,[F_calib(4,:).' ones(length(F_calib(1,:)),1)], [F_calib(5,:).' ones(length(F_calib(1,:)),1)], [F_calib(6,:).' ones(length(F_calib(1,:)),1)]);
        A_pinv(:,i,ii) = pinv(F)*F_ref;
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constr

load 'C_dist_constr'
load 'C_opt_tot_constr'

height = 0.16;
T = eye(6,6);
T(4,2) = -height;  T(5,1) = height;

sS = 0;



% R2 on set i from matrix ii

for i = 1:s_sets
    
    
    Sens = load(sets(i,:));
    Sens(:,2:13) = filter(b,a,Sens(:,2:13));
    
    F_sample_ref = T*Sens(:,8:13).';
    F_ref = reshape(F_sample_ref,6*length(Sens),1);
    
    for ii = 1:sc
        
        
        F_calib = C_dist(:,:,ii)*Sens(:,2:7).';
        R2_dist_cnstr(:,i,ii) = R_sqr2(F_sample_ref.',F_calib);
        
        F = blkdiag([F_calib(1,:).' ones(length(F_calib(1,:)),1)], [F_calib(2,:).' ones(length(F_calib(1,:)),1)], [F_calib(3,:).' ones(length(F_calib(1,:)),1)] ...
            ,[F_calib(4,:).' ones(length(F_calib(1,:)),1)], [F_calib(5,:).' ones(length(F_calib(1,:)),1)], [F_calib(6,:).' ones(length(F_calib(1,:)),1)]);
        A_dist_cnstr(:,i,ii) = pinv(F)*F_ref;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F_calib = C_opt_tot(:,:,ii)*Sens(:,2:7).';
        R2_tot_cnstr(:,i,ii) = R_sqr2(F_sample_ref.',F_calib);
        
        F = blkdiag([F_calib(1,:).' ones(length(F_calib(1,:)),1)], [F_calib(2,:).' ones(length(F_calib(1,:)),1)], [F_calib(3,:).' ones(length(F_calib(1,:)),1)] ...
            ,[F_calib(4,:).' ones(length(F_calib(1,:)),1)], [F_calib(5,:).' ones(length(F_calib(1,:)),1)], [F_calib(6,:).' ones(length(F_calib(1,:)),1)]);
        A_tot_cnstr(:,i,ii) = pinv(F)*F_ref;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F_calib = C_Pinv(:,:,ii)*Sens(:,2:7).';
        R2_pinv_cnstr(:,i,ii) = R_sqr2(F_sample_ref.',F_calib);
        
        F = blkdiag([F_calib(1,:).' ones(length(F_calib(1,:)),1)], [F_calib(2,:).' ones(length(F_calib(1,:)),1)], [F_calib(3,:).' ones(length(F_calib(1,:)),1)] ...
            ,[F_calib(4,:).' ones(length(F_calib(1,:)),1)], [F_calib(5,:).' ones(length(F_calib(1,:)),1)], [F_calib(6,:).' ones(length(F_calib(1,:)),1)]);
        A_pinv_cnstr(:,i,ii) = pinv(F)*F_ref;
        
    end
end

end

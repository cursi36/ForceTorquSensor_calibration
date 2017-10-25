function plot_figure
% 
% %while( ii== 0 )
% 
% 
    Sensors{:,1} = load ('sens2.txt');
Sensors{:,2} = load ('sens3.txt');
Sensors{:,3} = load ('sens4.txt');
Sensors{:,4} = load ('sens5.txt');

S(:,:,1) = load('Mat1.txt'); 
S(:,:,2) = load('Mat2.txt');
S(:,:,3) = load('Mat3.txt');
S(:,:,4) = load('Mat4.txt');
S(:,:,5) = load('Mat5.txt');


S_mine = zeros(6,6,5);

S_mine(:,:,1) = load('S_mine_1.txt'); 
S_mine(:,:,2) = load('S_mine_2.txt');
S_mine(:,:,3) = load('S_mine_3.txt');
S_mine(:,:,4) = load('S_mine_4.txt');


height = 0.05;
T = eye(6,6);
T(4,2) = -height;  T(5,1) = height;
S_pinv = zeros(6,6,5);




for ii = 1:4
Sens = Sensors{:,ii};
    F_ref = T*Sens(:,8:13).';
    
    S_pinv(:,:,ii) = Sens(:,2:7).'*pinv(F_ref);
    
end
for ii = 1:4
Sens = Sensors{:,ii};
    F_ref = T*Sens(:,8:13).';
    
for i = 1:5
    C(:,:,i) = inv(S(:,:,i));
    F_calib = C(:,:,i)*Sens(:,2:7).';
    [R2(:,i,ii)] = R_sqr2(F_ref.',F_calib);
    
    F_calib = inv(S_pinv(:,:,i))*Sens(:,2:7).';
    [R2_pinv(:,i,ii)] = R_sqr2(F_ref.',F_calib);
    
    F_calib = inv(S_mine(:,:,i))*Sens(:,2:7).';
    [R2_mine(:,i,ii)] = R_sqr2(F_ref.',F_calib);
    
end




end




% L = 0.02;
% 
% S(:,4:6,1) = S(:,4:6,1);
% 
% m = max(abs(S(:,:,1)),[],2);
% 
% Gamma = S(:,:,1)./m
% 
% Gamma_ref = zeros(6);
% Gamma_ref(1,3:
 
figure(1)
for i = 1:6
subplot(3,2,i)
hold on
plot(Sens(:,1),F_calib(i,:),'b');
plot(Sens(:,1),F_ref(i,:),'r');
hold off
end

% figure(2)
% for i = 1:6
% subplot(3,2,i)
% plot(Sensors(:,1),Sensors(:,7+i),'r');
% end
% pause(0.01);

%end



end

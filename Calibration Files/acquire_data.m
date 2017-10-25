%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Figures and start either LiveCalibration or Simulation
%
% - Figure(1) = Real Time Forces (red = reference, blue = calibrated);
% - Figure(2) = Real Time Moments (red = reference, blue = calibrated);
% - Figure(3) = minimum R2 over time and Stop Indicator;
% - Figure(4) = Forces Population;
% - Figure(5) = Matrices' Fingerprints;

% LiveCalibrate(f,hs,hati,ax_pop,limits(:,2),res,h_in,T);  Simulate(f,hs,hati,ax_pop,limits(:,2),res,h_in,T)
%
% - f = figures' handels;
% - hs = sensor components' handles;
% - hati = ATI components' handles;
% - ax_pop = populations' handles;
% - limits = components' limits;
% - res = resolution;
% - h_in = Inliers' handels;
% - T = transformation matrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function acquire_data

close all


R = input('Do you wish to LiveCalibrate [C] or Simulate [S]? \n','s');

lim = 1;

dir_sens=['Fx';'Fy';'Fz';'Mx';'My';'Mz'];

% Transformation Matrix (F_sample = T*F_ATI);

height = 0.05;
T = eye(6,6);
T(4,2) = -height;  T(5,1) = height;

% Forces and Moments Limits and Resolution (used for population)

limits(:,2) = [300;300;300;25; 25; 25];
limits(:,1) = -limits(:,2);

res = ones(6,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Real time forces/torques
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f(1) = figure(1); % Forces
set(f(1),'position',[0 0 500 450])
ax = zeros(1,6);
hs = zeros(1,6);
hati = zeros(1,6);

for j=1:3
    ax(j) = subplot(3,1,j);
    hold on
    hs(j) = plot(ax(j),0,0,'b');
    hati(j) = plot(ax(j),0,0,'r');
    h_in(j) = plot(ax(j),0,0,'.y');  %Inliers Label
    hold off
    ylabel(dir_sens(j,:));
    xlabel('t(s)');
    ylim(limits(j,:));
end

f(2) = figure(2); % Moments
set(f(2),'position',[600 0 500 450])

for j = 1:3
    ax(j+3) = subplot(3,1,j);
    hold on
    hs(j+3) = plot(ax(j+3),0,0,'b');
    hati(j+3) = plot(ax(j+3),0,0,'r');
    h_in(j+3) = plot(ax(j+3),0,0,'.y');
    hold off
    ylabel(dir_sens(j+3,:));
    xlabel('t(s)');
    ylim(limits(j+3,:));
end

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %% R2
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f(3) = figure(3);
set(f(3),'position',[0 550 500 450])

p(2,:) = [0.13 0.11 0.5 0.775];
p(1,:) = [0.5+0.15 (0.11+0.775)/2 0.2 0.2];  % subplot position

ax_R2(1) = subplot(2,1,1);
axis equal

ax_R2(2) = subplot(2,1,2);

set(ax_R2(1),'position',p(1,:));
set(ax_R2(2),'position',p(2,:));


%% Stop Indicator

plot(ax_R2(1),0,0,'.r','MarkerSize',50)


hold on
plot(ax_R2(2),0,0,'.')
hold off
title('R2')
xlabel('t(s)');
ylim([0 1.2])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax_pop = zeros(1,6);

f(4) = figure(4);
set(f(4),'position',[600 550 500 450])
xlim ([0 7 ]);
ylim ([-0.1 1.1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f(5) = figure(5);
set(f(5),'position',[1200 0 500 450])
subplot(2,1,1);

colormap(jet)
h_mat(1) = imagesc(0);
colorbar;
axis equal;
title('C_sample');


subplot(2,1,2);
colormap(jet)
h_mat(2) = imagesc(0);
colorbar;
axis equal;
title('C_ref');


if R == 'C'
    
    LiveCalibrate(f,hs,hati,ax_pop,ax_R2,limits(:,2),res,h_in,T,h_mat)
    
    
elseif R == 'S'
    
    Simulate(f,hs,hati,ax_pop,ax_R2,limits(:,2),res,h_in,T,h_mat)
    
    
end
end
  

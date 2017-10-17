%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  plotFigures_shape(f,hs,hati,h_in,F_calib,F_ref,Sensors,ax_R2,R2,R2_dot,mse,lab,Inliers,h_mat,S_ref,S_sample)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Live forces and moments in time

figure(f(1))
for j=1:3
    
    hold on
    %sensor
    set(hs(1,j),'xdata',Sensors(:,1),'ydata',F_calib(j,:));
    %ati
    set(hati(1,j),'xdata',Sensors(:,1),'ydata',F_ref(j,:));
    
    set(h_in(1,j),'xdata',Inliers{4,j},'ydata',Inliers{3,j});
    
    
    hold off
end
figure(f(2))
for jj=1:3
    
    hold on
    %sensor
    set(hs(1,j+jj),'xdata',Sensors(:,1),'ydata',F_calib(j+jj,:));
    %ati
    set(hati(1,j+jj),'xdata',Sensors(:,1),'ydata',F_ref(j+jj,:));
    
    
    set(h_in(1,j+jj),'xdata',Inliers{4,j+jj},'ydata',Inliers{3,j+jj});
    
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% R2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(f(3))

hold on
plot(ax_R2(2),Sensors(end,1),R2,'.')
ylabel(lab(2,:));
xlabel('t(s)');
ylim([0 1.2])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrices' Fingerprints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (S_sample == zeros(6))
    
else
    L = 20/1000;
    
    S_sample(:,4:6) = S_sample(:,4:6)*L;
    
    m = max(abs(S_sample),[],2);
    S_sample = S_sample./repmat(m,1,6);
    
    figure(f(5))
    subplot(2,1,1);
    imagesc(abs(S_sample));
    colorbar;
    
    subplot(2,1,2);
    imagesc(abs(S_ref));
    colorbar;
end

end
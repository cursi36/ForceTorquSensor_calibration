  
  
  function  plotFigures(f,hs,hati,h_in,F_calib,F_ref,Sensors,ax_R2,R2,R2_dot,mse,lab,Inliers,h_mat,C_ref,C_sample)
  
  [s,~]=size(Sensors);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Live Plot in time
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      
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
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%% R2_dot
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure(f(4))
%       
%        hold on
%         plot(Sensors(end,1),R2_dot,'.')
%          ylabel(lab(2,:));
%         hold off
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Live C matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (C_sample == zeros(6)) 

else

m = max(abs(C_sample),[],2);
C_sample = C_sample./repmat(m,1,6);

figure(f(5))
subplot(2,1,1);
imagesc(abs(C_sample));
colorbar;

subplot(2,1,2);
imagesc(abs(C_ref));
colorbar;
end


  
  end
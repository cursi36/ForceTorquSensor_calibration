%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize data from Sensor
 function sensor_data_plot
 clear all
 close all
 
 i = 0;
 dir = [ 'Fx'; 'Fy'; 'Fz'; 'Mx'; 'My'; 'Mz'];
 
 while i >= 0
     
  try
     Sensors = load ('sens.txt');
 
     figure(1)
     for ii = 1:6
     subplot(3,2,ii)
    hold on
     plot(Sensors(:,1)/1e+09,Sensors(:,ii+1),'r');
     plot(Sensors(:,1)/1e+09,Sensors(:,10),'b');
     ylabel(dir(ii,:));
     hold off
     end
     
     pause(0.005);
     i = i+1;
     
  catch 
%      
  end

 end
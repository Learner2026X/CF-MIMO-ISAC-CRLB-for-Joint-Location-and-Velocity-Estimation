%% plot histgram for the MLE
close all
SNRdB = 0;%dB
NSR = 1/10^(SNRdB/10);
p_tar = [4000,5000].';
v_tar = [20,30].';
x_trust = 3.2*sqrt(27.6724*NSR);
y_trust = 3.2*sqrt(25.5569*NSR);
vx_trust = 3.2*sqrt(0.0045*NSR);
vy_trust = 3.2*sqrt(0.0065*NSR);
Num_search = 20;
grid_search = zeros(2*Num_search+1,4);
grid_search(:,1) = -x_trust+p_tar(1):x_trust/Num_search:x_trust+p_tar(1);
grid_search(:,2) = -y_trust+p_tar(2):y_trust/Num_search:y_trust+p_tar(2);
grid_search(:,3) = -vx_trust+v_tar(1):vx_trust/Num_search:vx_trust+v_tar(1);
grid_search(:,4) = -vy_trust+v_tar(2):vy_trust/Num_search:vy_trust+v_tar(2);
load('data\Value_estimate0dB_1kHz.mat','Value_estimate')
figure
for i=1:4
    
    subplot(2,2,i)
    histogram(Value_estimate(:,i),grid_search(:,i))
%     xlabel(label{i})
    ax=gca;
    set(ax,'FontSize',10);
end
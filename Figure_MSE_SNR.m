close all
SNR = -30:5:5;%dB
len_SNR = length(SNR);
MSE_location = zeros(len_SNR,2);
MSE_velocity = zeros(len_SNR,2);
MSE_location_MT = zeros(len_SNR,2);
MSE_velocity_MT = zeros(len_SNR,2);
MSE_location_DF = zeros(len_SNR,2);
MSE_velocity_DF = zeros(len_SNR,2);
for i_SNR=1:len_SNR
    if SNR(i_SNR)<0
        load(['data/MSE_', num2str(-SNR(i_SNR)),'dB_1kHz.mat'],'MSE');
        MSE_location(i_SNR,:) = MSE(1:2);
        MSE_velocity(i_SNR,:) = MSE(3:4);
        load(['data/MSE_', num2str(-SNR(i_SNR)),'dB_1kHz_MT.mat'],'MSE');
        MSE_location_MT(i_SNR,:) = MSE(1:2);
        MSE_velocity_MT(i_SNR,:) = MSE(3:4);
        load(['data/MSE_', num2str(-SNR(i_SNR)),'dB_1kHz_DFMT.mat'],'MSE');
        MSE_location_DF(i_SNR,:) = MSE(1:2);
        MSE_velocity_DF(i_SNR,:) = MSE(3:4);
    else
        load(['data/MSE', num2str(SNR(i_SNR)),'dB_1kHz.mat'],'MSE');
        MSE_location(i_SNR,:) = MSE(1:2);
        MSE_velocity(i_SNR,:) = MSE(3:4);
        load(['data/MSE', num2str(SNR(i_SNR)),'dB_1kHz_MT.mat'],'MSE');
        MSE_location_MT(i_SNR,:) = MSE(1:2);
        MSE_velocity_MT(i_SNR,:) = MSE(3:4);
        load(['data/MSE', num2str(SNR(i_SNR)),'dB_1kHz_DFMT.mat'],'MSE');
        MSE_location_DF(i_SNR,:) = MSE(1:2);
        MSE_velocity_DF(i_SNR,:) = MSE(3:4);
    end
end
load('data/CRB_loc_accurate.mat','CRB_loc_accurate')
load('data/CRB_loc.mat','CRB_loc')
load('data/CRB_vel_accurate.mat','CRB_vel_accurate')
load('data/CRB_vel.mat','CRB_vel')
%%
figure
subplot(1,2,1)
semilogy(SNR,MSE_location(:,1),'o-',SNR,MSE_location(:,2),'>-',...
    SNR,MSE_location_MT(:,1),'ch-',SNR,MSE_location_MT(:,2),'m^-',...
    SNR,MSE_location_DF(:,1),'kd-',SNR,MSE_location_DF(:,2),'gv-',...
    SNR,CRB_loc_accurate(:,1),'bs-',SNR,CRB_loc_accurate(:,2),'r<-',...
    'MarkerSize',10,'LineWidth',1.2)
xlabel('SNR (dB)')
ylabel('MSE (m^2)')
% xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('x, single','y, single','x, Delay','y, Delay','x, Dopp.',...
    'y, Dopp.','CRLB, x','CRLB, y')
grid on
subplot(1,2,2)
semilogy(SNR,MSE_velocity(:,1),'o-',SNR,MSE_velocity(:,2),'>-',...
    SNR,MSE_velocity_MT(:,1),'ch-',SNR,MSE_velocity_MT(:,2),'m^-',...
    SNR,MSE_velocity_DF(:,1),'kd-',SNR,MSE_velocity_DF(:,2),'gv-',...
    SNR,CRB_vel_accurate(:,1),'bs-',SNR,CRB_vel_accurate(:,2),'r<-',...
    'MarkerSize',10,'LineWidth',1.2)
xlabel('SNR (dB)')
ylabel('MSE (m^2/s^2)')
% xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('v_x, single','v_y, single','v_{x}, Delay','v_{y}, Delay',...
    'v_{x}, Dopp.','v_{y}, Dopp.','CRLB v_x','CRLB v_y')
grid on
%%
figure
subplot(1,2,1)
semilogy(SNR,CRB_loc_accurate(:,1),'bs-',SNR,CRB_loc_accurate(:,2),'r<-',...
    SNR,CRB_loc(:,1),'cd-',SNR,CRB_loc(:,2),'mo-',...
    'MarkerSize',10,'LineWidth',1.2)
xlabel('SNR (dB)')
ylabel('CRLB (m^2)')
% xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('CRLB x','CRLB y','App. CRLB x','App. CRLB y')
grid on
subplot(1,2,2)
semilogy(SNR,CRB_vel_accurate(:,1),'bs-',SNR,CRB_vel_accurate(:,2),'r<-',...
    SNR,CRB_vel(:,1),'cd-',SNR,CRB_vel(:,2),'mo-',...
    'MarkerSize',10,'LineWidth',1.2)
xlabel('SNR (dB)')
ylabel('CRLB (m^2/s^2)')
% xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('CRLB v_x','CRLB v_y','App. CRLB v_x','App. CRLB v_y')
grid on

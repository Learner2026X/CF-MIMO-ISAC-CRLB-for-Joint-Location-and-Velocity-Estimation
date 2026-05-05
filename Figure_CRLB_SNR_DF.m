close all
SNR = -20:5:20;%dB
load('data/CRB_loc33_accurate.mat','CRB_loc_accurate')
CRB_loc33_accurate = CRB_loc_accurate;
load('data/CRB_loc77_accurate.mat','CRB_loc_accurate')
CRB_loc77_accurate = CRB_loc_accurate;
load('data/CRB_loc33_accurate_DF.mat','CRB_loc_accurate')
CRB_loc33_accurate_DF = CRB_loc_accurate;
load('data/CRB_loc77_accurate_DF.mat','CRB_loc_accurate')
CRB_loc77_accurate_DF = CRB_loc_accurate;
load('data/CRB_vel33_accurate.mat','CRB_vel_accurate')
CRB_vel33_accurate = CRB_vel_accurate;
load('data/CRB_vel77_accurate.mat','CRB_vel_accurate')
CRB_vel77_accurate = CRB_vel_accurate;
load('data/CRB_vel33_accurate_DF.mat','CRB_vel_accurate')
CRB_vel33_accurate_DF = CRB_vel_accurate;
load('data/CRB_vel77_accurate_DF.mat','CRB_vel_accurate')
CRB_vel77_accurate_DF = CRB_vel_accurate;
%%
figure
semilogy( SNR,CRB_loc33_accurate_DF(:,1),'b>--',SNR,CRB_loc33_accurate_DF(:,2),'r<--',...
    SNR,CRB_loc77_accurate_DF(:,1),'c^--',SNR,CRB_loc77_accurate_DF(:,2),'mv--',...
    SNR,CRB_loc33_accurate(:,1),'bp-',SNR,CRB_loc33_accurate(:,2),'rh-',...
    SNR,CRB_loc77_accurate(:,1),'co-',SNR,CRB_loc77_accurate(:,2),'ms-',...
    'MarkerSize',10,'LineWidth',1.2)
xlabel('SNR (dB)')
ylabel('CRLB (m^2)')
% xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('x, 3/times 3, DF','y, 3/times 3, DF','x, 7/times 7, DF','y, 7/times 7, DF',...
    'x, 3/times 3','y, 3/times 3','x, 7/times 7','y, 7/times 7')
grid on
%%
figure
semilogy( SNR,CRB_vel33_accurate_DF(:,1),'b>--',SNR,CRB_vel33_accurate_DF(:,2),'r<--',...
    SNR,CRB_vel77_accurate_DF(:,1),'c^--',SNR,CRB_vel77_accurate_DF(:,2),'mv--',...
    SNR,CRB_vel33_accurate(:,1),'bp-',SNR,CRB_vel33_accurate(:,2),'rh-',...
    SNR,CRB_vel77_accurate(:,1),'co-',SNR,CRB_vel77_accurate(:,2),'ms-',...
    'MarkerSize',10,'LineWidth',1.2)
xlabel('SNR (dB)')
ylabel('CRLB (m^2)')
% xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('v_x, 3/times 3, DF','v_y, 3/times 3, DF','v_x, 7/times 7, DF','v_y, 7/times 7, DF',...
    'v_x, 3/times 3','v_y, 3/times 3','v_x, 7/times 7','v_y, 7/times 7')
grid on
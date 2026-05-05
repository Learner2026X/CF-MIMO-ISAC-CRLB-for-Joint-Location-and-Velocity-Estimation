close all
SNRdB = -20:5:20;
load('data\CRB_loc77_accurate_T2.mat','CRB_loc_accurate')
CRB_loc77_accurate_T2 = CRB_loc_accurate;
load('data\CRB_vel77_accurate_T2.mat','CRB_vel_accurate')
CRB_vel77_accurate_T2 = CRB_vel_accurate;
load('data\CRB_loc77_T2.mat','CRB_loc')
CRB_loc77_T2 = CRB_loc;
load('data\CRB_vel77_T2.mat','CRB_vel')
CRB_vel77_T2 = CRB_vel;
load('data\CRB_loc77_accurate_T3.mat','CRB_loc_accurate')
CRB_loc77_accurate_T3 = CRB_loc_accurate;
load('data\CRB_vel77_accurate_T3.mat','CRB_vel_accurate')
CRB_vel77_accurate_T3 = CRB_vel_accurate;
load('data\CRB_loc77_T3.mat','CRB_loc')
CRB_loc77_T3 = CRB_loc;
load('data\CRB_vel77_T3.mat','CRB_vel')
CRB_vel77_T3 = CRB_vel;
%%
figure
% subplot(2,1,1)
semilogy( SNRdB,CRB_loc77_accurate_T3(:,1),'b>--',SNRdB,CRB_loc77_accurate_T3(:,2),'r<--',...
   SNRdB,CRB_loc77_T3(:,1),'c^--',SNRdB,CRB_loc77_T3(:,2),'mv--',...
   SNRdB,CRB_loc77_accurate_T2(:,1),'bp-',SNRdB,CRB_loc77_accurate_T2(:,2),'rh-',...
   SNRdB,CRB_loc77_T2(:,1),'cs-',SNRdB,CRB_loc77_T2(:,2),'md-',...
    'MarkerSize',10,'LineWidth',1.2)
xlabel('SNR (dB)')
ylabel('CRLB')
% xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('x, T=1e-3','y, T=1e-3','x, T=1e-3, App.','y, T=1e-3, App.',...
    'x, T=1e-2','y, T=1e-2','x, T=1e-2, App.','y, T=1e-2, App.')
grid on
%
figure
% subplot(2,1,2)
semilogy( SNRdB,CRB_vel77_accurate_T3(:,1),'b>--',SNRdB,CRB_vel77_accurate_T3(:,2),'r<--',...
   SNRdB,CRB_vel77_T3(:,1),'c^--',SNRdB,CRB_vel77_T3(:,2),'mv--',...
   SNRdB,CRB_vel77_accurate_T2(:,1),'bp-',SNRdB,CRB_vel77_accurate_T2(:,2),'rh-',...
   SNRdB,CRB_vel77_T2(:,1),'cs-',SNRdB,CRB_vel77_T2(:,2),'md-',...
    'MarkerSize',10,'LineWidth',1.2)
xlabel('SNR (dB)')
ylabel('CRLB')
% xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('v_x, T=1e-3','v_y, T=1e-3','v_x, T=1e-3, App.','v_y, T=1e-3, App.',...
    'v_x, T=1e-2','v_y, T=1e-2','v_x, T=1e-2, App.','v_y, T=1e-2, App.')
grid on
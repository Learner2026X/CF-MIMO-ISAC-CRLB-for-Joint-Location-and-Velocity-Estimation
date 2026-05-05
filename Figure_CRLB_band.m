close all
SNRdB = -20:5:20;
load('data/CRB_loc77_accurate_B2.mat','CRB_loc_accurate')
CRB_loc77_accurate_B2 = CRB_loc_accurate;
load('data/CRB_vel77_accurate_B2.mat','CRB_vel_accurate')
CRB_vel77_accurate_B2 = CRB_vel_accurate;
load('data/CRB_loc77_accurate_B3.mat','CRB_loc_accurate')
CRB_loc77_accurate_B3 = CRB_loc_accurate;
load('data/CRB_vel77_accurate_B3.mat','CRB_vel_accurate')
CRB_vel77_accurate_B3 = CRB_vel_accurate;
load('data/CRB_loc77_accurate_B1.mat','CRB_loc_accurate')
CRB_loc77_accurate_B1 = CRB_loc_accurate;
load('data/CRB_vel77_accurate_B1.mat','CRB_vel_accurate')
CRB_vel77_accurate_B1 = CRB_vel_accurate;
%%
figure
% subplot(2,1,1)
semilogy( SNRdB,CRB_loc77_accurate_B1(:,1),'b>--',SNRdB,CRB_loc77_accurate_B1(:,2),'r<--',...
   SNRdB,CRB_loc77_accurate_B2(:,1),'bp-.',SNRdB,CRB_loc77_accurate_B2(:,2),'rh-.',...
   SNRdB,CRB_loc77_accurate_B3(:,1),'bs:',SNRdB,CRB_loc77_accurate_B3(:,2),'ro:',...
   SNRdB,CRB_vel77_accurate_B1(:,1),'c^--',SNRdB,CRB_vel77_accurate_B1(:,2),'mv--',...
   SNRdB,CRB_vel77_accurate_B2(:,1),'cd-.',SNRdB,CRB_vel77_accurate_B2(:,2),'mx-.',...
   SNRdB,CRB_vel77_accurate_B3(:,1),'c+:',SNRdB,CRB_vel77_accurate_B3(:,2),'m*:',...
    'MarkerSize',10,'LineWidth',1.2)
xlabel('SNR (dB)')
ylabel('CRLB (m^2)')
% xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('x, M=16','y, M=16','x, M=128','y, M=128','x, M=1e3','y, M=1e3',...
    'v_x, M=16','v_y, M=16','v_x, M=128','v_y, M=128','v_x, M=1e3','v_y, M=1e3')
grid on
% %
% figure
% % subplot(2,1,2)
% semilogy(SNRdB,CRB_vel77_accurate_B1(:,1),'b>--',SNRdB,CRB_vel77_accurate_B1(:,2),'r<--',...
%    SNRdB,CRB_vel77_accurate_B2(:,1),'cp-.',SNRdB,CRB_vel77_accurate_B2(:,2),'mh-.',...
%    SNRdB,CRB_vel77_accurate_B3(:,1),'ks:',SNRdB,CRB_vel77_accurate_B3(:,2),'go:',...
%     'MarkerSize',10,'LineWidth',1.5)
% xlabel('SNR (dB)')
% ylabel('CRLB (m^2/s^2)')
% % xlim([-30 5])
% % ylim([1e-6 1e5+1000])
% ax=gca;
% set(ax,'FontSize',14);
% legend('v_x, M=16','v_y, M=16','v_x, M=128','v_y, M=128','v_x, M=1e3','v_y, M=1e3')
% grid on
close all
fs_db = 0:5:40;
fs = 10.^(fs_db/10);
load('data\CRB_loc77_accurate_fs.mat','CRB_loc_accurate')
CRB_loc77_accurate_fs = CRB_loc_accurate;
load('data\CRB_vel77_accurate_fs.mat','CRB_vel_accurate')
CRB_vel77_accurate_fs = CRB_vel_accurate;
%%
figure
semilogy( fs_db,CRB_loc77_accurate_fs(:,1),'b>--',fs_db,CRB_loc77_accurate_fs(:,2),'r<--',...
   fs_db,CRB_vel77_accurate_fs(:,1),'c^--',fs_db,CRB_vel77_accurate_fs(:,2),'mv--',...
    'MarkerSize',10,'LineWidth',1.2)
xlabel('fs (KHz \cdot dB)')
ylabel('CRLB')
% xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('x','y','v_x','v_y')
grid on

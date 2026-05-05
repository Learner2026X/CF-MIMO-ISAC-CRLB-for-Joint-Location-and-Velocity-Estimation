close all
SNR = -20:5:20;%dB
load('data/CRB_range_accurate.mat','CRB_loc_accurate')
CRB_loc_accurate_128 = CRB_loc_accurate;
load('data/CRB_range_OFDM_accurate.mat','CRB_loc_OFDM_accurate')
CRB_loc_OFDM_accurate_100k = CRB_loc_OFDM_accurate;
load('data/CRB_range_accurate_12.mat','CRB_loc_accurate')
CRB_loc_accurate_12 = CRB_loc_accurate;
%%
figure
semilogy( SNR,CRB_loc_OFDM_accurate_100k(:,1),'bo--',SNR,CRB_loc_OFDM_accurate_100k(:,2),'r>--',...
    SNR,CRB_loc_accurate_12(:,1),'cp-',SNR,CRB_loc_accurate_12(:,2),'md-',...
    SNR,CRB_loc_accurate_128(:,1),'ch-',SNR,CRB_loc_accurate_128(:,2),'ms-',...
    'MarkerSize',10,'LineWidth',1.2)
xlabel('SNR (dB)')
ylabel('CRLB (m^2)')
% xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('x, OFDM','y, OFDM','x, OCDM, M=12','y, OCDM, M=12',...
    'x, OCDM, M=128','y, OCDM, M=128')
grid on

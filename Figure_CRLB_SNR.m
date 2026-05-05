close all
SNR = -20:5:20;%dB
load('data/CRB_range_accurate.mat','CRB_loc_accurate')
CRB_loc_accurate_100k = CRB_loc_accurate;
load('data/CRB_range_OFDM_accurate.mat','CRB_loc_OFDM_accurate')
CRB_loc_OFDM_accurate_100k = CRB_loc_OFDM_accurate;
load('data/CRB_range_accurate_10MHz.mat','CRB_loc_accurate')
load('data/CRB_range_OFDM_accurate_10MHz.mat','CRB_loc_OFDM_accurate')
%%
figure
semilogy( SNR,CRB_loc_OFDM_accurate_100k(:,1),'bo--',SNR,CRB_loc_OFDM_accurate_100k(:,2),'r>--',...
    SNR,CRB_loc_accurate_100k(:,1),'ch--',SNR,CRB_loc_accurate_100k(:,2),'ms--',...
    SNR,CRB_loc_OFDM_accurate(:,1),'bo-',SNR,CRB_loc_OFDM_accurate(:,2),'r>-',...
    SNR,CRB_loc_accurate(:,1),'ch-',SNR,CRB_loc_accurate(:,2),'ms-',...
    'MarkerSize',10,'LineWidth',1.2)
xlabel('SNR (dB)')
ylabel('CRLB (m^2)')
% xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('x, OFDM, 100KHz','y, OFDM, 100KHz','x, OCDM, 100KHz','y, OCDM, 100KHz',...
    'x, OFDM, 10MHz','y, OFDM, 10MHz','x, OCDM, 10MHz','y, OCDM, 10MHz')
grid on

close all
SNR = -30:5:5;%dB
len_SNR = length(SNR);
MSE_location = zeros(len_SNR,2);
MSE_velocity = zeros(len_SNR,2);
for i_SNR=1:len_SNR
    if SNR(i_SNR)<0
        load(['data\MSE_', num2str(-SNR(i_SNR)),'dB_1kHz.mat'],'MSE');
        MSE_location(i_SNR,:) = MSE(1:2);
        MSE_velocity(i_SNR,:) = MSE(3:4);
    else
        load(['data\MSE', num2str(SNR(i_SNR)),'dB_1kHz.mat'],'MSE');
        MSE_location(i_SNR,:) = MSE(1:2);
        MSE_velocity(i_SNR,:) = MSE(3:4);
    end
end
load('data\CRB_loc.mat','CRB_loc')
load('data\CRB_vel.mat','CRB_vel')
%%
figure
% subplot(1,2,1)
semilogy(SNR,MSE_location(:,1),'o--',SNR,MSE_location(:,2),'>--',...
    SNR,CRB_loc_accurate(:,1),'ks-',SNR,CRB_loc_accurate(:,2),'g<-',...
    SNR,MSE_velocity(:,1),'ch--',SNR,MSE_velocity(:,2),'m^--',...
    SNR,CRB_vel_accurate(:,1),'bd-',SNR,CRB_vel_accurate(:,2),'rv-',...
    'MarkerSize',10,'LineWidth',1.5)
xlabel('SNR (dB)')
ylabel('MSE (m^2)')
xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('x, MLE','y, MLE','x, CRLB','y, CRLB',...
    'v_x, MLE','v_y, MLE','v_x, CRLB','v_y, CRLB')
grid on

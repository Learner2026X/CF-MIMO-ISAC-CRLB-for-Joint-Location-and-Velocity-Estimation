clearvars
close all
%%
numTr = 1;% Number of transmitters
numRe = 1;% Number of transmitters
radius = 1e4;
azimuth_Tr = [45 35 25 15]; 
azimuth_Re = [-45 -35 -25 -5];
p_tar = [0,0;0,0].';
% v_tar = [4,5].';
v_tar = [5,5;4,5].';
p_tx = [1e5,0].';
p_rx = p_tx;
% load('p_tx_wide.mat','p_tx')
% load('p_rx_wide.mat','p_rx')
% p_tx = p_tx(:,1:numTr);
% p_rx = p_rx(:,1:numRe);
% p_tx = radius*[cosd(azimuth_Tr);sind(azimuth_Tr)];
% p_rx = radius*[cosd(azimuth_Re);sind(azimuth_Re)];
%%
figure
plot(p_tx(1,:),p_tx(2,:),'bs',p_rx(1,:),p_rx(2,:),'r>',p_tar(1,1),p_tar(2,1),'ko',...
   'MarkerSize',8,'LineWidth',1.2)
xlabel('x/m')
ylabel('y/m')
ax=gca;
set(ax,'FontSize',14);
legend('Tx','Rx','Target 1')
grid on
grid on
%%
rng(2)
Num_Mont = 1200;
c = 3e8;
f_c = 3e9; % carrier frequency
lambda = c/f_c;
transmission_energy = 1; %w
T = 1e-3; % duration of transmitted waveform
M = 10; % number of durations of sampling
Num_chirp = 64; % number of chirps
sampling_time = 2*M*T;
SNRdB = 0;% dB
NSR = 1/10^(SNRdB/10);
load('RCS.mat','RCS')
RCS = RCS(1:numTr,1:numRe);
% RCS = sqrt(1/2)*(randn(numTr,numRe)+1j*randn(numTr,numRe));
p_noise = transmission_energy*NSR; % NO T?
delta_f = 1/T;
fs =2e3/T;
%%
Num_Tar = 1;
RCS_MT = zeros(numTr,numRe,Num_Tar);
Delay = zeros(numTr,numRe,Num_Tar);
Doppler = zeros(numTr,numRe,Num_Tar);
% parpool(16)
% c = parcluster('local'); % build the 'local' cluster object
% nw = c.NumWorkers        % get the number of workers
for i_q = 1:Num_Tar
    dis_tx_tar=sqrt(sum((p_tx-p_tar(:,i_q)).^2));
    dis_rx_tar=sqrt(sum((p_rx-p_tar(:,i_q)).^2));
    delay_tx_tar = dis_tx_tar/c;
    delay_rx_tar = dis_rx_tar/c;
    Delay(:,:,i_q) = delay_tx_tar.' + delay_rx_tar;
%     Delay = repmat(delay_tx_tar.',1, numRe)+repmat(delay_rx_tar,  numTr,1);
    Doppler_tx = -1/lambda*(v_tar(:,i_q).'*(p_tar(:,i_q)-p_tx)./dis_tx_tar);
    Doppler_rx = -1/lambda*(v_tar(:,i_q).'*(p_tar(:,i_q)-p_rx)./dis_rx_tar);
    Doppler(:,:,i_q) = Doppler_tx.' + Doppler_rx;
%     Doppler = repmat(Doppler_tx.',1, numRe)+repmat(Doppler_rx,  numTr,1);
end

t = (-sampling_time/2:1/fs:sampling_time/2-1/fs).';
% t = (0:1/fs:sampling_time).';
Num_sample = length(t);
% transmit_signal = zeros(Num_sample,numTr);
% for n=1:numTr
%     transmit_signal(:,n) = (2/T^2)^(1/4)*exp(-pi*t.^2/T^2).*...
%         exp(1j*pi*Num_chirp/T^2*(t-(n-1)/Num_chirp*T).^2);
% end
%% MLE
tic
rng('shuffle')
received_signal_unit = zeros(Num_sample,  numTr, numRe);
NumSamp_delay = round(Delay*fs);
% received signal
for n = 1:  numTr
    for k = 1: numRe
        for q = 1:Num_Tar
%             received_signal_unit(:,n,k) = ...
%                 (2/T^2)^(1/4)*exp(-pi*(t-Delay(n,k,q)).^2/T^2)...
%                 .*exp(1j*pi*Num_chirp/T^2*...
%                 (t-(n-1)/Num_chirp*T-Delay(n,k,q)).^2).*...
%                 exp(1j*2*pi*Doppler(n,k,q)*t); % OCDM
            received_signal_unit(:,n,k) = received_signal_unit(:,n,k) ...
                + (2/T^2)^(1/4)*exp(-pi*(t-Delay(n,k,q)).^2/T^2)...
                .*exp(1j*2*pi*n*delta_f*(t-Delay(n,k,q))).*...
                exp(1j*2*pi*Doppler(n,k,q)*t); % OFDM

%             received_signal_unit(:,n,k) = ...
%                 (2/T^2)^(1/4).*cosd(pi*Num_chirp/T^2*...
%                 (t-(n-1)/Num_chirp*T-Delay(n,k,q)).^2); % OCDM
%             received_signal_unit(:,n,k) = received_signal_unit(:,n,k) ...
%                 + (2/T^2)^(1/4).*cosd(2*pi*n*delta_f*(t-Delay(n,k,q))); % OFDM


%             received_signal(:,n,k) = awgn(received_signal(:,n,k),SNRdB,'measured');
        end
%         received_signal(:,n,k) = received_signal(:,n,k)+noise_mat(:,n,k);
    end
end
%%
Num_search = 30;
p_tar = [0,0;0,0].';
% v_tar = [4,5].';
v_tar = [5,5;4,5].';
%     x_trust = round(2*sqrt(8.8985e+05*NSR));
%     y_trust = round(2*sqrt(8.3939e+05*NSR));
%     vx_trust = 2*sqrt(3.5244*NSR);
%     vy_trust = 2*sqrt(3.7522*NSR);
% x_trust = 3.2*sqrt(4.5228e3*NSR);
% y_trust = 3.2*sqrt(2.5106e4*NSR);
% vx_trust = 3.2*sqrt(0.6108*NSR);
% vy_trust = 3.2*sqrt(0.5901*NSR);
% x_trust = 3.2*sqrt(27.6724*NSR);
% y_trust = 3.2*sqrt(25.5569*NSR);
% vx_trust = 3.2*sqrt(0.0045*NSR);
% vy_trust = 3.2*sqrt(0.0065*NSR);
x_trust = 5e4;
y_trust = 1e5;
vx_trust = 10;
vy_trust = 10;
x_search = -x_trust+p_tar(1,1):x_trust/Num_search:x_trust+p_tar(1,1);
y_search = -y_trust+p_tar(2,1):y_trust/Num_search:y_trust+p_tar(2,1);
vx_search = -vx_trust+v_tar(1,1):vx_trust/Num_search:vx_trust+v_tar(1,1);
vy_search = -vy_trust+v_tar(2,1):vy_trust/Num_search:vy_trust+v_tar(2,1);
len_x = length(x_search);
len_y = length(y_search);
len_vx = length(vx_search);
len_vy = length(vy_search);
%%
AF = zeros(len_x,len_y,len_vx,len_vy);
for i_x = 1:len_x
        i_x
    for i_y = 1:len_y
            i_y
        for i_vx = 31
%                 i_vx
            for i_vy = 31
                dis_tx_tar=sqrt(sum((p_tx-[x_search(i_x);y_search(i_y)]).^2));
                dis_rx_tar=sqrt(sum((p_rx-[x_search(i_x);y_search(i_y)]).^2));
                delay_tx_tar = dis_tx_tar/c;
                delay_rx_tar = dis_rx_tar/c;
                Delay_search = delay_tx_tar.' + delay_rx_tar;
%                     Delay_search = repmat(delay_tx_tar.',1, numRe)+repmat(delay_rx_tar,  numTr,1);
%                     NumSamp_delay = round(Delay_search*fs);
                Doppler_tx = -1/lambda*([vx_search(i_vx);vy_search(i_vy)].'...
                    *([x_search(i_x);y_search(i_y)]-p_tx)...
                    ./dis_tx_tar);
                Doppler_rx = -1/lambda*([vx_search(i_vx);vy_search(i_vy)].'...
                    *([x_search(i_x);y_search(i_y)]-p_rx)...
                    ./dis_rx_tar);
                Doppler_search = Doppler_tx.' + Doppler_rx;
%                     Doppler_search = repmat(Doppler_tx.',1, numRe)+...
%                         repmat(Doppler_rx,  numTr,1);
                sig_search = zeros(Num_sample,  numTr, numRe);
                for n=1:  numTr
                    for k=1: numRe
%                         sig_search(:,n,k) = ...
%                             (2/T^2)^(1/4)*exp(-pi*(t-Delay_search(n,k)).^2/T^2).*exp(1j*pi*Num_chirp/T^2*...
%                             (t-(n-1)/Num_chirp*T-Delay_search(n,k)).^2).*...
%                             exp(1j*2*pi*Doppler_search(n,k)*t); % OCDM
                        sig_search(:,n,k) = ...
                            (2/T^2)^(1/4)*exp(-pi*(t-Delay_search(n,k)).^2/T^2)...
                            .*exp(1j*2*pi*n*delta_f*(t-Delay_search(n,k))).*...
                            exp(1j*2*pi*Doppler_search(n,k)*t); % OFDM
%                         sig_search(:,n,k) = ...
%                             (2/T^2)^(1/4).*cosd(pi*Num_chirp/T^2*...
%                             (t-(n-1)/Num_chirp*T-Delay_search(n,k)).^2); % OCDM
%                         sig_search(:,n,k) = ...
%                             (2/T^2)^(1/4).*cosd(2*pi*n*delta_f*(t-Delay_search(n,k))); % OFDM
                        AF(i_x,i_y,i_vx,i_vy) = AF(i_x,i_y,i_vx,i_vy)+...
                            abs(RCS(n,k)*received_signal_unit(:,n,k)'*sig_search(:,n,k))^2/fs^2;
                    end
                end
            end
        end
    end
end
[AF_max, index] = max(AF(:));
[ind_x,ind_y,ind_vx,ind_vy] = ind2sub(size(AF), index);
x_est = x_search(ind_x);
y_est = y_search(ind_y);
vx_est = vx_search(ind_vx);
vy_est = vy_search(ind_vy);
%     disp(['complete',num2str(i_Mont)])
toc
%%
Value_estimate = [x_est,y_est,vx_est,vy_est];
AF_OFDM = AF;
save('data\AF_OFDM_Monostatic.mat','AF_OFDM')
% save('data\MSE0dB_1kHz_MT.mat','MSE')
%%
load('data/AF_OFDM_bistatic.mat','AF_OFDM') 
load('data/AF_OCDM_bistatic.mat','AF_OCDM')
% 
figure
subplot(1, 2, 1);
plot(x_search,AF_OFDM(:,31,31,31),'r',x_search,AF_OCDM(:,31,31,31),'b','LineWidth',1.2)
xlabel('x/m')
ylabel('AF')
% xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('OFDM','OCDM')
subplot(1, 2, 2);
plot(y_search,AF_OFDM(31,:,31,31),'r',y_search,AF_OCDM(31,:,31,31),'b','LineWidth',1.2)
xlabel('y/m')
ylabel('AF')
% xlim([-30 5])
% ylim([1e-6 1e5+1000])
ax=gca;
set(ax,'FontSize',14);
legend('OFDM','OCDM')
% subplot(2, 2, 3);
% plot(vx_search,squeeze(AF_OFDM(16,16,:,16)),'r',vx_search,squeeze(AF_OCDM(16,16,:,16)),'b')
% subplot(2, 2, 4);
% plot(vy_search,squeeze(AF_OFDM(16,16,16,:)),'r',vy_search,squeeze(AF_OCDM(16,16,16,:)),'b')

% figure
% subplot(2, 2, 1);
% plot(x_search,AF_OFDM_t(:,16,16,16),'r',x_search,AF_OCDM_t(:,16,16,16),'b')
% subplot(2, 2, 2);
% plot(y_search,AF_OFDM_t(16,:,16,16),'r',y_search,AF_OCDM_t(16,:,16,16),'b')
% subplot(2, 2, 3);
% plot(vx_search,squeeze(AF_OFDM_t(16,16,:,16)),'r',vx_search,squeeze(AF_OCDM_t(16,16,:,16)),'b')
% subplot(2, 2, 4);
% plot(vy_search,squeeze(AF_OFDM_t(16,16,16,:)),'r',vy_search,squeeze(AF_OCDM_t(16,16,16,:)),'b')
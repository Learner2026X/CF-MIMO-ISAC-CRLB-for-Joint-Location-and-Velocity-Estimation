clearvars
close all
%%
Num_Tr = 4;% Number of transmitters
Num_Re = 3;% Number of transmitters
% rng(2)
% p_tx = 10000*rand(2,Num_Tr);
% rng(0)
% p_rx = 10000*rand(2,Num_Re);
p_tar = [4000,5000;0,0].';
% v_tar = [4,5].';
v_tar = [20,30;80,120].';
load('p_tx_wide.mat','p_tx')
load('p_rx_wide.mat','p_rx')
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
T = 1e-2; % duration of transmitted waveform
M = 1.5; % number of durations of sampling
Num_chirp = 16; % number of chirps
sampling_time = 2*M*T;
SNRdB = -30;%dB
NSR = 1/10^(SNRdB/10);
load('RCS.mat','RCS')
RCS1 = sqrt(1/2)*(randn(Num_Tr,Num_Re)+1j*randn(Num_Tr,Num_Re));
p_noise = transmission_energy*NSR; % NO T?
fs =1e1/T;
%%
Num_Tar = 2;
RCS_MT = zeros(Num_Tr,Num_Re,Num_Tar);
RCS_MT(:,:,1) = RCS;
RCS_MT(:,:,2) = RCS1;
Delay = zeros(Num_Tr,Num_Re,Num_Tar);
Doppler = zeros(Num_Tr,Num_Re,Num_Tar);
% parpool(16)
% c = parcluster('local'); % build the 'local' cluster object
% nw = c.NumWorkers        % get the number of workers
for i_q = 1:Num_Tar
    dis_tx_tar=sqrt(sum((p_tx-p_tar(:,i_q)).^2));
    dis_rx_tar=sqrt(sum((p_rx-p_tar(:,i_q)).^2));
    delay_tx_tar = dis_tx_tar/c;
    delay_rx_tar = dis_rx_tar/c;
    Delay(:,:,i_q) = delay_tx_tar.' + delay_rx_tar;
%     Delay = repmat(delay_tx_tar.',1, Num_Re)+repmat(delay_rx_tar,  Num_Tr,1);
    Doppler_tx = -1/lambda*(v_tar(:,i_q).'*(p_tar(:,i_q)-p_tx)./dis_tx_tar);
    Doppler_rx = -1/lambda*(v_tar(:,i_q).'*(p_tar(:,i_q)-p_rx)./dis_rx_tar);
    Doppler(:,:,i_q) = Doppler_tx.' + Doppler_rx;
%     Doppler = repmat(Doppler_tx.',1, Num_Re)+repmat(Doppler_rx,  Num_Tr,1);
end

t = (-sampling_time/2:1/fs:sampling_time/2-1/fs).';
Num_sample = length(t);
% transmit_signal = zeros(Num_sample,Num_Tr);
% for n=1:Num_Tr
%     transmit_signal(:,n) = (2/T^2)^(1/4)*exp(-pi*t.^2/T^2).*...
%         exp(1j*pi*Num_chirp/T^2*(t-(n-1)/Num_chirp*T).^2);
% end
%% MLE
tic
x_est = zeros(Num_Mont,1);
y_est = zeros(Num_Mont,1);
vx_est = zeros(Num_Mont,1);
vy_est = zeros(Num_Mont,1);
MSE_x_MC = zeros(Num_Mont,1);
MSE_y_MC = zeros(Num_Mont,1);
MSE_vx_MC = zeros(Num_Mont,1);
MSE_vy_MC = zeros(Num_Mont,1);
% parpool(8)
parfor (i_Mont = 1:Num_Mont,8)
    i_Mont
    rng('shuffle')
    noise_mat = sqrt(p_noise/2)*(randn(Num_sample,  Num_Tr, Num_Re)+...
        1j*randn(Num_sample,  Num_Tr, Num_Re));
    received_signal = zeros(Num_sample,  Num_Tr, Num_Re);
    NumSamp_delay = round(Delay*fs);
    % received signal
    for n = 1:  Num_Tr
        for k = 1: Num_Re
            for q = 1:Num_Tar
                received_signal(:,n,k) = received_signal(:,n,k) ...
                    + sqrt(transmission_energy)*RCS_MT(n,k,q)*...
                    (2/T^2)^(1/4)*exp(-pi*(t-(n-1)/Num_chirp*T-Delay(n,k,q)).^2/T^2).*...
                    exp(1j*pi*Num_chirp/T^2*(t-(n-1)/Num_chirp*T-Delay(n,k,q)).^2).*...
                    exp(1j*2*pi*Doppler(n,k,q)*t); % Herein Gaussian pluse 
                % should be delayed too with (n-1)/Num_chirp*T?
%             received_signal(:,n,k) = awgn(received_signal(:,n,k),SNRdB,'measured');
            end
            received_signal(:,n,k) = received_signal(:,n,k)+noise_mat(:,n,k);
        end
    end
    %%
    Num_search = 15;
%     x_trust = round(2*sqrt(8.8985e+05*NSR));
%     y_trust = round(2*sqrt(8.3939e+05*NSR));
%     vx_trust = 2*sqrt(3.5244*NSR);
%     vy_trust = 2*sqrt(3.7522*NSR);
    x_trust = 3.2*sqrt(27.6724*NSR);
    y_trust = 3.2*sqrt(25.5569*NSR);
    vx_trust = 3.2*sqrt(0.0045*NSR);
    vy_trust = 3.2*sqrt(0.0065*NSR);
%     x_trust = 10;
%     y_trust = 10;
%     vx_trust = 0.5;
%     vy_trust = 0.5;
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
%         i_x
        for i_y = 1:len_y
%             i_y
            for i_vx = 1:len_vx
%                 i_vx
                for i_vy = 1:len_vy
                    dis_tx_tar=sqrt(sum((p_tx-[x_search(i_x);y_search(i_y)]).^2));
                    dis_rx_tar=sqrt(sum((p_rx-[x_search(i_x);y_search(i_y)]).^2));
                    delay_tx_tar = dis_tx_tar/c;
                    delay_rx_tar = dis_rx_tar/c;
                    Delay_search = delay_tx_tar.' + delay_rx_tar;
%                     Delay_search = repmat(delay_tx_tar.',1, Num_Re)+repmat(delay_rx_tar,  Num_Tr,1);
%                     NumSamp_delay = round(Delay_search*fs);
                    Doppler_tx = -1/lambda*([vx_search(i_vx);vy_search(i_vy)].'...
                        *([x_search(i_x);y_search(i_y)]-p_tx)...
                        ./dis_tx_tar);
                    Doppler_rx = -1/lambda*([vx_search(i_vx);vy_search(i_vy)].'...
                        *([x_search(i_x);y_search(i_y)]-p_rx)...
                        ./dis_rx_tar);
                    Doppler_search = Doppler_tx.' + Doppler_rx;
%                     Doppler_search = repmat(Doppler_tx.',1, Num_Re)+...
%                         repmat(Doppler_rx,  Num_Tr,1);
                    sig_search = zeros(Num_sample,  Num_Tr, Num_Re);
                    for n=1:  Num_Tr
                        for k=1: Num_Re
                            sig_search(:,n,k) = ...
                                (2/T^2)^(1/4)*exp(-pi*(t-(n-1)/Num_chirp*T ...
                                -Delay_search(n,k)).^2/T^2).*exp(1j*pi*Num_chirp/T^2*...
                                (t-(n-1)/Num_chirp*T-Delay_search(n,k)).^2).*...
                                exp(1j*2*pi*Doppler_search(n,k)*t);
                            AF(i_x,i_y,i_vx,i_vy) = AF(i_x,i_y,i_vx,i_vy)+...
                                abs(received_signal(:,n,k)'*sig_search(:,n,k))^2;
                        end
                    end
                end
            end
        end
    end
    [AF_max, index] = max(AF(:));
    [ind_x,ind_y,ind_vx,ind_vy] = ind2sub(size(AF), index);
    x_est(i_Mont) = x_search(ind_x);
    y_est(i_Mont) = y_search(ind_y);
    vx_est(i_Mont) = vx_search(ind_vx);
    vy_est(i_Mont) = vy_search(ind_vy);
    MSE_x_MC(i_Mont) = abs(x_est(i_Mont)-p_tar(1,1))^2;
    MSE_y_MC(i_Mont) = abs(y_est(i_Mont)-p_tar(2,1))^2;
    MSE_vx_MC(i_Mont) = abs(vx_est(i_Mont)-v_tar(1,1))^2;
    MSE_vy_MC(i_Mont) = abs(vy_est(i_Mont)-v_tar(2,1))^2;
%     disp(['complete',num2str(i_Mont)])
end
toc
%%
MSE_x = mean(MSE_x_MC);
MSE_y = mean(MSE_y_MC);
MSE_vx = mean(MSE_vx_MC);
MSE_vy = mean(MSE_vy_MC);
Value_estimate = [x_est,y_est,vx_est,vy_est];
MSE = [MSE_x,MSE_y,MSE_vx,MSE_vy];
save('data\Value_estimate_30dB_1kHz_DFMT.mat','Value_estimate')
save('data\MSE_30dB_1kHz_DFMT.mat','MSE')
%%
% space = 1e1;
% figure
% subplot(3, 1, 1);
% plot(1:space:Num_sample,real(received_signal(1:space:end,4,3)))
% subplot(3, 1, 2);
% plot(1:space:Num_sample,real(transmit_signal(1:space:end,4)))
% subplot(3, 1, 3);
% plot(1:space:Num_sample,real(sig_search(1:space:end,4,3)))
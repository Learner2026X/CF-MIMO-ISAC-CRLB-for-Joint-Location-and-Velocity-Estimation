
%%
rng('shuffle')
M = 10;
N = 12; %number of subcarriers
T = 1e-3;
width = M*T;
fs = 1e5;
fs_integral = 1e5;
len = 2*width*fs_integral;
t = linspace(-width,width,len);
dt = t(2) - t(1);
f = (-len/2:len/2-1)/(len*dt); % Frequency vector
%%
numTr = 7;% Number of transmitters
numRe = 7;% Number of transmitters
radius = 5e3;
% azimuth_Tr = [85 65 45 25 5]; 
% azimuth_Re = [-85 -65 -45 -25 -5];
azimuth_Tr = [105 85 65 45 25 5 -15]+45; 
azimuth_Re = [-105 -85 -65 -45 -25 -5 15]-45;
p_tar = [-0e3,0].';
% len_num = length(num);
% v_tar = [4,5].';
% load('p_tx_wide.mat','p_tx')
% load('p_rx_wide.mat','p_rx')
p_tx = radius*[cosd(azimuth_Tr(1:numTr));sind(azimuth_Tr(1:numTr))];
p_rx = radius*[cosd(azimuth_Re(1:numRe));sind(azimuth_Re(1:numRe))];
%%
close all
figure
plot(p_tx(1,:),p_tx(2,:),'b>',p_rx(1,:),p_rx(2,:),'rs',p_tar(1),p_tar(2),'ko',...
    'MarkerSize',8,'LineWidth',1.2)
xlabel('X/m')
ylabel('Y/m')
ax=gca;
set(ax,'FontSize',14);
legend('Tx','Rx','Target')
grid on
%%
rng(1)
c = 3e8;
f_c = 3e9; % carrier frequency
lambda = c/f_c;
transmission_energy = 1; %w
SNR = -20:5:20;%dB
% SNR = 0;%dB
len_SNR = length(SNR);
NSR = 1./10.^(SNR./10);
% Eff_bandw = 2.045141018530318e+07;% for T=1e-3
% Eff_time = 7.957747154036620e-08;
% TB = 8;
% load('RCS_var.mat','RCS_var')
% load('RCS.mat','RCS')
noise_power = transmission_energy*NSR;
RCS = sqrt(1/2)*(ones(numTr,numRe)+1j*ones(numTr,numRe));
% RCS = sqrt(1/2)*(randn(numTr,numRe)+1j*randn(numTr,numRe));
RCS_var = abs(RCS).^2;
RCS = RCS(1:numTr,1:numRe);
RCS_var = RCS_var(1:numTr,1:numRe);
CRB_loc = zeros(len_SNR,2);
CRB_vel = zeros(len_SNR,2);
CRB_loc_accurate = zeros(len_SNR,2);
CRB_vel_accurate = zeros(len_SNR,2);
for i_SNR = 1:len_SNR
    cons_f = 2*transmission_energy/noise_power(i_SNR)*fs;
%%
Num_Tar = 1;
for i_q = 1:  Num_Tar
    dis_tx_tar=sqrt(sum((p_tx-p_tar(:,i_q)).^2));
    dis_rx_tar=sqrt(sum((p_rx-p_tar(:,i_q)).^2));
    delay_tx_tar = dis_tx_tar/c;
    delay_rx_tar = dis_rx_tar/c;
    Delay = repmat(delay_tx_tar.',1, numRe)+repmat(delay_rx_tar,  numTr,1);
end
%% intermediate parameter
zero_mat = zeros(numTr*numRe);
cons = 2*fs/noise_power(i_SNR);
cons_A = zeros(numTr,numRe);
cons_B = zeros(numTr,numRe);
cons_D = zeros(numTr,numRe);
G_Re = zeros(numTr*numRe);
G_Im = zeros(numTr*numRe);
E_Re = zeros(numTr*numRe);
E_Im = zeros(numTr*numRe);
PA = ones(numTr,1);
weight_P = zeros(3,numTr,numRe);
weight_V = zeros(4,numTr,numRe);
weight_Y = zeros(3,numTr,numRe);
for i_tar = 1:1
    beta = zeros(2,numTr*numRe);
    eta = zeros(2,numTr*numRe);
    xi = zeros(2,numTr*numRe);
    for n = 1:numTr
        for k = 1:numRe
            temp = (n-1)*numRe+k;
            t_delay = t-Delay(n,k);
            x = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t_delay-(n-1)*T/N).^2/T^2);
%             x = (2/T^2)^(1/4)*exp(-pi*t_delay.^2/T^2).*...
%                 exp(1j*pi*N/T^2*(t_delay-(n-1)*T/N).^2);
            c1 = -2*pi*(1j*N-1)*(t_delay-(n-1)*T/N)/T^2;
%             c1 = -2*pi/T^2*((1j*N-1)*t_delay-1j*(n-1)*T);
            y_gRe = transmission_energy*PA(n).*c1.*x.*conj(x);
            G_Re(temp,temp) = -cons*real(conj(RCS(n,k))*trapz(t, y_gRe));
            G_Im(temp,temp) = cons*imag(conj(RCS(n,k))*trapz(t, y_gRe));
            y = t.*conj(x).*x.*c1;
            % Calculate the Fourier Transform of the signal
            X = fftshift(fft(x));
            % Calculate effective pulse width (EPW)
            Eff_time = trapz(t, t.^2 .* abs(x).^2);
            % Calculate effective bandwidth (EBW)
            Eff_bandw = trapz(dt, f.^2 .* abs(X).^2)/len;
%             sum(f.^2.*abs(X).^2) *dt/len
            TB = imag(sum(y)*dt);
            c_a = 4*pi^2/T^4*(1-1j*N)^2*(t-Delay(n,k)-(n-1)*T/N).^2-2*pi/T^2*(1-1j*N);
%             y1 = t.*conj(x).*x.*c1;
            y_a_sing = transmission_energy*PA(n).*c_a.*x.*conj(x);
            cons_A(n,k) = -cons*real(trapz(t, y_a_sing)*RCS_var(n,k));
%             cons_A(n,k) = 8*pi^2*transmission_energy*RCS_var(n,k)*Eff_bandw/noise_power(i_SNR)*fs; % OCDM
            beta(:,temp) = 1/c*((p_tar(:,i_tar)-p_tx(:,n))/norm(p_tar(:,i_tar)-p_tx(:,n))+...
                (p_tar(:,i_tar)-p_rx(:,k))/norm(p_tar(:,i_tar)-p_rx(:,k)));            
        end
    end
    aleph = blkdiag(beta,eye(2*numTr*numRe));
end
A_mat = cons_A.';
A = diag(A_mat(:));
% E_Im_mat = (cons_e*real(RCS).*Delay.*repmat(PA,1,numRe)).';
G = [G_Re,G_Im];
F = diag(repmat(cons_f*PA,2*numRe,1));
J = [A,G;G.',F];
FIM = aleph*J*aleph.';
CRB_mat = inv(FIM);
CRB_loc_accurate(i_SNR,:) = [CRB_mat(1,1),CRB_mat(2,2)];
end
%%
% save('data\CRB_range.mat','CRB_loc')
save('data\CRB_range_accurate_12.mat','CRB_loc_accurate')
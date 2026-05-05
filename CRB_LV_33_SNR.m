
%%
rng('shuffle')
M = 10;
N = 128;
T = 1e-2;
width = M*T;
fs = 1e3;
fs_integral = 1e5;
len = 2*width*fs_integral;
t = linspace(-width,width,len);
dt = t(2) - t(1);
f = (-len/2:len/2-1)/(len*dt); % Frequency vector
%%
numTr = 3;% Number of transmitters
numRe = 3;% Number of transmitters
radius = 5e3;
% azimuth_Tr = [85 65 45 25 5]; 
% azimuth_Re = [-85 -65 -45 -25 -5];
azimuth_Tr = [65 45 25]+45; 
azimuth_Re = [-65 -45 -25]-45;
p_tar = [0,0;0,0].';
% v_tar = [4,5].';
v_tar = [-15,0;4,5].';
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
noise_power = transmission_energy*NSR;
RCS = sqrt(1/2)*(ones(numTr,numRe)+1j*ones(numTr,numRe));
% RCS = 0.1*(ones(numTr,numRe)+1j*ones(numTr,numRe));  % Deep fading (DF)
RCS_var = abs(RCS).^2;
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
%             x = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t_delay-(n-1)*T/N).^2/T^2);
            x = (2/T^2)^(1/4)*exp(-pi*t_delay.^2/T^2).*...
                exp(1j*pi*N/T^2*(t_delay-(n-1)*T/N).^2);
%             c1 = -2*pi*(1j*N-1)*(t_delay-(n-1)*T/N)/T^2;
            c1 = -2*pi/T^2*((1j*N-1)*t_delay-1j*(n-1)*T);
            y_gRe = transmission_energy*PA(n).*c1.*x.*conj(x);
            G_Re(temp,temp) = -cons*real(conj(RCS(n,k))*trapz(t, y_gRe));
            G_Im(temp,temp) = cons*imag(conj(RCS(n,k))*trapz(t, y_gRe));
            y_eRe = transmission_energy*PA(n)*1j*2*pi.*t.*x.*conj(x);
            E_Re(temp,temp) = -cons*real(conj(RCS(n,k))*trapz(t, y_eRe));
            E_Im(temp,temp) = cons*imag(conj(RCS(n,k))*trapz(t, y_eRe));
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
            y_b_sing = transmission_energy*PA(n)*1j*2*pi.*t.*c1.*x.*conj(x);
            y_d_sing = -transmission_energy*PA(n)*4*pi^2.*t.^2.*x.*conj(x);
            cons_A(n,k) = -cons*real(trapz(t, y_a_sing)*RCS_var(n,k));
            cons_B(n,k) = -cons*real(trapz(t, y_b_sing)*RCS_var(n,k));
            cons_D(n,k) = -cons*real(trapz(t, y_d_sing)*RCS_var(n,k));
%             cons_A(n,k) = 8*pi^2*transmission_energy*RCS_var(n,k)*Eff_bandw/noise_power(i_SNR)*fs; % OCDM
%             cons_B(n,k) = 4*pi*transmission_energy*RCS_var(n,k)*TB/noise_power(i_SNR)*fs; 
% %             cons_B = 0;
%             cons_D(n,k) = 8*pi^2*transmission_energy*RCS_var(n,k)*Eff_time/noise_power(i_SNR)*fs; 
            beta(:,temp) = 1/c*((p_tar(:,i_tar)-p_tx(:,n))/norm(p_tar(:,i_tar)-p_tx(:,n))+...
                (p_tar(:,i_tar)-p_rx(:,k))/norm(p_tar(:,i_tar)-p_rx(:,k)));
            eta(1,temp) = 1/lambda*...
                ((p_tar(1,i_tar)-p_tx(1,n))*v_tar(:,i_tar).'*...
               (p_tar(:,i_tar)-p_tx(:,n))/norm(p_tar(:,i_tar)-p_tx(:,n))^3-...
               v_tar(1,i_tar)/norm(p_tar(:,i_tar)-p_tx(:,n))+...
               (p_tar(1,i_tar)-p_rx(1,k))*v_tar(:,i_tar).'*...
               (p_tar(:,i_tar)-p_rx(:,k))/norm(p_tar(:,i_tar)-p_rx(:,k))^3-...
               v_tar(1,i_tar)/norm(p_tar(:,i_tar)-p_rx(:,k)));
            eta(2,temp) = 1/lambda*...
                ((p_tar(2,i_tar)-p_tx(2,n))*v_tar(:,i_tar).'*...
               (p_tar(:,i_tar)-p_tx(:,n))/norm(p_tar(:,i_tar)-p_tx(:,n))^3-...
               v_tar(2,i_tar)/norm(p_tar(:,i_tar)-p_tx(:,n))+...
               (p_tar(2,i_tar)-p_rx(2,k))*v_tar(:,i_tar).'*...
               (p_tar(:,i_tar)-p_rx(:,k))/norm(p_tar(:,i_tar)-p_rx(:,k))^3-...
               v_tar(2,i_tar)/norm(p_tar(:,i_tar)-p_rx(:,k)));
            xi(1,temp) = -1/lambda*((p_tar(1,i_tar)-p_tx(1,n))/...
               norm(p_tar(:,i_tar)-p_tx(:,n))+(p_tar(1,i_tar)-p_rx(1,k))/...
               norm(p_tar(:,i_tar)-p_rx(:,k)));
            xi(2,temp) = -1/lambda*((p_tar(2,i_tar)-p_tx(2,n))/...
               norm(p_tar(:,i_tar)-p_tx(:,n))+(p_tar(2,i_tar)-p_rx(2,k))/...
               norm(p_tar(:,i_tar)-p_rx(:,k)));
            weight_P(:,n,k) = [beta(1,temp)^2,beta(1,temp)*beta(2,temp),beta(2,temp)^2].'...
                *cons_A(n,k)+[2*beta(1,temp)*eta(1,temp),...
                beta(1,temp)*eta(2,temp)+eta(1,temp)*beta(2,temp),2*beta(2,temp)*eta(2,temp)].'...
                *cons_B(n,k)+[eta(1,temp)^2,eta(1,temp)*eta(2,temp),eta(2,temp)^2].'...
                *cons_D(n,k);
            weight_V(:,n,k) = [beta(1,temp)*xi(1,temp),beta(2,temp)*xi(1,temp)...
                beta(1,temp)*xi(2,temp),beta(2,temp)*xi(2,temp)].'*cons_B(n,k)...
                +[eta(1,temp)*xi(1,temp),eta(2,temp)*xi(1,temp),eta(1,temp)*xi(2,temp),...
                eta(2,temp)*xi(2,temp)].'*cons_D(n,k);
            weight_Y(:,n,k) = [xi(1,temp)^2,xi(1,temp)*xi(2,temp),xi(2,temp)^2].'...
                *cons_D(n,k);
        end
    end
    aleph = blkdiag([beta,eta;zeros(2,numTr*numRe),xi],eye(2*numTr*numRe));
end
A_mat = cons_A.';
A = diag(A_mat(:));
B_mat = cons_B.';
B = diag(B_mat(:));
D_mat = cons_D.';
D = diag(D_mat(:));
% E_Im_mat = (cons_e*real(RCS).*Delay.*repmat(PA,1,numRe)).';
G = [G_Re,G_Im];
E = [E_Re,E_Im];
F = diag(repmat(cons_f*PA,2*numRe,1));
J = [A,B,G;B,D,E;G.', E.',F];
J_appr = blkdiag([A,B;B,D],F);
FIM = aleph*J*aleph.';
FIM_appr = aleph*J_appr*aleph.';
CRB_mat = inv(FIM);
CRB_mat_appr = inv(FIM_appr);
CRB_loc_accurate(i_SNR,:) = [CRB_mat(1,1),CRB_mat(2,2)];
CRB_vel_accurate(i_SNR,:) = [CRB_mat(3,3),CRB_mat(4,4)];
CRB_loc(i_SNR,:) = [CRB_mat_appr(1,1),CRB_mat_appr(2,2)];
CRB_vel(i_SNR,:) = [CRB_mat_appr(3,3),CRB_mat_appr(4,4)];
% rank(FIM(1:4,5:end))
% AA = FIM(1:4,1:4);
% BB = FIM(1:4,5:end)*inv(FIM(5:end,5:end))*FIM(1:4,5:end).';
% bb_eig = eig(BB);
% aa_eig = eig(AA);
% bb_norm = norm(BB);
% aa_norm = norm(AA);
end
%%
save('data\CRB_loc33.mat','CRB_loc')
save('data\CRB_vel33.mat','CRB_vel')
save('data\CRB_loc33_accurate.mat','CRB_loc_accurate')
save('data\CRB_vel33_accurate.mat','CRB_vel_accurate')
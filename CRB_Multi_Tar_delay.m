%%
close all
M = 10;
N = 128;
T = 1e-3;
width = M*T;
fs_integral = 1e5;
len = 2*width*fs_integral;
fs = 1e5;
t = linspace(-width,width,len);
dt = t(2) - t(1);
f = (-len/2:len/2-1)/(len*dt); % Frequency vector
df = f(2)-f(1);

%%
numTr = 7;% Number of transmitters
numRe = 7;% Number of transmitters
radius = 5e3;
% azimuth_Tr = [85 65 45 25 5]; 
% azimuth_Re = [-85 -65 -45 -25 -5];
azimuth_Tr = [105 85 65 45 25 5 -15]+45; 
azimuth_Re = [-105 -85 -65 -45 -25 -5 15]-45;
% num = [3:2:29,29.05:0.05:30.95,31:2:51];
num = -9:2:9;
% num = 20;
len_num = length(num);
% v_tar = [4,5].';
% load('p_tx_wide.mat','p_tx')
% load('p_rx_wide.mat','p_rx')
p_tx = radius*[cosd(azimuth_Tr(1:numTr));sind(azimuth_Tr(1:numTr))];
p_rx = radius*[cosd(azimuth_Re(1:numRe));sind(azimuth_Re(1:numRe))];
%%
p_tar = [0e3,0].';
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
SNR = -0;%dB
NSR = 1/10^(SNR/10);
% Eff_bandw = 2.045141018530318e+07;% for T=1e-3
% Eff_time = 7.957747154036620e-08;
% TB = 8;
noise_power = transmission_energy*NSR;
% RCS = sqrt(1/2)*(randn(numTr,numRe)+1j*randn(numTr,numRe));
RCS = sqrt(1/2)*(ones(numTr,numRe)+1j*ones(numTr,numRe));
RCS_var = abs(RCS).^2;
% load('RCS_var.mat','RCS_var')
% load('RCS.mat','RCS') 
RCS = RCS(1:numTr,1:numRe);
RCS_var = RCS_var(1:numTr,1:numRe);
rng(2)
% RCS1 = sqrt(1/2)*(randn(numTr,numRe)+1j*randn(numTr,numRe));
RCS1 = 0.2*RCS;
%%
Num_Tar = 2;
CRB_loc = zeros(4,len_num);
CRB_vel = zeros(4,len_num);
CRB_loc_single = zeros(4,len_num);
CRB_vel_single = zeros(4,len_num);
CRB_loc_accurate = zeros(4,len_num);
CRB_vel_accurate = zeros(4,len_num);
RCS_sq_MT = zeros(numTr,numRe,Num_Tar);
RCS_sq_MT(:,:,1) = RCS_var;
RCS_sq_MT(:,:,2) = abs(RCS1).^2;
RCS_MT = zeros(numTr,numRe,Num_Tar);
RCS_MT(:,:,1) = RCS;
RCS_MT(:,:,2) = RCS1;
RCS_product = RCS.*conj(RCS1);
tic
% space = 3e3;
space = 120;
for i_dis = 1:len_num
    p_tar = [-0e3,0;0,0].';
%     p_tar = [-0e3,0;-0e3+num(i_dis)*space,0].';
    v_tar = [-15,0;-15+num(i_dis)*space,0].';
%     v_tar = [-15,0;-15,0].';
    Delay = zeros(numTr,numRe,Num_Tar);
    Doppler = zeros(numTr,numRe,Num_Tar);
    for i_tar = 1:Num_Tar
        dis_tx_tar=sqrt(sum((p_tx-p_tar(:,i_tar)).^2));
        dis_rx_tar=sqrt(sum((p_rx-p_tar(:,i_tar)).^2));
        delay_tx_tar = dis_tx_tar/c;
        delay_rx_tar = dis_rx_tar/c;
        Delay(:,:,i_tar) = delay_tx_tar.' + delay_rx_tar;
    %     Delay = repmat(delay_tx_tar.',1, numRe)+repmat(delay_rx_tar,  numTr,1);
        Doppler_tx = -1/lambda*(v_tar(:,i_tar).'*(p_tar(:,i_tar)-p_tx)./dis_tx_tar);
        Doppler_rx = -1/lambda*(v_tar(:,i_tar).'*(p_tar(:,i_tar)-p_rx)./dis_rx_tar);
        Doppler(:,:,i_tar) = Doppler_tx.' + Doppler_rx;
    %     Doppler = repmat(Doppler_tx.',1, numRe)+repmat(Doppler_rx,  numTr,1);
    end
    Delay_diff = Delay(:,:,1)-Delay(:,:,end);
    Doppler_diff = Doppler(:,:,1)-Doppler(:,:,end);
    %% Jql
    PA = ones(numTr,1);
    const12 = 2*fs/noise_power;
    A_12 = zeros(numTr*numRe);
    B1 = zeros(numTr*numRe);
    B2 = zeros(numTr*numRe);
    D_12 = zeros(numTr*numRe);
    G1_Re = zeros(numTr*numRe);
    G1_Im = zeros(numTr*numRe);
    G2_Re = zeros(numTr*numRe);
    G2_Im = zeros(numTr*numRe);
    E1_Re = zeros(numTr*numRe);
    E1_Im = zeros(numTr*numRe);
    E2_Re = zeros(numTr*numRe);
    E2_Im = zeros(numTr*numRe);
    F_12 = zeros(2*numTr*numRe);
    for n = 1:numTr
        for k = 1:numRe
            temp = (n-1)*numRe+k;
            x1 = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t-Delay(n,k,1)-(n-1)*T/N).^2/T^2);
            x2 = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t-Delay(n,k,2)-(n-1)*T/N).^2/T^2);
            c1 = -2*pi*(1j*N-1)*(t-Delay(n,k,1)-(n-1)*T/N)/T^2;
            c2 = 2*pi*(1j*N+1)*(t-Delay(n,k,2)-(n-1)*T/N)/T^2;
            y_a = transmission_energy*PA(n).*c1.*c2.*x1.*conj(x2)...
                .*exp(1j*2*pi*Doppler_diff(n,k)*t);
            A_12(temp,temp) = const12*real(trapz(t, y_a)*RCS_product(n,k));
            y_b1 = -transmission_energy*PA(n)*1j*2*pi.*t.*c1.*x1.*conj(x2)...
                .*exp(1j*2*pi*Doppler_diff(n,k)*t);
            B1(temp,temp) = const12*real(trapz(t, y_b1)*RCS_product(n,k));
            y_b2 = transmission_energy*PA(n)*1j*2*pi.*t.*c2.*x1.*conj(x2)...
                .*exp(1j*2*pi*Doppler_diff(n,k)*t);
            B2(temp,temp) = const12*real(trapz(t, y_b2)*RCS_product(n,k));
            y_d = transmission_energy*PA(n)*4*pi^2.*t.^2.*x1.*conj(x2)...
                .*exp(1j*2*pi*Doppler_diff(n,k)*t);
            D_12(temp,temp) = const12*real(trapz(t, y_d)*RCS_product(n,k));
            y_g1Re = transmission_energy*PA(n).*c1.*x1.*conj(x2)...
                .*exp(1j*2*pi*Doppler_diff(n,k)*t);
            G1_Re(temp,temp) = const12*real(trapz(t, y_g1Re)*RCS_MT(n,k,1));
            G1_Im(temp,temp) = const12*imag(trapz(t, y_g1Re)*RCS_MT(n,k,1));
            y_g2Re = transmission_energy*PA(n).*c2.*x1.*conj(x2)...
                .*exp(1j*2*pi*Doppler_diff(n,k)*t);
            G2_Re(temp,temp) = const12*real(trapz(t, y_g2Re)*conj(RCS_MT(n,k,2)));
            G2_Im(temp,temp) = -const12*imag(trapz(t, y_g2Re)*conj(RCS_MT(n,k,2)));
            y_e1Re = transmission_energy*PA(n)*1j*2*pi.*t.*x1.*conj(x2)...
                .*exp(1j*2*pi*Doppler_diff(n,k)*t);
            E1_Re(temp,temp) = const12*real(trapz(t, y_e1Re)*RCS_MT(n,k,1));
            E1_Im(temp,temp) = const12*imag(trapz(t, y_e1Re)*RCS_MT(n,k,1));
            y_e2Re = -transmission_energy*PA(n)*1j*2*pi.*t.*x1.*conj(x2)...
                .*exp(1j*2*pi*Doppler_diff(n,k)*t);
            E2_Re(temp,temp) = const12*real(trapz(t, y_e2Re)*conj(RCS_MT(n,k,2)));
            E2_Im(temp,temp) = -const12*imag(trapz(t, y_e2Re)*conj(RCS_MT(n,k,2)));
            y_12 = transmission_energy*PA(n).*x1.*conj(x2)...
                .*exp(1j*2*pi*Doppler_diff(n,k)*t);
            Y_12 = trapz(t, y_12);
            F_12(temp,temp) = const12*real(Y_12);
            F_12(temp,numTr*numRe+temp) = const12*imag(Y_12);
            F_12(numTr*numRe+temp,temp) = -const12*imag(Y_12);
            F_12(numTr*numRe+temp,numTr*numRe+temp) = const12*real(Y_12);
        end
    end
    J_12 = [A_12,B1,G1_Re,G1_Im;B2,D_12,E1_Re,E1_Im;...
        [G2_Re,G2_Im].', [E2_Re,E2_Im].',F_12];
    %% Alpha
    aleph = zeros(4+2*numTr*numRe,4*numTr*numRe,Num_Tar); 
    for i_tar = 1:Num_Tar
        beta = zeros(2,numTr*numRe);
        eta = zeros(2,numTr*numRe);
        xi = zeros(2,numTr*numRe);
        for n = 1:numTr
            for k = 1:numRe
                temp = (n-1)*numRe+k;
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
            end
        end
        aleph(:,:,i_tar) = blkdiag([beta,eta;zeros(2,numTr*numRe),xi],eye(2*numTr*numRe));
    end
    %% Diagonal
    FIM = zeros(4+2*numTr*numRe,4+2*numTr*numRe,Num_Tar);
    FIM_single = zeros(4+2*numTr*numRe,4+2*numTr*numRe,Num_Tar);
    cons_f = 2*transmission_energy/noise_power*fs;
    cons = 2*fs/noise_power;
    cons_A = zeros(numTr,numRe);
    cons_B = zeros(numTr,numRe);
    cons_D = zeros(numTr,numRe);
    cons_A_single = zeros(numTr,numRe);
    cons_B_single = zeros(numTr,numRe);
    cons_D_single = zeros(numTr,numRe);
    G_Re = zeros(numTr*numRe);
    G_Im = zeros(numTr*numRe);
    E_Re = zeros(numTr*numRe);
    E_Im = zeros(numTr*numRe);
    G_Re_single = zeros(numTr*numRe);
    G_Im_single = zeros(numTr*numRe);
    E_Re_single = zeros(numTr*numRe);
    E_Im_single = zeros(numTr*numRe);
    int_A = zeros(numTr,numRe);
    int_B = zeros(numTr,numRe);
    int_D = zeros(numTr,numRe);
    int_G_Re = zeros(numTr*numRe);
    int_G_Im = zeros(numTr*numRe);
    int_E_Re = zeros(numTr*numRe);
    int_E_Im = zeros(numTr*numRe);
    %% target 1
    for n = 1:numTr
        for k = 1:numRe
            temp = (n-1)*numRe+k;
    %             t_delay = t-Delay(n,k,i_tar);
    %             x = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t_delay-(n-1)*T/N).^2/T^2);
    %             x = (2/T^2)^(1/4)*exp(-pi*t_delay.^2/T^2).*...
    %                 exp(1j*pi*N/T^2*(t_delay-(n-1)*T/N).^2);
    %             c1 = -2*pi*(1j*N-1)*(t_delay-(n-1)*T/N)/T^2;
    %             c1 = -2*pi/T^2*((1j*N-1)*t_delay-1j*(n-1)*T);
            x1 = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t-Delay(n,k,1)-(n-1)*T/N).^2/T^2);
            x2 = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t-Delay(n,k,2)-(n-1)*T/N).^2/T^2);
            c1 = -2*pi*(1j*N-1)*(t-Delay(n,k,1)-(n-1)*T/N)/T^2;
            y_gRe = transmission_energy*PA(n).*c1.*x1.*conj(x1);
            i_gRe = transmission_energy*PA(n).*c1.*x1.*conj(x2)...
                .*exp(1j*2*pi*Doppler_diff(n,k)*t);
            int_G_Re(temp,temp) = -cons*real(conj(RCS_MT(n,k,2))*trapz(t, i_gRe));
            int_G_Im(temp,temp) = cons*imag(conj(RCS_MT(n,k,2))*trapz(t, i_gRe));
            G_Re_single(temp,temp) = -cons*real(conj(RCS_MT(n,k,1))*trapz(t, y_gRe));
            G_Im_single(temp,temp) = cons*imag(conj(RCS_MT(n,k,1))*trapz(t, y_gRe));
            G_Re(temp,temp) = int_G_Re(temp,temp)+G_Re_single(temp,temp);
            G_Im(temp,temp) = int_G_Im(temp,temp)+G_Im_single(temp,temp);
            y_eRe = transmission_energy*PA(n)*1j*2*pi.*t.*x1.*conj(x1);
            i_eRe = transmission_energy*PA(n)*1j*2*pi.*t.*x1.*conj(x2)...
                .*exp(1j*2*pi*Doppler_diff(n,k)*t);
            int_E_Re(temp,temp) = -cons*real(conj(RCS_MT(n,k,2))*trapz(t, i_eRe));
            int_E_Im(temp,temp) = cons*imag(conj(RCS_MT(n,k,2))*trapz(t, i_eRe));
            E_Re_single(temp,temp) = -cons*real(conj(RCS_MT(n,k,1))*trapz(t, y_eRe));
            E_Im_single(temp,temp) = cons*imag(conj(RCS_MT(n,k,1))*trapz(t, y_eRe));
            E_Re(temp,temp) = int_E_Re(temp,temp) + E_Re_single(temp,temp);
            E_Im(temp,temp) = int_E_Im(temp,temp) + E_Im_single(temp,temp);
            c_a = 4*pi^2/T^4*(1-1j*N)^2*(t-Delay(n,k,1)-(n-1)*T/N).^2-2*pi/T^2*(1-1j*N);
            y_a = transmission_energy*PA(n).*c_a.*x1.*conj(x2)...
                .*exp(1j*2*pi*Doppler_diff(n,k)*t);
            int_A(n,k) = -cons*real(trapz(t, y_a)*RCS_product(n,k));
            y_b = transmission_energy*PA(n)*1j*2*pi.*t.*c1.*x1.*conj(x2)...
                .*exp(1j*2*pi*Doppler_diff(n,k)*t);
            int_B(n,k) = -cons*real(trapz(t, y_b)*RCS_product(n,k));
            y_d = -transmission_energy*PA(n)*4*pi^2.*t.^2.*x1.*conj(x2)...
                .*exp(1j*2*pi*Doppler_diff(n,k)*t);
            int_D(n,k) = -cons*real(trapz(t, y_d)*RCS_product(n,k));
    
            y1 = t.*conj(x1).*x1.*c1;
            y_a_sing = transmission_energy*PA(n).*c_a.*x1.*conj(x1);
            y_b_sing = transmission_energy*PA(n)*1j*2*pi.*t.*c1.*x1.*conj(x1);
            y_d_sing = -transmission_energy*PA(n)*4*pi^2.*t.^2.*x1.*conj(x1);
            % Calculate the Fourier Transform of the signal
            X1 = fftshift(fft(x1));
            % Calculate effective pulse width (EPW)
            Eff_time = trapz(t, t.^2 .* abs(x1).^2);
            % Calculate effective bandwidth (EBW)
            Eff_bandw = trapz(dt, f.^2 .* abs(X1).^2)/len;
            cons_A_single(n,k) = -cons*real(trapz(t, y_a_sing)*RCS_sq_MT(n,k,1));
            cons_A(n,k) = int_A(n,k)+cons_A_single(n,k); % OCDM
            cons_B_single(n,k) = -cons*real(trapz(t, y_b_sing)*RCS_sq_MT(n,k,1));
            cons_B(n,k) = int_B(n,k)+cons_B_single(n,k);
            cons_D_single(n,k) = -cons*real(trapz(t, y_d_sing)*RCS_sq_MT(n,k,1));
            % cons_B = 0;
            cons_D(n,k) = int_D(n,k)+cons_D_single(n,k);
                
        end
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
    FIM(:,:,1) = aleph(:,:,1)*J*aleph(:,:,1).';
    %
    A_mat_single = cons_A_single.';
    A_single = diag(A_mat_single(:));
    B_mat_single = cons_B_single.';
    B_single = diag(B_mat_single(:));
    D_mat_single = cons_D_single.';
    D_single = diag(D_mat_single(:));
    % E_Im_mat = (cons_e*real(RCS).*Delay.*repmat(PA,1,numRe)).';
    G_single = [G_Re_single,G_Im_single];
    E_single = [E_Re_single,E_Im_single];
    J_single = [A_single,B_single,G_single;B_single,D_single,E_single;...
        G_single.', E_single.',F];
    FIM_single(:,:,1) = aleph(:,:,1)*J_single*aleph(:,:,1).';
    %% target 2
    for n = 1:numTr
        for k = 1:numRe
            temp = (n-1)*numRe+k;
    %             t_delay = t-Delay(n,k,i_tar);
    %             x = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t_delay-(n-1)*T/N).^2/T^2);
    %             x = (2/T^2)^(1/4)*exp(-pi*t_delay.^2/T^2).*...
    %                 exp(1j*pi*N/T^2*(t_delay-(n-1)*T/N).^2);
    %             c1 = -2*pi*(1j*N-1)*(t_delay-(n-1)*T/N)/T^2;
    %             c1 = -2*pi/T^2*((1j*N-1)*t_delay-1j*(n-1)*T);
            x1 = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t-Delay(n,k,1)-(n-1)*T/N).^2/T^2);
            x2 = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t-Delay(n,k,2)-(n-1)*T/N).^2/T^2);
            c2 = -2*pi*(1j*N-1)*(t-Delay(n,k,2)-(n-1)*T/N)/T^2;
            y_gRe = transmission_energy*PA(n).*c2.*x2.*conj(x2);
            i_gRe = transmission_energy*PA(n).*c2.*x2.*conj(x1)...
                .*exp(-1j*2*pi*Doppler_diff(n,k)*t);
            int_G_Re(temp,temp) = -cons*real(conj(RCS_MT(n,k,1))*trapz(t, i_gRe));
            int_G_Im(temp,temp) = cons*imag(conj(RCS_MT(n,k,1))*trapz(t, i_gRe));
            G_Re_single(temp,temp) = -cons*real(conj(RCS_MT(n,k,2))*trapz(t, y_gRe));
            G_Im_single(temp,temp) = cons*imag(conj(RCS_MT(n,k,2))*trapz(t, y_gRe));
            G_Re(temp,temp) = int_G_Re(temp,temp)+G_Re_single(temp,temp);
            G_Im(temp,temp) = int_G_Im(temp,temp)+G_Im_single(temp,temp);
            y_eRe = transmission_energy*PA(n)*1j*2*pi.*t.*x2.*conj(x2);
            i_eRe = transmission_energy*PA(n)*1j*2*pi.*t.*x2.*conj(x1)...
                .*exp(-1j*2*pi*Doppler_diff(n,k)*t);
            int_E_Re(temp,temp) = -cons*real(conj(RCS_MT(n,k,1))*trapz(t, i_eRe));
            int_E_Im(temp,temp) = cons*imag(conj(RCS_MT(n,k,1))*trapz(t, i_eRe));
            E_Re_single(temp,temp) = -cons*real(conj(RCS_MT(n,k,2))*trapz(t, y_eRe));
            E_Im_single(temp,temp) = cons*imag(conj(RCS_MT(n,k,2))*trapz(t, y_eRe));
            E_Re(temp,temp) = int_E_Re(temp,temp) + E_Re_single(temp,temp);
            E_Im(temp,temp) = int_E_Im(temp,temp) + E_Im_single(temp,temp);
            c_a = 4*pi^2/T^4*(1-1j*N)^2*(t-Delay(n,k,2)-(n-1)*T/N).^2-2*pi/T^2*(1-1j*N);
            y_a = transmission_energy*PA(n).*c_a.*x2.*conj(x1)...
                .*exp(-1j*2*pi*Doppler_diff(n,k)*t);
            int_A(n,k) = -cons*real(trapz(t, y_a)*conj(RCS_product(n,k)));
            y_b = transmission_energy*PA(n)*1j*2*pi.*t.*c2.*x2.*conj(x1)...
                .*exp(-1j*2*pi*Doppler_diff(n,k)*t);
            int_B(n,k) = -cons*real(trapz(t, y_b)*conj(RCS_product(n,k)));
            y_d = -transmission_energy*PA(n)*4*pi^2.*t.^2.*x2.*conj(x1)...
                .*exp(-1j*2*pi*Doppler_diff(n,k)*t);
            int_D(n,k) = -cons*real(trapz(t, y_d)*conj(RCS_product(n,k)));
    
            y2 = t.*conj(x2).*x2.*c2;
            y_a_sing = transmission_energy*PA(n).*c_a.*x2.*conj(x2);
            y_b_sing = transmission_energy*PA(n)*1j*2*pi.*t.*c2.*x2.*conj(x2);
            y_d_sing = -transmission_energy*PA(n)*4*pi^2.*t.^2.*x2.*conj(x2);
            % Calculate the Fourier Transform of the signal
            X2 = fftshift(fft(x2));
            % Calculate effective pulse width (EPW)
            Eff_time = trapz(t, t.^2 .* abs(x2).^2);
            % Calculate effective bandwidth (EBW)
            Eff_bandw = trapz(dt, f.^2 .* abs(X2).^2)/len;
            TB = imag(sum(y2)*dt);
            cons_A_single(n,k) = -cons*real(trapz(t, y_a_sing)*RCS_sq_MT(n,k,2));
            cons_A(n,k) = int_A(n,k)+cons_A_single(n,k); % OCDM
            cons_B_single(n,k) = -cons*real(trapz(t, y_b_sing)*RCS_sq_MT(n,k,2));
            cons_B(n,k) = int_B(n,k)+cons_B_single(n,k);
            cons_D_single(n,k) = -cons*real(trapz(t, y_d_sing)*RCS_sq_MT(n,k,2));
            % cons_B = 0;
            cons_D(n,k) = int_D(n,k)+cons_D_single(n,k);
        end
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
    FIM(:,:,2) = aleph(:,:,2)*J*aleph(:,:,2).';
    %
    A_mat_single = cons_A_single.';
    A_single = diag(A_mat_single(:));
    B_mat_single = cons_B_single.';
    B_single = diag(B_mat_single(:));
    D_mat_single = cons_D_single.';
    D_single = diag(D_mat_single(:));
    % E_Im_mat = (cons_e*real(RCS).*Delay.*repmat(PA,1,numRe)).';
    G_single = [G_Re_single,G_Im_single];
    E_single = [E_Re_single,E_Im_single];
    J_single = [A_single,B_single,G_single;B_single,D_single,E_single;...
        G_single.', E_single.',F];
    FIM_single(:,:,2) = aleph(:,:,2)*J_single*aleph(:,:,2).';
    %%
    leng = 4+2*numTr*numRe;
    
    CRB_mat_single= inv(blkdiag(FIM_single(:,:,1),FIM_single(:,:,2)));
    CRB_loc_single(:,i_dis) = [CRB_mat_single(1,1),CRB_mat_single(2,2),...
        CRB_mat_single(1+leng,1+leng),CRB_mat_single(2+leng,2+leng)].';
    CRB_vel_single(:,i_dis) = [CRB_mat_single(3,3),CRB_mat_single(4,4),...
        CRB_mat_single(3+leng,3+leng),CRB_mat_single(4+leng,4+leng)].';
    %
    CRB_mat= inv(blkdiag(FIM(:,:,1),FIM(:,:,2)));
    CRB_loc(:,i_dis) = [CRB_mat(1,1),CRB_mat(2,2),...
        CRB_mat(1+leng,1+leng),CRB_mat(2+leng,2+leng)].';
    CRB_vel(:,i_dis) = [CRB_mat(3,3),CRB_mat(4,4),...
        CRB_mat(3+leng,3+leng),CRB_mat(4+leng,4+leng)].';
    %
    FIM_12 = aleph(:,:,1)*J_12*aleph(:,:,2).';
    FIM_21 = FIM_12.';
    FIM_accurate = [FIM(:,:,1),FIM_12;FIM_21,FIM(:,:,2)];
    CRB_mat_accurate = pinv(FIM_accurate);
    CRB_loc_accurate(:,i_dis) = [CRB_mat_accurate(1,1),CRB_mat_accurate(2,2),...
        CRB_mat_accurate(1+leng,1+leng),CRB_mat_accurate(2+leng,2+leng)].';
    CRB_vel_accurate(:,i_dis) = [CRB_mat_accurate(3,3),CRB_mat_accurate(4,4),...
        CRB_mat_accurate(3+leng,3+leng),CRB_mat_accurate(4+leng,4+leng)].';
end
toc
%%
figure
x_index = num*space;
% x_index = num*space;
subplot(2, 2, 1);
plot(x_index,CRB_loc_accurate(1,:),'ro-',x_index,CRB_loc(1,:),'bs-',...
    x_index,CRB_loc_single(1,:),'g>-',...
    'MarkerSize',8,'LineWidth',1)
% xlim([-77e3 83e3])
xlabel('x/m')
% xlabel('v_x (m/s)')
ylabel('CRLB x (m^2)')
ax=gca;
set(ax,'FontSize',13);
legend('CRLB','Approx CRLB','Single CRLB','FontSize',11)
grid on
subplot(2, 2, 2);
plot(x_index,CRB_loc_accurate(2,:),'ro-',x_index,CRB_loc(2,:),'bs-',...
    x_index,CRB_loc_single(2,:),'g>-',...
    'MarkerSize',8,'LineWidth',1)
% xlim([-77e3 83e3])
xlabel('x/m')
% xlabel('v_x (m/s)')
ylabel('CRLB y (m^2)')
ax=gca;
set(ax,'FontSize',13);
legend('CRLB','Approx CRLB','Single CRLB','FontSize',11)
grid on
subplot(2, 2, 3);
plot(x_index,CRB_vel_accurate(1,:),'ro-',x_index,CRB_vel(1,:),'bs-',...
    x_index,CRB_vel_single(1,:),'g>-',...
    'MarkerSize',8,'LineWidth',1)
% xlim([-77e3 83e3])
xlabel('x/m')
% xlabel('v_x (m/s)')
ylabel('CRLB v_x (m/s)^2')
ax=gca;
set(ax,'FontSize',13);
legend('CRLB','Approx CRLB','Single CRLB','FontSize',11)
grid on
subplot(2, 2, 4);
plot(x_index,CRB_vel_accurate(2,:),'ro-',x_index,CRB_vel(2,:),'bs-',...
    x_index,CRB_vel_single(2,:),'g>-',...
    'MarkerSize',8,'LineWidth',1)
% xlim([-77e3 83e3])
xlabel('x/m')
% xlabel('v_x (m/s)')
ylabel('CRLB v_y (m/s)^2')
ax=gca;
set(ax,'FontSize',13);
legend('CRLB','Approx CRLB','Single CRLB','FontSize',11)
grid on
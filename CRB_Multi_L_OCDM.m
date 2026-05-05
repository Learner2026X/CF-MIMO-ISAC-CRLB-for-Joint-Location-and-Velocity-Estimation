clearvars
% clc
close all
%%
M = 10;
N = 128;
T = 1e-3;
delta_f = 1/T;
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
num = -11:2:11;
% num = 1;
len_num = length(num);
% v_tar = [4,5].';
% load('p_tx_wide.mat','p_tx')
% load('p_rx_wide.mat','p_rx')
p_tx = radius*[cosd(azimuth_Tr(1:numTr));sind(azimuth_Tr(1:numTr))];
p_rx = radius*[cosd(azimuth_Re(1:numRe));sind(azimuth_Re(1:numRe))];
%%
p_tar = [0,0].';
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
rng('shuffle')
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
% load('RCS_multi.mat','RCS') 
% load('RCS1_multi.mat','RCS1') 
RCS = sqrt(1/2)*(ones(numTr,numRe)+1j*ones(numTr,numRe));
% RCS = sqrt(1/2)*(randn(numTr,numRe)+1j*randn(numTr,numRe));
RCS = RCS(1:numTr,1:numRe);
RCS_var = abs(RCS).^2;
% rng(1)
% RCS1 = 0.3*ones(numTr,numRe);
% RCS1 = sqrt(1/2)*(randn(numTr,numRe)+1j*randn(numTr,numRe));
RCS1 = 0.2*RCS;
RCS1 = RCS1(1:numTr,1:numRe);
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
space = 3.15e3;
for i_dis = 1:len_num
%     p_tar = [-3e3,0;-0e3,0].';
    p_tar = [-0e3,0e3;-0e3+num(i_dis)*space,0].';
    Delay = zeros(numTr,numRe,Num_Tar);
    Doppler = zeros(numTr,numRe,Num_Tar);
    for i_tar = 1:Num_Tar
        dis_tx_tar=sqrt(sum((p_tx-p_tar(:,i_tar)).^2));
        dis_rx_tar=sqrt(sum((p_rx-p_tar(:,i_tar)).^2));
        delay_tx_tar = dis_tx_tar/c;
        delay_rx_tar = dis_rx_tar/c;
        Delay(:,:,i_tar) = delay_tx_tar.' + delay_rx_tar;
    %     Delay = repmat(delay_tx_tar.',1, numRe)+repmat(delay_rx_tar,  numTr,1);
    end
    %% Jql
    PA = ones(numTr,1);
    const12 = 2*fs/noise_power;
    A_12 = zeros(numTr*numRe);
    G1_Re = zeros(numTr*numRe);
    G1_Im = zeros(numTr*numRe);
    G2_Re = zeros(numTr*numRe);
    G2_Im = zeros(numTr*numRe);
    F_12 = zeros(2*numTr*numRe);
    for n = 1:numTr
        for k = 1:numRe
            temp = (n-1)*numRe+k;
            x1 = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t-Delay(n,k,1)-(n-1)*T/N).^2/T^2);
            x2 = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t-Delay(n,k,2)-(n-1)*T/N).^2/T^2);
            c1 = -2*pi*(1j*N-1)*(t-Delay(n,k,1)-(n-1)*T/N)/T^2;
            c2 = 2*pi*(1j*N+1)*(t-Delay(n,k,2)-(n-1)*T/N)/T^2; % taking the conjugate
%             x1 = (2/T^2)^(1/4)*exp(-pi*(t-Delay(n,k,1)).^2/T^2).*...
%                 exp(1j*2*pi*n*delta_f*(t-Delay(n,k,1)));
%             x2 = (2/T^2)^(1/4)*exp(-pi*(t-Delay(n,k,2)).^2/T^2).*...
%                 exp(1j*2*pi*n*delta_f*(t-Delay(n,k,2)));
%             c1 = 2*pi*((t-Delay(n,k,1))/T^2-1j*n*delta_f);
%             c2 = 2*pi*((t-Delay(n,k,2))/T^2+1j*n*delta_f);% taking the conjugate
            y_a = transmission_energy*PA(n).*c1.*c2.*x1.*conj(x2);
            A_12(temp,temp) = const12*real(trapz(t, y_a)*RCS_product(n,k));
            y_g1Re = transmission_energy*PA(n).*c1.*x1.*conj(x2);
            G1_Re(temp,temp) = const12*real(trapz(t, y_g1Re)*RCS_MT(n,k,1));
            G1_Im(temp,temp) = const12*imag(trapz(t, y_g1Re)*RCS_MT(n,k,1));
%             int_G_Re(temp,temp) = -cons*real(conj(RCS_MT(n,k,2))*trapz(t, i_gRe));
%             int_G_Im(temp,temp) = cons*imag(conj(RCS_MT(n,k,2))*trapz(t, i_gRe));
            y_g2Re = transmission_energy*PA(n).*c2.*x1.*conj(x2);
            G2_Re(temp,temp) = const12*real(trapz(t, y_g2Re)*conj(RCS_MT(n,k,2)));
            G2_Im(temp,temp) = -const12*imag(trapz(t, y_g2Re)*conj(RCS_MT(n,k,2)));
            y_12 = transmission_energy*PA(n).*x1.*conj(x2);
            Y_12 = trapz(t, y_12);
            F_12(temp,temp) = const12*real(Y_12);
            F_12(temp,numTr*numRe+temp) = const12*imag(Y_12);
            F_12(numTr*numRe+temp,temp) = -const12*imag(Y_12);
            F_12(numTr*numRe+temp,numTr*numRe+temp) = const12*real(Y_12);
        end
    end
    J_12 = [A_12,G1_Re,G1_Im;[G2_Re,G2_Im].',F_12];
    %% Alpha
    aleph = zeros(2+2*numTr*numRe,3*numTr*numRe,Num_Tar); 
    for i_tar = 1:Num_Tar
        beta = zeros(2,numTr*numRe);
        eta = zeros(2,numTr*numRe);
        xi = zeros(2,numTr*numRe);
        for n = 1:numTr
            for k = 1:numRe
                temp = (n-1)*numRe+k;
                beta(:,temp) = 1/c*((p_tar(:,i_tar)-p_tx(:,n))/norm(p_tar(:,i_tar)-p_tx(:,n))+...
                    (p_tar(:,i_tar)-p_rx(:,k))/norm(p_tar(:,i_tar)-p_rx(:,k)));
            end
        end
        aleph(:,:,i_tar) = blkdiag(beta,eye(2*numTr*numRe));
    end
    %% Diagonal
    FIM = zeros(2+2*numTr*numRe,2+2*numTr*numRe,Num_Tar);
    FIM_single = zeros(2+2*numTr*numRe,2+2*numTr*numRe,Num_Tar);
    cons_f = 2*transmission_energy/noise_power*fs;
    cons = 2*fs/noise_power;
    cons_A = zeros(numTr,numRe);
    cons_A_single = zeros(numTr,numRe);
    G_Re = zeros(numTr*numRe);
    G_Im = zeros(numTr*numRe);
    G_Re_single = zeros(numTr*numRe);
    G_Im_single = zeros(numTr*numRe);
    int_A = zeros(numTr,numRe);
    int_G_Re = zeros(numTr*numRe);
    int_G_Im = zeros(numTr*numRe);
    %% target 1
    for n = 1:numTr
        for k = 1:numRe
            temp = (n-1)*numRe+k;
            x1 = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t-Delay(n,k,1)-(n-1)*T/N).^2/T^2);
            x2 = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t-Delay(n,k,2)-(n-1)*T/N).^2/T^2);
            c1 = -2*pi*(1j*N-1)*(t-Delay(n,k,1)-(n-1)*T/N)/T^2;
            c_a = 4*pi^2/T^4*(1-1j*N)^2*(t-Delay(n,k,1)-(n-1)*T/N).^2-2*pi/T^2*(1-1j*N);         
%             x1 = (2/T^2)^(1/4)*exp(-pi*(t-Delay(n,k,1)).^2/T^2).*...
%                 exp(1j*2*pi*n*delta_f*(t-Delay(n,k,1)));
%             x2 = (2/T^2)^(1/4)*exp(-pi*(t-Delay(n,k,2)).^2/T^2).*...
%                 exp(1j*2*pi*n*delta_f*(t-Delay(n,k,2)));
%             c1 = 2*pi*((t-Delay(n,k,1))/T^2-1j*n*delta_f);
%             c_a = 4*pi^2*((t-Delay(n,k,1))/T^2-1j*n*delta_f).^2-2*pi/T^2;
            y_gRe = transmission_energy*PA(n).*c1.*x1.*conj(x1);
            i_gRe = transmission_energy*PA(n).*c1.*x1.*conj(x2);
            G_Re_single(temp,temp) = -cons*real(conj(RCS_MT(n,k,1))*trapz(t, y_gRe));
            G_Im_single(temp,temp) = cons*imag(conj(RCS_MT(n,k,1))*trapz(t, y_gRe));
            int_G_Re(temp,temp) = -cons*real(conj(RCS_MT(n,k,2))*trapz(t, i_gRe));
            int_G_Im(temp,temp) = cons*imag(conj(RCS_MT(n,k,2))*trapz(t, i_gRe));
            G_Re(temp,temp) = int_G_Re(temp,temp)+G_Re_single(temp,temp);
            G_Im(temp,temp) = int_G_Im(temp,temp)+G_Im_single(temp,temp);
            
            y_a = transmission_energy*PA(n).*c_a.*x1.*conj(x2);
            int_A(n,k) = -cons*real(trapz(t, y_a)*RCS_product(n,k));
    
%             y1 = t.*conj(x1).*x1.*c1;
            y_a_sing = transmission_energy*PA(n).*c_a.*x1.*conj(x1);
            % Calculate the Fourier Transform of the signal
            X1 = fftshift(fft(x1));
            % Calculate effective pulse width (EPW)
            Eff_time = trapz(t, t.^2 .* abs(x1).^2);
            % Calculate effective bandwidth (EBW)
            Eff_bandw = trapz(dt, f.^2 .* abs(X1).^2)/len;
            cons_A_single(n,k) = -cons*real(trapz(t, y_a_sing)*RCS_sq_MT(n,k,1));
            cons_A(n,k) = int_A(n,k)+cons_A_single(n,k); % OCDM                
        end
    end
    A_mat = cons_A.';
    A = diag(A_mat(:));
    % E_Im_mat = (cons_e*real(RCS).*Delay.*repmat(PA,1,numRe)).';
    G = [G_Re,G_Im];
    F = diag(repmat(cons_f*PA,2*numRe,1));
    J = [A,G;G.',F];
    FIM(:,:,1) = aleph(:,:,1)*J*aleph(:,:,1).';
    %
    A_mat_single = cons_A_single.';
    A_single = diag(A_mat_single(:));
   
    G_single = [G_Re_single,G_Im_single];
    J_single = [A_single,G_single;G_single.',F];
    FIM_single(:,:,1) = aleph(:,:,1)*J_single*aleph(:,:,1).';
    %% target 2
    for n = 1:numTr
        for k = 1:numRe
            temp = (n-1)*numRe+k;
            x1 = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t-Delay(n,k,1)-(n-1)*T/N).^2/T^2);
            x2 = (2/T^2)^(1/4)*exp(pi*(1j*N-1)*(t-Delay(n,k,2)-(n-1)*T/N).^2/T^2);
            c2 = -2*pi*(1j*N-1)*(t-Delay(n,k,2)-(n-1)*T/N)/T^2;
            c_a = 4*pi^2/T^4*(1-1j*N)^2*(t-Delay(n,k,2)-(n-1)*T/N).^2-2*pi/T^2*(1-1j*N);
%             x1 = (2/T^2)^(1/4)*exp(-pi*(t-Delay(n,k,1)).^2/T^2).*...
%                 exp(1j*2*pi*n*delta_f*(t-Delay(n,k,1)));
%             x2 = (2/T^2)^(1/4)*exp(-pi*(t-Delay(n,k,2)).^2/T^2).*...
%                 exp(1j*2*pi*n*delta_f*(t-Delay(n,k,2)));
%             c2 = 2*pi*((t-Delay(n,k,2))/T^2-1j*n*delta_f);
%             c_a = 4*pi^2*((t-Delay(n,k,2))/T^2-1j*n*delta_f).^2-2*pi/T^2;
            y_gRe = transmission_energy*PA(n).*c2.*x2.*conj(x2);
            i_gRe = transmission_energy*PA(n).*c2.*x2.*conj(x1);
            G_Re_single(temp,temp) = -cons*real(conj(RCS_MT(n,k,2))*trapz(t, y_gRe));
            G_Im_single(temp,temp) = cons*imag(conj(RCS_MT(n,k,2))*trapz(t, y_gRe));
            int_G_Re(temp,temp) = -cons*real(conj(RCS_MT(n,k,1))*trapz(t, i_gRe));
            int_G_Im(temp,temp) = cons*imag(conj(RCS_MT(n,k,1))*trapz(t, i_gRe));
            G_Re(temp,temp) = int_G_Re(temp,temp)+G_Re_single(temp,temp);
            G_Im(temp,temp) = int_G_Im(temp,temp)+G_Im_single(temp,temp);
            
            y_a = transmission_energy*PA(n).*c_a.*x2.*conj(x1);
            int_A(n,k) = -cons*real(trapz(t, y_a)*conj(RCS_product(n,k)));
    
            y2 = t.*conj(x2).*x2.*c2;
            y_a_sing = transmission_energy*PA(n).*c_a.*x2.*conj(x2);
            cons_A_single(n,k) = -cons*real(trapz(t, y_a_sing)*RCS_sq_MT(n,k,2));
            cons_A(n,k) = int_A(n,k)+cons_A_single(n,k); % OCDM
            % Calculate the Fourier Transform of the signal
            X2 = fftshift(fft(x2));
            % Calculate effective pulse width (EPW)
            Eff_time = trapz(t, t.^2 .* abs(x2).^2);
            % Calculate effective bandwidth (EBW)
            Eff_bandw = trapz(dt, f.^2 .* abs(X2).^2)/len;
            TB = imag(sum(y2)*dt);
            
        end
    end
    A_mat = cons_A.';
    A = diag(A_mat(:));
    % E_Im_mat = (cons_e*real(RCS).*Delay.*repmat(PA,1,numRe)).';
    G = [G_Re,G_Im];
    F = diag(repmat(cons_f*PA,2*numRe,1));
    J = [A,G;G.',F];
    FIM(:,:,2) = aleph(:,:,2)*J*aleph(:,:,2).';
    %
    A_mat_single = cons_A_single.';
    A_single = diag(A_mat_single(:));
    % E_Im_mat = (cons_e*real(RCS).*Delay.*repmat(PA,1,numRe)).';
    G_single = [G_Re_single,G_Im_single];
    J_single = [A_single,G_single;G_single.',F];
    FIM_single(:,:,2) = aleph(:,:,2)*J_single*aleph(:,:,2).';
    %%
    leng = 2+2*numTr*numRe;
    CRB_mat_single= pinv(blkdiag(FIM_single(:,:,1),FIM_single(:,:,2)));
    CRB_loc_single(:,i_dis) = [CRB_mat_single(1,1),CRB_mat_single(2,2),...
        CRB_mat_single(1+leng,1+leng),CRB_mat_single(2+leng,2+leng)].';
    %
    CRB_mat= pinv(blkdiag(FIM(:,:,1),FIM(:,:,2)));
    CRB_loc(:,i_dis) = [CRB_mat(1,1),CRB_mat(2,2),...
        CRB_mat(1+leng,1+leng),CRB_mat(2+leng,2+leng)].';
    %
    FIM_12 = aleph(:,:,1)*J_12*aleph(:,:,2).';
    FIM_21 = FIM_12.';
    FIM_accurate = [FIM(:,:,1),FIM_12;FIM_21,FIM(:,:,2)];
    CRB_mat_accurate = pinv(FIM_accurate);
    CRB_loc_accurate(:,i_dis) = [CRB_mat_accurate(1,1),CRB_mat_accurate(2,2),...
        CRB_mat_accurate(1+leng,1+leng),CRB_mat_accurate(2+leng,2+leng)].';
end
CRB_OCDM = CRB_loc;
CRB_OCDM_accurate = CRB_loc_accurate;
CRB_OCDM_single = CRB_loc_single;
save('data\CRB_MT_L_OCDM.mat','CRB_OCDM')
save('data\CRB_MT_L_OCDM_accurate.mat','CRB_OCDM_accurate')
save('data\CRB_MT_L_OCDM_single.mat','CRB_OCDM_single')
toc
%%
x_index = num*space;
figure
subplot(2, 1, 1);
plot(x_index,CRB_loc_accurate(1,:),'ro-',x_index,CRB_loc(1,:),'bs-',...
    x_index,CRB_loc_single(1,:),'g>-',...
    'MarkerSize',7,'LineWidth',1)
% ylim([150 850])
xlabel('x/m')
ylabel('CRLB x/m^2')
ax=gca;
set(ax,'FontSize',13);
legend('CRLB','Approx CRLB','Single CRLB')
grid on
subplot(2, 1, 2);
plot(x_index,CRB_loc_accurate(2,:),'ro-',x_index,CRB_loc(2,:),'bs-',...
    x_index,CRB_loc_single(2,:),'g>-',...
    'MarkerSize',7,'LineWidth',1)
% xlim([-77e3 83e3])
xlabel('x/m')
ylabel('CRLB y/m^2')
ax=gca;
set(ax,'FontSize',13);
grid on
legend('CRLB','Approx CRLB','Single CRLB')
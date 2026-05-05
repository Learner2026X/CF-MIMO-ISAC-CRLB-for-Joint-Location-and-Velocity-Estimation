clearvars
rng(0)
numTr = 10;% Number of transmitters
numRe = 5;% Number of transmitters
% p_tx = [0,0;1500,1000;500,1000;300,500;2000,500].';
% p_rx = [0,1500;990,0;].';
% p_tx = 2000*rand(2,numTr);
% p_rx = 2000*rand(2,numRe);
load('p_tx.mat','p_tx')
load('p_rx.mat','p_rx')
gain_channel = abs(randn(numTr,1)+1j*randn(numTr,1)).^2;
p_tx = p_tx(:,1:numTr);
p_rx = p_rx(:,1:numRe);
p_tar = [400,500].';
v_tar = [4,5].';
c = 3e8;
f_c = 3e9; % carrier frequency
lambda = c/f_c;
transmission_power = 1; %w
T = 1e-2; % duration of transmitted waveform
SNR = 30;%dB
NSR = 1/10^(SNR/10);
Eff_bandw = 2.0451e+05;
Eff_time = 7.9577e-06;
TB = 8;
noise_power = transmission_power*T*NSR;
% fading_random = 0.5*randn(numTr,1)+0.5j*randn(numTr,1);%fading
% RCS_var = abs(randn(numTr,numRe)+1j*randn(numTr,numRe)).^2;
% RCS_var = ones(numTr,numRe);
load('RCS_var.mat','RCS_var')                                                                                                     
cons_a = 8*pi^2*transmission_power*T*Eff_bandw/noise_power; % OCDM
cons_b = 4*pi*transmission_power*T*TB/noise_power; 
% cons_b = 0;
cons_d = 8*pi^2*transmission_power*T*Eff_time/noise_power; 
%% intermediate parameter
beta = zeros(2,numTr,numRe);
eta = zeros(2,numTr,numRe);
xi = zeros(2,numTr,numRe);
weight_P = zeros(3,numTr,numRe);
weight_V = zeros(4,numTr,numRe);
weight_Y = zeros(3,numTr,numRe);
i_tar = 1;
for n = 1:numTr
    for k = 1:numRe
        beta(:,n,k) = 1/c*((p_tar(:,i_tar)-p_tx(:,n))/norm(p_tar(:,i_tar)-p_tx(:,n))+...
            (p_tar(:,i_tar)-p_rx(:,k))/norm(p_tar(:,i_tar)-p_rx(:,k)));
        eta(1,n,k) = 1/lambda*((p_tar(1,i_tar)-p_tx(1,n))*v_tar(:,i_tar).'*...
           (p_tar(:,i_tar)-p_tx(:,n))/norm(p_tar(:,i_tar)-p_tx(:,n))^3-...
           v_tar(1,i_tar)/norm(p_tar(:,i_tar)-p_tx(:,n))+...
           (p_tar(1,i_tar)-p_rx(1,k))*v_tar(:,i_tar).'*...
           (p_tar(:,i_tar)-p_rx(:,k))/norm(p_tar(:,i_tar)-p_rx(:,k))^3-...
           v_tar(1,i_tar)/norm(p_tar(:,i_tar)-p_rx(:,k)));
        eta(2,n,k) = 1/lambda*((p_tar(2,i_tar)-p_tx(2,n))*v_tar(:,i_tar).'*...
           (p_tar(:,i_tar)-p_tx(:,n))/norm(p_tar(:,i_tar)-p_tx(:,n))^3-...
           v_tar(2,i_tar)/norm(p_tar(:,i_tar)-p_tx(:,n))+...
           (p_tar(2,i_tar)-p_rx(2,k))*v_tar(:,i_tar).'*...
           (p_tar(:,i_tar)-p_rx(:,k))/norm(p_tar(:,i_tar)-p_rx(:,k))^3-...
           v_tar(2,i_tar)/norm(p_tar(:,i_tar)-p_rx(:,k)));
        xi(1,n,k) = -1/lambda*((p_tar(1,i_tar)-p_tx(1,n))/...
           norm(p_tar(:,i_tar)-p_tx(:,n))+(p_tar(1,i_tar)-p_rx(1,k))/...
           norm(p_tar(:,i_tar)-p_rx(:,k)));
        xi(2,n,k) = -1/lambda*((p_tar(2,i_tar)-p_tx(2,n))/...
           norm(p_tar(:,i_tar)-p_tx(:,n))+(p_tar(2,i_tar)-p_rx(2,k))/...
           norm(p_tar(:,i_tar)-p_rx(:,k)));
        weight_P(:,n,k) = [beta(1,n,k)^2,beta(1,n,k)*beta(2,n,k),beta(2,n,k)^2].'...
            *RCS_var(n,k)*cons_a+[2*beta(1,n,k)*eta(1,n,k),...
            beta(1,n,k)*eta(2,n,k)+eta(1,n,k)*beta(2,n,k),2*beta(2,n,k)*eta(2,n,k)].'...
            *RCS_var(n,k)*cons_b+[eta(1,n,k)^2,eta(1,n,k)*eta(2,n,k),eta(2,n,k)^2].'...
            *RCS_var(n,k)*cons_d;
        weight_V(:,n,k) = [beta(1,n,k)*xi(1,n,k),beta(2,n,k)*xi(1,n,k)...
            beta(1,n,k)*xi(2,n,k),beta(2,n,k)*xi(2,n,k)].'*RCS_var(n,k)*cons_b...
            +[eta(1,n,k)*xi(1,n,k),eta(2,n,k)*xi(1,n,k),eta(1,n,k)*xi(2,n,k),...
            eta(2,n,k)*xi(2,n,k)].'*RCS_var(n,k)*cons_d;
        weight_Y(:,n,k) = [xi(1,n,k)^2,xi(1,n,k)*xi(2,n,k),xi(2,n,k)^2].'...
            *RCS_var(n,k)*cons_d;
    end
end
%% CRB matrix
PA = ones(numTr,1)/numTr;
W_P = zeros(numTr,4);
W_V = zeros(numTr,4);
W_Y = zeros(numTr,4);                                                                                                             
P_mat = zeros(2,2);
V_mat = zeros(2,2);
Y_mat = zeros(2,2);
for iter = 1:3
    W_P(:,iter) = sum(squeeze(weight_P(iter,:,:)),2);
    W_V(:,iter) = sum(squeeze(weight_V(iter,:,:)),2);
    W_Y(:,iter) = sum(squeeze(weight_Y(iter,:,:)),2);
end
W_P(:,4) = W_P(:,3);
W_P(:,3) = W_P(:,2);
W_Y(:,4) = W_Y(:,3);
W_Y(:,3) = W_Y(:,2);
W_V(:,4) = sum(squeeze(weight_V(4,:,:)),2);
W_y0 = W_Y(:,1)*W_Y(:,4).'-W_Y(:,3)*W_Y(:,2).';
W_y1 = W_V(:,3)*W_Y(:,1).'-W_V(:,1)*W_Y(:,3).';
W_y2 = W_V(:,4)*W_Y(:,1).'-W_V(:,2)*W_Y(:,3).';
W_y3 = W_V(:,1)*W_Y(:,4).'-W_V(:,3)*W_Y(:,2).';
W_y4 = W_V(:,2)*W_Y(:,4).'-W_V(:,4)*W_Y(:,2).';
W_p0 = W_P(:,1)*W_P(:,4).'-W_P(:,3)*W_P(:,2).';
W_p1 = W_V(:,2)*W_P(:,1).'-W_V(:,1)*W_P(:,3).';
W_p2 = W_V(:,4)*W_P(:,1).'-W_V(:,3)*W_P(:,3).';
W_p3 = W_V(:,1)*W_P(:,4).'-W_V(:,2)*W_P(:,2).';
W_p4 = W_V(:,3)*W_P(:,4).'-W_V(:,4)*W_P(:,2).';
%%
tic
rng('shuffle')
large_num = 200000;
lip_constant = zeros(1,large_num);
rho_min = 0.01*ones(numTr,1);
rho_max = 0.3*ones(numTr,1);
mat_projection = eye(numTr)-ones(numTr,1)*ones(numTr,1).'/numTr;
% rho(:,1) = gain_channel/sum(gain_channel);
delta_l = 120;
delta_v = 0.02;
penalty_factor = 1e6;
for i = 1:large_num
    rho1 = rand(numTr,1);
    rho1= rho1/sum(rho1);
    rho2 = rand(numTr,1);
    rho2= rho2/sum(rho2);
    Diff_min1 = rho_min-rho1;
    Diff_max1 = rho1-rho_max;
    Deriv_min1 = -2*max([Diff_min1,zeros(numTr,1)].').';
    Deriv_max1 = 2*max([Diff_max1,zeros(numTr,1)].').';
    Diff_min2 = rho_min-rho2;
    Diff_max2 = rho2-rho_max;
    Deriv_min2 = -2*max([Diff_min2,zeros(numTr,1)].').';
    Deriv_max2 = 2*max([Diff_max2,zeros(numTr,1)].').';
    grad1 = -gain_channel + penalty_factor*...
        (Derivative_l_square(W_P,W_V,W_y0,W_y1,W_y2,W_y3,W_y4,rho1,delta_l)+...
        Derivative_v_square(W_Y,W_V,W_p0,W_p1,W_p2,W_p3,W_p4,rho1,delta_v)+...
        Deriv_min1+Deriv_max1);
    grad2 = -gain_channel + penalty_factor*...
        (Derivative_l_square(W_P,W_V,W_y0,W_y1,W_y2,W_y3,W_y4,rho2,delta_l)+...
        Derivative_v_square(W_Y,W_V,W_p0,W_p1,W_p2,W_p3,W_p4,rho2,delta_v)+...
        Deriv_min2+Deriv_max2);
%     grad1 = Derivative_l_square(W_P,W_V,W_y0,W_y1,W_y2,W_y3,W_y4,rho1,delta_l);
%     grad2 = Derivative_l_square(W_P,W_V,W_y0,W_y1,W_y2,W_y3,W_y4,rho2,delta_l);
    lip_constant(i) = norm(grad2-grad1)/norm(rho2-rho1);
end
max_lip_constant = max(lip_constant);
toc
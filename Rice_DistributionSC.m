% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

%% Rice Selection Combining channels simulation
clear; close all; clc;

addpath('Functions/')                           % It adds the path where the functions are located

N = 2;                                          % Number of receive branches
Nsim = 1e6;                                     % Number of random variable realizations
R = 1.7;                                        % Threshold rate for Outage metrics (R)
W_th = 2^(R) - 1;                               % Minimum power required to receive the packet sent at a transmission rate of R

% The axis of average SNRs on which it operates is defined
Omega_dB = -0:0.5:25;                           % Average SNR (in dB)
Omega = 10.^(Omega_dB/10);                      % Average SNR

% Parameter for Figures
markers_ind = ceil(length(Omega_dB)/100*4);

%%  Outage probability
%   oP  --> Simulation of outage probability from averaged realizations
K = 10;

% Nsim realizations other than the channel are generated
h = sqrt(K/(K + 1)).*ones(N, Nsim) + 1/sqrt(K + 1).*(randn(N, Nsim) + 1i.*randn(N, Nsim))/sqrt(2);
gSC =max(abs(h),[],1);

SNR = ((abs(gSC).^2)'*Omega)';                  % Realizations of instantaneous SNR

% (a) Simulated outage probabiilty -->  Pr(log2(1 + SNR < Rth));
% the number of times the channel capacity is less than the threshold rate is counted
oP = sum(log2(1 + SNR) < R, 2)/Nsim;

% (b) Theoretical outage probability Pr(Rth < log2(1 + SNR)) --> Pr(SNR < 2^Rth - 1)) = CDF(2^Rth - 1);
% %   oPtheo --> Theoretical expression of outage probability from the CDF
oPtheo = cdfSNRRiceThSC(W_th./Omega, K, 1,N);


%% Theoretical linear approximation for several SNR values
Omega0dB_var = [10 20];
Omega0_var = 10.^(Omega0dB_var/10);

linearAppox_log_Theo_var = zeros(length(Omega0dB_var), length(Omega_dB));
oPtheo_Om0_var = zeros(1, length(Omega0dB_var));
slope_temp = zeros(1, length(Omega0dB_var));
legendInfo{1} = 'Theoretical';
legendInfo{2} = 'Simulated';
for i = 1:length(Omega0dB_var)
    [slope_temp(i), oPtheo_Om0_var(i)] = Theoretical_Slope(Omega0_var(i), W_th, 'RiceSC', [K,N]);
    linearAppox_log_Theo_var(i,:) = oPtheo_Om0_var(i)*(Omega0_var(i)./Omega).^(slope_temp(i));
    legendInfo{i+2} = ['$\Omega_0$ = ' num2str(Omega0dB_var(i)),' dB'];
end

% %   oPasym --> approximated OP for high SNR (asymptotic)
oPasym = acdfSNRRiceTh(W_th./Omega, K, 1).^N;
legendInfo{end+1} = '$\Omega_0$ $\rightarrow$ $\infty$';

f_OP_Rice = figure;
set(f_OP_Rice, 'Position',  [1260 360, 560, 420])
set(f_OP_Rice, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_OP_Rice, 'defaultLegendInterpreter','latex');
set(f_OP_Rice, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_OP_Rice, 'defaultLineLineWidth',1);
set(f_OP_Rice, 'color','w');
semilogy(Omega_dB, oPtheo, 'b', 'LineWidth',2)
hold on
semilogy(Omega_dB, oP, 'o','MarkerFaceColor', 'g', 'MarkerEdgeColor', 'b', 'MarkerIndices',1:markers_ind:length(Omega_dB))
semilogy(Omega_dB,linearAppox_log_Theo_var,'--','LineWidth', 2)
semilogy(Omega_dB, oPasym,':k', 'LineWidth',2)
plot(Omega0dB_var, oPtheo_Om0_var,'xk','LineWidth', 2 ,MarkerSize = 10)
xline(Omega0dB_var,'k--', 'LineWidth',1)
grid on;
legend(legendInfo,'Location','northeast');
xlabel('$\Omega$ (dB)'); 
ylabel('OP');
axis([-inf inf 1e-10 2])
hold off;





%% Variation of parameters a and b (intercept and slope)
K = [0.05 1 2 5 10 15];

Omega_dB = -10:0.1:40;                      % Average SNR (in dB)
Omega = 10.^(Omega_dB/10);                  % Average SNR
 
% Parameter for Figures
markers_ind = ceil(length(Omega_dB)/100*4);

m_var = zeros(length(K), length(Omega));
a_var = zeros(length(K), length(Omega));
for i = 1:length(K)
    [m_var(i,:), a_var(i,:)] = Theoretical_Slope(Omega, W_th, 'Rice', K(i));
end


%% Numerical and theoretical log10(CDF)'s derivative comparison
oPtheo = zeros(length(K), length(Omega));
G = zeros(length(K), length(Omega));
Gprime = zeros(length(K), length(Omega) - 1);
delta_num = zeros(length(K), length(Omega) - 1);
m_var = zeros(length(K), length(Omega));
delta_theo = zeros(length(K), length(Omega));
a_var = zeros(length(K), length(Omega));
for k = 1:length(K)
    % oPtheo --> Theoretical expression of outage probability from the CDF
    oPtheo(k,:) = cdfSNRRiceThSC(W_th./Omega, K(k), 1,N);

    % Numerical and theoretical log10(CDF)'s derivative
    G(k,:) = log10(oPtheo(k,:));                                                                        % log10 of the CDF

    Gprime(k,:) = -diff(G(k,:))./diff(Omega_dB);                                                        % Numerical derivative
    delta_num(k,:) = 10*Gprime(k,:);

    [m_var(k,:), a_var(k,:)] = Theoretical_Slope(Omega, W_th, 'RiceSC', [K(k),N]);                      % Theoretical derivative
    delta_theo(k,:) = m_var(k,:);
end



h1 = zeros(1,length(K));
legendInfo_slope = cell(1, length(K));

f_ODO_Rice = figure;
set(f_ODO_Rice, 'Position',  [40 360, 560, 420])
set(f_ODO_Rice, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_ODO_Rice, 'defaultLegendInterpreter','latex');
set(f_ODO_Rice, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_ODO_Rice,'defaultLineLineWidth',1.5);
set(f_ODO_Rice,'color','w');
hold on
for i = 1:length(K)
    h1(i) = plot(Omega_dB, delta_theo(i, :));
end

for i = 1:length(K)
    plot(Omega_dB(2:end), delta_num(i, :),'s--','Color',get(h1(i), 'color'),'MarkerIndices', 2*i:markers_ind:length(Omega_dB)-1)
    legendInfo_slope{i} = ['K = ' num2str(K(i))];
end
yline(N,':k','LineWidth',2);
xlabel('$\Omega_0$ (dB)');
ylabel('$\delta_{\rm Rice}$');
legend(legendInfo_slope,'Location','northeast','AutoUpdate','off');
grid on
hold off




h1 = zeros(1,length(K));
c_theo = 10./delta_theo;
c_num = 10./delta_num;

f_c = figure;
set(f_c, 'Position',  [640, 360, 560, 420])
set(f_c, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_c, 'defaultLegendInterpreter','latex');
set(f_c, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_c, 'defaultLineLineWidth',1.5);
set(f_c, 'color','w');
hold on
for i = 1:length(K)
    h1(i) = plot(Omega_dB, c_theo(i,:));
end
for i = 1:length(K)
    plot(Omega_dB(2:end), c_num(i,:), 's--','Color',get(h1(i), 'color'),'MarkerIndices', 2*i:markers_ind:length(Omega_dB)-1)    
end
yline(10/N,':k','LineWidth',2);
grid on
xlabel('$\Omega_0$ (dB)');
ylabel('$c_{\rm Rice}$ (dB)');
legend(legendInfo_slope)
hold off;
axis([0 Omega_dB(end) 0 20])
% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

%% Two-Wave with Difusse Power (TWDP) channels simulation
clear; close all; clc;

addpath('Functions/')                                       % It adds the path where the functions are located

Nsim = 1e6;                                                 % Number of random variable realizations
R = 1.7;                                                  % Threshold rate for Outage metrics (R)
W_th = 2^(R) - 1;                                           % Minimum power required to receive the packet sent at a transmission rate of R

%% Numerical and theoretical log10(CDF)'s derivative comparison
% Axis of average SNRs on which it operates
Omega_dB = -10:0.1:40;                                      % Average SNR (in dB)
Omega = 10.^(Omega_dB/10);                                 % Average SNR

% Parameter for Figures
markers_ind = ceil(length(Omega_dB)/100*4);

%% TWDP parameters
K = 12;
Delta_var = 0:0.1:1;

type = 0;                                                   % type: 0 - amplitude  / 1- power

%% K constant and Delta variable
legendInfo_dVar = cell(1, length(Delta_var));
Param_dVar = zeros(length(Delta_var), 2);
oPtheo_dVar = zeros(length(Delta_var), length(Omega_dB));
G_dVar = zeros(length(Delta_var), length(Omega_dB));
Gprime_dVar = zeros(length(Delta_var), length(Omega_dB) - 1);
delta_num_dVar = zeros(length(Delta_var), length(Omega_dB) - 1);
GprimeTheo_dVar = zeros(length(Delta_var), length(Omega_dB));
delta_theo_dVar = zeros(length(Delta_var), length(Omega_dB));
oPtheo_Om0_var_dVar = zeros(length(Delta_var), length(Omega_dB));
for d = 1:length(Delta_var)
    % oPtheo --> Theoretical expression of outage probability from the CDF
    Param_dVar(d,:) = [K, Delta_var(d)];
    oPtheo_dVar(d,:) = cdfSNRTWDPTh(Param_dVar(d,:), 1, W_th./Omega);

    %% Numerical and theoretical log10(CDF)'s derivative
    G_dVar(d,:) = log10(oPtheo_dVar(d,:));                                                                      % log10 of the CDF

    Gprime_dVar(d,:) = -diff(G_dVar(d,:))./diff(Omega_dB);                                                      % Numerical derivative
    delta_num_dVar(d,:) = 10*Gprime_dVar(d,:);

    [GprimeTheo_dVar(d,:), oPtheo_Om0_var_dVar(d,:)] = Theoretical_Slope(1, W_th./Omega, 'TWDP', Param_dVar(d,:));
    delta_theo_dVar(d,:) = GprimeTheo_dVar(d,:);

    legendInfo_dVar{d} = ['$\Delta = $ ',num2str(Delta_var(d))];        
end

h1 = zeros(1,length(Delta_var));
f_ODO_TWDP = figure;
set(f_ODO_TWDP, 'Position',  [520, 360, 560, 420])
set(f_ODO_TWDP, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_ODO_TWDP, 'defaultLegendInterpreter','latex');
set(f_ODO_TWDP, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_ODO_TWDP, 'defaultLineLineWidth',1);
set(f_ODO_TWDP, 'color','w');
hold on
for d = 1:length(Delta_var)
    h1(d) = plot(Omega_dB, delta_theo_dVar(d,:));
end
for d = 1:length(Delta_var)
    plot(Omega_dB(2:end), delta_num_dVar(d,:), 's--','Color',get(h1(d), 'color'),'MarkerIndices', 2*d:markers_ind:length(Omega_dB)-1)    
end
yline(1,':k','LineWidth',2);
grid on
xlabel('$\Omega_0$ (dB)');
ylabel('$\delta_{\rm TWDP}$');
legend(legendInfo_dVar)
hold off;


c_theo = 10./delta_theo_dVar;
c_num = 10./delta_num_dVar;

f_c = figure;
set(f_c, 'Position',  [1100, 360, 560, 420])
set(f_c, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_c, 'defaultLegendInterpreter','latex');
set(f_c, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_c, 'defaultLineLineWidth',1);
set(f_c, 'color','w');
hold on
for d = 1:length(Delta_var)
    h1(d) = plot(Omega_dB, c_theo(d,:));
end
for d = 1:length(Delta_var)
    plot(Omega_dB(2:end), c_num(d,:), 's--','Color',get(h1(d), 'color'),'MarkerIndices', 2*d:markers_ind:length(Omega_dB)-1)    
end
yline(10,':k','LineWidth',2);
grid on
xlabel('$\Omega_0$ (dB)');
ylabel('$c_{\rm TWDP}$ (dB)');
legend(legendInfo_dVar)
hold off;
axis([0 Omega_dB(end) 2 20])


disp('Solid line --> Theoretical derivative - Dashed line with markers --> Numerical derivative')


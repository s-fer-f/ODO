% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336


%% Cascaded Rayleigh channels simulation
clear; close all; clc;

addpath('Functions/')                                       % It adds the path where the functions are located

Nsim = 1e6;                                                 % Number of random variable realizations
R = 1.7;                                                  % Threshold rate for Outage metrics (R)
W_th = 2^(R) - 1;                                           % Minimum power required to receive the packet sent at a transmission rate of R

% The axis of average SNRs on which it operates is defined
Omega_dB = -10:0.5:40;                                        % Average SNR (in dB)
Omega = 10.^(Omega_dB/10);                                 % Average SNR

% Parameter for Figures
markers_ind = ceil(length(Omega_dB)/100*4);


%% Cascaded Rayleigh
h1 = (randn(1, Nsim) + 1i.*randn(1, Nsim))/sqrt(2);         % Nsim realizations other than the channel are generated.
h2 = (randn(1, Nsim) + 1i.*randn(1, Nsim))/sqrt(2);         % Nsim realizations other than the channel are generated.
h_casc = h1.*h2;

SNR_casc = ((abs(h_casc).^2)'*Omega)';                             % Realizations of instantaneous SNR

%%  Outage probability
% (a) Simulated outage probabiilty -->  Pr(log2(1 + SNR < Rth));
% the number of times the channel capacity is less than the threshold rate is counted
%   oP  --> Simulation of outage probability from averaged realizations
oP_casc = sum(log2(1 + SNR_casc) < R, 2)/Nsim;

% (b) Theoretical outage probability Pr(Rth < log2(1 + SNR)) --> Pr(SNR < 2^Rth - 1)) = CDF(2^Rth - 1);
% %   oPtheo --> Theoretical expression of outage probability from the CDF
oPtheo_casc = cdfSNRProductTh(W_th./Omega, 1);


%% Rayleigh
h_Ray = (randn(1, Nsim) + 1i.*randn(1, Nsim))/sqrt(2);         % Nsim realizations other than the channel are generated.

SNR_Ray = ((abs(h_Ray).^2)'*Omega)';                             % Realizations of instantaneous SNR

%%  Outage probability
% (a) Simulated outage probabiilty -->  Pr(log2(1 + SNR < Rth));
% the number of times the channel capacity is less than the threshold rate is counted
%   oP  --> Simulation of outage probability from averaged realizations
oP_Ray = sum(log2(1 + SNR_Ray) < R, 2)/Nsim;

% (b) Theoretical outage probability Pr(Rth < log2(1 + SNR)) --> Pr(SNR < 2^Rth - 1)) = CDF(2^Rth - 1);
% %   oPtheo --> Theoretical expression of outage probability from the CDF
oPtheo_Ray = cdfSNRRayleighTh(W_th./Omega, 1);



%% Theoretical linear approximation for several SNR values
Omega0dB_var = 5:5:30;
Omega0_var = 10.^(Omega0dB_var/10);

linearAppox_log_Theo_var_Cascade = zeros(length(Omega0dB_var), length(Omega_dB));
linearAppox_log_Theo_var_Ray = zeros(length(Omega0dB_var), length(Omega_dB));
oPtheo_Om0_var_Cascade = zeros(1, length(Omega0dB_var));
oPtheo_Om0_var_Ray = zeros(1, length(Omega0dB_var));
slope_temp_cascade = zeros(1, length(Omega0dB_var));
slope_temp_Ray = zeros(1, length(Omega0dB_var));
legendInfo{1} = 'Simulated';
legendInfo{2} = 'Theoretical';
for i = 1:length(Omega0dB_var)
    %% Cascaded Rayleigh
    [slope_temp_cascade(i), oPtheo_Om0_var_Cascade(i)] = Theoretical_Slope(Omega0_var(i), W_th, 'Cascaded');
    linearAppox_log_Theo_var_Cascade(i,:) = oPtheo_Om0_var_Cascade(i)*(Omega0_var(i)./Omega).^(slope_temp_cascade(i));
    
    %% Rayleigh
    [slope_temp_Ray(i), oPtheo_Om0_var_Ray(i)] = Theoretical_Slope(Omega0_var(i), W_th, 'Rayleigh');
    linearAppox_log_Theo_var_Ray(i,:) = oPtheo_Om0_var_Ray(i)*(Omega0_var(i)./Omega).^(slope_temp_Ray(i));

    legendInfo{i+2} = ['$\Omega_0$ = ' num2str(Omega0dB_var(i)),' dB'];
end

%% Cascaded Rayleigh
% %   oPasym --> approximated OP for high SNR (asymptotic)
% DOES NOT EXIST FOR THE CASCADED RAYLEIGH DISTRIBUTION

%% Rayleigh
% %   oPasym --> approximated OP for high SNR (asymptotic)
oPasym_Ray = acdfSNRRayleighTh(W_th./Omega, 1);

legendInfo{end+1} = '$\Omega_0$ $\rightarrow$ $\infty$';

f_OP_Rayleigh = figure;
set(f_OP_Rayleigh, 'Position',  [1160 360, 560, 420])
set(f_OP_Rayleigh, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_OP_Rayleigh, 'defaultLegendInterpreter','latex');
set(f_OP_Rayleigh, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_OP_Rayleigh, 'defaultLineLineWidth',1);
set(f_OP_Rayleigh, 'color','w');
semilogy(Omega_dB, oP_casc, 'o','MarkerFaceColor', 'g', 'MarkerEdgeColor', 'c', 'MarkerIndices',1:markers_ind:length(Omega_dB))
hold on
semilogy(Omega_dB, oPtheo_casc, 'b','LineWidth',2)
semilogy(Omega_dB, linearAppox_log_Theo_var_Cascade,'-.','LineWidth', 2)

semilogy(Omega_dB, oPasym_Ray,'k:','LineWidth',2)
semilogy(Omega_dB, oP_Ray, 'o','MarkerFaceColor', 'm', 'MarkerEdgeColor', 'y', 'MarkerIndices',3:markers_ind:length(Omega_dB))
semilogy(Omega_dB, oPtheo_Ray, 'b', 'LineWidth',2)
semilogy(Omega_dB, linearAppox_log_Theo_var_Ray,'--','LineWidth', 2)

plot(Omega0dB_var, oPtheo_Om0_var_Cascade,'xk','LineWidth', 2 ,MarkerSize = 10)
plot(Omega0dB_var, oPtheo_Om0_var_Ray,'xk','LineWidth', 2 ,MarkerSize = 10)
xline(Omega0dB_var,'k--', 'LineWidth',1)
grid on;
legend(legendInfo,'Location','southwest');
xlabel('$\Omega$ (dB)'); ylabel('OP');
axis([-inf inf 1e-4 2])
hold off;










%% Numerical and theoretical log10(CDF)'s derivative comparison
% The axis of average SNRs on which it operates is defined
Omega_dB = -10:0.1:40;                                        % Average SNR (in dB)
Omega = 10.^(Omega_dB/10);                                 % Average SNR

% Parameter for Figures
markers_ind = ceil(length(Omega_dB)/100*4);


%% Cascaded Rayleigh
% oPtheo --> Theoretical expression of outage probability from the CDF
oPtheo_casc = cdfSNRProductTh(W_th./Omega, 1);

%% Rayleigh
% oPtheo --> Theoretical expression of outage probability from the CDF
oPtheo_Ray = cdfSNRRayleighTh(W_th./Omega, 1);


%% Numerical and theoretical log10(CDF)'s derivative
%% Cascaded Rayleigh
G_casc = log10(oPtheo_casc);                                                                  % log10 of the CDF

Gprime_casc = -diff(G_casc)./diff(Omega_dB);                                                  % Numerical derivative
delta_num_casc = 10*Gprime_casc;

[GprimeTheo_casc, a_var_casc] = Theoretical_Slope(Omega, W_th, 'Cascaded');
delta_theo_casc = GprimeTheo_casc;

%% Rayleigh
G_Ray = log10(oPtheo_Ray);                                                                  % log10 of the CDF

Gprime_Ray = -diff(G_Ray)./diff(Omega_dB);                                                  % Numerical derivative
delta_num_Ray = 10*Gprime_Ray;

[GprimeTheo_Ray, a_var_Ray] = Theoretical_Slope(Omega, W_th, 'Rayleigh');
delta_theo_Ray = GprimeTheo_Ray;


f_derivative_logCDF = figure;
set(f_derivative_logCDF, 'Position',  [460, 360, 560, 420])
set(f_derivative_logCDF, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_derivative_logCDF, 'defaultLegendInterpreter','latex');
set(f_derivative_logCDF, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_derivative_logCDF, 'defaultLineLineWidth',1);
set(f_derivative_logCDF, 'color','w');
plot(Omega_dB, delta_theo_Ray, 'b')
hold on;
plot(Omega_dB, delta_theo_casc, 'Color',[255/255 140/255 0/255])
plot(Omega_dB(2:end), delta_num_Ray,'b--s','MarkerIndices', 1:markers_ind:length(Omega_dB))
plot(Omega_dB(2:end), delta_num_casc,'--s','Color',[255/255 140/255 0/255],'MarkerIndices', 1:markers_ind:length(Omega_dB))
grid on
xlabel('$\Omega_0$ (dB)');
ylabel('$\delta_{\rm cas}$');
legend('Cascaded Rayleigh','Rayleigh','Location','southeast')
hold off;



disp('Solid line --> Theoretical derivative - Dashed line with markers --> Numerical derivative')
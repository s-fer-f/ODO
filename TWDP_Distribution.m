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
Omega_dB = -10:0.5:40;                                      % Average SNR (in dB)
Omega = 10.^(Omega_dB/10);                                 % Average SNR

%% Specific SNR values
Omega0dB = 5:5:30;
Omega0 = 10.^(Omega0dB/10);

% Parameter for Figures
markers_ind = ceil(length(Omega_dB)/100*4);

%% TWDP parameters
K = 12;
Delta_var = 0:0.1:1;

type = 0;                                                   % type: 0 - amplitude  / 1- power

%% K constant and Delta variable
legendInfo_dVar = cell(1, length(Delta_var));
oP_dVar = zeros(length(Delta_var), length(Omega_dB));
Param_dVar = zeros(length(Delta_var), 2);
oPtheo_dVar = zeros(length(Delta_var), length(Omega_dB));
G_dVar = zeros(length(Delta_var), length(Omega_dB));
Gprime_dVar = zeros(length(Delta_var), length(Omega_dB) - 1);
delta_num_dVar = zeros(length(Delta_var), length(Omega_dB) - 1);
GprimeTheo_dVar = zeros(length(Delta_var), length(Omega_dB));
delta_theo_dVar = zeros(length(Delta_var), length(Omega_dB));
oPasym_dVar = zeros(length(Delta_var), length(Omega_dB));
oPtheo_Om0_var_dVar = zeros(length(Delta_var), length(Omega_dB));
slope_temp_dVar = zeros(length(Delta_var), length(Omega0dB));
linearApprox_log_Theo_var_dVar = zeros(length(Delta_var), length(Omega0dB), length(Omega_dB));
for d = 1:length(Delta_var)
    h_dVar = genTWDPsim(K, Delta_var(d), Nsim, type).';                 % Nsim realizations other than the channel are generated.

    SNR_dVar = ((abs(h_dVar).^2)'*Omega)';                             % Realizations of instantaneous SNR

    %%  Outage probability
    % (a) Simulated outage probabiilty -->  Pr(log2(1 + SNR < Rth));
    % the number of times the channel capacity is less than the threshold rate is counted
    %   oP  --> Simulation of outage probability from averaged realizations
    oP_dVar(d,:) = sum(log2(1 + SNR_dVar) < R, 2)/Nsim;

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

    % %   oPasym --> approximated OP for high SNR (asymptotic)
    oPasym_dVar(d,:) = acdfSNR_TWDP_Th(W_th, K, Delta_var(d), Omega);
    % For high SNR, the slope of the CDF can be approximated by the scaled CDF of the Rice --> "the tail can be lower-bounded through a scaled Rician tail"

    %% Theoretical linear approximation for several SNR values
    for i = 1:length(Omega0)
        [slope_temp_dVar(d,i), oPtheo_Om0_var_dVar(d,i)] = Theoretical_Slope(Omega0(i), W_th, 'TWDP', Param_dVar(d,:));

        linearApprox_log_Theo_var_dVar(d,i,:) = oPtheo_Om0_var_dVar(d,i)*(Omega0(i)./Omega).^(slope_temp_dVar(d,i));
    end
end






h1 = zeros(1,length(Delta_var));

f_ODO_TWDP = figure;
set(f_ODO_TWDP, 'Position',  [40, 360, 560, 420])
set(f_ODO_TWDP, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_ODO_TWDP, 'defaultLegendInterpreter','latex');
set(f_ODO_TWDP, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_ODO_TWDP, 'defaultLineLineWidth',1.5);
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







h1 = zeros(1,length(Delta_var));
c_theo = 10./delta_theo_dVar;
c_num = 10./delta_num_dVar;

f_c = figure;
set(f_c, 'Position',  [640, 360, 560, 420])
set(f_c, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_c, 'defaultLegendInterpreter','latex');
set(f_c, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_c, 'defaultLineLineWidth',1.5);
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
legend(legendInfo_dVar, 'NumColumns',2)
hold off;
axis([0 Omega_dB(end) 2 20])











pp = zeros(1, length(Omega0));

f_snr_var = figure;
set(f_snr_var, 'Position',  [1260 360, 560, 420])
set(f_snr_var, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_snr_var, 'defaultLegendInterpreter','latex');
set(f_snr_var, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_snr_var, 'defaultLineLineWidth',1);
set(f_snr_var, 'color','w');
semilogy(Omega_dB, oPtheo_dVar(1,:), 'b', 'LineWidth',2)
hold on
semilogy(Omega_dB, oP_dVar(1,:), 'o','MarkerFaceColor', 'g', 'MarkerEdgeColor', 'b', 'MarkerIndices',1:markers_ind:length(Omega_dB))

for i = 1:length(Omega0)
    pp(i) = semilogy(Omega_dB,squeeze(linearApprox_log_Theo_var_dVar(1,i,:)),'--','LineWidth', 1.5);
end
semilogy(Omega_dB, oPasym_dVar(1,:), 'k:', 'LineWidth',2)

semilogy(Omega_dB, oPtheo_dVar(ceil(end/2),:), 'b', 'LineWidth',2)
semilogy(Omega_dB, oP_dVar(ceil(end/2),:), 'o','MarkerFaceColor', 'g', 'MarkerEdgeColor', 'b', 'MarkerIndices',1:markers_ind:length(Omega_dB))
semilogy(Omega_dB, oPtheo_dVar(end-2,:), 'b', 'LineWidth',2)
semilogy(Omega_dB, oP_dVar(end-2,:), 'o','MarkerFaceColor', 'g', 'MarkerEdgeColor', 'b', 'MarkerIndices',1:markers_ind:length(Omega_dB))
semilogy(Omega_dB, oPtheo_dVar(end,:), 'b', 'LineWidth',2)
semilogy(Omega_dB, oP_dVar(end,:), 'o','MarkerFaceColor', 'g', 'MarkerEdgeColor', 'b', 'MarkerIndices',1:markers_ind:length(Omega_dB))

for i = 1:length(Omega0)
    semilogy(Omega_dB,squeeze(linearApprox_log_Theo_var_dVar(ceil(end/2),i,:)),'-.','Color',get(pp(i), 'color'),'LineWidth', 1.5)    
    semilogy(Omega_dB,squeeze(linearApprox_log_Theo_var_dVar(end-2,i,:)),'-^','Color',get(pp(i), 'color'),'LineWidth', 1.5, 'MarkerSize',3,'MarkerFaceColor', 'b', 'MarkerIndices',1:markers_ind:length(Omega_dB))        
    semilogy(Omega_dB,squeeze(linearApprox_log_Theo_var_dVar(end,i,:)),'-x','Color',get(pp(i), 'color'),'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerIndices',1:markers_ind:length(Omega_dB))    
end
semilogy(Omega_dB, oPasym_dVar(ceil(end/2),:), 'k:', 'LineWidth',2)
semilogy(Omega_dB, oPasym_dVar(end-2,:), 'k:', 'LineWidth',2)
semilogy(Omega_dB, oPasym_dVar(end,:), 'k:', 'LineWidth',2)

for i = 1:length(Omega0)
    semilogy(10*log10(Omega0(i)),oPtheo_Om0_var_dVar(1,i),'xk','LineWidth', 2, MarkerSize = 10)
    semilogy(10*log10(Omega0(i)),oPtheo_Om0_var_dVar(ceil(end/2),i),'xk','LineWidth', 2, MarkerSize = 10)
    semilogy(10*log10(Omega0(i)),oPtheo_Om0_var_dVar(end-2,i),'xk','LineWidth', 2, MarkerSize = 10)
    semilogy(10*log10(Omega0(i)),oPtheo_Om0_var_dVar(end,i),'xk','LineWidth', 2, MarkerSize = 10)
    xline(10*log10(Omega0(i)),'k--')
end
axis([-inf inf 1e-6 2])
xlabel('$\Omega$ (dB)');
ylabel('OP');
grid on;

an_x = [.42 .8];
an_y = [.2 .8];

annotation('arrow',an_x,an_y)

dim = [.73 .57 .3 .3];
str = {'$0 \leq \Delta \leq 1$'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Color','k','BackgroundColor','w','Interpreter','latex');


% Generating legends
qw{1} = plot(nan, 'b', 'LineWidth',2);
qw{2} = plot(nan, 'o','MarkerFaceColor', 'g', 'MarkerEdgeColor', 'b');
for i = 1:length(Omega0)
    qw{2 + i} = plot(nan, 'Color',get(pp(i), 'color'), 'LineWidth',2);
end
qw{i+2+1} = plot(nan, 'k:','LineWidth',2);

legend([qw{:}], {'Theoretical','Simulated',...
    ['$\Omega_0$ = ' num2str(Omega0dB(1)),' dB'], ...
    ['$\Omega_0$ = ' num2str(Omega0dB(2)),' dB'] ...
    ['$\Omega_0$ = ' num2str(Omega0dB(3)),' dB'] ...
    ['$\Omega_0$ = ' num2str(Omega0dB(4)),' dB'] ...
    ['$\Omega_0$ = ' num2str(Omega0dB(5)),' dB'] ...
    ['$\Omega_0$ = ' num2str(Omega0dB(6)),' dB'] ...
    '$\Omega_0$ $\rightarrow$ $\infty$'}, 'location', 'southwest')
% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

%% Rice channels simulation. Comparison SISO, SC and MRC
clear; close all; clc;

addpath('Functions/')                               % It adds the path where the functions are located

N = 4;                                              % Number of receive branches
Nsim = 1e6;                                         % Number of random variable realizations
R = 1.7;                                            % Threshold rate for Outage metrics (R)
W_th = 2^(R) - 1;                                   % Minimum power required to receive the packet sent at a transmission rate of R

% The axis of average SNRs per branch on which it operates is defined
Omega_dB = -0:0.5:25;                               % Average SNR (in dB) per receive bracn
Omega = 10.^(Omega_dB/10);                          % Average SNR per receive branch

% Parameter for Figures
markers_ind = ceil(length(Omega_dB)/100*4);

%%  Outage probability
%   oP  --> Simulation of outage probability from averaged realizations
K = 5;

% Nsim realizations other than the channel are generated
h = sqrt(K/(K + 1)).*ones(N, Nsim) + 1/sqrt(K + 1).*(randn(N, Nsim) + 1i.*randn(N, Nsim))/sqrt(2);

g = abs(h(1,:));                                    % SISO
gSC =max(abs(h),[],1);                              % SC
gMRC=sqrt(sum(abs(h).^2,1));                        % MRC

SNR = ((abs(g).^2)'*Omega)';                        % Realizations of instantaneous SNR
SNRsc = ((abs(gSC).^2)'*Omega)';                    % Realizations of instantaneous SNR SC
SNRmrc = ((abs(gMRC).^2)'*Omega)';                  % Realizations of instantaneous SNR MRC

% (a) Simulated outage probabiilty -->  Pr(log2(1 + SNR < Rth));
% the number of times the channel capacity is less than the threshold rate is counted
oP = sum(log2(1 + SNR) < R, 2)/Nsim;
oPsc = sum(log2(1 + SNRsc) < R, 2)/Nsim;
oPmrc = sum(log2(1 + SNRmrc) < R, 2)/Nsim;

% (b) Theoretical outage probability Pr(Rth < log2(1 + SNR)) --> Pr(SNR < 2^Rth - 1)) = CDF(2^Rth - 1);
% %   oPtheo --> Theoretical expression of outage probability from the CDF
oPtheo = cdfSNRRiceThSC(W_th./Omega, K, 1,1);
oPtheoSC = cdfSNRRiceThSC(W_th./Omega, K, 1,N);
oPtheoMRC = cdfSNRRiceThMRC(W_th./Omega, K, 1,N);

eje = 0.01:1/100:2*N;

% (c) Sanity check for the PDF expressions
pdftheoSISO = pdfSNRRiceThSC(eje, K, 1,1);
pdftheoSISO2 = pdfSNRRiceThMRC(eje, K, 1,1);
pdftheoSC = pdfSNRRiceThSC(eje, K, 1,N);
pdftheoMRC = pdfSNRRiceThMRC(eje, K, 1,N);

[pdf_sim,x] = ksdensity((abs(g).^2),(0.01:1/100:2*N),'Support','positive');
[pdf_simsc,xsc] = ksdensity((abs(gSC).^2),(0.01:1/100:2*N),'Support','positive');
[pdf_simmrc,xmrc] = ksdensity((abs(gMRC).^2),(0.01:1/100:2*N),'Support','positive');



f_check = figure;
set(f_check, 'Position',  [40 560, 560, 420])
set(f_check, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_check, 'defaultLegendInterpreter','latex');
set(f_check, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_check,'defaultLineLineWidth',1.5);
set(f_check,'color','w');
plot(eje,pdftheoSISO,'b')
hold on
plot(eje,pdftheoSISO2,'x', 'MarkerIndices',1:length(eje)/80:length(eje));
plot(eje,pdftheoSC,'r')
plot(eje,pdftheoMRC,'k')
plot(x,pdf_sim,'bo', 'MarkerIndices',5:length(eje)/80:length(eje));
plot(xsc,pdf_simsc,'ro', 'MarkerIndices',1:length(eje)/80:length(eje));
plot(xmrc,pdf_simmrc,'ko', 'MarkerIndices',1:length(eje)/80:length(eje));
grid on
ylabel('PDF')
xlabel('x')



%% Theoretical linear approximation for several SNR values
Omega0dB_var = 15;
Omega0_var = 10.^(Omega0dB_var/10);

linearAppox_log_Theo_var = zeros(length(Omega0dB_var), length(Omega_dB));
linearAppox_log_Theo_varsc = zeros(length(Omega0dB_var), length(Omega_dB));
linearAppox_log_Theo_varmrc = zeros(length(Omega0dB_var), length(Omega_dB));
oPtheo_Om0_var = zeros(1, length(Omega0dB_var));
slope_temp = zeros(1, length(Omega0dB_var));
oPtheo_Om0_varsc = zeros(1, length(Omega0dB_var));
slope_tempsc = zeros(1, length(Omega0dB_var));
oPtheo_Om0_varmrc = zeros(1, length(Omega0dB_var));
slope_tempmrc = zeros(1, length(Omega0dB_var));
legendInfo{1} = 'Theoretical SISO';
legendInfo{2} = 'Simulated SISO';
legendInfo{3} = 'Theoretical SC';
legendInfo{4} = 'Simulated SC';
legendInfo{5} = 'Theoretical MRC';
legendInfo{6} = 'Simulated MRC';
for i = 1:length(Omega0dB_var)
    [slope_temp(i), oPtheo_Om0_var(i)] = Theoretical_Slope(Omega0_var(i), W_th, 'Rice', K);
    linearAppox_log_Theo_var(i,:) = oPtheo_Om0_var(i)*(Omega0_var(i)./Omega).^(slope_temp(i));
    [slope_tempsc(i), oPtheo_Om0_varsc(i)] = Theoretical_Slope(Omega0_var(i), W_th, 'RiceSC', [K,N]);
    linearAppox_log_Theo_varsc(i,:) = oPtheo_Om0_varsc(i)*(Omega0_var(i)./Omega).^(slope_tempsc(i));
    [slope_tempmrc(i), oPtheo_Om0_varmrc(i)] = Theoretical_Slope(Omega0_var(i), W_th, 'RiceMRC', [K,N]);
    linearAppox_log_Theo_varmrc(i,:) = oPtheo_Om0_varmrc(i)*(Omega0_var(i)./Omega).^(slope_tempmrc(i));
    legendInfo{i+6} = ['$\Omega_0$ = ' num2str(Omega0dB_var(i)),' dB'];
end
% 
% % %   oPasym --> approximated OP for high SNR (asymptotic)
oPasym = acdfSNRRiceTh(W_th./Omega, K, 1);
oPasymsc = acdfSNRRiceTh(W_th./Omega, K, 1).^N;

legendInfo{i+6+1} = '$\Omega_0$ $\rightarrow$ $\infty$';




f_OP_Rice = figure;
set(f_OP_Rice, 'Position',  [1260 560, 560, 420])
set(f_OP_Rice, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_OP_Rice, 'defaultLegendInterpreter','latex');
set(f_OP_Rice, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_OP_Rice, 'defaultLineLineWidth',1);
set(f_OP_Rice, 'color','w');
pp{1} = semilogy(Omega_dB, oPtheo, 'b', 'LineWidth',2);
hold on
pp{2} = semilogy(Omega_dB, oP, 'o','MarkerFaceColor', 'g', 'MarkerEdgeColor', 'b', 'MarkerIndices',1:markers_ind:length(Omega_dB));
pp{3} = semilogy(Omega_dB, oPtheoSC, 'r', 'LineWidth',2);
pp{4} = semilogy(Omega_dB, oPsc, 'o','MarkerFaceColor', 'g', 'MarkerEdgeColor', 'r', 'MarkerIndices',1:markers_ind:length(Omega_dB));
pp{5} = semilogy(Omega_dB, oPtheoMRC, 'k', 'LineWidth',2);
pp{6} = semilogy(Omega_dB, oPmrc, 'o','MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'MarkerIndices',1:markers_ind:length(Omega_dB));
pp{7} = semilogy(Omega_dB,linearAppox_log_Theo_var,'m--','LineWidth', 2);
pp{8} = semilogy(Omega_dB,linearAppox_log_Theo_varsc,'m--','LineWidth', 2);
pp{9} = semilogy(Omega_dB,linearAppox_log_Theo_varmrc,'m--','LineWidth', 2);
pp{10} = semilogy(Omega_dB, oPasym,':k', 'LineWidth',2);
pp{11} = semilogy(Omega_dB, oPasymsc,':k', 'LineWidth',2);
grid on;
legend([pp{1:7}, pp{10}], {legendInfo{1:7}, legendInfo{8}},'Location','southeast');
xlabel('$\Omega$ (dB)'); 
ylabel('OP');
axis([-inf inf 1e-10 2])
hold off;


%% Variation of parameters a and b (intercept and slope)
K = [1 5 10 15];

Omega_dB = -10:0.1:35;                    % Average SNR (in dB)
Omega = 10.^(Omega_dB/10);             % Average SNR
 
% Parameter for Figures
markers_ind = ceil(length(Omega_dB)/100*4);

m_var = zeros(length(K), length(Omega));
a_var = zeros(length(K), length(Omega));
for i = 1:length(K)
    [m_var(i,:), a_var(i,:)] = Theoretical_Slope(Omega, W_th, 'Rice', K(i));
end


%% Numerical and theoretical log10(CDF)'s derivative comparison
oPtheo = zeros(length(K), length(Omega));
G_SISO = zeros(length(K), length(Omega));
Gprime_SISO = zeros(length(K), length(Omega) - 1);
delta_num_SISO = zeros(length(K), length(Omega) - 1);
m_var = zeros(length(K), length(Omega));
delta_theo_SISO = zeros(length(K), length(Omega));
a_var = zeros(length(K), length(Omega));
oPtheomrc = zeros(length(K), length(Omega));
Gmrc = zeros(length(K), length(Omega));
Gprimemrc = zeros(length(K), length(Omega) - 1);
delta_nummrc = zeros(length(K), length(Omega) - 1);
m_varmrc = zeros(length(K), length(Omega));
delta_theomrc = zeros(length(K), length(Omega));
a_varmrc = zeros(length(K), length(Omega));
for k = 1:length(K)
    % oPtheo --> Theoretical expression of outage probability from the CDF
    oPtheo(k,:) = cdfSNRRiceThSC(W_th./Omega, K(k), 1,N);
    oPtheomrc(k,:) = cdfSNRRiceThMRC(W_th./Omega, K(k), 1,N);

    % Numerical and theoretical log10(CDF)'s derivative
    G_SISO(k,:) = log10(oPtheo(k,:));                                                                       % log10 of the CDF
    Gmrc(k,:) = log10(oPtheomrc(k,:));                                                                      % log10 of the CDF

    Gprime_SISO(k,:) = -diff(G_SISO(k,:))./diff(Omega_dB);                                                  % Numerical derivative
    delta_num_SISO(k,:) = 10*Gprime_SISO(k,:);
    Gprimemrc(k,:) = -diff(Gmrc(k,:))./diff(Omega_dB);                                                      % Numerical derivative
    delta_nummrc(k,:) = 10*Gprimemrc(k,:);

    [m_var(k,:), a_var(k,:)] = Theoretical_Slope(Omega, W_th, 'RiceSC', [K(k),N]);                          % Theoretical derivative
    delta_theo_SISO(k,:) = m_var(k,:);
    [m_varmrc(k,:), a_varmrc(k,:)] = Theoretical_Slope(Omega, W_th, 'RiceMRC', [K(k),N]);                   % Theoretical derivative
    delta_theomrc(k,:) = m_varmrc(k,:);
end




h1 = zeros(1,length(K));
h1mrc = zeros(1,length(K));
legendInfo_slope = cell(1, length(K));

f_ODO_Rice = figure;
set(f_ODO_Rice, 'Position',  [640 560, 560, 420])
set(f_ODO_Rice, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_ODO_Rice, 'defaultLegendInterpreter','latex');
set(f_ODO_Rice, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_ODO_Rice,'defaultLineLineWidth',1.5);
set(f_ODO_Rice,'color','w');
hold on
for i = 1:length(K)
    h1(i) = plot(Omega_dB(1:10:end), delta_theo_SISO(i, (1:10:end)),'--');
    h1mrc(i) = plot(Omega_dB(1:10:end), delta_theomrc(i, (1:10:end)),'Color',get(h1(i), 'color'));
end

for i = 1:length(K)
    plot(Omega_dB(2:end), delta_num_SISO(i, :),'s--','Color',get(h1(i), 'color'),'MarkerIndices', 2*i:markers_ind:length(Omega_dB)-1)
    plot(Omega_dB(2:end), delta_nummrc(i, :),'o--','Color',get(h1mrc(i), 'color'),'MarkerIndices', 3*i:markers_ind:length(Omega_dB)-1)
    legendInfo_slope{i} = ['K = ' num2str(K(i)), ' (both)'];
end
    

yline(N,':k','LineWidth',2);
xlabel('$\Omega_0$ (dB)');
ylabel('$\delta_{\rm Rice}$');
legend([h1mrc(:)], legendInfo_slope,'Location','northeast','AutoUpdate','off');
grid on

pos = [1 1 1.2 2];
rectangle('Position',pos,'Curvature',[1 1],'LineStyle','-','LineWidth',1)

an_1_x = [.36 .32];
an_1_y = [.2 .2];
annotation('textarrow',an_1_x,an_1_y,'String','SC')


pos = [-3.5 1 1.2 2];
rectangle('Position',pos,'Curvature',[1 1],'LineStyle','-','LineWidth',1)

an_2_x = [.2 .23];
an_2_y = [.23 .2];
annotation('textarrow',an_2_x,an_2_y,'String','MRC')

hold off





c_theo_SISO = 10./delta_theo_SISO;
c_num_SISO = 10./delta_num_SISO;

c_theo_MRC = 10./delta_theomrc;
c_num_MRC = 10./delta_nummrc;

h1 = zeros(1,length(K));
h2 = zeros(1,length(K));

f_c = figure;
set(f_c, 'Position',  [640, 50, 560, 420])
set(f_c, 'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',12);
set(f_c, 'defaultLegendInterpreter','latex');
set(f_c, 'defaultTextInterpreter','latex','defaultTextFontSize',14);
set(f_c, 'defaultLineLineWidth',1.5);
set(f_c, 'color','w');
hold on
for i = 1:length(K)
    h1(i) = plot(Omega_dB, c_theo_SISO(i,:),'--');
    h2(i) = plot(Omega_dB, c_theo_MRC(i,:),'Color',get(h1(i), 'color'));
end
for i = 1:length(K)
    plot(Omega_dB(2:end), c_num_SISO(i,:), 's--','Color',get(h1(i), 'color'),'MarkerIndices', 2*i:markers_ind:length(Omega_dB)-1)    
    plot(Omega_dB(2:end), c_num_MRC(i,:), 'o--','Color',get(h2(i), 'color'),'MarkerIndices', 2*i:markers_ind:length(Omega_dB)-1)    
end
yline(10/N,':k','LineWidth',2);
grid on
xlabel('$\Omega_0$ (dB)');
ylabel('$c_{\rm Rice}$ (dB)');

legend([h2(:)], legendInfo_slope,'Location','northeast','AutoUpdate','off');

pos = [0.5 5.5 1.5 1];
rectangle('Position',pos,'Curvature',[1 1],'LineStyle','-','LineWidth',1)

an_1_x = [.39 .27];
an_1_y = [.61 .61];
annotation('textarrow',an_1_x,an_1_y,'String','SC')

pos = [-4 7.5 1.5 1];
rectangle('Position',pos,'Curvature',[1 1],'LineStyle','-','LineWidth',1)

an_2_x = [.3 .18];
an_2_y = [.77 .77];
annotation('textarrow',an_2_x,an_2_y,'String','MRC')

hold off
axis([-5 Omega_dB(end) 0 10])
% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

function [m,a] = Theoretical_Slope(Omega0, PR, distribution, parameters)
% Function to obtain the analytical slope for an specific value of SNR
% Possible distributions:
% --------------------------------------
% INPUT --> DISTRIBUTION
% --------------------------------------
% Rayleigh --> Rayleigh
% Rice --> Rice
% Cascaded --> Cascaded Rayleigh
% TWDP --> Two-Wave with Diffuse Power
% RiceSC --> Rice (selection combining)
% RiceMRC --> Rice (Maximal ratio combining)
% -------------------------------------

if strcmp(distribution, 'Rayleigh')
    a = cdfSNRRayleighTh(PR./Omega0, 1);
    pdf = pdfSNRRayleighTh(PR./Omega0, 1);
    m = (PR./Omega0).*pdf./a;
elseif strcmp(distribution, 'Cascaded')
    a = cdfSNRProductTh(PR./Omega0, 1);
    pdf = pdfSNRProductTh(PR./Omega0, 1);
    m = (PR./Omega0).*pdf./a;
elseif strcmp(distribution, 'Rice')
    if nargin < 4
        disp('Parameter K is necessary for Rice distribution')
    end
    K = parameters;
    a = cdfSNRRiceTh(PR./Omega0, K, 1);
    pdf = pdfSNRRiceTh(PR./Omega0, K, 1);
    m = (PR./Omega0).*pdf./a;
elseif strcmp(distribution, 'TWDP')
    if nargin < 4
        disp('Parameter K and Delta are necessary for TWDP distribution')
    else
        if length(parameters) < 2
            disp('Parameter K and Delta are necessary for TWDP distribution')
        end
    end
    a = cdfSNRTWDPTh(parameters, 1, PR./Omega0);
    pdf = pdfSNRTWDPTh(PR./Omega0, parameters, 1);
    m = (PR./Omega0).*pdf./a;
elseif strcmp(distribution, 'RiceSC')
    if nargin < 4
        disp('Parameters K and N are necessary for Rice SC distribution')
    else
        if length(parameters) < 2
            disp('Parameter K and N are necessary for Rice SC distribution')
        end
    end
    K = parameters(1);
    N = parameters(2);
    a = cdfSNRRiceTh(PR./Omega0, K, 1).^N;
    pdf = pdfSNRRiceThSC(PR./Omega0,K,1,N);
    m = (PR./Omega0).*pdf./a;
elseif strcmp(distribution, 'RiceMRC')
    if nargin < 4
        disp('Parameters K and N are necessary for Rice MRC distribution')
    else
        if length(parameters) < 2
            disp('Parameter K and N are necessary for Rice MRC distribution')
        end
    end
    K = parameters(1);
    N = parameters(2);
    a = cdfSNRRiceThMRC(PR./Omega0, K, 1,N);
    pdf = pdfSNRRiceThMRC(PR./Omega0, K, 1,N);
    m = (PR./Omega0).*pdf./a;
else
    disp('Unsupported distribution')
    return
end

end
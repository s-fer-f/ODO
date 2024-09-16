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
% Two-Wave with Diffuse Power (TWDP)
% -------------------------------------

if strcmp(distribution, 'Rayleigh')
    a = cdfSNRRayleighTh(PR./Omega0, 1);
    m = (PR./Omega0).*pdfSNRRayleighTh(PR./Omega0, 1)./a;
elseif strcmp(distribution, 'Cascaded')
    a = cdfSNRProductTh(PR./Omega0, 1);
    m = (PR./Omega0).*pdfSNRProductTh(PR./Omega0, 1)./a;
elseif strcmp(distribution, 'Rice')
    if nargin < 4
        disp('Parameter K is necessary for Rice distribution')
    end
    K = parameters;
    a = cdfSNRRiceTh(PR./Omega0, K, 1);
    m = (PR./Omega0).*pdfSNRRiceTh(PR./Omega0, K, 1)./a;
elseif strcmp(distribution, 'TWDP')
    if nargin < 4
        disp('Parameter K and Delta are necessary for TWDP distribution')
    else
        if length(parameters) < 2
            disp('Parameter K and Delta are necessary for TWDP distribution')
        end
    end
    a = cdfSNRTWDPTh(parameters, 1, PR./Omega0);
    m = (PR./Omega0).*pdfSNRTWDPTh(PR./Omega0, parameters, 1)./a;
else
    disp('Unsupported distribution')
    return
end

end
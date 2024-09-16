% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

function CDF = cdfSNRRayleighTh(PR, SNRmedia)
% Theoretical power CDF (SNR) under Rayleigh channels
% Eq. 16 - Wireless Channel Modeling Perspectives for Ultra-Reliable Communications
% https://ieeexplore.ieee.org/document/8660712

CDF = 1 - exp(-PR./SNRmedia);
end

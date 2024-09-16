% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

function asint_CDF = acdfSNR_TWDP_Th(PR, K, Delta, SNR_av)
% Asymptotic CDF of the TWDP distribution for high SNR
% Eq. 40 - Wireless Channel Modeling Perspectives for Ultra-Reliable Communications
% https://ieeexplore.ieee.org/document/8660712
% "the tail can be lower-bounded through a scaled Rician tail"

asint_CDF = (acdfSNRRiceTh(PR, K, SNR_av)).'.*besseli(0, K*Delta);

end
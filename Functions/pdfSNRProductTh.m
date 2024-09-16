% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

function PDF = pdfSNRProductTh(PR, SNR_av)
% Theoretical power PDF (SNR) under Cascade Rayleigh channels
% https://ieeexplore.ieee.org/abstract/document/8950282
% Eq. 3 with  m = \hat{m} = 1 y \hat\omega = 1
% or
% Eq. 5 - Keyholes, Correlations, and Capacities of Multielement Transmit and Receive Antennas
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=994830

PDF = 2./SNR_av.*besselk(0, 2*sqrt(PR./SNR_av));
end
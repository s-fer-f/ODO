% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

function PDF = pdfSNRRayleighTh(PR, SNR_av)
% Theoretical power PDF (SNR) under Rayleigh channels

PDF = 1./SNR_av.*exp(-PR./(SNR_av));
end
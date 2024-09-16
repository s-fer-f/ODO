% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

function PDF = pdfSNRRiceTh(PR, K, SNR_av)
% Theoretical power PDF (SNR) under Rice channels

PDF = (1 + K)./(SNR_av).*exp(- K - (1 + K).*PR./(SNR_av) + 2*sqrt(PR.*K.*(1 + K)./(SNR_av))).*besseli(0, 2*sqrt(PR.*K.*(1 + K)./(SNR_av)), 1);
end
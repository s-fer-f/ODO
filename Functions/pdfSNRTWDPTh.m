% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

function pdf = pdfSNRTWDPTh(g, Param, gmed)
% Theoretical power PDF (SNR) under TWDP channels
% Eq. 23 - Statistics and system performance metrics for the Two Wave with Diffuse Power fading model
% https://ieeexplore.ieee.org/document/6814106

K = Param(1);
Delta = Param(2);

term_1 = (1 + K)./gmed;

f = @(theta) exp(-K*(1 + Delta*cos(theta)) - g*(1 + K)./gmed + ...
    (2*sqrt(g*K*(1 + Delta*cos(theta))*(K + 1)./gmed))).*...
    besseli(0, 2*sqrt(g*K*(1 + Delta*cos(theta))*(K+1)./gmed),1);

pdf = term_1.*integral(f, 0, 2*pi, 'ArrayValued', true)/(2*pi);
end
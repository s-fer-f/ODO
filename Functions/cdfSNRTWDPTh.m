% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

function CDF = cdfSNRTWDPTh(Param, SRN_av, PR)
% Eq. 23 - Statistics and System Performance Metrics for the Two Wave With Diffuse Power Fading Model
% https://ieeexplore.ieee.org/document/6814106

K = Param(1);
delta = Param(2);

CDF = 1 - integral(@(theta)(marcumq(sqrt(2*K*(1 + delta*cos(theta))), sqrt(PR)*sqrt(2*(K + 1)./SRN_av))), 0, pi, 'ArrayValued', true)./(pi);
end
% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

% PDF SC
function PDF = pdfSNRRiceThSC(g, K, gmed,N)
PDF = N.*pdfSNRRiceTh(g./gmed, K, 1).*cdfSNRRiceTh(g./gmed, K, 1).^(N-1);
end

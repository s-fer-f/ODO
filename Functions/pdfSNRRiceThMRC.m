% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

% PDF MRC
function PDF = pdfSNRRiceThMRC(g, K, gmed,N)
PDF = ((1+K)/gmed).^((N+1)./2).*(g./(N*K)).^((N-1)./2).*besseli(N-1,2*sqrt((1+K).*K.*N.*g./gmed),1).*exp(-K.*N-(1+K).*g./gmed+2*sqrt((1+K).*K.*N.*g./gmed));
end

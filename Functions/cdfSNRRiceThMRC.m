% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

% CDF MRC
function CDF = cdfSNRRiceThMRC(g, K, gmed,N)
CDF = zeros(1,length(g));
for ind = 1:length(g)
    CDF(ind) = integral(@(t)pdfSNRRiceThMRC(t, K, gmed,N),0,g(ind),'ArrayValued',true);
end
%yy = 1 - marcumq(sqrt(2*K*N), sqrt(2*(K + 1).*g./gmed),N);
end

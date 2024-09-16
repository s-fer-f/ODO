% Authors: F. Javier López-Martínez & Santiago Fernández
% Departamento de Teoría de la Señal, Telemática y Comunicaciones (TSTC)
% Universidad de Granada (UGR) - Granada, España
% Centro de Investigación en Tecnologías de la Información y las Comunicaciones CITIC-UGR - Granada, España
% 2024
%
% If you want to use these scripts, please reference the following article: https://arxiv.org/abs/2405.09336

function rvsim = genTWDPsim(K, Delta, Nsim, type)
%% GENERATING NORMALIZED TWDP VARIATES FROM ITS PHYSICAL MODEL
% Parameters:
%   K: Rician-like K parameter
%   Delta:  [0,1] amplitude imbalances between the two LoS components
%   Nsim: number of IFTR variates to be generated
%   type: 0 - amplitude  / 1 - power

gb = 1;
Sig = sqrt(gb/sqrt(2*(1 + K)));

b = Sig*Sig*K;
% Specular component calculation for special channel parameters (K,Delta)

V1 = sqrt(b*(1 + sqrt(1 - Delta.^2)));
V2 = sqrt(b*(1 - sqrt(1 - Delta.^2)));

% Parameters Check
% Delta should be (V1*V2*2)/(V1^2+V2^2);
% K should be (V1^2+V2^2)/(2*mu*Sig^2);

F1 = 2*pi*rand(Nsim,1);         % Phase 1
F2 = 2*pi*rand(Nsim,1);         % Phase 2

X = sqrt(Sig^2)*randn(Nsim,1); % Normal random
Y = sqrt(Sig^2)*randn(Nsim,1); %Normal random

% TWDP power Channel
h = (V1.*cos(F1) + V2.*cos(F2) + X).^2 + (V1.*sin(F1) + V2.*sin(F2) + Y).^2;

% Normalization to unit power
h = h./(V1^2 + V2^2 + 2*Sig^2);

if type==0
    rvsim = sqrt(h);
else
    rvsim = abs(h);
end

end
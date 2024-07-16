%% Upper boundary calculations El amin
% Calculate the upper boundary (i.e. z=0) conditions for a given S and c
%% Sw, c and normalized saturation
S = 1; c = 0.004; Sro = 0.001; Siw = 0.001;
Sw = S*(1 - Sro - Siw) + Siw;
%% Capillary pressure, pc
bW = 500; bO = -500; aW = 0.5; aO = 0.5;
% method 1
pc = bW*(S + eps(10^15))^-aW + bO*(1-S + eps(10^15))^-aO;  % eps(10^15) = 0.1250
% % method 2
% pc = bW*(S + 1/eps(10))^-aW + bO*(1-S + 1/eps(10))^-aO;
%% Water Velocity (Phase Flux)
krW = 1*S^4; krO = 1*(1-S)^4; % relperm
muW = 0.001; muO =0.001; rhoW = 1000; rhoO = 600; K0 = 0.002e-12;
H = 0.2; Nc = 100; g = - 9.81; 
dz = 0.2/100;
% method 1
vW = K0*(krW/muW)*((krO/muO)/(krW/muW + krO/muO))*(pc/dz - (rhoW - rhoO)*g); 
%% Water & oil pressure
pW = dz*((muW/krW)*(1/K0)*vW + rhoW*g);
pO = pW + pc;
% based on all examples, the pressure given in initial state and/or
% boundary conditions, represent the phase pressure of the oil phase. When
% Capillary pressure is negligible --> pW = p0. Contrary, when pc is
% significant --> pW = pO + pc (N.B. numbers sign)
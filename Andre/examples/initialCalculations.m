%% Initial calculations El amin
% Calculate the initial values and the values of the first iteration for a
% selected concentration and saturation
%% Sw, c and normalized saturation
Sw = 0.27; c =0.004;
S = (Sw-0.001)/(1-0.001-0.001);
%% Algebraic Equations
krW = 1*S^4; krO = 1*(1-S)^4; % relperm
bW = 500; bO = -500; aW = 0.5; aO = 0.5;
% method 1
pc = bW*(S + eps(10^15))^-aW + bO*(1-S + eps(10^15))^-aO; % eps(10^15) = 0.1250
% % method 2
% pc = bW*(S + 1/eps(10))^-aW + bO*(1-S + 1/eps(10))^-aO;
%% Velocity (Phase Flux)
muW = 0.001; muO =0.001; rhoW = 1000; rhoO = 600; K0 = 0.002e-12;
H = 0.2; Nc = 100; g = 9.81; 
dz = 0.2/100;
p = 0; pW = p + pc; pO = p;
% method 1
vW = K0*(krW/muW)*((krO/muO)/(krW/muW + krO/muO))*(pc/dz - (rhoW - rhoO)*g); 
% method 2
vW1 = (krW/muW)*K0*(pW/dz - rhoW*g);
vO = (krO/muO)*K0*(pO/dz - rhoO*g);
%%  Differential equations
dt = 43746.8354430380; D = 5.6e-8; phi0 = 0.3;


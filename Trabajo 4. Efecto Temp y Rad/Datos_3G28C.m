clear all; clc; close all

%% DATOS DEL PANEL
% Orden de los datos (FichaTecnica:FT) -> FT1BOL, FT1_2.5E14, FT1_2.5E14, FT1_2.5E14, Experimental
load('z.mat')
dat_exp = z;

Imp = [0.478, 0.4821, 0.4724, 0.4578, 0.4783];   %[A]
Isc = [0.506, 0.5009, 0.5009, 0.4858, 0.502925]; %[A]
Voc = [2.667, 2.560, 2.534, 2.480, 19.0442];     %[V]
Vmp = [2.371, 2.276, 2.229, 2.205, 17.3681];     %[V]
n     = Voc(5)/Voc(1);
Vmp(1:4) = n*Vmp(1:4);
Voc(1:4) = n*Voc(1:4);
alpha =  Vmp./Voc;
beta  =  Imp./Isc;
dat   = [Isc;Imp;Vmp;Voc;beta;alpha];

V = zeros(size(dat,2),200);                      % Inicializar Voltaje
for i = 1: size(dat,2)
    V(i,:) = linspace(0,dat(4,i),200);
end

%% DATOS DEL SATÉLITE

h = 450e3;                      % [m]
RT = 6378e3;                    % [m]
omega = 0.052;                  % [rad/s]
mu_T = 5.986e14;                % [m^3/s^2]
G = 1360;                       % [W/m2]
% PARAMETROS ORBITALES
% Órbita SS
a = h + RT;
n = sqrt(mu_T./a.^3);                   % [rad/s]
e = 0;
w = 0;
T_periodo = 2*pi*sqrt(a.^3/mu_T);       % [s]
pasoT=1e4;                              % paso temporal
T_max = 80+273.15;                      % [K]
T_min = -20+273.15;                     % [K]

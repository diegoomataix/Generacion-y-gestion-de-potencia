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

alpha_Imp = 1e-3*[0.28, 0.36, 0.20, 0.29,0.28 ];
alpha_Isc = 1e-3*[0.32, 0.33, 0.31, 0.39, 0.32 ];
alpha_Vmp = 1e-3*[-6.1, -6.8, -6.3, -6.4, -6.1 ];
alpha_Voc = 1e-3*[-6.0, -6.4, -6.2, -6.3, -6.0 ];

T_0  = [20, 20, 20, 20, 28] + 273.15;             % Temperatura nominal [K]

k = 1.3806503e-23;   % Boltzmann [J/K]
q = 1.60217646e-19;  % Carga electron [C]

Vt = n.*k.*T_0./q;                        % Voltaje termico
a_n = 1.3;
%% DATOS DEL SATÉLITE

h = 450e3;                      % [m]
RT = 6378e3;                    % [m]
omega = 0.052;                  % [rad/s]
mu_T = 5.986e14;                % [m^3/s^2]
G_0 = 1360;                       % [W/m2]
% PARAMETROS ORBITALES
% Órbita SS
a = h + RT;
e = 0;
w = 0;
T_periodo = (2*pi*sqrt(a.^3/mu_T))/100; % [s]
pasoT=1e2;                              % paso temporal
T_max = 80+273.15;                      % [K]
T_min = -20+273.15;                     % [K]

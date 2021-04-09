clear all; clc; close all

%% DATOS DEL PANEL
% Orden de los datos (FichaTecnica:FT) -> FT1BOL, FT1_2.5E14, FT1_2.5E14, FT1_2.5E14, Experimental
load('z.mat')
dat_exp = z;

Imp = [0.478, 0.4821, 0.4724, 0.4578, 0.4783];   %[A]
Isc = [0.506, 0.5009, 0.5009, 0.4858, 0.502925]; %[A]
Voc = [2.667, 2.560, 2.534, 2.480, 19.0442];     %[V]
Vmp = [2.371, 2.276, 2.229, 2.205, 17.3681];     %[V]

% Efecto temperatura
alpha_Imp = 1e-3*[0.28, 0.36, 0.20, 0.29,0.28 ];
alpha_Isc = 1e-3*[0.32, 0.33, 0.31, 0.39, 0.32 ];
alpha_Vmp = 1e-3*[-6.1, -6.8, -6.3, -6.4, -6.1 ];
alpha_Voc = 1e-3*[-6.0, -6.4, -6.2, -6.3, -6.0 ];

for i = 1:4
Isc_amb(i) =  (Isc(i) + alpha_Isc(i) * (8));
Imp_amb(i) =  (Imp(i) + alpha_Imp(i) * (8));
Voc_amb(i) = Voc(i) + (alpha_Voc(i) * (8));
Vmp_amb(i) = Vmp(i) + (alpha_Vmp(i) * (8));

Imp(i) = Imp_amb(i);     %[A]
Isc(i) = Isc_amb(i);     %[A]
Voc(i) = Voc_amb(i);     %[V]
Vmp(i) = Vmp_amb(i);     %[V]
end

n     = round(Voc(5)/Voc(1));
Vmp(1:4) = n*Vmp(1:4);
Voc(1:4) = n*Voc(1:4);
alpha =  Vmp./Voc;
beta  =  Imp./Isc;
dat   = [Isc;Imp;Vmp;Voc;beta;alpha];

V = zeros(size(dat,2),200);                       % Inicializar Voltaje
for i = 1: size(dat,2)
    V(i,:) = linspace(0,dat(4,i),200);
end

% Datos para determinar los efectos de la variación de Irradiancia y Temp
alpha_Imp = 1e-3*[0.28, 0.36, 0.20, 0.29,0.28 ];
alpha_Isc = 1e-3*[0.32, 0.33, 0.31, 0.39, 0.32 ];
alpha_Vmp = 1e-3*[-6.1, -6.8, -6.3, -6.4, -6.1 ];
alpha_Voc = 1e-3*[-6.0, -6.4, -6.2, -6.3, -6.0 ];

T_0  = [20, 20, 20, 20, 28] + 273.15;             % Temperatura nominal [K]

%% DATOS DEL SATÉLITE
h = 450e3;                                        % [m]
RT = 6378e3;                                      % [m]
omega = 0.052;                                    % [rad/s]
mu_T = 5.986e14;                                  % [m^3/s^2]
G_0 = 1367;                                       % [W/m2]
% PARAMETROS ORBITALES
% Órbita  
aa = h + RT;                                       % [m]
T_periodo = 4*pi/omega; % (2*pi*sqrt(aa.^3/mu_T))/40;            % [s]
pasoT = 1e2;                                      % paso temporal
T_max = 80+273.15;                                % [K]
T_min = -20+273.15;                               % [K]

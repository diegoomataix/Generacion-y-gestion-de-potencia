clear all; clc; close all

%% DATOS

h = [450e3 500e3 600e3];        % [m]
RT = 6378e3;                    % [m]
omega = [0.05 0.1 0.5];         % [rad/s]
eta = 0.3;                      % eficiencia paneles
fo = 0.8;                       % factor de ocupación
J2 =  1.0827e-3;                % J2 tierra
mu_T = 5.986e14; 
beta = 30;                      % Ángulo incidencia solar
G = 1360;                       % [W/m2]
Ap = 0.1*0.3;                   % [m^2]


%% PARAMETROS ORBITALES

% Órbita SS
a = h(1) + RT;                  % [m]
n = sqrt(mu_T/a^3);             % [rad/s]
e = 0;
i = rad2deg(acos((2*pi/(365.25*3600*24))/(-3/2*J2*RT^2/(a^2)*sqrt(mu_T/a^3)))); % [deg]
w = 0;
T = 2*pi*sqrt(a^3/mu_T);        % [s]

%% COSENO THETA

t = 1000; % linspace(0,T,100);  % [s]
RAAN = 15/3600 * t;             % [deg]
AV = n * t;                     % [rad/s]

rs_ECI = [1 0 0];               % vector sol (equinocio de primavera), en ejes ECI
nA_B = [cos(omega(1)*t) 0 0];   % vector normal a la cara del satelite, en ejes Body

Ci = [1 0        0;             % matriz de cambio de coordenadas (inclinacion), equivalente a matrix Cx en las diapos
      0 cosd(i)  sind(i);
      0 -sind(i) cosd(i)];
  
CRAAN = [cosd(RAAN)  sind(RAAN) 0;  % matriz de cambio de coordenadas (RAAN), equivalente a matrix Cy en las diapos
         -sind(RAAN) cosd(RAAN) 0;
         0           0         1 ];

CAV = [cos(AV+w)  sin(AV+w) 0;      % matriz de cambio de coordenadas (AV+w), equivalente a matrix Cy en las diapos
       -sin(AV+w) cos(AV+w) 0;
       0          0         1 ];
   
CECItoB = CAV*Ci*CRAAN;             % matriz de cambio de coordenadas
rs_B = inv(CECItoB)*rs_ECI';        % vector sol (equinocio de primavera), en ejes Body

THETA = acos(dot(rs_B,nA_B)/((norm(rs_B)*(norm(nA_B))))); % [rad/s]

%% ECLIPSES

nOP_B = [0 0 1];

rho = asind(RT/(RT+h(1)));              % [deg]
beta_s = 90 - rad2deg(acos(dot(rs_B,nOP_B)/((norm(rs_B)*(norm(nOP_B)))))); % [deg]
phi = 2*acosd(cosd(rho)/cosd(beta_s));  % [deg]

%% POTENCIA

Pi = G*eta*Ap*fo*sind(beta)*cos(THETA); % [W]
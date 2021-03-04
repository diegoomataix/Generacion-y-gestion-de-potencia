clear all; clc; close all

%% DATOS

h = [450e3 500e3 600e3];        % [m]
RT = 6378e3;                    % [m]
omega = [0.05 0.1 0.5];         % [rad/s]
eta = 0.3;                      % eficiencia paneles
fo = 0.8;                       % factor de ocupación
J2 =  1.0827e-3;                % J2 tierra
mu_T = 5.986e14;                % [m^3/s^2]
beta = 30;                      % Ángulo incidencia solar [deg]
G = 1360;                       % [W/m2]
Ap = 0.1*0.3;                   % [m^2]
fg = (1 + sqrt(2)) / 2;         % factor geometrico

%% PARAMETROS ORBITALES
% Órbita SS
a = h + RT;
n = sqrt(mu_T./a.^3);           % [rad/s]
e = 0;
w = 0;
T = 2*pi*sqrt(a.^3/mu_T);       % [s]

%% ACTITUD DEL SATELITE

for i = 1:3
    t(i,:) = linspace(0,T(i),100);
    alfa(i,:) = t(i,:)*n(i);             % [rad]
    for j = 1:3
        roll(i,j,:) = t(i,:) * omega(j); % [rad]
    end
end
   
%% ECLIPSES

rho = asind(RT./(RT+h));                % [deg]
alfa_eclipse_in = 180 - rho;            % [deg]
alfa_eclipse_out = 180 + rho;           % [deg]

%% POTENCIA

% Componente Perpendicular
P_my = G*eta*Ap*fo*fg*sind(beta)*(2/3); % [W]

% Componente Horizontal


% Pi = G*eta*Ap*fo*sind(beta)*cos(THETA); % [W]
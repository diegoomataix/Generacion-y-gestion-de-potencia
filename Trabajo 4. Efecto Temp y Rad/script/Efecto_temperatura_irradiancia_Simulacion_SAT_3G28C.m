clear all; clc; close all

Datos_3G28C
%% SIMULACION DE LA ACTITUD DEL SATELITE
[t, G, roll,roll_raw, roll_Tmax, T] = ...
    attitude(T_periodo, pasoT, omega,G_0, T_max, T_min); % [s], [deg], [deg], [K]

figure();plot(t,roll_raw)
figure();hold on; plot(cos(roll_raw - roll_Tmax));plot(cos(roll_raw));plot(G/G_0);legend('1','2','3');hold off
figure();plot(G);
%% EFECTO COND. AMBIENTALES
k_b = 1.3806503e-23;                                % Boltzmann [J/K]
q = 1.60217646e-19;                               % Carga electron [C]
for k = 1:size(T,2)
    Vt(k) =  n*k_b*T(k)/q;
    for j = 1:size(G,2)
        for i = 1:size(dat,2)
            %%% en algunos casos no tienen muy buena pinta los órdenes de
            %%% magnitud
            Isc_amb(i,j,k) = (G(j)/G_0) * (Isc(i) + alpha_Isc(i) * (T(k) - T_0(i)));
            Imp_amb(i,j,k) = (G(j)/G_0) * (Imp(i) + alpha_Imp(i) * (T(k) - T_0(i)));
            Voc_amb(i,j,k) = Voc(i) + a(i)*Vt(k) * log(G(j)/G_0) + (alpha_Voc(i) * (T(k) - T_0(i)));
            Vmp_amb(i,j,k) = Vmp(i) + a(i)*Vt(k) * log(G(j)/G_0) + (alpha_Vmp(i) * (T(k) - T_0(i)));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Funciones %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, G, roll, roll_raw, roll_Tmax, T] = attitude(T_periodo, pasoT, omega, G_0,T_max, T_min) % [s], [deg], [deg], [K]

t = linspace(0,T_periodo,pasoT); % [s]
% ANGULO INDICENCIA SOLAR
roll_raw = t * omega;            % [rad]
roll = wrapTo2Pi(roll_raw);      % [rad]
roll = rad2deg(roll);            % [deg]
roll_Tmax = (15 * omega);        % [rad] roll angle for 15s T delay
G = G_0 * cos(roll_raw);
T = ((T_max-T_min)/2) * cos(roll_raw -roll_Tmax) +T_min+((T_max-T_min)/2);
for i = 1:size(G,2)
    if G(i) <= G_0*cosd(75)
        G(i) = 1e-12;
    end
end
end

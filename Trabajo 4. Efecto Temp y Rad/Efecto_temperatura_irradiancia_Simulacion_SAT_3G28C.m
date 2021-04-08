clear all; clc; close all

Datos_3G28C
%% SIMULACION DE LA ACTITUD DEL SATELITE
[t, G, roll,roll_raw, roll_Tmax, T] = ...
    attitude(T_periodo, pasoT, omega,G_0, T_max, T_min); % [s], [deg], [deg], [K]

figure();plot(t,roll_raw)
figure();hold on; plot(cos(roll_raw - roll_Tmax));plot(cos(roll_raw));plot(G/G_0);legend('1','2','3');hold off
figure();plot(G);
%% EFECTO COND. AMBIENTALES
for k = 1:size(T,2)
    for j = 1:size(G,2)
        for i = 1:size(dat,2)
            %%% en algunos casos no tienen muy buena pinta los órdenes de
            %%% magnitud
            Isc_amb(i,j,k) = (G(j)/G_0) * (Isc(i) + alpha_Isc(i) * (T(k) - T_0(i)));
            Imp_amb(i,j,k) = (G(j)/G_0) * (Imp(i) + alpha_Imp(i) * (T(k) - T_0(i)));
            Voc_amb(i,j,k) = Voc(i) + a(i)*Vt(i) * log(G(j)/G_0) + (alpha_Voc(i) * (T(k) - T_0(i)));
            Vmp_amb(i,j,k) = Vmp(i) + a(i)*Vt(i) * log(G(j)/G_0) + (alpha_Vmp(i) * (T(k) - T_0(i)));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Funciones %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ipv,I0,Rs,Rsh] = param_1D_2R_Lap(Isc,Voc,Imp,Vmp,a,Vt)

A=-(2*Vmp-Voc)/(a*Vt)+(Vmp*Isc-Voc*Imp)/(Vmp*Isc+Voc*(Imp-Isc));
B=-Vmp*(2*Imp-Isc)/(Vmp*Isc+Voc*(Imp-Isc));
C=a*Vt/Imp;
D=(Vmp-Voc)/(a*Vt);

M1=0.3361;
M2=-0.0042;
M3=-0.0201;
sigma = -1-log(-B)-A;
Wn =-1-sigma -2/M1* (1-1./(1+M1*sqrt(sigma/2)./(1+M2*sigma.*exp(M3*sqrt(sigma)))) );

Rs=C*(Wn-(D+A));
Rsh=(Vmp-Imp*Rs)*(Vmp-Rs*(Isc-Imp)-a*Vt)/((Vmp-Imp*Rs)*(Isc-Imp)-a*Vt*Imp);
Ipv=(Rsh+Rs)/Rsh*Isc;
I0=((Rsh+Rs)/Rsh*Isc-Voc/Rsh)/(exp((Voc)/(a*Vt)));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ipv, I0, Rs, Rsh] = UND2R(Isc,Voc,Imp,Vmp,n, T)
k = 1.3806503e-23;   % Boltzmann [J/K]
q = 1.60217646e-19;  % Carga electron [C]

Vt = n.*k.*T./q;     % Voltaje termico
a_n = 1.3;
a = a_n;
[Ipv,I0,Rs,Rsh] = param_1D_2R_Lap(Isc,Voc,Imp,Vmp,a,Vt);


% Calculo I-V mediante funcion de Lambert W
% for i = 1:size(dat,2)     % bucle para los distintos casos
%     for j = 1:size(V,2)   % bucle para el voltaje
%         I(i,j) = (Rsh(i)*(Ipv(i) + I0(i)) - V(i,j))/(Rs(i) + Rsh(i)) - ((a(i)*Vt(i))/Rs(i))*...
%             lambertw(0,(((Rs(i)*Rsh(i)*I0(i))/((a(i)*Vt(i))*(Rs(i) + Rsh(i))))*exp((Rsh(i)*...
%             (Rs(i)*Ipv(i) + Rs(i)*I0(i) + V(i,j)))/(a(i)*Vt(i)*(Rs(i) + Rsh(i))))));
%     end
% end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C, m_func, m, gamma] = KH_model(betha, alpha)
%[ Isc, Imp, Vmp, Voc, betha, alpha, k]
%%% Determine coefficients %%%

% for i = 1:size(dat,2)
C = (1- betha - alpha)/ ((2*betha)-1);
m_func = lambertw(-1, (-(alpha^ (-1/C) * log(alpha))/C));   %
m= (m_func/log(alpha)) + (C^-1) + 1 ; % iterar para cada uno de los paneles
% m=real(m);
gamma = (2*betha -1)/ ((alpha^m *(m-1)));
% gamma=real(gamma);

%%% I y V %%%
% I_Kar = zeros(size(dat,2),size(V,2));
% for i = 1:size(dat,2)
%     for j = 1: size(V,2)
%         I_Kar(i,j) = Isc(1,i) * ( 1 - (1- gamma(i))*(V(i,j)/Voc(1,i)) - gamma(i)* ( (V(i,j)/Voc(1,i))^m(i) ));
%     end
% end
end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myplot(I, V)
for i = 1: size(dat,2)
    % Plot I-V
    figure()
    hold on
    grid on
    box on
    plot(  V(i,:) , I(i,:),'-k','LineWidth',2)
%     plot(dat_exp(:,1), dat_exp(:,2), '--k')
    axis tight
%     axis([0 dat(4,i)*1.2 0 dat(1,i)*1.2])
    xlabel('{\it I} [A]')
    ylabel('{\it t} [s]');
%     legend({'Modelo completo','Resultados experimentales'},'Location','northeast','NumColumns',2)
    box on
    set(gca,'FontSize',18)
    hold off
end
end

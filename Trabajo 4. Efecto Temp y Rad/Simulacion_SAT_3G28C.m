clear all; clc; close all

Datos_3G28C
%% SIMULACION DE LA ACTITUD DEL SATELITE
% La temperatura del panel varia entre -80 y -20 grados, estando su máximo
% retrasado 15 s respecto al máximo de la radiación sobre el panel.
[t, alfa, roll, T] = attitude(T_periodo, pasoT, n, omega); % [s], [deg], [deg], [K]       % roll es la actitud en radianes

figure();plot(t,roll)
figure();plot(T);
%% CALCULO CON LAS 3 RESISTENCIAS
R = [35, 37.2 , 42];           % Ohm

%% MODELOS DEL PANEL (emplean los datos que se obtienen con la simulación)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caso = 1;              %   1: Modelo explícito K&H         2: Modelo 1D2R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch(caso)
    %%% Modelo de Kalmalkar-Haneefa %%%
    case 1
        %%% Determine coefficients %%%
        [C, m_func, m, gamma, I_Kar] = KH_model(dat,V, dat(1,:), dat(4,:), dat(5,:), dat(6,:));
        %%% PLOT %%%
%         myplot(I_Kar, V, dat, dat_exp)
    %%% Modelo 1D2R %%%
    case 2
        %%% Determine coefficients %%%
        [Vt, I] = UND2R(dat, V, Isc,Voc,Imp,Vmp,n,T);
        %%% PLOT %%%
%         myplot(I, V, dat, dat_exp)
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
function [Vt, I] = UND2R(dat, V, Isc,Voc,Imp,Vmp,n, T)
k = 1.3806503e-23;   % Boltzmann [J/K]
q = 1.60217646e-19;  % Carga electron [C]

Vt = n.*k.*T./q;                        % Voltaje termico
a_n = 1.3;
for i = 1:size(dat,2)
    a(i) = a_n;
    [Ipv(i),I0(i),Rs(i),Rsh(i)] = param_1D_2R_Lap(Isc(i),Voc(i),Imp(i),Vmp(i),a(i),Vt(i));
end

% Calculo I-V mediante funcion de Lambert W
for i = 1:size(dat,2)     % bucle para los distintos casos
    for j = 1:size(V,2)   % bucle para el voltaje
        I(i,j) = (Rsh(i)*(Ipv(i) + I0(i)) - V(i,j))/(Rs(i) + Rsh(i)) - ((a(i)*Vt(i))/Rs(i))*...
            lambertw(0,(((Rs(i)*Rsh(i)*I0(i))/((a(i)*Vt(i))*(Rs(i) + Rsh(i))))*exp((Rsh(i)*...
            (Rs(i)*Ipv(i) + Rs(i)*I0(i) + V(i,j)))/(a(i)*Vt(i)*(Rs(i) + Rsh(i))))));
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C, m_func, m, gamma, I_Kar] = KH_model(dat,V, Isc, Voc, betha, alpha)
%[ Isc, Imp, Vmp, Voc, betha, alpha, k]
%%% Determine coefficients %%%
for i = 1:size(dat,2)

    C(i) = (1- betha(1,i) - alpha(1,i))/ ((2*betha(1,i))-1);
    m_func(i) = lambertw(-1, (-(alpha(1,i)^ (-1/C(i)) * log(alpha(1,i)))/C(i)));   %
    m(i)= (m_func(i)/log(alpha(1,i))) + (C(i)^-1) + 1 ; % iterar para cada uno de los paneles
    gamma(i) = (2*betha(1,i) -1)/ ((alpha(1,i)^m(i) *(m(i)-1)));

end

%%% I y V %%%
I_Kar = zeros(size(dat,2),size(V,2));
for i = 1:size(dat,2)
    for j = 1: size(V,2)
        I_Kar(i,j) = Isc(1,i) * ( 1 - (1- gamma(i))*(V(i,j)/Voc(1,i)) - gamma(i)* ( (V(i,j)/Voc(1,i))^m(i) ));
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, alfa, roll, T] = attitude(T_periodo, pasoT, n, omega) % [s], [deg], [deg], [K]

    t = linspace(0,T_periodo,pasoT); % [s]
    alfa = rad2deg(t*n);             % [deg]
    % ANGULO INDICENCIA SOLAR
    roll_raw = t * omega;            % [rad]
    roll = wrapTo2Pi(roll_raw);      % [rad] 
    roll = rad2deg(roll);            % [deg]
    % roll es la actitud en grados
    % TEMPERATURA
    roll_Tmax = rad2deg(15 * omega); % [deg] roll angle for 15s T delay
    T = ones(1,size(roll,2));
    % from roll = 0 to roll=roll_Tmax: y = 2.2376020x - 20.0000000
    % from roll = roll_Tmax to roll=360: y = -0.3171489x + 94.1736067
    for i = 1:size(roll,2)
    if roll(i) <= roll_Tmax
        T(i) = 2.2376020*(roll(i)) - 20;
    else
        T(i) = -0.3171489*(roll(i)) + 94.1736067;
    end
    if T(i) < 0
%         T(i)=0;
    end
    end
    % La temperatura del panel varia entre -80 y -20 grados, estando su máximo
    % retrasado 15 s respecto al máximo de la radiación sobre el panel.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myplot(I, V, dat, dat_exp)
for i = 1: size(dat,2)
    % Plot I-V
    figure()
    hold on
    grid on
    box on
    plot(  V(i,:) , I(i,:),'-k','LineWidth',2)
    plot(dat_exp(:,1), dat_exp(:,2), '--k')
    axis tight
    axis([0 dat(4,i)*1.2 0 dat(1,i)*1.2])
    xlabel('{\it V} [V]')
    ylabel('{\it I} [A]');
    legend({'Modelo completo','Resultados experimentales'},'Location','northeast','NumColumns',2)
    box on
    set(gca,'FontSize',18)
    hold off
end
end

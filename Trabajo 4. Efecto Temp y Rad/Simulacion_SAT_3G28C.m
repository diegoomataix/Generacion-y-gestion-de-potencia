clear all; clc; close all

Datos_3G28C
%% SIMULACION DE LA ACTITUD DEL SATELITE
% La temperatura del panel varia entre -80 y -20 grados, estando su máximo
% retrasado 15 s respecto al máximo de la radiación sobre el panel.
[t, G, roll,roll_raw, roll_Tmax, T] = attitude(T_periodo, pasoT, omega,G_0, T_max, T_min); % [s], [deg], [deg], [K]

% figure();plot(t,roll_raw)
% figure();hold on; plot(cos(roll_raw + roll_Tmax));plot(cos(roll_raw));hold off
% figure();plot(G);

%% EFECTO COND. AMBIENTALES
for j = 1:size(G,2)
for i = 1:size(dat,2)
Isc_amb(i,j) = (G(j)/G_0) * (Isc(i) + alpha_Isc(i) * (T(j) - T_0(i)));
Imp_amb(i,j) = (G(j)/G_0) * (Imp(i) + alpha_Imp(i) * (T(j) - T_0(i)));
Voc_amb(i,j) = Voc(i) + a*Vt(i) * log(G(j)/G_0) + (alpha_Voc(i) * (T(j) - T_0(i)));
Vmp_amb(i,j) = Vmp(i) + a*Vt(i) * log(G(j)/G_0) + (alpha_Vmp(i) * (T(j) - T_0(i)));
end
end
%% CALCULO CON LAS 3 RESISTENCIAS
R = [35, 37.2 , 42];           % Ohm
% V = IR

%% MODELOS DEL PANEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caso = 1;              %   1: Modelo explícito K&H         2: Modelo 1D2R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch(caso)
    %%% Modelo de Kalmalkar-Haneefa %%%
    case 1
        %%% Determine coefficients %%%
        syms I_Kar_sym
        for k = 1:size(G,2)
%             if G(k) <=0
                        
            for j = 1:size(dat,2)
                    [C, m_func, m, gamma] = KH_model((Imp_amb(j,k)/Isc_amb(j,k)), (Vmp_amb(j,k)/Voc_amb(j,k)));
                for i = 1:size(R)
                    I_Kar2(j,i,k) = double( vpasolve( I_Kar_sym == Isc_amb(j,k) * ( 1 - (1- gamma)*(I_Kar_sym*R(i)/Voc_amb(j,k)) - gamma* ( (I_Kar_sym*R(i)/Voc_amb(j,k))^m )) ) );
                end
            end
        end
        %%% PLOT %%%
         myplot(I_Kar2, V, dat, dat_exp)
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
function [C, m_func, m, gamma] = KH_model(betha, alpha)
%[ Isc, Imp, Vmp, Voc, betha, alpha, k]
%%% Determine coefficients %%%

% for i = 1:size(dat,2)
    C = (1- betha - alpha)/ ((2*betha)-1);
    m_func = lambertw(-1, (-(alpha^ (-1/C) * log(alpha))/C));   %
    m= (m_func/log(alpha)) + (C^-1) + 1 ; % iterar para cada uno de los paneles
    gamma = (2*betha -1)/ ((alpha^m *(m-1)));
    
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
    T = ((T_max-T_min)/2) * cos(roll_raw + roll_Tmax)+T_min;
    for i = 1:size(G,2)
    if G(i) <= G_0*cosd(75)
        G(i) = 0;
    end
    end
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
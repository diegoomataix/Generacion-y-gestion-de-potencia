clear all; clc; close all

Datos_3G28C
%% SIMULACION DE LA ACTITUD DEL SATELITE
[t, G, roll,roll_raw, roll_Tmax, T] = ...
    attitude(T_periodo, pasoT, omega,G_0, T_max, T_min); % [s], [deg], [deg], [K]

figure();hold on;plot(t,G,'k', 'Linewidth',1);grid on;axis tight;ylabel('{\it G} [W/m^2]');xlabel('{\it t} [s]');box on;set(gca,'FontSize',18)
%         legend('{\it G} [W/m^2]'); hold off

figure();hold on;plot(t,T,'-k', 'Linewidth',1);grid on;axis tight;ylabel('{\it T} [K]');xlabel('{\it t} [s]');box on;set(gca,'FontSize',18)
%         legend('{\it T} [K]'); hold off
%% EFECTO COND. AMBIENTALES
% Determinar a_n y voltaje termico
k = 1.3806503e-23;                                % Boltzmann [J/K]
q = 1.60217646e-19;                               % Carga electron [C]
% Vt = n.*k.*T_0./q;                                % Voltaje termico

for j = 1:size(G,2)
    Vt(j) =  n*k*T(j)/q;
    for i = 1:size(dat,2)
        Isc_amb(i,j) = (G(j)/G_0) * (Isc(i) + alpha_Isc(i) * (T(j) - T_0(i)));
        Imp_amb(i,j) = (G(j)/G_0) * (Imp(i) + alpha_Imp(i) * (T(j) - T_0(i)));
        Voc_amb(i,j) = Voc(i) + a(i)*Vt(j) * log(G(j)/G_0) + (alpha_Voc(i) * (T(j) - T_0(i)));
        Vmp_amb(i,j) = Vmp(i) + a(i)*Vt(j) * log(G(j)/G_0) + (alpha_Vmp(i) * (T(j) - T_0(i)));
    end
end
%% CALCULO CON LAS 3 RESISTENCIAS
R = [35, 37.2 , 42];           % Ohm
% V = IR
%% MODELOS DEL PANEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caso = 2;              %   1: Modelo explícito K&H         2: Modelo 1D2R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch(caso)
    %%% Modelo de Kalmalkar-Haneefa %%%
    case 1
        %%% Determine coefficients %%%
        syms I_Kar_sym
        % Inicializar vectores
        for k = 1:size(G,2)             % Size Irradiance vector
            for j = 1:size(dat,2)       % Size Data vector
                alpha_amb(j,k) = (Vmp_amb(j,k)/Voc_amb(j,k));
                betha_amb(j,k) = (Imp_amb(j,k)/Isc_amb(j,k));
                [C(j,k), m_func(j,k), m(j,k), gamma(j,k)] =...
                    KH_model( betha_amb(j,k), alpha_amb(j,k) );
                
                for i = 1:size(R,2)     % Size Resistance vector
                    % I_Kar(Irradiancia,CasoPanel,Resist)
                    I_Kar2(k,j,i) = double( vpasolve( I_Kar_sym == ...
                        Isc_amb(j,k) * ( 1 - (1- gamma(j,k))*(I_Kar_sym*R(i)/Voc_amb(j,k)) ...
                        - gamma(j,k)* ( (I_Kar_sym*R(i)/Voc_amb(j,k))^m(j,k) )), I_Kar_sym ) );
                end
            end
        end
        %%% PLOT %%%
        figure()
        hold on
        grid on
        plot(t, I_Kar2(:,5,1), 'k', 'Linewidth',1)
        plot(t, I_Kar2(:,5,2), '--k', 'Linewidth',1)
        plot(t, I_Kar2(:,5,3), '-.k', 'Linewidth',1)
        axis tight
        ylabel('{\it I} [V]')
        xlabel('{\it t} [s]');
        box on
        set(gca,'FontSize',18)
        %         legend('Modelo Kalmalkar & Haneefa', 'Modelo 1D2R')
        legend('Modelo Kalmalkar & Haneefa. {\it R} = 35.0 \Omega', ...
            'Modelo Kalmalkar & Haneefa. {\it R} = 37.2 \Omega',...
            'Modelo Kalmalkar & Haneefa. {\it R} = 42.0 \Omega')
        hold off
        
        %%%%%%%%%%%%% Modelo 1D2R %%%%%%%%%%%%%%%%%
    case 2
        %%% Determine coefficients %%%
        syms I_sym
        for k = 1:size(G,2)             % Size Irradiance vector
            for j = 1:size(dat,2)       % Size Data vector
                [Ipv(j,k), I0(j,k), Rs(j,k), Rsh(j,k)] = ...
                    UND2R(Isc_amb(j,k),Voc_amb(j,k),Imp_amb(j,k),Vmp_amb(j,k),n, T_0(j));
                
                for i = 1:size(R,2)     % Size Resistance vector
                    I(k,j,i) = double( vpasolve( I_sym ==(Rsh(j,k)*(Ipv(j,k) ...
                        + I0(j,k)) - I_sym*R(i))/(Rs(j,k) + Rsh(j,k))...
                        - ((a(j)*Vt(j))/Rs(j,k))*...
                        lambertw(0,(((Rs(j,k)*Rsh(j,k)*...
                        I0(j,k))/((a(j)*Vt(j))*(Rs(j,k) + Rsh(j,k))))*exp((Rsh(j,k)*...
                        (Rs(j,k)*Ipv(j,k) + Rs(j,k)*I0(j,k) + I_sym*R(i)))/(a(j)*Vt(j)*...
                        (Rs(j,k) + Rsh(j,k)))))),I_sym ) );
                end
            end
        end
        
        %%% PLOT %%%
        figure()
        hold on
        grid on
        plot(t, I(:,5,1), 'k', 'Linewidth',1)
        plot(t, I(:,5,2), '--k', 'Linewidth',1)
        plot(t, I(:,5,3), '-.k', 'Linewidth',1)
        axis tight
        ylabel('{\it I} [V]')
        xlabel('{\it t} [s]');
        box on
        set(gca,'FontSize',18)
        legend('Modelo 1D2R. {\it R} = 35.0 \Omega', ...
            'Modelo 1D2R. {\it R} = 37.2 \Omega',...
            'Modelo 1D2R. {\it R} = 42.0 \Omega')
        hold off
        
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
    if roll(i) >= 75 && roll(i) <= 285
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

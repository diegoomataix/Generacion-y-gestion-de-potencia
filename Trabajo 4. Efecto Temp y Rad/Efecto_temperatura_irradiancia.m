clear all; clc; close all

%% DATOS
% Orden de los datos (FichaTecnica:FT) -> FT1BOL, FT1_2.5E14, FT1_2.5E14, FT1_2.5E14, Experimental

load('z.mat')
dat_exp = z;
load('Isc_amb.mat')
load('Imp_amb.mat')
load('Vmp_amb.mat')
load('Voc_amb.mat')
load('Temp')
load('G')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pintar = 2;         %   1: Pintar efecto T        2: Pintar efecto G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch(pintar)
    case 1
        Irr = 1;
        Case = 5;
        Imp = [Imp_amb(Case,Irr,4), Imp_amb(Case,Irr,8), Imp_amb(Case,Irr,11), Imp_amb(Case,Irr,15), Imp_amb(Case,Irr,17)];   %[A]
        Isc = [Isc_amb(Case,Irr,4), Isc_amb(Case,Irr,8), Isc_amb(Case,Irr,11), Isc_amb(Case,Irr,15), Isc_amb(Case,Irr,17)]; %[A]
        Voc = [Voc_amb(Case,Irr,4), Voc_amb(Case,Irr,8), Voc_amb(Case,Irr,11), Voc_amb(Case,Irr,15), Voc_amb(Case,Irr,17)]; %[V]
        Vmp = [Vmp_amb(Case,Irr,4), Vmp_amb(Case,Irr,8), Vmp_amb(Case,Irr,11), Vmp_amb(Case,Irr,15), Vmp_amb(Case,Irr,17)]; %[V]
        n     = Voc(5)/Voc(1);
        alpha =  Vmp./Voc;
        beta  =  Imp./Isc;
        dat   = [Isc;Imp;Vmp;Voc;beta;alpha];
        
    case 2
        Temp = 11;
        Case = 5;
        Imp = [Imp_amb(Case,1,Temp), Imp_amb(Case,2,Temp), Imp_amb(Case,3,Temp), Imp_amb(Case,4,Temp), Imp_amb(Case,5,Temp)];   %[A]
        Isc = [Isc_amb(Case,1,Temp), Isc_amb(Case,2,Temp), Isc_amb(Case,3,Temp), Isc_amb(Case,4,Temp), Isc_amb(Case,5,Temp)]; %[A]
        Voc = [Voc_amb(Case,1,Temp), Voc_amb(Case,2,Temp), Voc_amb(Case,3,Temp), Voc_amb(Case,4,Temp), Voc_amb(Case,5,Temp)]; %[V]
        Vmp = [Vmp_amb(Case,1,Temp), Vmp_amb(Case,2,Temp), Vmp_amb(Case,3,Temp), Vmp_amb(Case,4,Temp), Vmp_amb(Case,5,Temp)]; %[V]
        
        n     = Voc(5)/Voc(1);
        alpha =  Vmp./Voc;
        beta  =  Imp./Isc;
        dat   = [Isc;Imp;Vmp;Voc;beta;alpha];
        
        V = zeros(size(dat,2),200);     % Inicializar Voltaje
        for i = 1: size(dat,2)
            V(i,:) = linspace(0,dat(4,i),200);
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caso = 2;              %   1: Modelo explícito K&H         2: Modelo 1D2R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(caso)
    %%% Modelo de Kalmalkar-Haneefa %%%
    case 1
        %%% Determine coefficients %%%
        [C, m_func, m, gamma, I_Kar] = KH_model(dat,V, dat(1,:), dat(4,:), dat(5,:), dat(6,:));
        %%% PLOT %%%
        myplot(I_Kar, V, dat, dat_exp)
    %%% Modelo 1D2R %%%
    case 2
        %%% Determine coefficients %%%
        T_0  = [20, 20, 20, 20, 28] + 273.15;     % Temperatura nominal [K]
        [Vt, I] = UND2R(dat, V, Isc,Voc,Imp,Vmp,n,T_0);
        %%% PLOT %%%
        myplot(I, V, dat, dat_exp)
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
function myplot(I, V, dat, dat_exp)
    figure()
    plot(dat_exp(:,1), dat_exp(:,2), '--k')
    legend({'1','2','3','4','5'},'Location','northeast','NumColumns',2)
for i = 1: size(dat,2)
    % Plot I-V
    hold on
    grid on
    box on
    plot(  V(i,:) , I(i,:),'LineWidth',1)
    axis tight
    axis([0 dat(4,i)*1.2 0 dat(1,i)*2])
    xlabel('{\it V} [V]')
    ylabel('{\it I} [A]');
    box on
    set(gca,'FontSize',18)
end
end

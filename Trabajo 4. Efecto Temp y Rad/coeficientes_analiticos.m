clear all; clc; close all

%% DATOS
%Orden de los datos (FichaTecnica:FT) -> FT1BOL, FT1_2.5E14, FT1_2.5E14, FT1_2.5E14, Experimental
%                                        [ Isc, Imp, Vmp, Voc, betha, alpha, k]

load('z.mat')
dat_exp = z;

Imp = [0.478, 0.4821, 0.4724, 0.4578, 0.4783];   %[A]
Isc = [0.506, 0.5009, 0.5009, 0.4858, 0.502925]; %[A]
Voc = [2.667, 2.560, 2.534, 2.480, 19.0442];     %[V]
Vmp = [2.371, 2.276, 2.229, 2.205, 17.3681];     %[V]
n     = round(Voc(5)/Voc(1));
Vmp(1:4) = n*Vmp(1:4);
Voc(1:4) = n*Voc(1:4);
alpha =  Vmp./Voc;
beta  =  Imp./Isc;
dat   = [Isc;Imp;Vmp;Voc;beta;alpha];

V = zeros(size(dat,2),200);     % Inicializar Voltaje
for i = 1: size(dat,2)
    V(i,:) = linspace(0,dat(4,i),200);
end


%% Resultados experimentales

caso = 1;
switch(caso)
    
    case 1

%% Modelo de Kalmalkar-Haneefa

%Ajuste mediante dos parametros gamma_k y m
% i = 1 - (1- gamma)* v - gamma * v^m
% i y v son la intensidad y el voltaje en forma adimensional.
gamma = zeros(1,size(dat,2));           % matriz de 1x(n de datos) para el coeficiente k del modelo de Das
m = zeros(1,size(dat,2));
m_simp = zeros(1,size(dat,2));
m_func = zeros(1,size(dat,2));
C = zeros(1,size(dat,2));

%%% Determine coefficients %%%
for i = 1:size(dat,2)
    
    C(i) = (1- dat(5,i) - dat(6,i))/ ((2*dat(5,i))-1);
    
    m_func(i) = lambertw(-1, ((dat(6,i)* (C(i)^-1) * log(dat(6,i)))/C(i)));   %
    m(i)= (m_func(i)/log(dat(6,i))) + (C(i)^-1) + 1 ; % iterar para cada uno de los paneles
    
    gamma(i) = (2*dat(5,i) -1)/ ((dat(6,i)^m(i) *(m(i)-1)));
    
end

%%% I y V %%%
I_Kar = zeros(size(dat,2),size(V,2));
for i = 1:size(dat,2)
    for j = 1: size(V,2)
        I_Kar(i,j) = dat(1,i) * ( 1 - (1- gamma(i))*(V(i,j)/dat(4,i)) - gamma(i)* ( (V(i,j)/dat(4,i))^m(i) ));
    end
end

%%% PLOT %%%
for i = 1: size(dat,2)
    
    % Plot I-V
    figure()
    hold on
    grid on
    box on
    % axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
    plot(  V(i,:) , I_Kar(i,:),'-k','LineWidth',2)
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

%% Modelo 1D2R

    case 2
        
% Constantes

k = 1.3806503e-23;   %Boltzmann [J/K]
q = 1.60217646e-19;  %Carga electron [C]

T  = 20 + 273.15;    %Temperatura nominal [K]
Vt = n*k*T/q;        %Voltaje termico  

a_n = 1.3; 

for i = 1:size(dat,2)
    a(i) = a_n;
    [Ipv(i),I0(i),Rs(i),Rsh(i)] = param_1D_2R_Lap(Isc(i),Voc(i),Imp(i),Vmp(i),a(i),Vt);   
end

% Calculo I-V mediante funcion de Lambert W

for i = 1:size(dat,2)     % bucle para los distintos casos
    for j = 1:size(V,2)   % bucle para el voltaje      
    I(i,j) = (Rsh(i)*(Ipv(i) + I0(i)) - V(i,j))/(Rs(i) + Rsh(i)) - ((a(i)*Vt)/Rs(i))*...
        lambertw(0,(((Rs(i)*Rsh(i)*I0(i))/((a(i)*Vt)*(Rs(i) + Rsh(i))))*exp((Rsh(i)*...
        (Rs(i)*Ipv(i) + Rs(i)*I0(i) + V(i,j)))/(a(i)*Vt*(Rs(i) + Rsh(i))))));  
    end
end

%%% PLOT %%%

for i = 1: size(dat,2)  
    % Plot I-V
    figure()
    hold on
    grid on
    box on
    % axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
    plot(V(i,:) , I(i,:),'-k','LineWidth',2)
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

%% Funciones

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
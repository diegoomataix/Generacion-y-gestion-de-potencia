clear all; clc; close all

%% DATOS
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

V = [dat_exp(:,1), dat_exp(:,1), dat_exp(:,1), dat_exp(:,1), dat_exp(:,1)]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caso = 1;              %   1: Modelo explícito K&H         2: Modelo 1D2R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(caso)
    %%% Modelo de Kalmalkar-Haneefa %%%
    case 1
        %%% Determine coefficients %%%
        [C, m_func, m, gamma, I_Kar] = KH_model(dat,V, dat(1,:), dat(4,:), dat(5,:), dat(6,:));
        I_Kar = real(I_Kar);
        I_result = I_Kar(5,:)';
        %%% PLOT %%%
        myplot(I_Kar, V, dat, dat_exp)
        
        
    %%% Modelo 1D2R %%%
    case 2
        %%% Determine coefficients %%%
        T  = [20, 20, 20, 20, 28] + 273.15;     % Temperatura nominal [K]
        [Vt, I] = UND2R(dat, V, Isc,Voc,Imp,Vmp,n,T);
        I_result = I(5,:)';
        %%% PLOT %%%
        myplot(I, V, dat, dat_exp)
end


%%%%%% CALCULO DE ERRORES %%%%%%%%%
[RMSE_PbP_Dim,RMSE_PbP_NonDim] = error_PbPfun(I_result, dat_exp(:,2));
RMSE_dim = RMSE(I_result, dat_exp(:,2));
RMSE_aDim = RMSE_NonDim(I_result, dat_exp(:,2));

plot_fun(dat_exp(:,1),RMSE_PbP_NonDim,'{\it V} [V]' ,'{\it \xi}',' ')


%% 
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

    % Plot I-V
    figure()
    hold on
    grid on
    box on
    plot(  V(1,:) , I(1,:),'k','LineWidth',2)
    plot(  V(2,:) , I(2,:),'--k','LineWidth',2)
    plot(  V(3,:) , I(3,:),':k','LineWidth',2)
    plot(  V(4,:) , I(4,:),'-.k','LineWidth',2)
    plot(  V(5,:) , I(5,:),'Color','#797d7f','LineWidth',2)
    axis tight
    axis([0 dat(4,5)*1.2 0 dat(1,5)*1.2])
    xlabel('{\it V} [V]')
    ylabel('{\it I} [A]');
    legend({'Caso 1', 'Caso 2', 'Caso 3', 'Caso 4', 'BOL'},'Location','northeast','NumColumns',2)
    box on
    set(gca,'FontSize',18)
    hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_fun(axisx,axisy,labelx,labely,yourtitle)
figure()

axis tight
% axis([min(axisx) max(axisx)  min(axisy) max(axisy)])
xlabel(labelx)
ylabel(labely)
title(yourtitle)
box on
set(gca,'FontSize',18)
grid on
grid minor
hold on
plot_result = plot(axisx,axisy,'-k','LineWidth',2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FORMULAS ERRORES 

function [RMSE_PbP_Dim,RMSE_PbP_NonDim] = error_PbPfun(I_modelo, I_exp)
I_sc = max(I_exp);
    RMSE_PbP_Dim(:) = abs(I_modelo(:) - I_exp(:));
    RMSE_PbP_NonDim(:) = (1/I_sc)*abs(I_modelo(:) - I_exp(:));
end 

function RMSE_dim = RMSE(I_modelo, I_exp)

RMSE_dim = ((sum((I_modelo - I_exp).^2)/size(I_exp,1))^0.5);
end 

function RMSE_adim = RMSE_NonDim(I_modelo, I_exp)
 
I_sc = max(I_exp);

RMSE_adim = (1/I_sc)*((sum((I_modelo - I_exp).^2)/size(I_exp,1))^0.5); %/sum(I_exp);
end 

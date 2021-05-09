clc;close all;clear all

%% Importacion y carga de datos
global limites
global n_dat
global caso
global eta_max
%CASO 3,3 V

% V1_33  = xlsread('efficiencies_dc_dc.xlsx',1,'A2:A26');
% I1_33  = xlsread('efficiencies_dc_dc.xlsx',1,'B2:B26');
% V2_33  = xlsread('efficiencies_dc_dc.xlsx',1,'D2:D26');
% I2_33  = xlsread('efficiencies_dc_dc.xlsx',1,'E2:E26');

load('I1_33');load('V1_33');load('I2_33');load('V2_33');
P1_33  = I1_33.*V1_33;
P2_33  = I2_33.*V2_33;
eta_33 = P2_33./P1_33;
eta_33(1) = 0;

%CASO 5 V

% V1_5  = xlsread('efficiencies_dc_dc.xlsx',2,'A2:A26');
% I1_5  = xlsread('efficiencies_dc_dc.xlsx',2,'B2:B26');
% V2_5  = xlsread('efficiencies_dc_dc.xlsx',2,'D2:D26');
% I2_5  = xlsread('efficiencies_dc_dc.xlsx',2,'E2:E26');

load('I1_5');load('V1_5');load('I2_5');load('V2_5');
P1_5  = I1_5.*V1_5;
P2_5  = I2_5.*V2_5;
eta_5 = P2_5./P1_5;
eta_5(1) = 0;

%CASO 15 V

% V1_15  = xlsread('efficiencies_dc_dc.xlsx',3,'A3:A9');
% I1_15  = xlsread('efficiencies_dc_dc.xlsx',3,'C3:C9');
% V2_15  = xlsread('efficiencies_dc_dc.xlsx',3,'D3:D9');
% I2_15  = xlsread('efficiencies_dc_dc.xlsx',3,'E3:E9');

load('I1_15');load('V1_15');load('I2_15');load('V2_15');
P1_15  = I1_15.*V1_15;
P2_15  = I2_15.*V2_15;
eta_15 = P2_15./P1_15;
V2_15 = [15; V2_15];
P2_15 = [0;P2_15 ];
eta_15 = [0;eta_15];
%% Estructuracion de datos
% Columnas:
limites = [size(I1_33,1), size(I1_5,1), size(I1_15,1)];
n_dat = limites


% Esto a tomar por culo de aqui
eta_exp = [eta_33; eta_5 ; eta_15];
eta_max = [max(eta_33), max(eta_5), max(eta_15)];
V_exp = [V2_33  ; V2_5  ; V2_15];
P_exp =  [P2_33  ; P2_5  ; P2_15];

%% Escoger el metodo de ajuste
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caso = 2;          % 1. Lineal 2. Polinomico 2º orden 3. Polinomico 3er orden
%dispositivo = 1;   % 1. 3,3 V  2. 5 V 3.15V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for caso = 1:3

for dispositivo = 1:3
switch (caso)
case 1

	if dispositivo ==1
U0 = -0.5;
[umin,fval]=fminsearch(@(u) RMSE_eta(u, eta_33, P2_33, V2_33), U0);

Param_1_33 = umin
RMSE_1_33 = fval

elseif dispositivo ==2
U0 = -0.5;
[umin,fval]=fminsearch(@(u) RMSE_eta(u, eta_5, P2_5, V2_5), U0);

Param_1_5 = umin
RMSE_1_5 = fval

elseif dispositivo ==3
U0 = -0.5;
[umin,fval]=fminsearch(@(u) RMSE_eta(u, eta_15, P2_15, V2_15), U0);

Param_1_15 = umin
RMSE_1_15 = fval
end

case 2

if dispositivo == 1
U_2 = [-0.5, 0.1];
[umin,fval]=fminsearch(@(u) RMSE_eta(u, eta_33, P2_33, V2_33), U_2);

Param_2_33 = umin
RMSE_2_33 = fval

elseif dispositivo == 2
U_2 = [-0.5, 0.1];
[umin,fval]=fminsearch(@(u) RMSE_eta(u, eta_5, P2_5, V2_5), U_2);

Param_2_5 = umin
RMSE_2_5 = fval

elseif dispositivo == 3
U_2 = [-0.5, 0.1];
[umin,fval]=fminsearch(@(u) RMSE_eta(u, eta_15, P2_15, V2_15), U_2);

Param_2_15 = umin
RMSE_2_15 = fval
end


case 3

if dispositivo == 1
U_3 = [-0.5, 0.1, 0.1];
[umin,fval]=fminsearch(@(u) RMSE_eta(u, eta_33, P2_33, V2_33), U_3);

Param_3_33 = umin
RMSE_3_33 = fval

elseif dispositivo == 2
	U_3 = [-0.5, 0.1, 0.1];
	[umin,fval]=fminsearch(@(u) RMSE_eta(u, eta_5, P2_5, V2_5), U_3);

	Param_3_5 = umin
	RMSE_3_5 = fval

elseif dispositivo == 3

	U_3 = [-0.5, 0.1, 0.01];
	[umin,fval]=fminsearch(@(u) RMSE_eta(u, eta_15, P2_15, V2_15), U_3);

	Param_3_15 = umin
	RMSE_3_15 = fval
end

end
end

%% PLOTS

if caso == 1

eta_modelo_33_c1 = DCDC(Param_1_33, V2_33, P2_5, max(eta_33));
eta_modelo_5_c1 = DCDC(Param_1_5, V2_5, P2_5, max(eta_5));
eta_modelo_15_c1 = DCDC(Param_1_15, V2_15, P2_15, max(eta_15));

elseif caso == 2

eta_modelo_33_c2 = DCDC2(Param_2_33, V2_33, P2_33, max(eta_33));
eta_modelo_5_c2  = DCDC2(Param_2_5, V2_5, P2_5, max(eta_5));
eta_modelo_15_c2 = DCDC2(Param_2_15, V2_15, P2_15, max(eta_15));

elseif caso == 3

eta_modelo_33_c3 = DCDC3(Param_3_33, V2_33, P2_33, max(eta_33));
eta_modelo_5_c3  = DCDC3(Param_3_5, V2_5, P2_5, max(eta_5));
eta_modelo_15_c3 = DCDC3(Param_3_15, V2_15, P2_15, max(eta_15));

end
end
%end CASO


figure(1)
hold on
grid on
plot(P2_33,eta_33, '-k', 'LineWidth' ,1.5)
plot(P2_33,eta_modelo_33_c1, '--*k')
plot(P2_33,eta_modelo_33_c2, '--*k')
plot(P2_33,eta_modelo_33_c3, '--dk')
set(gca,'FontSize',18)
xlabel('{\it{P_{out}}} [W]')
ylabel('\eta')
legend('Datos experimentales', 'Ajuste 1' ,'Ajuste 2', 'Ajuste 3')
hold off

figure(2)
hold on
grid on
plot(P2_5,eta_5, '-k','LineWidth', 1.5)
plot(P2_5,eta_modelo_5_c1,'--*k')
plot(P2_5,eta_modelo_5_c2, '--+k')
plot(P2_5,eta_modelo_5_c3, '--dk')
set(gca,'FontSize',18)
xlabel('{\it{P_{out}}} [W]')
ylabel('\eta')
legend('Datos experimentales', 'Ajuste 1', 'Ajuste 2', 'Ajuste 3')
hold off


figure(3)
hold on
grid on
plot(P2_15,eta_15, '-k' , 'LineWidth', 1.5)
plot(P2_15,eta_modelo_15_c1, '--*k')
plot(P2_15,eta_modelo_15_c2,'--+k')
plot(P2_15,eta_modelo_15_c3, '--dk')
set(gca,'FontSize',18)
xlabel('{\it{P_{out}}} [W]')
ylabel('\eta')
legend('Datos experimentales', 'Ajuste 1', 'Ajuste 2', 'Ajuste 3' )
hold off

%end

%Funciones

function eta_ajust_1 = DCDC(u, V2, P2, eta_max)
f  = u(1).*P2;
eta_ajust_1 = eta_max*(1-exp(f));
end

function eta_ajust_2 = DCDC2(u, V2, P2, eta_max)
f = u(1).*P2 + u(2).*P2.^2;
eta_ajust_2 = eta_max*(1-exp(f));
end

function eta_ajust_3 = DCDC3(u, V2, P2, eta_max)
f = u(1).*P2 + u(2).*P2.^2 + u(3).*P2.^3;
eta_ajust_3 = eta_max*(1-exp(f));
end

function error = RMSE_eta(u, eta_exp, P_exp, V_exp)
global n_dat; global limites; global caso
global eta_max

eta_modelo = zeros (1,size(eta_exp,1));

for i = 1:size(eta_exp,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if caso == 1
	eta_modelo(i) = DCDC(u, V_exp(i), P_exp(i), max(eta_exp))  ;
    elseif caso == 2
    eta_modelo(i) = DCDC2(u, V_exp(i), P_exp(i), max(eta_exp));
    elseif  caso == 3
    eta_modelo(i) = DCDC3(u, V_exp(i), P_exp(i), max(eta_exp));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
		eta_mod_err = eta_modelo;
		eta_exp_err = eta_exp;
		error_vect = ((sum (( eta_mod_err' - eta_exp_err).^2 ) / size(eta_exp,1))^0.5);

% Error resultante
error = error_vect;
end

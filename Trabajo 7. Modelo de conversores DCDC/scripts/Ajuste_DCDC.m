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

%% Estructuracion de datos
% Columnas:
limites = [size(I1_33,1), size(I1_5,1), size(I1_15,1)];
n_dat = limites


eta_exp = [eta_33; eta_5 ; eta_15];
eta_max = [max(eta_33), max(eta_5), max(eta_15)];
V_exp = [V2_33  ; V2_5  ; V2_15];
P_exp =  [P2_33  ; P2_5  ; P2_15];

%% Escoger el metodo de ajuste
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caso = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch (caso)
case 1
U0 = -0.5;
[umin,fval]=fminsearch(@(u) RMSE_eta(u, eta_exp, P_exp, V_exp), U0);

Param_1 = umin
RMSE_1 = fval

case 2

U_2 = [-0.5, 0.1];
[umin,fval]=fminsearch(@(u) RMSE_eta(u, eta_exp, P_exp, V_exp), U_2);

Param_2 = umin
RMSE_2 = fval

case 3

U_3 = [-0.5, 0.1, 0.1];
[umin,fval]=fminsearch(@(u) RMSE_eta(u, eta_exp, P_exp, V_exp), U_3);

Param_3 = umin
RMSE_3= fval
end

%% PLOTS

eta_modelo = zeros (1,sum(limites));
for i = 1:size(eta_exp,1)
	if i <= limites(1)
		j = 1;
	elseif i <= limites(2) + limites(1)
		j = 2;
	else
		j = 3;
	end
if caso == 1
eta_modelo(i) = DCDC(Param_1, V_exp(i), P_exp(i), eta_max(j));
elseif caso == 2
eta_modelo(i) = DCDC2(Param_2, V_exp(i), P_exp(i), eta_max(j));
elseif  caso == 3
eta_modelo(i) = DCDC3(Param_3, V_exp(i), P_exp(i), eta_max(j));
end
end

figure()
hold on
plot(P2_33,eta_33,  '--k')
plot(P2_33,eta_modelo(1:limites(1)))
hold off

figure()
hold on
plot(P2_5,eta_5, '--k')
plot(P2_5,eta_modelo(limites(1)+1:limites(1)+limites(2)))
hold off

figure()
hold on
plot(P2_5,eta_5,  '--k')
plot(P2_15,eta_modelo(limites(1)+limites(2)+1:limites(1)+limites(2)+limites(3)))
hold off

%Funciones

function eta_ajust_1 = DCDC(u, V2, P2, eta_max)
f  = (u(1)/V2)*P2;
eta_ajust_1 = eta_max*(1-exp(f));
end

function eta_ajust_2 = DCDC2(u, V2, P2, eta_max)
f = (u(1)/V)*P2 + (u(2)/V2)*P2^2;
eta_ajust_2 = eta_max*(1-exp(f));
end

function eta_ajust_3 = DCDC3(u, V2, P2, eta_max)
f = (u(1)/V2)*P2 + (u(2)/V2)*P2^2 + (u(3)/V2)*P2^3;
eta_ajust_3 = eta_max*(1-exp(f));
end

function error = RMSE_eta(u, eta_exp, P_exp, V_exp)
global n_dat; global limites; global caso
global eta_max

eta_modelo = zeros (1,sum(limites));
for i = 1:size(eta_exp,1)
	if i <= limites(1)
		j = 1;
	elseif i <= limites(2) + limites(1)
		j = 2;
	else
		j = 3;
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if caso == 1
	eta_modelo(i) = DCDC(u, V_exp(i), P_exp(i), eta_max(j))  ;
    elseif caso == 2
    eta_modelo(i) = DCDC2(u, V_exp(i), P_exp(i), eta_max(j));
    elseif  caso == 3
    eta_modelo(i) = DCDC3(u, V_exp(i), P_exp(i), eta_max(j));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if i == limites(1)
		eta_mod_err = eta_modelo(1:limites(1));
		eta_exp_err = eta_exp(1:limites(1));
		error_vect(j) = ((sum (( eta_mod_err' - eta_exp_err).^2 ) / n_dat(j))^0.5);

	elseif i == limites(1) + limites(2)
		eta_mod_err = eta_modelo( limites(1)+1:limites(1)+limites(2));
		eta_exp_err = eta_exp( limites(1)+1:limites(1)+limites(2));
		error_vect(j) = ((sum (( eta_mod_err' - eta_exp_err).^2 ) / n_dat(j))^0.5);

	elseif  i == sum(limites)
		eta_mod_err = eta_modelo(limites(1)+limites(2)+1:limites(1)+limites(2)+limites(3));
		eta_exp_err = eta_exp(limites(1)+limites(2)+1:limites(1)+limites(2)+limites(3));
		error_vect(j) = ((sum (( eta_mod_err' - eta_exp_err).^2 ) / n_dat(j))^0.5);
     else
	end

end

error = error_vect(1) + error_vect(2)+ error_vect(3);
end

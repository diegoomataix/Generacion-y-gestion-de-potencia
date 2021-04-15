clc; close all; clear all;

%% CARGAR DATOS
load('descarga5A')
load('carga5A')
load('descarga2_5A')
load('carga2_5A')
load('descarga1_5A')
load('carga1_5A')

%% PROCESAR DATOS
descarga5A(:,4) = descarga5A(:,1).*descarga5A(:,2);
carga5A(:,4) = carga5A(:,1).*carga5A(:,2);
descarga2_5A(:,4) = descarga2_5A(:,1).*descarga2_5A(:,2);
carga2_5A(:,4) = carga2_5A(:,1).*carga2_5A(:,2);
descarga1_5A(:,4) = descarga1_5A(:,1).*descarga1_5A(:,2);
carga1_5A(:,4) = carga1_5A(:,1).*carga1_5A(:,2);

%% APARTADO 1
V_nom_i = 4.2;
C_nom_cell = 2750;
n_serie = round([(max(max(descarga5A(:,3))))/V_nom_i ; (max(max(descarga1_5A(:,3))))/V_nom_i ;  (max(max(descarga2_5A(:,3))))/V_nom_i ])

t = [10615,21365,36528];    % [s]
I = [5, 2.5, 1.5];

k = log(t(3)/t(1))/ log(I(1)/I(3));

C_p = I.^k.*t;              % [A·s]
n_par = round(C_p ./ C_nom_cell);

%% APARTADO 2



%% BACKUP
% {Tiempo , Intensidad , Voltaje}
% [descarga5A,~,~]=(xlsread('../ensayos_bateria.xlsx','descarga 5A'));
% [carga5A,~,~]=(xlsread('../ensayos_bateria.xlsx','carga 5A'));
% [descarga2_5A,~,~]=(xlsread('../ensayos_bateria.xlsx','descarga 2.5A'));
% [carga2_5A,~,~]=(xlsread('../ensayos_bateria.xlsx','carga 2.5A'));
% [descarga1_5A,~,~]=(xlsread('../ensayos_bateria.xlsx','descarga 1.5A'));
% [carga1_5,~,~]=(xlsread('../ensayos_bateria.xlsx','carga 1.5A'));

% Ley de Peukert
%I = 1.5;                                         % [A]
%C_espec = 2.750*3600*[1, 0.97, 0.95, 0.92];      % [A·s]
%I_espec = [0.55, 2.75, 5.5, 8,25]
% K = -(log(t)-log(C))/log(I);        % Constante de Peukert
%alpha_espec = log(C_espec(4)/C_espec(1))/log(I_espec(4)/I_espec(1))
%k_espec = (2*alpha_espec+1)/(1+alpha_espec)

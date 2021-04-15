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
n_serie = round([(max(max(descarga5A(:,3))))/V_nom_i ; (max(max(descarga1_5A(:,3))))/V_nom_i ;  (max(max(descarga2_5A(:,3))))/V_nom_i ])

%Ley de Peukert
I = 1.5;           %[A]
t = 36528;         %[s]
C = [2.750*3600, 2.750*0.97*3600];     %[A·s]
I = [0.55, 2.75, 5.5, 8,25]

% K = -(log(t)-log(C))/log(I);        % Constante de Peukert

alpha = log(C(2)/C(1))/log(I(2)/I(1))

k = (2*alpha+1)/(1+alpha)



%% BACKUP
% {Tiempo , Intensidad , Voltaje}
% [descarga5A,~,~]=(xlsread('../ensayos_bateria.xlsx','descarga 5A'));
% [carga5A,~,~]=(xlsread('../ensayos_bateria.xlsx','carga 5A'));
% [descarga2_5A,~,~]=(xlsread('../ensayos_bateria.xlsx','descarga 2.5A'));
% [carga2_5A,~,~]=(xlsread('../ensayos_bateria.xlsx','carga 2.5A'));
% [descarga1_5A,~,~]=(xlsread('../ensayos_bateria.xlsx','descarga 1.5A'));
% [carga1_5,~,~]=(xlsread('../ensayos_bateria.xlsx','carga 1.5A'));

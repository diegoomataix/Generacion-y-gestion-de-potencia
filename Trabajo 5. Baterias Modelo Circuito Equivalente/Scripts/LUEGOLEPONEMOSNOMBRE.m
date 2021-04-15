clc; close all; clear all;

%% CARGAR DATOS

% {Tiempo , Intensidad , Voltaje}
[descarga5A,~,~]=(xlsread('../ensayos_bateria.xlsx','descarga 5A'));
[carga5A,~,~]=(xlsread('../ensayos_bateria.xlsx','carga 5A'));
[descarga2_5A,~,~]=(xlsread('../ensayos_bateria.xlsx','descarga 2.5A'));
[carga2_5A,~,~]=(xlsread('../ensayos_bateria.xlsx','carga 2.5A'));
[descarga1_5A,~,~]=(xlsread('../ensayos_bateria.xlsx','descarga 1.5A'));
[carga1_5,~,~]=(xlsread('../ensayos_bateria.xlsx','carga 1.5A'));


%% PROCESAR DATOS

descarga5A(:,4) = descarga5A(:,1).*descarga5A(:,2);
carga5A(:,4) = carga5A(:,1).*carga5A(:,2);
descarga2_5A(:,4) = descarga2_5A(:,1).*descarga2_5A(:,2);
carga2_5A(:,4) = carga2_5A(:,1).*carga2_5A(:,2);
descarga1_5A(:,4) = descarga1_5A(:,1).*descarga1_5A(:,2);
carga1_5(:,4) = carga1_5(:,1).*carga1_5(:,2);

%% APARTADO 1
V_nom_i = 4.2;
n_serie = round([(max(max(descarga5A(:,3))))/V_nom_i ; (max(max(descarga1_5A(:,3))))/V_nom_i ;  (max(max(descarga2_5A(:,3))))/V_nom_i ])

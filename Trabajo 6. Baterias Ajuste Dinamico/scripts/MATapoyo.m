clc;close all;clear all

%% DESCARGAS DINAMICAS

% Importación datos

I_D  = xlsread('medidas_bateria.xlsx',1,'B2:B7213');
V_D  = xlsread('medidas_bateria.xlsx',1,'C2:C7213');
t_D  = xlsread('medidas_bateria.xlsx',1,'A2:A7213');
phi  = xlsread('medidas_bateria.xlsx',1,'D2:D7213');
phi2 = xlsread('medidas_bateria.xlsx',1,'E2:E7213');
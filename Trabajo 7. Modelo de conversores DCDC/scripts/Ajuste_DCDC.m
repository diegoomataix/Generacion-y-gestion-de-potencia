clc;close all;clear all

%% Importacion y carga de datos

%CASO 3,3 V

% V1_33  = xlsread('efficiencies_dc_dc.xlsx',1,'A2:A26');
% I1_33  = xlsread('efficiencies_dc_dc.xlsx',1,'B2:B26');
% V2_33  = xlsread('efficiencies_dc_dc.xlsx',1,'D2:D26');
% I2_33  = xlsread('efficiencies_dc_dc.xlsx',1,'E2:E26');

load('I1_33');load('V1_33');load('I2_33');load('V2_33');

%CASO 5 V

% V1_5  = xlsread('efficiencies_dc_dc.xlsx',2,'A2:A26');
% I1_5  = xlsread('efficiencies_dc_dc.xlsx',2,'B2:B26');
% V2_5  = xlsread('efficiencies_dc_dc.xlsx',2,'D2:D26');
% I2_5  = xlsread('efficiencies_dc_dc.xlsx',2,'E2:E26');

load('I1_5');load('V1_5');load('I2_5');load('V2_5');

%CASO 15 V

% V1_15  = xlsread('efficiencies_dc_dc.xlsx',3,'A3:A9');
% I1_15  = xlsread('efficiencies_dc_dc.xlsx',3,'C3:C9');
% V2_15  = xlsread('efficiencies_dc_dc.xlsx',3,'D3:D9');
% I2_15  = xlsread('efficiencies_dc_dc.xlsx',3,'E3:E9');

load('I1_15');load('V1_15');load('I2_15');load('V2_15');

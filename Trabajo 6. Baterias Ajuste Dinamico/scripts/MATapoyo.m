clc;close all;clear all

%% DESCARGAS DINAMICAS

% Importación y carga de datos

%I_D  = xlsread('medidas_bateria.xlsx',1,'B2:B7213');
%V_D  = xlsread('medidas_bateria.xlsx',1,'C2:C7213');
%t_D  = xlsread('medidas_bateria.xlsx',1,'A2:A7213');
%phi  = xlsread('medidas_bateria.xlsx',1,'D2:D7213');
%phi2 = xlsread('medidas_bateria.xlsx',1,'E2:E7213');

load('I_D') ; load('V_D') ; load('t_D') ; load('phi'); load('phi2'); 
load('resultados_modelo_1C1R');load('resultados_modelo_2c2r');load('V_analitico');

% MODELO ANALITICO

% Division por tramos

I_1 = I_D(1:1775);
I_2 = I_D(1776:3535);
I_3 = I_D(3551:4734);
I_4 = I_D(4735:5326);
I_5 = I_D(5327:6509);
I_6 = I_D(6510:end);

V_1 = V_D(1:1775);
V_2 = V_D(1776:3535);
V_2 = - V_2 + 2*max(V_2);
V_3 = V_D(3551:4734);
V_4 = V_D(4735:5326);
V_5 = V_D(5327:6509);
V_5 = - V_5 + 2*max(V_5);
V_6 = V_D(6510:end);

t_1 = t_D(1:1775);
t_2 = t_D(1776:3535) - t_1(end);
t_3 = t_D(3551:4734) - t_1(end) - t_2(end);
t_4 = t_D(4735:5326) - t_1(end) - t_2(end) - t_3(end) ;
t_5 = t_D(5327:6509) - t_4(end) -  t_3(end) - t_2(end)  - t_1(end);
t_6 = t_D(6510:end)  - t_5(end) - t_4(end) - t_3(end) - t_2(end) - t_1(end);

% Filtrado

V2_rect = -0.0000701410.*t_2 + 23.9962177741;
V3_rect = -0.0000285504.*t_3 + 23.7911572283;
V4_rect = -0.0004534444.*t_4 + 23.2003273594;
V5_rect = -0.0000977002.*t_5 + 23.8667771768;
V6_rect = -0.0002121135.*t_6 + 23.3572374802;

V_2 = V_2 - V2_rect;
V_3 = V_3 - V3_rect;
V_4 = V_4 - V4_rect;
V_5 = V_5 - V5_rect;
V_6 = V_6 - V6_rect;

% Suavizado 

V_1 = smooth(V_1,180);
V_2 = smooth(V_2,180);
V_3 = smooth(V_3,180);
V_4 = smooth(V_4,180);
V_5 = smooth(V_5,180);
V_6 = smooth(V_6,180);

Dv_Di_2 = V_2./abs(I_1(300)-I_2(300));
Dv_Di_3 = V_3./abs(I_2(300)-I_3(300));
Dv_Di_4 = V_4./abs(I_3(300)-I_4(300));
Dv_Di_5 = V_5./abs(I_4(300)-I_5(300));
Dv_Di_6 = V_6./abs(I_5(300)-I_6(300));

% Calculo de las resistencias

R12 = Dv_Di_2(1);
R23 = Dv_Di_3(1);
R34 = Dv_Di_4(1);
R45 = Dv_Di_5(1);
R56 = Dv_Di_6(1);

R_1 = (R12 + R34 + R45 + R56) /4 ;

[a,Pos2] = min(abs(Dv_Di_2-0.368*max(Dv_Di_2)));
[b,Pos3] = min(abs(Dv_Di_3-0.368*max(Dv_Di_3)));
[c,Pos4] = min(abs(Dv_Di_4-0.368*max(Dv_Di_4)));
[d,Pos5] = min(abs(Dv_Di_5-0.368*max(Dv_Di_5)));
[e,Pos6] = min(abs(Dv_Di_6-0.368*max(Dv_Di_6)));

% Calculo de los condensadores

C12 = t_2(Pos2)/R12;
C23 = t_3(Pos3)/R23;
C34 = t_4(Pos4)/R34;
C45 = t_5(Pos5)/R45;
C56 = t_6(Pos6)/R56;

C_1 = (C12 + C34 + C45 + C56) /4 ;

% Expresión modelo simple 1C1R para cada tramo

deltaVdeltaI_12 = R12*exp(-(t_2./(R12*C12)));
deltaVdeltaI_23 = R23*exp(-(t_3./(R23*C23)));
deltaVdeltaI_34 = R34*exp(-(t_4./(R34*C34)));
deltaVdeltaI_45 = R45*exp(-(t_5./(R45*C45)));
deltaVdeltaI_56 = R56*exp(-(t_6./(R56*C56)));

% Calculo de las K modelo extendido

K12 = -1-log((1/R12)*Dv_Di_2(Pos2));
K23 = -1-log((1/R23)*Dv_Di_3(Pos3));
K34 = -1-log((1/R34)*Dv_Di_4(Pos4));
K45 = -1-log((1/R45)*Dv_Di_5(Pos5));
K56 = -1-log((1/R56)*Dv_Di_6(Pos6));

% Expresión modelo extendido 1C1R para cada tramo

deltaVdeltaI_12_exp = R12*exp(-(t_2./(R12*C12))+K12.*(t_2/(C12*R12)));
deltaVdeltaI_23_exp = R23*exp(-(t_3./(R23*C23))+K23.*(t_3/(C23*R23)));
deltaVdeltaI_34_exp = R34*exp(-(t_4./(R34*C34))+K34.*(t_4/(C34*R34)));
deltaVdeltaI_45_exp = R45*exp(-(t_5./(R45*C45))+K45.*(t_5/(C45*R45)));
deltaVdeltaI_56_exp = R56*exp(-(t_6./(R56*C56))+K56.*(t_6/(C56*R56)));

%PLOTS

figure()
box on
grid on
hold on
plot(t_2,deltaVdeltaI_12,'k','LineWidth',1.2)
plot(t_2,deltaVdeltaI_12_exp,'--k','LineWidth',1.2)
plot(t_2, Dv_Di_2, '--ok','MarkerIndices',1:100:length(Dv_Di_2),'LineWidth',1.2)
set(gca,'FontSize',18)
xlabel('{\it t} [s]')
%ylabel('\Delta{\it  V} /_\Delta{\it I} [V/A]')
ylabel('$\frac{\Delta{\it  V}}{\Delta{\it I}}$ [V/A]','Interpreter','latex')
legend('Ajuste 1','Ajuste 2','Datos experimentales filtrados - 2º Tramo')
hold off

figure()
box on
grid on
hold on
plot(t_3,deltaVdeltaI_23,'k','LineWidth',1.2)
plot(t_3,deltaVdeltaI_23_exp,'--k','LineWidth',1.2)
plot(t_3, Dv_Di_3,'--ok','MarkerIndices',1:100:length(Dv_Di_3),'LineWidth',1.2)
set(gca,'FontSize',18)
xlabel('{\it t} [s]')
% ylabel('{\it \Delta V/ \Delta I} [V/A]')
ylabel('$\frac{\Delta{\it  V}}{\Delta{\it I}}$ [V/A]','Interpreter','latex')
legend('Ajuste 1','Ajuste 2','Datos experimentales filtrados - 3^{er} Tramo')
hold off

figure()
box on
grid on
hold on
plot(t_4,deltaVdeltaI_34,'k','LineWidth',1.2)
plot(t_4,deltaVdeltaI_34_exp,'--k','LineWidth',1.2)
plot(t_4, Dv_Di_4,'--ok','MarkerIndices',1:100:length(Dv_Di_4),'LineWidth',1.2)
set(gca,'FontSize',18)
xlabel('{\it t} [s]')
% ylabel('{\it \Delta V/ \Delta I} [V/A]')
ylabel('$\frac{\Delta{\it  V}}{\Delta{\it I}}$ [V/A]','Interpreter','latex')
legend('Ajuste 1','Ajuste 2','Datos experimentales filtrados - 4º Tramo')
hold off

figure()
box on
grid on
hold on
plot(t_5,deltaVdeltaI_45,'k','LineWidth',1.2)
plot(t_5,deltaVdeltaI_45_exp,'--k','LineWidth',1.2)
plot(t_5, Dv_Di_5,'--ok','MarkerIndices',1:100:length(Dv_Di_5),'LineWidth',1.2)
set(gca,'FontSize',18)
xlabel('{\it t} [s]')
% ylabel('{\it \Delta V/ \Delta I} [V/A]')
ylabel('$\frac{\Delta{\it  V}}{\Delta{\it I}}$ [V/A]','Interpreter','latex')
legend('Ajuste 1','Ajuste 2','Datos experimentales filtrados - 5º Tramo')
hold off

figure()
box on
grid on
hold on
plot(t_6,deltaVdeltaI_56,'k','LineWidth',1.2)
plot(t_6,deltaVdeltaI_56_exp,'--k','LineWidth',1.2)
plot(t_6, Dv_Di_6, '--ok','MarkerIndices',1:100:length(Dv_Di_6),'LineWidth',1.2)
set(gca,'FontSize',18)
xlabel('{\it t} [s]')
% ylabel('{\it \Delta V/ \Delta I} [V/A]')
ylabel('$\frac{\Delta{\it  V}}{\Delta{\it I}}$ [V/A]','Interpreter','latex')
legend('Ajuste 1','Ajuste 2','Datos experimentales filtrados - 6º Tramo')
hold off

figure()
box on
grid on
hold on
plot(t_D,V_D,'k','LineWidth',1.2)
plot(V_analitico, 'd--k', 'MarkerIndices',1:120:7210, 'LineWidth',1.2)
xlabel('{\it t} [s]')
ylabel('{\it V} [V]')
legend('Datos experimentales' , 'Modelo 1RC analítico')
set(gca,'FontSize',18)
hold off

figure()
box on
grid on
hold on
plot(t_D,V_D,'k','LineWidth',1.2)
plot(Resultados_Modelo_1C1R,'d--k','MarkerIndices',1:120:7210, 'LineWidth',1.2)
xlabel('{\it t} [s]')
ylabel('{\it V} [V]')
legend('Datos experimentales' , 'Modelo 1RC')
set(gca,'FontSize',18)
hold off

figure()
box on
grid on
hold on
plot(t_D,V_D,'k','LineWidth',1.2)
plot(Modelo_2c2r, 'd--k', 'MarkerIndices',1:120:7210, 'LineWidth',1.2)
xlabel('{\it t} [s]')
ylabel('{\it V} [V]')
legend('Datos experimentales' , 'Modelo 2RC')
set(gca,'FontSize',18)
hold off

%% MODELO SIMULINK

% Descarga

R_d = 0.1315;
E_d0 = 24.3;
E_d1 = -3.547e-6;
E_d2_0 = 1.439e-9;
E_d2_1 = -1.603e-9;
E_d2_2 = 2.1092e-10;
E_d3_0 = 1.834e-5;
E_d3_1 = 2.019e-15;

%Carga
R_c = 0.0875;
E_c0 = 24.3;
E_c1 = 3.265e-6;
E_c2_0 = 5.1084e-14;
E_c3_1 = 2.458e-5;
E_c3_2 = 1.8e-7;
phi_d = phi + phi2 * R_d;
phi_c = phi - phi2 * R_c;


for i=1:size(I_D,1)

    if I_D(i)>0
    E(i)= E_d0 + E_d1*phi_d(i) + (E_d2_0 + E_d2_1*abs(I_D(i)) + E_d2_2*I_D(i)^2)*exp((E_d3_0 + E_d3_1*abs(I_D(i)))*phi_d(i));
    V(i) = E(i) - R_d*(I_D(i));

    else
    E(i)= E_c0 - E_c1*phi_c(i)-E_c2_0*exp((E_c3_1+E_c3_2*abs(I_D(i)))*phi_c(i));
    V(i) = E(i) + R_c*(abs(I_D(i)));

    end
end

E = E';

figure()
box on
grid on
hold on
plot(t_D,V,'--k','LineWidth',1.2)
plot(t_D,V_D,'k','LineWidth',1.2)
xlabel('{\it t} [s]')
ylabel('{\it V} [V]')
legend('Ajuste estático' , 'Datos experimentales')
set(gca,'FontSize',18)
hold off

%% VALORES PARA SIMULINK

I_D = - I_D;

E_D = E_d0 + E_d1.*phi_d + (E_d2_0 + E_d2_1.*I_D + E_d2_2.*I_D.^2).*exp((E_d3_0 + E_d3_1.*I_D).*phi_d);
E_C = E_c0 - E_c1.*phi_c - E_c2_0.*exp((E_c3_1+E_c3_2.*I_D).*phi_c);

E_CD = E_C - E_D;

Rint = 0.1315; %Rd
Rcd  = 0.0875; %Rc
R1   = 0.01;
R2   = 0.01;
C1   = 1000;
C2   = 1000;

figure()
box on
grid on
hold on
plot(t_D,E_C)
plot(t_D,E)
xlabel('{\it t} [s]')
ylabel('{\it V} [V]')
hold off



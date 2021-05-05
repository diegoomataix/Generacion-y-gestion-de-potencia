clc;close all;clear all

%% DESCARGAS DINAMICAS

% Importación y carga de datos

%I_D  = xlsread('medidas_bateria.xlsx',1,'B2:B7213');
%V_D  = xlsread('medidas_bateria.xlsx',1,'C2:C7213');
%t_D  = xlsread('medidas_bateria.xlsx',1,'A2:A7213');
%phi  = xlsread('medidas_bateria.xlsx',1,'D2:D7213');
%phi2 = xlsread('medidas_bateria.xlsx',1,'E2:E7213');

load('I_D') ; load('V_D') ; load('t_D') ; load('phi'); load('phi2');

% MODELO ANALITICO

I_1 = I_D(1:1774);
I_2 =  I_D(1775:3532);
I_3 = I_D(3533:4733);
I_4 = I_D(4734:5325);
I_5 = I_D(5326:6508);
I_6 = I_D(6509:end);

V_1 = V_D(1:1774);
V_2 = V_D(1775:3532);
V_3 = V_D(3533:4733);
V_4 = V_D(4734:5325);
V_5 = V_D(5326:6508);
V_6 = V_D(6509:end);

t_1 = t_D(1:1774);
t_2 = t_D(1775:3532) - t_1(end);
t_3 = t_D(3533:4733) - t_1(end) - t_2(end);
t_4 = t_D(4734:5325) - t_1(end) - t_2(end) - t_3(end) ;
t_5 = t_D(5326:6508) - t_4(end) -  t_3(end) - t_2(end)  - t_1(end);
t_6 = t_D(6509:end) -  t_5(end) - t_4(end) - t_3(end) - t_2(end) - t_1(end);

V1_rect = -0.00048441.*t_1 + 23.97875522;
V2_rect = 0.0000701410.*t_2 + 23.6394221899;
V3_rect = -0.00003619.*t_3 + 23.92519957;
V4_rect = -0.000461*t_4 + 25.384720;
V5_rect = 0.0000879957.*t_5 + 23.1738203106;
V6_rect = -0.00021592.*t_6 + 24.76426337

V_1 = (V_1 - V1_rect);
V_2 = (V_2 - V2_rect);
V_3 = (V_3 - V3_rect);
V_4 = (V_4 - V4_rect);
V_5 = V_5 - V5_rect;
V_6 = V_6 - V6_rect;

V_1 = V_1 + abs(min(V_1));
V_2 = abs(V_2 -max(V_2));
V_3 = V_3 + abs(min(V_3));
V_4 = V_4 + abs(min(V_4));
V_5 = abs(V_5 - max(V_5));
V_6 = V_6 + abs(min(V_6));

V_1 = smooth(V_1);
V_2 = smooth(V_2);
V_3 = smooth(V_3);
V_4 = smooth(V_4);
V_5 = smooth(V_5);
V_6 = smooth(V_6);

R12 = max(V_2)/abs(I_1(300)-I_2(300));
R13 = max(V_3)/abs(I_2(300)-I_3(300));
R14 = max(V_4)/abs(I_3(300)-I_4(300));
R15 = max(V_5)/abs(I_4(300)-I_5(300));
R16 = max(V_6)/abs(I_5(300)-I_6(300));

R_inc = (R12 + R13 + R14 + R15 + R16) /5 ;

Dv_Di_2 = V_2./(I_1(300)-I_2(300));
Dv_Di_3 = V_3./(I_2(300)-I_3(300));
Dv_Di_4 = V_4./(I_3(300)-I_4(300));
Dv_Di_5 = V_5./(I_4(300)-I_5(300));
Dv_Di_6 = V_6./(I_5(300)-I_6(300));

[a,Pos2] = min(abs(Dv_Di_2-0.368*max(Dv_Di_2)));
[b,Pos3] = min(abs(Dv_Di_3-0.368*max(Dv_Di_3)));
[c,Pos4] = min(abs(Dv_Di_4-0.368*max(Dv_Di_4)));
[d,Pos5] = min(abs(Dv_Di_5-0.368*max(Dv_Di_5)));
[e,Pos6] = min(abs(Dv_Di_6-0.368*max(Dv_Di_6)));

C12 = t_2(Pos2)/R12;
C13 = t_3(Pos3)/R13;
C14 = t_4(Pos4)/R14;
C15 = t_5(Pos5)/R15;
C16 = t_6(Pos6)/R16;

%plot_fun(t_2,Dv_Di_2,'','','')
%plot_fun(t_3,Dv_Di_3,'','','')
%plot_fun(t_4,Dv_Di_4,'','','')
%plot_fun(t_5,Dv_Di_5,'','','')
%plot_fun(t_6, Dv_Di_6,'','','')

%% MODELO SIMULINK

% Descarga

%R_d = 0.1315;
%E_d0 = 24.3;
%E_d1 = -3.547e-6;
%E_d2_0 = 1.439e-9;
%E_d2_1 = -1.603e-9;
%E_d2_2 = 2.1092e-10;
%E_d3_0 = 1.834e-5;
%E_d3_1 = 2.019e-15;
%Con phi2 correcto
R_d = [0.118667455470931];
E_d0 = [24.1112203383152];
E_d1 = [-3.51909965949880e-06];
E_d2_0 = [5.03123121278738e-10];
E_d2_1 = [-6.84502902638074e-10];
E_d2_2 = [9.60096310808297e-11];
E_d3_0 = [1.87750230665809e-05];
E_d3_1 = [-1.68795045318578e-13];


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
    %I_D(i)
    %E(i) = E_d0 + E_d1*phi(i);
    %V(i) = E(i) - R_d*abs(I_D(i));

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
plot(t_D,V)
plot(t_D,V_D)
xlabel('{\it t} [s]')
ylabel('{\it V} [V]')
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



%% BACKUP
% figure()
% box on
% grid on
% plot(t_1,V_1)
% figure()
% box on
% grid on
% plot(t_2,V_2)
% figure()
% box on
% grid on
% plot(t_3,V_3)
% figure()
% box on
% grid on
% plot(t_4,V_4)
% figure()
% box on
% grid on
% plot(t_5,V_5)
% figure()
% box on
% grid on
% plot(t_6,V_6)
%
 function plot_result = plot_fun(axisx,axisy,labelx,labely,yourtitle)
 figure()

 axis([min(axisx) max(axisx)  min(axisy)*0.95 max(axisy)*1.05])
 xlabel(labelx)
 ylabel(labely)
 title(yourtitle)
 box on
 set(gca,'FontSize',18)
 grid on
 grid minor
 hold on
 plot_result = plot(axisx,axisy,'k','LineWidth',2);
 end

clc; close all; clear all;

% CARGAR DATOS
% % {Tiempo , Intensidad , Voltaje} (en columnas)
%[descarga5A,~,~]=(xlsread('../ensayos_bateria.xlsx','descarga 5A'));
%[carga5A,~,~]=(xlsread('../ensayos_bateria.xlsx','carga 5A'));
%[descarga2_5A,~,~]=(xlsread('../ensayos_bateria.xlsx','descarga 2.5A'));
%[carga2_5A,~,~]=(xlsread('../ensayos_bateria.xlsx','carga 2.5A'));
%[descarga1_5A,~,~]=(xlsread('../ensayos_bateria.xlsx','descarga 1.5A'));
%[carga1_5A,~,~]=(xlsread('../ensayos_bateria.xlsx','carga 1.5A'));
 load('descarga5A'); load('carga5A');
 load('descarga2_5A'); load('carga2_5A');
 load('descarga1_5A');load('carga1_5A');

%% APARTADO 1
V_nom_i = 4.2;
C_nom_cell = 2750;
n_serie = round([(max(max(descarga5A(:,3))))/V_nom_i ; (max(max(descarga1_5A(:,3))))/V_nom_i ;  (max(max(descarga2_5A(:,3))))/V_nom_i ]);

t = [10615,21365,36528];    % [s]
I = [5, 2.5, 1.5];

k = log(t(3)/t(1))/ log(I(1)/I(3));

C_p = I.^k.*t;              % [A·s]
n_par = round(C_p ./ C_nom_cell);

%% APARTADO 2 - DESCARGA
R_d = 0.1472;               % [Ohm]

% %syms E_d_sym
% %E_d = double( vpasolve( V = E_d_sym - R_d .*I ) )
% U0 = [24, R_d];
% V_exp = [descarga5A(:, 3); descarga2_5A(:, 3); descarga1_5A(:, 3)];
% I = [descarga5A(:, 2); descarga2_5A(:, 2); descarga1_5A(:, 2)];
%
% [umin,fval]=fminsearch(@(u) RMSE_V(u,  V_exp, I, phi), U0);
%
% RMSE = fval

%%%%%%%%%%%%%%%%%%%%%%%Funciones
function error = RMSE_V(u, V_exp, I, phi)
global n_dat
for j=1:3
for i=1:size(phi,2)
	V_modelo(j, i) = lineal(u, I(j), phi(j, i));
end
error_vect(j) = (sum((V_modelo(j, :) - V_exp(j, :)).^2)/n_dat)^0.5;
end
error = error_vect(1)*17.1816 + error_vect(2)*8.5446 + error_vect(3)
end


function V_lineal = lineal(u, I, phi)
global Vt
%R_d = u(1) E_d1 = u(2) E_d2 = u(3)
E_d = E_d0 + E_d1*phi
V = E_d - R_d .*I
end
%% BACKUP
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


%C_test = [100, 1000, 3000] % Capacidades para sacar la primera Rd [A·s]

%for i = 1:3
%t(i, :) = (round(C_test(:)./I(i)))
%end
%t = round(t./5)
%t = t.*5

%ind5 = find(descarga5A(:, 1)==t(1, 1))
%ind25 = find(descarga2_5A(:, 1)==t(1, 1))

%R_d_test1 = (descarga5A(ind5,3)-descarga2_5A(ind25, 3))./(descarga5A(ind5,2)-descarga2_5A(ind25, 2))
%R_d_test2 = (descarga5A(t(3, :),3)-descarga1_5A(t(2, :), 3))./(descarga5A(t(3, :),2)-descarga1_5A(t(2, :), 2))

%Capacidad
% descarga5A(:,5) = descarga5A(:,1).*descarga5A(:,2);
% carga5A(:,5) = carga5A(:,1).*carga5A(:,2);
% descarga2_5A(:,5) = descarga2_5A(:,1).*descarga2_5A(:,2);
% carga2_5A(:,5) = carga2_5A(:,1).*carga2_5A(:,2);
% descarga1_5A(:,5) = descarga1_5A(:,1).*descarga1_5A(:,2);
% carga1_5A(:,5) = carga1_5A(:,1).*carga1_5A(:,2);

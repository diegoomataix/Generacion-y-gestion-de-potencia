clc; close all; clear all;

%% CARGAR DATOS
% {Tiempo , Intensidad , Voltaje} (en columnas)
load('descarga5A'); load('carga5A');
load('descarga2_5A'); load('carga2_5A');
load('descarga1_5A');load('carga1_5A');

global limites; global pesos; global n_dat; global caso;
global phi2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCARGA [1: lineal   % 2: exp. 1era aprox.   % 3: exp. completo ]
% CARGA    [4: lineal   % 5: exp. completo]
caso = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelo = 1;        % 1: descarga       % 2: carga
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PROCESAR DATOS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% APARTADO 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_nom_i = 4.2;
C_nom_cell = 2850*3600/1000; % [A.s]
n_serie = round([(max(max(descarga5A(:,3))))/V_nom_i ; (max(max(descarga1_5A(:,3))))/V_nom_i ;  (max(max(descarga2_5A(:,3))))/V_nom_i ]);

t = [10615,21365,36528];    % [s]
I = [5, 2.5, 1.5];

k = log(t(3)/t(1))/ log(I(1)/I(3));

C_p = I.^k.*t;              % [A·s]
n_par = round(C_p ./ C_nom_cell);

switch(modelo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% APARTADO 2 - DESCARGA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1

%% LINEAL
R_d_init = 0.1472;               % [Ohm]

% U0 = [R_d_init, 24, -1e-5 ]; % Primera Iteracion
U0 = [0.129898838224734,24.3041841517924,-3.99223729243347e-06]; % Segunda Iteracion
t_vect = [descarga5A(:, 1); descarga2_5A(:, 1); descarga1_5A(:, 1)];
I = [descarga5A(:, 2); descarga2_5A(:, 2); descarga1_5A(:, 2)];
V_exp = [descarga5A(:, 3); descarga2_5A(:, 3); descarga1_5A(:, 3)];
phi = [descarga5A(:,4);descarga2_5A(:,4);descarga1_5A(:,4)];
phi2 = [descarga5A(:,5);descarga2_5A(:,5);descarga1_5A(:,5)];


pesos = [17.1816, 8.5446, 1];
limites = [size(descarga5A,1), size(descarga2_5A,1), size(descarga1_5A,1)];
n_dat = limites;

% M =[ I1_5,phi1_5,phi2;
%	   I2_5,phi2_5,phi2;
%	   I5,phi5,phi2 ]
% M = [I(3),phi(3);]
switch(caso)
case 1
[umin,fval]=fminsearch(@(u) RMSE_V(u,  V_exp, I, phi), U0);

Params = umin
RMSE1 = fval

%% EXPONENCIAL MALA
case 2
U01 = [0.129898838224734,24.3041841517924,-3.99223729243347e-06, -1e-12, 3.35e-5];
[umin_exp1,fval_exp1]=fminsearch(@(u1) RMSE_V(u1,  V_exp, I, phi), U01);
Params_exp1 = umin_exp1
RMSE2 = fval_exp1

%% EXPONENCIAL BUENA
case 3
%U02 = [0.129898838224734,24.3041841517924,-3.99223729243347e-06, -1e-12, 3.35e-5 ,    ,    ,    ];
% valores del paper
U02 = [0.134, 24.369,-1.227e-2/3600, -1e-9, -4e-10 , 4.7e-11 , 5.8e-2/3600, -3.36e-10/3600];

[umin_exp2,fval_exp2]=fminsearch(@(u2) RMSE_V(u2,  V_exp, I, phi), U02);
Params_exp2 = umin_exp2
RMSE2 = fval_exp2


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% APARTADO 3 - CARGA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2
%% LINEAL
R_c_init = 0.1472;               % [Ohm]

U0 = [0.17547,20,3.6e-03/3600]; %
t_vect = [carga5A(:, 1); carga2_5A(:, 1); carga1_5A(:, 1)];
I = [carga5A(:, 2); carga2_5A(:, 2); carga1_5A(:, 2)];
V_exp = [carga5A(:, 3); carga2_5A(:, 3); carga1_5A(:, 3)];
phi = [carga5A(:,4);carga2_5A(:,4);carga1_5A(:,4)];

limites = [size(carga5A,1), size(carga2_5A,1), size(carga1_5A,1)];
lim = [1, limites(1), limites(1)+1, limites(1)+limites(2), limites(1)+limites(2)+1, limites(1)+limites(2)+limites(3)];

n_dat = limites;

switch(caso)
case 4
[umin,fval]=fminsearch(@(u) RMSE_V(u,  V_exp, I, phi), U0);

Params_cargal = umin
RMSE_cargal = fval

case 5
U1 = [-0.128,20.2091,3.3937e-6, 2.64e-14, 8.44e-2/3600, -3.89e-3/3600]; %
[umin,fval]=fminsearch(@(u) RMSE_V(u,  V_exp, I, phi), U1);

Params_cargaexp = umin
RMSE_carga2 = fval


end % END SWITCH CASE CASO CARGA

%%%%% END SWITCH CASE CARGA/DESCARGA
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lim = [1, limites(1), limites(1)+1, limites(1)+limites(2), limites(1)+limites(2)+1, limites(1)+limites(2)+limites(3)];
%phi_mat = [phi(lim(1):lim(2)), phi(lim(3):lim(4)), phi(lim(5):lim(6))];

figure()
hold on;grid on;box on
V_modelo = zeros(limites(1), 1);
for i=1:lim(2)
if caso == 1
   V_modelo(i) = lineal(umin,5, phi(i));                        % Parte Lineal DESCARGA
elseif caso == 2
   V_modelo(i) = V_exp1(Params_exp1,5,phi(i));             % Parte Exponencial Mala DESCARGA
elseif caso == 3
   V_modelo(i) = V_exp2(Params_exp2,5,phi(i),phi2(i));             % Parte Exponencial Buena  DESCARGA
elseif caso == 4
   V_modelo(i) = V_lineal_carga(Params_cargal,5,phi(i));             % Parte Lineal CARGA
elseif caso == 5                                                % exp Completo   CARGA
    V_modelo(i) = V_exp_carga(Params_cargaexp,5,phi(i));     %TBD
end
end
plot (phi(lim(1):lim(2))/3600, V_modelo, '--k');                    % Parte Lineal DESCARGA


V_modelo = zeros(limites(2), 1);
for i=lim(3):lim(4)
    if caso == 1
    V_modelo(i-lim(3)+1) = lineal(umin,2.5, phi(i));             % Parte Lineal DESCARGA
elseif caso == 2
    V_modelo(i-lim(3)+1) = V_exp1(Params_exp1,2.5,phi(i));  % Parte Exponencial Mala DESCARGA
elseif caso == 3
   V_modelo(i-lim(3)+1) = V_exp2(Params_exp2,2.5,phi(i),phi2(i));   % Parte Exponencial  Buena DESCARGA
elseif caso == 4
   V_modelo(i-lim(3)+1) = V_lineal_carga(Params_cargal,2.5,phi(i));
elseif caso == 5                                                    % exp Completo   CARGA
   V_modelo(i-lim(3)+1) = V_exp_carga(Params_cargaexp,2.5,phi(i));
end
end
plot (phi(lim(3):lim(4))/3600, V_modelo, '--k');



V_modelo = zeros(limites(3), 1);
for i=lim(5):(lim(6))
    if caso == 1
    V_modelo(i-lim(5)+1) = lineal(umin,1.5, phi(i));             % Parte Lineal
elseif caso == 2
   V_modelo(i-lim(5)+1) = V_exp1(Params_exp1,1.5,phi(i));
elseif caso == 3
    V_modelo(i-lim(5)+1) = V_exp2(Params_exp2,1.5,phi(i),phi2(i));    % Parte Exponencial
elseif caso == 4
    V_modelo(i-lim(5)+1) = V_lineal_carga(Params_cargal,1.5,phi(i));
elseif caso == 5
    V_modelo(i-lim(5)+1) = V_exp_carga(Params_cargaexp,1.5,phi(i)); %TBD
end
end
plot (phi(lim(5):lim(6))/3600, V_modelo, '--k');


if modelo == 1
    % Curvas Experimentales DESCARGA
plot (phi(1:limites(1))/3600,V_exp(1:limites(1)),'k');
plot (phi(limites(1)+1:limites(1)+limites(2))/3600,V_exp( limites(1)+1:limites(1)+limites(2)),'k');
plot (phi(limites(1)+limites(2)+1:limites(1)+limites(2)+limites(3))/3600,V_exp(limites(1)+limites(2)+1:limites(1)+limites(2)+limites(3)),'k');
elseif modelo == 2
    % Curvas experimentales de CARGA
plot (phi(1:limites(1))/3600,V_exp(1:limites(1)),'k');
plot (phi(limites(1)+1:limites(1)+limites(2))/3600,V_exp( limites(1)+1:limites(1)+limites(2)),'k');
plot (phi(limites(1)+limites(2)+1:limites(1)+limites(2)+limites(3))/3600,V_exp(limites(1)+limites(2)+1:limites(1)+limites(2)+limites(3)),'k');
end
ylabel('{\it V} [V]')
xlabel('{\it \phi} [W·h]');
%legend({'Aproximación numérica','Datos experimentales'},'Location','northeast','NumColumns',1)
axis([0 12e5/3600 15 25])
set(gca,'FontSize',18)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Funciones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function error = RMSE_V(u, V_exp, I, phi)
global n_dat; global pesos; global limites;
global caso
global phi2
V_modelo = zeros (1,sum(limites));
for i = 1:size(phi,1)
	if i <= limites(1)
		j = 1;
	elseif i <= limites(2) + limites(1)
		j = 2;
	else
		j = 3;
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if caso == 1
	V_modelo(i) = lineal(u, I(i), phi(i));
    elseif caso == 2
    V_modelo(i) = V_exp1(u, I(i), phi(i));
    elseif  caso == 3
    V_modelo(i) = V_exp2(u, I(i), phi(i),phi2(i));
    elseif caso == 4
    V_modelo(i) = V_lineal_carga(u, I(i), phi(i));
    elseif caso == 5
    V_modelo(i) = V_exp_carga(u, I(i), phi(i));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if i == limites(1)
		V_mod_err = V_modelo(1:limites(1));
		V_exp_err = V_exp(1:limites(1));
		error_vect(j) = ((sum (( V_mod_err' - V_exp_err).^2 ) / n_dat(j))^0.5);

	elseif i == limites(1) + limites(2)
		V_mod_err = V_modelo( limites(1)+1:limites(1)+limites(2));
		V_exp_err = V_exp( limites(1)+1:limites(1)+limites(2));
		error_vect(j) = ((sum (( V_mod_err' - V_exp_err).^2 ) / n_dat(j))^0.5);

	elseif  i == sum(limites)
		V_mod_err = V_modelo(limites(1)+limites(2)+1:limites(1)+limites(2)+limites(3));
		V_exp_err = V_exp(limites(1)+limites(2)+1:limites(1)+limites(2)+limites(3));
		error_vect(j) = (sum((V_mod_err' - V_exp_err).^2)/n_dat(j))^0.5;
     else
	end

end

error = error_vect(1) + error_vect(2)+ error_vect(3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V_lineal = lineal(u, I, phi)
global Vt
%R_d = u(1); E_d0 = u(2); E_d1 = u(3);
E_d = u(2) + u(3)*phi;
V_lineal = E_d - u(1)*I;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V_exp1 = V_exp1(u, I, phi)
global Vt
%R_d = u(1); E_d0 = u(2); E_d1 = u(3); E_d2=u(4); E_d3=u(5);
E_d = u(2) + u(3)*phi + u(4)*exp(u(5)*phi);
V_exp1 = E_d - u(1)*I;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V_exp2 = V_exp2(u, I, phi, phi2)
global Vt
%R_d = u(1); E_d0 = u(2); E_d1 = u(3); E_d2_0=u(4); E_d2_1=u(5); E_d2_2=u(6); E_d3_0=u(7); E_d3_1=u(8);
phi = phi + phi2 * u(1);
E_d = u(2) + u(3)*phi + (u(4) + u(5)*I + u(6)*I^2)*exp((u(7) + u(8)*I)*phi);
V_exp2 = E_d - u(1)*I;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V_lineal_carga = V_lineal_carga(u, I, phi)
global Vt
%R_c = u(1); E_c0 = u(2); E_c1 = u(3);
E_c = u(2) + u(3)*phi;
V_lineal_carga = E_c - u(1)*I;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V_exp_carga = V_exp_carga(u, I, phi)
global Vt
%R_c = u(1); E_c0 = u(2); E_c1 = u(3); E_c2_0=u(4); E_c3_1=u(5); E_3_2=u(6);
E_c = u(2) + u(3)*phi+u(4)*exp((u(5)+u(6)*I)*phi);
V_exp_carga = E_c - u(1)*I;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% BACKUP
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

% {Tiempo , Intensidad , Voltaje} (en columnas)
% [descarga5A,~,~]=(xlsread('ensayos_bateria.xlsx','descarga 5A'));
% [carga5A,~,~]=(xlsread('ensayos_bateria.xlsx','carga 5A'));
% [descarga2_5A,~,~]=(xlsread('ensayos_bateria.xlsx','descarga 2.5A'));
% [carga2_5A,~,~]=(xlsread('ensayos_bateria.xlsx','carga 2.5A'));
% [descarga1_5A,~,~]=(xlsread('ensayos_bateria.xlsx','descarga 1.5A'));
% [carga1_5A,~,~]=(xlsread('ensayos_bateria.xlsx','carga 1.5A'));

%test sin mover
% U_lin = [0.129898838224734,24.3041841517924,-3.99223729243347e-06];
% U01 = [-1e-17, 3.35e-5];
% [umin_exp1,fval_exp1]=fminsearch(@(u1) RMSE_exp1(u1, U_lin,  V_exp, I, phi), U01);

% Params_exp1 = [U_lin,umin_exp1]
% RMSE2 = fval_exp1

% function error = RMSE_exp1(u, V_exp, I, phi)
% global n_dat; global pesos; global limites;
% %V_modelo = zeros (1,sum(limites));
% for i = 1:size(phi,1)
% 	if i <= limites(1)
% 		j = 1;
% 	elseif i <= limites(2) + limites(1)
% 		j = 2;
% 	else
% 		j = 3;
% 	end
% 	V_modelo(i) = V_exp1(u, I(i), phi(i));
% 	if i == limites(1)
% 		V_mod_err = V_modelo(1:limites(1));
% 		V_exp_err = V_exp(1:limites(1));
% 		error_vect(j) = ((sum (( V_mod_err' - V_exp_err).^2 ) / n_dat(j))^0.5);
%
% 	elseif i == limites(1) + limites(2)
% 		V_mod_err = V_modelo( limites(1)+1:limites(1)+limites(2));
% 		V_exp_err = V_exp( limites(1)+1:limites(1)+limites(2));
% 		error_vect(j) = ((sum (( V_mod_err' - V_exp_err).^2 ) / n_dat(j))^0.5);
%
% 	elseif  i == sum(limites)
% 		V_mod_err = V_modelo(limites(1)+limites(2)+1:limites(1)+limites(2)+limites(3));
% 		V_exp_err = V_exp(limites(1)+limites(2)+1:limites(1)+limites(2)+limites(3));
% 		error_vect(j) = (sum((V_mod_err' - V_exp_err).^2)/n_dat(j))^0.5;
%      else
% 	end
%
% end
%
% error = error_vect(1) + error_vect(2)+ error_vect(3);
% end
% [descarga5A,~,~]=(xlsread('ensayos_bateria.xlsx','descarga 5A'));
% [carga5A,~,~]=(xlsread('ensayos_bateria.xlsx','carga 5A'));
% [descarga2_5A,~,~]=(xlsread('ensayos_bateria.xlsx','descarga 2.5A'));
% [carga2_5A,~,~]=(xlsread('ensayos_bateria.xlsx','carga 2.5A'));
% [descarga1_5A,~,~]=(xlsread('ensayos_bateria.xlsx','descarga 1.5A'));
% [carga1_5A,~,~]=(xlsread('ensayos_bateria.xlsx','carga 1.5A'));

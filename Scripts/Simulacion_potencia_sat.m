clear all; clc; close all

%% DATOS

h = [450e3 500e3 600e3];        % [m]
RT = 6378e3;                    % [m]
omega = [0.05 0.1 0.5];         % [rad/s]
eta = 0.3;                      % eficiencia paneles
fo = 0.8;                       % factor de ocupación
fg = (1 + sqrt(2)) / 2;         % factor geometrico
J2 =  1.0827e-3;                % J2 tierra
mu_T = 5.986e14;                % [m^3/s^2]
beta = 30;                      % Ángulo incidencia solar [deg]
G = 1360;                       % [W/m2]
Ap = 0.1*0.3;                   % [m^2]


%% PARAMETROS ORBITALES
% Órbita SS
a = h + RT;
n = sqrt(mu_T./a.^3);           % [rad/s]
e = 0;
w = 0;
T = 2*pi*sqrt(a.^3/mu_T);       % [s]

%% ACTITUD DEL SATELITE
pasoT=1000;

for i = 1:3
    t(i,:) = linspace(0,T(i),pasoT);
    alfa(i,:) = rad2deg(t(i,:)*n(i));             % [deg]
    for j = 1:3
        roll(i,j,:) = t(i,:) * omega(j);          % [rad]
    end
end

% Definición caras (desfase)

phi = deg2rad([0 90 180 270]);   %[rad]

%% ECLIPSES

rho = asind(RT./(RT+h));                % [deg]

alfa_eclipse_in = 180 - rho;            % [deg]
alfa_eclipse_out = 180 + rho;           % [deg]


for i = 1:3
    for j = 1:pasoT
        if alfa(i,j) >= alfa_eclipse_in(i) && alfa(i,j) <= alfa_eclipse_out(i)
           alfa(i,j) = 0; 
        end
    end
end

%% POTENCIA

% Componente Perpendicular
P_mP = G*eta*Ap*fo*fg*sind(beta)*(2/3); % [W]

% Componente Paralelo

for i = 1:4                    %caras
    for j = 1:3                %alturas
        for k = 1:3            %actitud
            for l = 1:pasoT    %tiempo
                P_mPa(i,j,k,l) = G*eta*Ap*fo*cosd(beta)*cos(roll(j,k,l)+phi(i))*sind(alfa(j,l)); % [W]
                if P_mPa(i,j,k,l)<0
                   P_mPa(i,j,k,l)=0;
                end
            end
        end
    end
end


P_mPa_total(:,:,:) = P_mPa(1,:,:,:)+P_mPa(2,:,:,:)+P_mPa(3,:,:,:)+P_mPa(4,:,:,:)+P_mP;

%Redefinición de alfa

for i=1:3
    alfa(i,:) = rad2deg(t(i,:)*n(i)); 
end

%Potencia media total nula durante eclipse
 
for i = 1:3
    for j = 1:pasoT
        if alfa(i,j) >= alfa_eclipse_in(i) && alfa(i,j) <= alfa_eclipse_out(i)
           P_mPa_total(i,:,j) = 0; 
        end
    end
end

%% FIGURAS

%Altura
height = 3;
switch(height)
    case 1
        alfaplot(:) = alfa(1,:);
        Pplot(:,i) = P_mPa(:,1,1,i);
        PtotalPLOT(:,:,i) = P_mPa_total(1,1,i);
    case 2
        alfaplot(:) = alfa(2,:);
        Pplot(:,i) = P_mPa(:,2,1,i);
        PtotalPLOT(:,:,i) = P_mPa_total(2,1,i);
    case 3
        alfaplot(:) = alfa(3,:);
        Pplot(:,i) = P_mPa(:,3,1,i);
        PtotalPLOT(:,:,i) = P_mPa_total(3,1,i);
end
        
%Actitud
act = 1;
switch(act)
    case 1
        Pplot(:,i) = P_mPa(:,height,1,i);
        PtotalPLOT(:,:,i) = P_mPa_total(height,1,i);
    case 2
        Pplot(:,i) = P_mPa(:,height,2,i);
        PtotalPLOT(:,:,i) = P_mPa_total(height,2,i);
    case 3
        Pplot(:,i) = P_mPa(:,height,3,i);
        PtotalPLOT(:,:,i) = P_mPa_total(height,3,i);
end

for i = 1:pasoT
    Pplot(:,i) = P_mPa(:,height,act,i);
end

for i = 1:pasoT
    PtotalPLOT(i) = P_mPa_total(height,act,i);
end

%Todas las caras

figure()
hold on
grid on
plot(alfaplot(:),Pplot(:,:))

%Cara X+ X-

figure()
hold on
grid on
plot(alfaplot(:),Pplot(1,:))
plot(alfaplot(:),Pplot(3,:))
hold off

%Cara Y+ Y-

figure()
hold on
grid on
plot(alfaplot(:),Pplot(2,:))
plot(alfaplot(:),Pplot(4,:))
hold off

%Potencia total

figure()
hold on
grid on
plot(alfaplot(:),PtotalPLOT(:))

%Potencia media generada por órbita

Pmedia = trapz(t(height,:),PtotalPLOT)/T(height);



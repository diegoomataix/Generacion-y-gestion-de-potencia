%% **********************************************************************************
%                      MÉTODOS EXPLÍCITOS - PANELES/CÉLULAS SOLARES
%____________________________________________________________________________________
clear all; clc
%% Escoger apartado
choose = 1;
%____________________________________________________________________________________
%% Puntos caracteristicos
% Orden de los datos de la matriz A:
% [ Isc, Imp, Vmp, Voc, betha, alpha, k ]   (cada uno es una fila de la matriz dat)

switch(choose)
    case 1  % Apartado 1
% cada variable en una fila en el orden indicado, cada columna corresponde a cada uno de las 8 células solares que se estudian
        dat = [0.760500000000000,0.523900000000000,0.462800000000000,0.520200000000000,1.03200000000000, 8.21000000000000, 0.503440000000000,7.55141000000000;
               0.689400000000000,0.496000000000000,0.438900000000000,0.504400000000000,0.925500000000000,7.61000000000000, 0.484760000000000,4.53786986400000;
               0.450700000000000,2.27000000000000, 2.41000000000000, 2.41100000000000, 12.4930000000000, 26.3000000000000, 12.0990000000000, 0.561760000000000;
               0.572700000000000,2.56500000000000, 2.72600000000000, 2.70000000000000, 16.7780000000000, 32.9000000000000, 13.5750000000000, 0.753649000000000;
               0.906508876000000,0.946745562000000,0.948357822000000,0.969627067000000,0.896802326000000,0.926918392000000,0.962895280000000,0.600930139000000;
               0.786973983000000,0.884990253000000,0.884079237000000,0.892962963000000,0.744606032000000,0.799392097000000,0.891270718000000,0.745386778000000;
               10.0367720400000, 27.6047703600000, 27.2574316300000, 30.4460176700000, 6.93744975500000, 11.0813335000000, 29.8260968800000, 9.33467154700000];
    case 2  % Apartado 2
% cada variable en una fila en el orden indicado, cada columna corresponde a cada uno de las 3 células solares que se estudian
        dat = ones(7); % tenemos que sacar los puntos críticos
end

% About lambertw:
%  For real x where −e^−1 < x < 0, the equation has exactly two real solutions. The larger
%  solution is represented by y = lambertW(x) and the smaller solution by y = lambertW(–1,x).

%____________________________________________________________________________________
%% Das model
% I / I_sc = ( 1 - ( V / V_oc ) ^k ) / ( 1 + h * ( V / V_oc ) )
k = zeros(1,size(dat,2));           % matriz de 1x(n de datos) para el coeficiente k del modelo de Das
h = zeros(1,size(dat,2));           % matriz de 1x(n de datos) para el coeficiente h del modelo de Das
k_func = zeros(1,size(dat,2));

% Determine coefficients
for i = 1:size(dat,2)               % iterar para cada uno de los paneles
    
    k_func(i) = lambertw(-1,dat(2,i)/dat(1,i)*log(dat(3,i)/dat(4,i)));
    k(i)=k_func(i) / log(dat(3,i)/dat(4,i) );  % W_func( Imp / Isc * log(Vmp/Voc) ) / log(Vmp/Voc)
    
    h(i) = (dat(4,i)/dat(3,i)) * (dat(1,i)/dat(2,i) - (1/k(i)) - 1);    % Voc / Vmp * (Isc/Imp) - 1/k -1)
end


V = zeros(size(dat,2),200);
for i = 1: size(dat,2)
V(i,:) = linspace(0,dat(4,i),200);
end 


I_das = zeros(size(dat,2),size(V,2));
for i = 1:size(dat,2)
    for j = 1: size(V,2)
I_das(i,j) = dat(1,i) .* ( ( 1 - ( V(i,j) ./ dat(4,i) ) .^k(i) ./ ( 1 + h(i) .* ( V(i,j) ./ dat(4,i) ) ) ) );
    end 
end 


for i = 1: size(dat,2)
figure(i)
hold on
grid on
box on
% axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
plot(  V(i,:) , I_das(i,:))

figure(i+size(dat,2))
hold on
grid on
box on
% axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
plot( V(i,:), I_das(i,:) .* V(i,:))

end 

%____________________________________________________________________________________
%% Karmalkar & Hannefa's model
%Ajuste mediante dos parametros gamma_k y m
% i = 1 - (1- gamma)* v - gamma * v^m         
% i y v son la intensidad y el voltaje en forma adimensional.  
gamma = zeros(1,size(dat,2));           % matriz de 1x(n de datos) para el coeficiente k del modelo de Das
gamma_simp = zeros(1,size(dat,2)); 
m = zeros(1,size(dat,2)); 
m_simp = zeros(1,size(dat,2)); 
m_func = zeros(1,size(dat,2)); 
C = zeros(1,size(dat,2)); 

% Determine coefficients
for i = 1:size(dat,2)   
    
    C(i) = (1- dat(5,i) - dat(6,i))/ ((2*dat(5,i))-1);
    
    m_func(i) = lambertw(-1, ((dat(6,i)* (C(i)^-1) * log(dat(6,i)))/C(i)));   %
    m(i)= (m_func(i)/log(dat(6,i))) + (C(i)^-1) + 1 ; % iterar para cada uno de los paneles
    
    gamma(i) = (2*dat(5,i) -1)/ ((dat(6,i)^m(i) *(m(i)-1)));    
    gamma_simp(i) = 1 + (1-dat(5,i))/dat(6,i);
    m_simp (i) = (log(1-dat(5,i))/log(dat(6,i)));
                       
end


I_Kar = zeros(size(dat,2),size(V,2));
I_Kar_simp = zeros(size(dat,2),size(V,2));
for i = 1:size(dat,2)
    for j = 1: size(V,2)
I_Kar(i,j) = dat(1,i) * ( 1 - (1- gamma(i))*(V(i,j)/dat(4,i)) - gamma(i)* ( (V(i,j)/dat(4,i))^m(i) ));
I_Kar_simp(i,j) = dat(1,i) * ( 1 - (1- gamma_simp(i))*(V(i,j)/dat(4,i)) - gamma_simp(i)* ( (V(i,j)/dat(4,i))^m_simp(i) ));
    end 
end 


for i = 1: size(dat,2)
figure()
hold on
grid on
box on
% axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
plot( V(i,:) , I_Kar(i,:))

figure()
hold on
grid on
box on
% axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
plot( V(i,:), I_Kar(i,:) .* V(i,:))


figure()
hold on
grid on
box on
% axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
plot( V(i,:) , I_Kar_simp(i,:))


figure()
hold on
grid on
box on
% axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
plot( V(i,:), I_Kar_simp(i,:) .* V(i,:))

end 


%% Pindado & Cubas's model
% Ajuste mediante un unico parametros eta
% I = Isc(1 -(1 - Imp/Isc)(V/Vmp)^(Imp/(Isc - Imp)))  para V<Vmp
% I = Imp*(Vmp/V)*(1 - ((V - Vmp)/(Voc - Vmp))^eta)   para V>Vmp
 
eta = zeros(1,size(dat,2));
Datos_paneles

for i=1:size(dat,2)
    eta(i)  = (dat(1,i)/dat(2,i))*(dat(1,i)/(dat(1,i) - dat(2,i)))*((dat(4,i) - dat(3,i))/dat(4,i));
    if i == 1 
        a = RTC(19,1);
        b = RTC(19,2);
    elseif i == 2
        a = TNJ(48,1);
        b = TNJ(48,2);
    elseif i == 3
        a = ZTJ(60,1);
        b = ZTJ(60,2);
    elseif i == 4
        a = G30C(940,1);
        b = G30C(940,2);
    elseif i == 5
        a = PWP(20,1);
        b = PWP(20,2);
    elseif i == 6
        a = KC2(70,1);
        b = KC2(70,2);
    elseif i == 7
        a = SPV(1130,1);
        b = SPV(1130,2);
    elseif i == 8
        a = PSC(17,1);
        b = PSC(17,2);
    end      
    eta_ba(i) = (log(dat(2,i)*dat(3,i)-a*b) - log(dat(2,i)*dat(3,i)))/(log(a - dat(3,i)) - log(dat(4,i) - dat(3,i)));
end


%____

% V = linspace(0,30,8)                                                    % V [V]

% Current
% I(:) = dat(1,:) .* ( ( 1 - ( V ./ dat(4,:) ) .^k ) ./ ( 1 + h .* ( V ./ dat(4,:) ) ) )

% for i = length(V)
%     I(:) = dat(1,:) * ( ( 1 - ( V(i) / dat(4,:) ) ^k(:) ) / ( 1 + h(:) * ( V(i) / dat(4,:) ) ) )
% end
%______________________________________________________________________________________
%% PLOTS
%%% Lambert W function %%%
% figure()
% syms x
% fplot(lambertw(x))
% hold on
% fplot(lambertw(-1,x))
% hold off
% axis([-0.5 4 -4 2])
% title('Lambert W function, two main branches')
% legend('k=0','k=1','Location','best')
%
% figure()
% syms x y
% f = lambertw(x + 1i*y);
% interval = [-100 100 -100 100];
% fmesh(real(f),interval,'ShowContours','On')

%%% Da's model %%%

%%% Karmalkar & Haneefa’s model %%%

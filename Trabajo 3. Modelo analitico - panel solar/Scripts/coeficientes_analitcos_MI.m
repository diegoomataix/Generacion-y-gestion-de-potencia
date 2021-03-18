%% **********************************************************************************
%                      MÉTODOS IMPLÍCITOS - PANELES/CÉLULAS SOLARES
%____________________________________________________________________________________
clear all; clc; close all
Datos_paneles

%% Escoger apartado
choose = 2;         % 1: Apartado 1      2: Apartado 2
%____________________________________________________________________________________
%% Puntos caracteristicos
% Orden de los datos de la matriz A:
% [ Isc, Imp, Vmp, Voc, betha, alpha, k ]   (cada uno es una fila de la matriz dat)
global V_oc
global I_sc
global V_mp
global I_mp
global n
global Vt 

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
        dat = [0.467533021000000,0.473000000000000,8.89000000000000;
               0.459346650000000,0.454000000000000,8.18000000000000;
               9.81821900000000,2.31000000000000,31.2000000000000;
               10.9998740000000,2.61000000000000,37.8000000000000;
               0.892576000000000,0.885057000000000,0.825397000000000;
               0.982490000000000,0.959831000000000,0.920135000000000]; % tenemos que sacar los puntos críticos
end

V = zeros(size(dat,2),200);     % Inicializar Voltaje
for i = 1: size(dat,2)
    V(i,:) = linspace(0,dat(4,i),200);
end

% About lambertw:
%  For real x where −e^−1 < x < 0, the equation has exactly two real solutions. The larger
%  solution is represented by y = lambertW(x) and the smaller solution by y = lambertW(–1,x).
%____________________________________________________________________________________

%% MODELO DE CIRCUITO EQUIVALENTE 1D2R

% Determinacion de los parametros [ Ipv, a, I0, Rs, Rsh ]

k = 1.3806503e-23;            % Boltzmann [J/K]
q = 1.60217646e-19;           % Electron charge [C]
T = 273.15+28;                % Temperature of the Essay [K]
n = round(dat(4,:)/2.5,0);    % Number of cells in the Panel (x3 if the cell is triple junction)
% n(1) = 1;
% n(8) = 1;
Vt = n*k*T/q;                 % Thermal Voltage   
a  = 1.35;

%Calculo de parametros 1D2R

for i = 1:size(dat,2)
    [Ipv(i),I0(i),Rs(i),Rsh(i)] = param_1D_2R(dat(1,i),dat(4,i),dat(2,i),dat(3,i),a,Vt(i));
end

%Calculo I-V mediante funcion de Lambert W

for i = 1:size(dat,2)     % bucle para los distintos paneles
    for j = 1:size(V,2)   % bucle para el voltaje      
    I(i,j) = (Rsh(i)*(Ipv(i) + I0(i)) - V(i,j))/(Rs(i) + Rsh(i)) - ((a*Vt(i))/Rs(i))*...
        lambertw(0,(((Rs(i)*Rsh(i)*I0(i))/((a*Vt(i))*(Rs(i) + Rsh(i))))*exp((Rsh(i)*...
        (Rs(i)*Ipv(i) + Rs(i)*I0(i) + V(i,j)))/(a*Vt(i)*(Rs(i) + Rsh(i))))));  
    end
end

% for i = 1: size(dat,2)
% 
%     % Plot I-V
%     figure()
%     hold on
%     grid on
%     box on
%     % axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
%     plot(  V(i,:) , I(i,:),'-k','LineWidth',2)
%     switch(choose)
%         case 1
%             switch(i)
%                 case 1
%                     plot(RTC(:,1), RTC(:,2), ':k','LineWidth',2)
%                 case 2
%                     plot(TNJ(:,1), TNJ(:,2), ':k','LineWidth',2)
%                 case 3
%                     plot(ZTJ(:,1), ZTJ(:,2), ':k','LineWidth',2)
%                 case 4
%                     plot(G30C(:,1), G30C(:,2), ':k','LineWidth',2)
%                 case 5
%                     plot(PWP(:,1), PWP(:,2), ':k','LineWidth',2)
%                 case 6
%                     plot(KC2(:,1), KC2(:,2), ':k','LineWidth',2)
%                 case 7
%                     plot(SPV(:,1), SPV(:,2), ':k','LineWidth',2)
%                 case 8
%                     plot(PSC(:,1), PSC(:,2), ':k','LineWidth',2)
%             end
%             
%         case 2
%             switch(i)
%                 case 1
%                     plot(SIP(:,1), SIP(:,2), ':k','LineWidth',2)
%                 case 2
%                     plot(CTJ30(:,1), CTJ30(:,2), ':k','LineWidth',2)
%                 case 3
%                     plot(MLU(:,1), MLU(:,2), ':k','LineWidth',2)
%             end
%     end
%     
%     axis tight
%     axis([0 dat(4,i)*1.2 0 dat(1,i)*1.2])
%     xlabel('{\it V} [V]')
%     ylabel('{\it I} [A]');
%     legend({'Modelo completo','Resultados experimentales'},'Location','northeast','NumColumns',2)
%     box on
%     set(gca,'FontSize',18)
%     hold off
%     
%     % Plot P-V
%             figure()
%             hold on
%             grid on
%             box on
%             % axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
%             plot( V(i,:), I(i,:) .* V(i,:), '-k','LineWidth',1)
%             switch(choose)
%                 case 1
%             switch(i)
%                 case 1
%                     plot(RTC(:,1), RTC(:,1).* RTC(:,2), ':k','LineWidth',2)      % 'Color', '#494949'
%                 case 2
%                     plot(TNJ(:,1), TNJ(:,1) .* TNJ(:,2), ':k','LineWidth',2)
%                 case 3
%                     plot(ZTJ(:,1), ZTJ(:,1) .* ZTJ(:,2), ':k','LineWidth',2)
%                 case 4
%                     plot(G30C(:,1), G30C(:,1) .* G30C(:,2), ':k','LineWidth',2)
%                 case 5
%                     plot(PWP(:,1), PWP(:,1).* PWP(:,2), ':k','LineWidth',2)
%                 case 6
%                     plot(KC2(:,1), KC2(:,1).* KC2(:,2), ':k','LineWidth',2)
%                 case 7
%                     plot(SPV(:,1), SPV(:,1) .* SPV(:,2), ':k','LineWidth',2)
%                 case 8
%                     plot(PSC(:,1), PSC(:,1) .* PSC(:,2), ':k','LineWidth',2)
%             end
%             
%             case 2
%             switch(i)
%                 case 1
%                     plot(SIP(:,1), SIP(:,1).* SIP(:,2), ':k','LineWidth',2)      % 'Color', '#494949'
%                 case 2
%                     plot(CTJ30(:,1), CTJ30(:,1) .* CTJ30(:,2), ':k','LineWidth',2)
%                 case 3
%                     plot(MLU(:,1), MLU(:,1) .* MLU(:,2), ':k','LineWidth',2)
%             end
%             end
%             
%             axis tight
%             axis([0 dat(4,i)*1.2 0 dat(2,i)*dat(3,i)*1.2])
%             xlabel('{\it V} [V]')
%             ylabel('{\it P} [W]');
%             legend({'Modelo explícito','Resultados experimentales'},'Location','northeast','NumColumns',2)
%             box on
%             set(gca,'FontSize',18)
%             hold off
%         
% end

%% MODELO DE CIRCUITO EQUIVALENTE 2D2R (queda por revisar, da intensidades de mierda)

clear Rs Rsh Ipv 

% Determinacion de los parametros [Ipv, a1, a2, I01, I02, Rs, Rs0, Rsh, Rsh0]

% Estimar a2

a2 = 2;

% Obtener Rs

syms Rs

switch(choose)
    case 1
        Rsh0 = [-0.0263,-0.00107,-0.00091,-0.000153,-0.00198,-0.0084,-0.002642,-4.6366];
        Rs0  = [-8.9199,-2.92969,-2.44245,-2.23729,-0.3419,-2.3299,-0.864,-0.7202];
    case 2
        Rsh0 = [-0.0004673,-0.0131853,-0.0046559];
        Rs0  = [-0.57789,-1.6902,-2.899855];
end

Isc(:) = dat(1,:);
Imp(:) = dat(2,:);
Vmp(:) = dat(3,:);
Voc(:) = dat(4,:);

Rsvec = zeros(1,size(dat,2));

for i = 1:size(dat,2)
    
equ1(i) = log( ((Rsh0(i)*(Isc(i)-Imp(i))-Vmp(i)) - a2*Vt(i)*((Rsh0(i) - (Vmp(i)/Imp(i)))/ ((Vmp(i)/Imp(i)) - Rs))) / ((Rsh0(i)*Isc(i) - Voc(i))-a2*Vt(i)*((Rsh0(i)-Rs0(i))/(Rs0(i)-Rs)))) ...
- ( (Vmp(i) + Imp(i)*Rs - Voc(i))*(((Rsh0(i) - (Vmp(i)/Imp(i)))/((Vmp(i)/Imp(i))-Rs)) - ((Rsh0(i) -Rs0(i))/(Rs0(i)-Rs))*exp((Vmp(i) + Imp(i)*Rs-Voc(i))/(a2*Vt(i))))) ...
/ (((Rsh0(i)*(Isc(i)-Imp(i))-Vmp(i))) - (Rsh0(i)*Isc(i) - Voc(i))*exp((Vmp(i)+Imp(i)*Rs-Voc(i))/a2*Vt(i)));

S1 = vpasolve(equ1(i) == 0, Rs);

Rsvec(i) = S1;

end


Rsvec    = double(real(Rsvec));

% Obtener a1
   
a1 = (1./Vt).*((Rsh0.*(Isc-Imp)-Vmp)-(Rs0.*Isc-Voc).*exp((Vmp+Imp.*Rsvec-Voc)./a2.*Vt))...
.*((Rsh0-Vmp./Imp)./(Vmp./Imp-Rsvec)-(Rsh0-Rs0)./(Rs0-Rsvec).*exp((Vmp+Imp.*Rsvec-Voc)./(a2.*Vt))).^(-1);

% Obtener I01 e I02

I01 = (a1./ (a2-a1) ) .* exp(-Voc./(a1.*Vt) ) .* ( (a2 .* Vt .* (Rsh0 - Rs0) - (Rs0 - Rsvec) .* (Rsh0.*Isc - Voc) ) ./ (Rsh0 - Rsvec) .* (Rs0 - Rsvec) );

I02 = (a2./ (a1-a2) ) .* exp(-Voc./(a2.*Vt) ) .* ( (a1 .* Vt .* (Rsh0 - Rs0) - (Rs0 - Rsvec) .* (Rsh0.*Isc - Voc) ) ./ (Rsh0 - Rsvec) .* (Rs0 - Rsvec) );

Rsh = Rsh0 - Rsvec;

Ipv = ( (Rsh + Rsvec) ./ Rsh ) .* Isc;

% Obtencion curvas I-V

for i =1:size(dat,2)
    for j = 1:size(V,2)
        
    I2D2R(i,j) = modelo2D2R(Ipv(i),Rsvec(i),Rsh(i),I01(i),I02(i),a1(i),a2,V(i,j),Vt(i));
    
    end
end

%PLOTS



%% Funciones

%Calculo de parametros 1D2R

function [Ipv,I0,Rs,Rsh] = param_1D_2R(Isc,Voc,Imp,Vmp,a,Vt)

    A=-(2*Vmp-Voc)./(a*Vt)+(Vmp*Isc-Voc*Imp)/(Vmp*Isc+Voc*(Imp-Isc));
    B=-Vmp*(2*Imp-Isc)/(Vmp*Isc+Voc*(Imp-Isc));
    C=a*Vt/Imp;
    D=(Vmp-Voc)./(a*Vt);
    
    M1=0.3361;
    M2=-0.0042;
    M3=-0.0201;
    sigma = -1-log(-B)-A;
    Wn =-1-sigma -2/M1* (1-1./(1+M1*sqrt(sigma/2)./(1+M2*sigma.*exp(M3*sqrt(sigma)))) );
    
    Rs=C*(Wn-(D+A));
    Rsh=(Vmp-Imp*Rs)*(Vmp-Rs*(Isc-Imp)-a*Vt)/((Vmp-Imp*Rs)*(Isc-Imp)-a*Vt*Imp);
    Ipv=(Rsh+Rs)/Rsh*Isc;
    I0=((Rsh+Rs)/Rsh*Isc-Voc/Rsh)/(exp((Voc)/(a*Vt)));
    
end

% I-V modelo 2D2R

function I_2D2R = modelo2D2R(Ipv,Rs,Rsh,I01,I02,a1,a2,V,Vt)


 I_2D2R = fzero(@(I) Ipv-I01*(exp((V+Rs*I)/(Vt*a1))-1)-I02*(exp((V+Rs*I)/(Vt*a2))-1)...
     -(V+Rs*I)/Rsh -I, 0);
 
end






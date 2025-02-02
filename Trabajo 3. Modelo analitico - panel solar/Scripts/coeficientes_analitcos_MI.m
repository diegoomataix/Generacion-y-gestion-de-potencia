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
if choose == 1
    n(1) = 1;
    n(8) = 1;
end
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

%Calculo de errores
switch(choose)
    case 1
N_exp = [size(RTC,1),size(TNJ,1),size(ZTJ,1),size(G30C,1),size(PWP,1),size(KC2,1),size(SPV,1),size(PSC,1)];
    case 2
N_exp = [size(SIP,1),size(CTJ30,1),size(MLU,1)];    
end 
I_modelo =zeros(size(dat,2),max(N_exp));
for i = 1:size(dat,2) % bucle para los distintos paneles
    N = N_exp(i);
    switch(choose)
        case 1
            switch(i)
                case 1
                    V_exp = RTC(:,1)';
                case 2
                    V_exp = TNJ(:,1)';
                case 3
                    V_exp = ZTJ(:,1)';
                case 4
                    V_exp = G30C(:,1)';
                case 5
                    V_exp = PWP(:,1)';
                case 6
                    V_exp = KC2(:,1)';
                case 7
                    V_exp = SPV(:,1)';
                case 8
                    V_exp = PSC(:,1)';
            end
        case 2
            
            switch(i)
                case 1
                    V_exp = SIP(:,1)';
                case 2
                    V_exp = CTJ30(:,1)';
                case 3
                    V_exp = MLU(:,1)';
            end
    end
%     V_exp(:) = linspace(0,dat(4,i),N);
    for j = 1:N   % bucle para el voltaje      
    I_modelo(i,j) = (Rsh(i)*(Ipv(i) + I0(i)) - V_exp(1,j))/(Rs(i) + Rsh(i)) - ((a*Vt(i))/Rs(i))*...
        lambertw(0,(((Rs(i)*Rsh(i)*I0(i))/((a*Vt(i))*(Rs(i) + Rsh(i))))*exp((Rsh(i)*...
        (Rs(i)*Ipv(i) + Rs(i)*I0(i) + V_exp(1,j)))/(a*Vt(i)*(Rs(i) + Rsh(i))))));  
    end
    clear V_exp
end

for i = 1:size(dat,2)
    N = N_exp(i);
    switch(choose)
        case 1
    switch(i)
        case 1
            I_exp = RTC(:,2)';
        case 2 
            I_exp = TNJ(:,2)';
        case 3
            I_exp = ZTJ(:,2)';
         case 4
            I_exp = G30C(:,2)';
        case 5 
            I_exp = PWP(:,2)';
        case 6
            I_exp = KC2(:,2)';
        case 7
            I_exp = SPV(:,2)';
        case 8 
            I_exp = PSC(:,2)';
    end 
        case 2
            switch(i)
        case 1
            I_exp = SIP(:,2)';
        case 2 
            I_exp = CTJ30(:,2)';
        case 3
            I_exp = MLU(:,2)';
            end
    end  
     
error_1D2R(i) =  RECT(I_modelo(i,1:N),I_exp);
RMSE_adim_1D2R(i) = RMSE_NonDim(I_modelo(i,1:N),I_exp);
[RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)] = error_PbPfun(I_modelo(i,1:N),I_exp);
switch(choose)
    case 1
    switch(i)
        case 1
            RTC_error_PbP_1D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 2 
            TNJ_error_PbP_1D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 3
            ZTJ_error_PbP_1D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
         case 4
            G30C_error_PbP_1D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 5 
            PWP_error_PbP_1D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 6
            KC2_error_PbP_1D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 7
            SPV_error_PbP_1D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 8 
            PSC_error_PbP_1D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
    end
      
    case 2
        switch(i)
        case 1
            SIP_error_PbP_1D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 2 
            CTJ30_error_PbP_1D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 3
            MLU_error_PbP_1D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        end
end
    clear RMSE_PbP_Dim
    clear RMSE_PbP_NonDim
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

a2 = 2*ones(1,size(dat,2));

% Obtener Rs

syms Rs

% Rsh0=[17.7010883542898,65.3807801517350,106.624885483321,382.909402620218,262.452140725135,46.4871041953477,...
%     971.407411855694,0.167179474494445,1631.70348889855,128.377472772625,31.9527297038283]
% Rs0 = [0.0937153988151533,0.367407671238610,0.403976141862510,0.279445832931270,2.26795829130260,0.570310624828062,...
%     1.22746947680144,0.0215896843209166,1.52675317332618,0.372264808649448,0.385924181841292]

switch(choose)
    case 1
        Rsh0 = 1./[0.0263,0.00107,0.00091,0.000153,1/262.452140725135,1/46.4871041953477,1/971.407411855694,1/0.167179474494445];
        Rs0  = 1./[8.9199,2.92969,2.44245,2.23729,1/2.26795829130260,1/0.570310624828062,1/1.22746947680144,1/0.0215896843209166];
    case 2
        Rsh0 = 1./[0.0004673,0.0131853,0.0046559];
        Rs0  = 1./[0.57789,1.6902,2.899855];
end

Isc(:) = dat(1,:);
Imp(:) = dat(2,:);
Vmp(:) = dat(3,:);
Voc(:) = dat(4,:);

Rsvec = zeros(1,size(dat,2));

if choose == 1
    rs0 = [0.1 0.29817 0.323727530 0.05 2.5 55e-2 65.115e-2 0.01];
elseif choose == 2
    rs0 = [1 -8.58 0.07];
end

for i = 1:size(dat,2)
    
equ1(i) = log( ((Rsh0(i)*(Isc(i)-Imp(i))-Vmp(i)) - a2(i)*Vt(i)*((Rsh0(i) - (Vmp(i)/Imp(i)))...
    / ((Vmp(i)/Imp(i)) - Rs))) / ((Rsh0(i)*Isc(i) - Voc(i))-a2(i)*Vt(i)*((Rsh0(i)-Rs0(i))/(Rs0(i)-Rs)))) ...
    - ( (Vmp(i) + Imp(i)*Rs - Voc(i))*(((Rsh0(i) - (Vmp(i)/Imp(i)))/((Vmp(i)/Imp(i))-Rs)) - ...
    ((Rsh0(i) -Rs0(i))/(Rs0(i)-Rs))*exp((Vmp(i) + Imp(i)*Rs-Voc(i))/(a2(i)*Vt(i))))) ...
    / (((Rsh0(i)*(Isc(i)-Imp(i))-Vmp(i))) - (Rsh0(i)*Isc(i) - Voc(i))*exp((Vmp(i)+...
    Imp(i)*Rs-Voc(i))/(a2(i)*Vt(i))));

S1 = vpasolve(equ1(i) == 0, Rs, rs0(i));

Rsvec(i) = S1;

end



% for i = 1:size(dat,2)
%     Rsvec(i) = CalculoRs(Isc(i),Imp(i),Vmp(i),Voc(i),Rsh0(i),Rs0(i),a2(i),Vt(i));
% end
% 
% Rsvec    = real(Rsvec);

Rsvec    = double(real(Rsvec));


% Obtener a1
  
a1 = (1./Vt).*(((Rsh0.*(Isc-Imp)-Vmp)-(Rsh0.*Isc-Voc).*exp((Vmp+Imp.*Rsvec-Voc)./(a2.*Vt)))...
.*((Rsh0-Vmp./Imp)./(Vmp./Imp-Rsvec)-(Rsh0-Rs0)./(Rs0-Rsvec).*exp((Vmp+Imp.*Rsvec-Voc)./(a2.*Vt))).^(-1));


% for i=1:size(dat,2)
% a1(i) = (1/Vt(i))*((Rsh0(i)*(Isc(i)-Imp(i))-Vmp(i))-(Rs0(i)*Isc(i)-Voc(i))*exp((Vmp(i)+Imp(i)*Rsvec(i)-Voc(i))/(a2(i)*Vt(i))))...
% *((Rsh0(i)-Vmp(i)/Imp(i))/(Vmp(i)/Imp(i)-Rsvec(i))-(Rsh0(i)-Rs0(i))/(Rs0(i)-Rsvec(i))*exp((Vmp(i)+Imp(i)*Rsvec(i)-Voc(i))/(a2(i)*Vt(i))))^(-1);
% end

% Obtener I01 e I02

I01 = (a1./ (a2-a1) ) .* exp(-Voc./(a1.*Vt) ) .* ( (a2 .* Vt .* (Rsh0 - Rs0) - (Rs0 - Rsvec) .* (Rsh0.*Isc - Voc) ) ./ ((Rsh0 - Rsvec) .* (Rs0 - Rsvec)) );

I02 = (a2./ (a1-a2) ) .* exp(-Voc./(a2.*Vt) ) .* ( (a1 .* Vt .* (Rsh0 - Rs0) - (Rs0 - Rsvec) .* (Rsh0.*Isc - Voc) ) ./ ((Rsh0 - Rsvec) .* (Rs0 - Rsvec)) );

Rsh = Rsh0 - Rsvec;

Ipv = ( (Rsh + Rsvec) ./ Rsh ) .* Isc;

% Obtencion curvas I-V

% syms I_2D2R
% 
% 
%     for j = 1:size(V,2)
%         
%     equ2(j) = Ipv(1)-I01(1)*(exp((V(1,j)+Rsvec(1)*I_2D2R)/(Vt(1)*a1(1)))-1)-I02(1)*(exp((V(1,j)+Rsvec(1)*I_2D2R)...
%         /(Vt(1)*a2(1)))-1)-(V(1,j)+Rsvec(1)*I_2D2R)/Rsh(1) -I_2D2R;
%     
%     S2 = vpasolve(equ2(j) == 0, I_2D2R);
%     
%     I_2D2Rvect(j) = S2;
%     
%     end


%I_2D2Rvect    = double(real(S2));

for i = 1:size(dat,2)
    for j = 1:size(V,2)
        I2D2R(i,j) = modelo2D2R(Ipv(i),Rsvec(i),Rsh(i),I01(i),I02(i),a1(i),a2(i),V(i,j),Vt(i));
    end
end

%Calculo de errores

I_modelo_2D2R =zeros(size(dat,2),max(N_exp));
for i = 1:size(dat,2) % bucle para los distintos paneles
    N = N_exp(i);
    switch(choose)
        case 1
            switch(i)
                case 1
                    V_exp = RTC(:,1)';
                case 2
                    V_exp = TNJ(:,1)';
                case 3
                    V_exp = ZTJ(:,1)';
                case 4
                    V_exp = G30C(:,1)';
                case 5
                    V_exp = PWP(:,1)';
                case 6
                    V_exp = KC2(:,1)';
                case 7
                    V_exp = SPV(:,1)';
                case 8
                    V_exp = PSC(:,1)';
            end
        case 2
            
            switch(i)
                case 1
                    V_exp = SIP(:,1)';
                case 2
                    V_exp = CTJ30(:,1)';
                case 3
                    V_exp = MLU(:,1)';
            end
    end
    %     V_exp(:) = linspace(0,dat(4,i),N);
    for j = 1:N   % bucle para el voltaje
        I_modelo_2D2R(i,j) = modelo2D2R(Ipv(i),Rsvec(i),Rsh(i),I01(i),I02(i),a1(i),a2(i),V_exp(1,j),Vt(i));
    end
    clear V_exp
end

for i = 1:size(dat,2)
    N = N_exp(i);
    switch(choose)
        case 1
    switch(i)
        case 1
            I_exp = RTC(:,2)';
        case 2 
            I_exp = TNJ(:,2)';
        case 3
            I_exp = ZTJ(:,2)';
         case 4
            I_exp = G30C(:,2)';
        case 5 
            I_exp = PWP(:,2)';
        case 6
            I_exp = KC2(:,2)';
        case 7
            I_exp = SPV(:,2)';
        case 8 
            I_exp = PSC(:,2)';
    end 
        case 2
            switch(i)
        case 1
            I_exp = SIP(:,2)';
        case 2 
            I_exp = CTJ30(:,2)';
        case 3
            I_exp = MLU(:,2)';
            end
    end 
            
     
error_2D2R(i) =  RECT(I_modelo_2D2R(i,1:N),I_exp);
RMSE_adim_2D2R(i) = RMSE_NonDim(I_modelo_2D2R(i,1:N),I_exp);
[RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)] = error_PbPfun(I_modelo_2D2R(i,1:N),I_exp);
    
switch(choose)
    case 1 
    switch(i)
        case 1
            RTC_error_PbP_2D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 2 
            TNJ_error_PbP_2D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 3
            ZTJ_error_PbP_2D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
         case 4
            G30C_error_PbP_2D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 5
            PWP_error_PbP_2D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 6
            KC2_error_PbP_2D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 7
            SPV_error_PbP_2D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 8 
            PSC_error_PbP_2D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
    end
    case 2
        switch(i)
        case 1
            SIP_error_PbP_2D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 2 
            CTJ30_error_PbP_2D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        case 3
            MLU_error_PbP_2D2R = [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)];
        end 
end 

        
    clear RMSE_PbP_Dim
    clear RMSE_PbP_NonDim
end

%PLOTS

for i = 1: size(dat,2)

    % Plot I-V
    figure()
    hold on
    grid on
    box on
    %axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
    plot(  V(i,:) , I2D2R(i,:),'-k','LineWidth',2)
    switch(choose)
        case 1
            switch(i)
                case 1
                    plot(RTC(:,1), RTC(:,2), ':k','LineWidth',2)
                case 2
                    plot(TNJ(:,1), TNJ(:,2), ':k','LineWidth',2)
                case 3
                    plot(ZTJ(:,1), ZTJ(:,2), ':k','LineWidth',2)
                case 4
                    plot(G30C(:,1), G30C(:,2), ':k','LineWidth',2)
                case 5
                    plot(PWP(:,1), PWP(:,2), ':k','LineWidth',2)
                case 6
                    plot(KC2(:,1), KC2(:,2), ':k','LineWidth',2)
                case 7
                    plot(SPV(:,1), SPV(:,2), ':k','LineWidth',2)
                case 8
                    plot(PSC(:,1), PSC(:,2), ':k','LineWidth',2)
            end
            
        case 2
            switch(i)
                case 1
                    plot(SIP(:,1), SIP(:,2), ':k','LineWidth',2)
                case 2
                    plot(CTJ30(:,1), CTJ30(:,2), ':k','LineWidth',2)
                case 3
                    plot(MLU(:,1), MLU(:,2), ':k','LineWidth',2)
            end
    end
    
    axis tight
    axis([0 dat(4,i)*1.2 0 dat(1,i)*1.2])
    xlabel('{\it V} [V]')
    ylabel('{\it I} [A]');
    legend({'Modelo completo','Resultados experimentales'},'Location','northeast','NumColumns',2)
    box on
    set(gca,'FontSize',18)
    hold off
    
%     % Plot P-V
%             figure()
%             hold on
%             grid on
%             box on
%             % axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
%             plot( V(i,:), I2D2R(i,:) .* V(i,:), '-k','LineWidth',1)
% %             switch(choose)
% %                 case 1
% %             switch(i)
% %                 case 1
% %                     plot(RTC(:,1), RTC(:,1).* RTC(:,2), ':k','LineWidth',2)      % 'Color', '#494949'
% %                 case 2
% %                     plot(TNJ(:,1), TNJ(:,1) .* TNJ(:,2), ':k','LineWidth',2)
% %                 case 3
% %                     plot(ZTJ(:,1), ZTJ(:,1) .* ZTJ(:,2), ':k','LineWidth',2)
% %                 case 4
% %                     plot(G30C(:,1), G30C(:,1) .* G30C(:,2), ':k','LineWidth',2)
% %                 case 5
% %                     plot(PWP(:,1), PWP(:,1).* PWP(:,2), ':k','LineWidth',2)
% %                 case 6
% %                     plot(KC2(:,1), KC2(:,1).* KC2(:,2), ':k','LineWidth',2)
% %                 case 7
% %                     plot(SPV(:,1), SPV(:,1) .* SPV(:,2), ':k','LineWidth',2)
% %                 case 8
% %                     plot(PSC(:,1), PSC(:,1) .* PSC(:,2), ':k','LineWidth',2)
% %             end
% %             
% %             case 2
% %             switch(i)
% %                 case 1
% %                     plot(SIP(:,1), SIP(:,1).* SIP(:,2), ':k','LineWidth',2)      % 'Color', '#494949'
% %                 case 2
% %                     plot(CTJ30(:,1), CTJ30(:,1) .* CTJ30(:,2), ':k','LineWidth',2)
% %                 case 3
% %                     plot(MLU(:,1), MLU(:,1) .* MLU(:,2), ':k','LineWidth',2)
% %             end
% %             end
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
end


%% Plot de errores 

switch(choose)
    case 1
    figure()
    axis tight
    xlabel('{\it V} [V]')
    ylabel('{\it | I-I_{exp} | / I_{sc}}');
    box on
    set(gca,'FontSize',18)
    hold on
    plot(G30C(:,1),G30C_error_PbP_1D2R(:,2), '-dk','LineWidth',2, 'MarkerSize', 12,'MarkerIndices',1:50:size(G30C,1))
    plot(G30C(:,1),G30C_error_PbP_2D2R(:,2), '-ok','LineWidth',2, 'MarkerSize', 12,'MarkerIndices',1:50:size(G30C,1))
    legend('Modelo 1D2R', 'Modelo 2D2R')
    
    
    figure()
    axis tight
    xlabel('{\it V} [V]')
    ylabel('{\it | I-I_{exp} | / I_{sc}}');
    box on
    set(gca,'FontSize',18)
    hold on
    plot(PSC(:,1),PSC_error_PbP_1D2R(:,2), '-dk','LineWidth',2, 'MarkerSize', 12,'MarkerIndices',1:2:size(PSC,1))
    plot(PSC(:,1),PSC_error_PbP_2D2R(:,2), '-ok','LineWidth',2, 'MarkerSize', 12,'MarkerIndices',1:2:size(PSC,1))
    legend('Modelo 1D2R', 'Modelo 2D2R')
    
    case 2
            figure()
    axis tight
    xlabel('{\it V} [V]')
    ylabel('{\it | I-I_{exp} | / I_{sc}}');
    box on
    set(gca,'FontSize',18)
    hold on
    plot(SIP(:,1),SIP_error_PbP_1D2R(:,2), '-dk','LineWidth',2, 'MarkerSize', 12,'MarkerIndices',1:2:size(SIP,1))
    plot(SIP(:,1),SIP_error_PbP_2D2R(:,2), '-ok','LineWidth',2, 'MarkerSize', 12,'MarkerIndices',1:2:size(SIP,1))
    legend('Modelo 1D2R', 'Modelo 2D2R')
    
    
    figure()
    axis tight
    xlabel('{\it V} [V]')
    ylabel('{\it | I-I_{exp} | / I_{sc}}');
    box on
    set(gca,'FontSize',18)
    hold on
    plot(MLU(:,1),MLU_error_PbP_1D2R(:,2), '-dk','LineWidth',1, 'MarkerSize', 8,'MarkerIndices',1:8:size(MLU,1))
    plot(MLU(:,1),MLU_error_PbP_2D2R(:,2), '-ok','LineWidth',1, 'MarkerSize', 8,'MarkerIndices',1:8:size(MLU,1))
    legend('Modelo 1D2R', 'Modelo 2D2R')
end
    
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
     -(V+Rs*I)/Rsh -I, 1);
 
end


% function cRss = CalculoRs(Isc,Imp,Vmp,Voc,Rsh0,Rs0,a2,Vt)
% 
% cRss = fzero(@(Rs) log( ((Rsh0*(Isc-Imp)-Vmp) - a2*Vt*((Rsh0 - (Vmp/Imp))...
%      / ((Vmp/Imp) - Rs))) / ((Rsh0*Isc - Voc)-a2*Vt*((Rsh0-Rs0)/(Rs0-Rs)))) ...
%      - ( (Vmp + Imp*Rs - Voc)*(((Rsh0 - (Vmp/Imp))/((Vmp/Imp)-Rs)) - ...
%      ((Rsh0 -Rs0)/(Rs0-Rs))*exp((Vmp + Imp*Rs-Voc)/(a2*Vt)))) ...
%      / (((Rsh0*(Isc-Imp)-Vmp)) - (Rsh0*Isc - Voc)*exp((Vmp+...
%      Imp*Rs-Voc)/(a2*Vt))),-0.5);
%  
% end


    
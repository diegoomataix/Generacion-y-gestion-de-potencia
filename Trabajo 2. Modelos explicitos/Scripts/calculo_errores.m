%% **********************************************************************************
%                      MÉTODOS EXPLÍCITOS - PANELES/CÉLULAS SOLARES
%____________________________________________________________________________________
clear all; clc; close all
Datos_paneles
global model
%% Escoger apartado
choose = 1;         % 1: Apartado 1      2: Apartado 2
model = 1;          % 1: Da's model      2: Karmalkar & Hannefa's model      3: Pindado & Cubas's model1
%____________________________________________________________________________________
%% Puntos caracteristicos
% Orden de los datos de la matriz A:
% [ Isc, Imp, Vmp, Voc, betha, alpha, k ]   (cada uno es una fila de la matriz dat)
global V_oc
global I_sc
global V_mp
global I_mp
global n


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
for model = 1:3
switch(model)

    %% Das model
    case 1
        % I / I_sc = ( 1 - ( V / V_oc ) ^k ) / ( 1 + h * ( V / V_oc ) )
        k = zeros(1,size(dat,2));           % matriz de 1x(n de datos) para el coeficiente k del modelo de Das
        h = zeros(1,size(dat,2));           % matriz de 1x(n de datos) para el coeficiente h del modelo de Das
        k_func = zeros(1,size(dat,2));

        %%% Determine coefficients %%%
        for i = 1:size(dat,2)               % iterar para cada uno de los paneles

            k_func(i) = lambertw(-1,dat(2,i)/dat(1,i)*log(dat(3,i)/dat(4,i)));
            k(i)=k_func(i) / log(dat(3,i)/dat(4,i) );  % W_func( Imp / Isc * log(Vmp/Voc) ) / log(Vmp/Voc)

            h(i) = (dat(4,i)/dat(3,i)) * (dat(1,i)/dat(2,i) - (1/k(i)) - 1);    % Voc / Vmp * (Isc/Imp) - 1/k -1)
        end

        %%% V y I %%%
        I_das = zeros(size(dat,2),size(V,2));
        for i = 1:size(dat,2)
            for j = 1: size(V,2)
                I_das(i,j) = dat(1,i) .* ( ( 1 - ( V(i,j) ./ dat(4,i) ) .^k(i) ./ ( 1 + h(i) .* ( V(i,j) ./ dat(4,i) ) ) ) );
            end
        end

        %%% PLOT %%%
%         for i = 1: size(dat,2)
% 
%             % Plot I-V
%             figure()
%             hold on
%             grid on
%             box on
%             % axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
%             plot(  V(i,:) , I_das(i,:),'-k','LineWidth',2)
%             switch(choose)
%                 case 1
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
%             case 2
%             switch(i)
%                 case 1
%                     plot(SIP(:,1), SIP(:,2), ':k','LineWidth',2)
%                 case 2
%                     plot(CTJ30(:,1), CTJ30(:,2), ':k','LineWidth',2)
%                 case 3
%                     plot(MLU(:,1), MLU(:,2), ':k','LineWidth',2)
%             end
%             end
%             
%             axis tight
%             axis([0 dat(4,i)*1.2 0 dat(1,i)*1.2])
%             xlabel('{\it V} [V]')
%             ylabel('{\it I} [A]');
%             legend({'Modelo completo','Resultados experimentales'},'Location','northeast','NumColumns',2)
%             box on
%             set(gca,'FontSize',18)
%             hold off
%             
%             % Plot P-V
%             figure()
%             hold on
%             grid on
%             box on
%             % axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
%             plot( V(i,:), I_das(i,:) .* V(i,:), '-k','LineWidth',1)
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
%         end
        %____________________________________________________________________________________
        %% Karmalkar & Hannefa's model
    case 2
        %Ajuste mediante dos parametros gamma_k y m
        % i = 1 - (1- gamma)* v - gamma * v^m
        % i y v son la intensidad y el voltaje en forma adimensional.
        gamma = zeros(1,size(dat,2));           % matriz de 1x(n de datos) para el coeficiente k del modelo de Das
        gamma_simp = zeros(1,size(dat,2));
        m = zeros(1,size(dat,2));
        m_simp = zeros(1,size(dat,2));
        m_func = zeros(1,size(dat,2));
        C = zeros(1,size(dat,2));

        %%% Determine coefficients %%%
        for i = 1:size(dat,2)

            C(i) = (1- dat(5,i) - dat(6,i))/ ((2*dat(5,i))-1);

            m_func(i) = lambertw(-1, ((dat(6,i)* (C(i)^-1) * log(dat(6,i)))/C(i)));   %
            m(i)= (m_func(i)/log(dat(6,i))) + (C(i)^-1) + 1 ; % iterar para cada uno de los paneles

            gamma(i) = (2*dat(5,i) -1)/ ((dat(6,i)^m(i) *(m(i)-1)));
            gamma_simp(i) = 1 + (1-dat(5,i))/dat(6,i);
            m_simp (i) = (log(1-dat(5,i))/log(dat(6,i)));

        end

        %%% I y V %%%
        I_Kar = zeros(size(dat,2),size(V,2));
        I_Kar_simp = zeros(size(dat,2),size(V,2));
        for i = 1:size(dat,2)
            for j = 1: size(V,2)
                I_Kar(i,j) = dat(1,i) * ( 1 - (1- gamma(i))*(V(i,j)/dat(4,i)) - gamma(i)* ( (V(i,j)/dat(4,i))^m(i) ));
                I_Kar_simp(i,j) = dat(1,i) * ( 1 - (1- gamma_simp(i))*(V(i,j)/dat(4,i)) - gamma_simp(i)* ( (V(i,j)/dat(4,i))^m_simp(i) ));
            end
        end

        %%% PLOT %%%
%         for i = 1: size(dat,2)
%             %%%% MODELO COMPLETO %%%%
%             % Plot I-V
%             figure()
%             hold on
%             grid on
%             box on
%             % axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
%             plot( V(i,:) , I_Kar(i,:), '--k','LineWidth',1)
%             plot( V(i,:) , I_Kar_simp(i,:), '-k','LineWidth',1)
%             switch(choose)
%                 case 1
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
%             case 2
%             switch(i)
%                 case 1
%                     plot(SIP(:,1), SIP(:,2), ':k','LineWidth',2)
%                 case 2
%                     plot(CTJ30(:,1), CTJ30(:,2), ':k','LineWidth',2)
%                 case 3
%                     plot(MLU(:,1), MLU(:,2), ':k','LineWidth',2)
%             end
%             end
%             
%             axis tight
%             axis([0 dat(4,i)*1.2 0 dat(1,i)*1.3])
%             xlabel('{\it V} [V]')
%             ylabel('{\it I} [A]');
%             legend({'Modelo completo','Modelo simplificado','Resultados experimentales'},'Location','northeast','NumColumns',2)
%             box on
%             set(gca,'FontSize',18)
%             hold off
% 
% %           Plot P-V
%             figure()
%             hold on
%             grid on
%             box on
%             % axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
%             plot( V(i,:), I_Kar(i,:) .* V(i,:), '-k','LineWidth',1)
%             plot( V(i,:), I_Kar_simp(i,:) .* V(i,:),'--k','LineWidth',1)
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
%             axis([0 dat(4,i)*1.2 0 dat(2,i)*dat(3,i)*1.5])
%             xlabel('{\it V} [V]')
%             ylabel('{\it P} [W]');
%             legend({'Modelo completo','Modelo simplificado','Resultados experimentales'},'Location','northeast','NumColumns',2)
%             box on
%             set(gca,'FontSize',18)
%             hold off
%         end

        %% Pindado & Cubas's model
    case 3
        % Ajuste mediante un unico parametros eta
        % I = Isc(1 -(1 - Imp/Isc)(V/Vmp)^(Imp/(Isc - Imp)))  para V<Vmp
        % I = Imp*(Vmp/V)*(1 - ((V - Vmp)/(Voc - Vmp))^eta)   para V>Vmp

        eta = zeros(1,size(dat,2));

        for i=1:size(dat,2)
            eta(i)  = (dat(1,i)/dat(2,i))*(dat(1,i)/(dat(1,i) - dat(2,i)))*((dat(4,i) - dat(3,i))/dat(4,i));
            
            switch(choose)
                case 1
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
                case 2
            if i == 1
                a = SIP(19,1);
                b = SIP(19,2);
            elseif i == 2
                a = CTJ30(100,1);
                b = CTJ30(100,2);
            elseif i == 3
                a = MLU(80,1);
                b = MLU(80,2);
            end
            end
            eta_ba(i) = (log(dat(2,i)*dat(3,i)-a*b) - log(dat(2,i)*dat(3,i)))/(log(a - dat(3,i)) - log(dat(4,i) - dat(3,i)));
        end

        %%% V y I %%%
        I_das = zeros(size(dat,2),size(V,2));
        for i = 1:size(dat,2)
            for j = 1: size(V,2)
                if V(i,j) < dat(3,i)
                    I_PC(i,j) = dat(1,i) * ( ( 1 - ( 1 - (dat(2,i) / dat(1,i) ) ) * (V(i,j)/dat(3,i))^((dat(2,i)/(dat(1,i) - dat(2,i))))));
                    I_PC_simp(i,j) = I_PC(i,j); % I = Isc(1 -(1 - Imp/Isc)(V/Vmp)^(Imp/(Isc - Imp)))  para V<Vmp
                    
                else
                    
                    I_PC(i,j) = dat(2,i) * (dat(3,i)/V(i,j)) * (1 - (( V(i,j) - dat(3,i)) / (dat(4,i) - dat(3,i) ))^eta_ba(i) );
                    I_PC_simp(i,j) = dat(2,i) * (dat(3,i)/V(i,j)) * (1 - (( V(i,j) - dat(3,i)) / (dat(4,i) - dat(3,i) ))^eta(i) );
                    % I = Imp*(Vmp/V)*(1 - ((V - Vmp)/(Voc - Vmp))^eta)   para V>Vmp
                end
            end
        end

        %%% PLOT %%%
%         for i = 1: size(dat,2)
%             %%%% MODELO COMPLETO %%%%
%             % Plot I-V
%             figure()
%             hold on
%             grid on
%             box on
%             % axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
%             plot( V(i,:) , I_PC(i,:), '--k','LineWidth',1)
%             plot( V(i,:) , I_PC_simp(i,:), '-k','LineWidth',1)
%             switch(choose)
%                 case 1
%             switch(i)
%                 case 1
%                     plot(RTC(:,1), RTC(:,2), ':k','LineWidth',2)      % 'Color', '#494949'
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
%                 case 2
%             switch(i)
%                 case 1
%                     plot(SIP(:,1), SIP(:,2), ':k','LineWidth',2)      % 'Color', '#494949'
%                 case 2
%                     plot(CTJ30(:,1), CTJ30(:,2), ':k','LineWidth',2)
%                 case 3
%                     plot(MLU(:,1), MLU(:,2), ':k','LineWidth',2)
%             end
%             end
%             
%             axis tight
%             axis([0 dat(4,i)*1.2 0 dat(1,i)*1.2])
%             xlabel('{\it V} [V]')
%             ylabel('{\it I} [A]');
%             legend({'Modelo completo','Modelo simplificado','Resultados experimentales'},'Location','northeast','NumColumns',2)
%             box on
%             set(gca,'FontSize',18)
%             hold off
% 
%             %Plot P-V
%             figure()
%             hold on
%             grid on
%             box on
%             % axis([V(i,1) V(i,end)  I_das(i,1)+1  0])
%             plot( V(i,:), I_PC(i,:) .* V(i,:), '-k','LineWidth',2)
%             plot( V(i,:), I_PC_simp(i,:) .* V(i,:), '--k','LineWidth',2)
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
%             legend({'Modelo completo','Modelo simplificado','Resultados experimentales'},'Location','northeast','NumColumns',2)
%             box on
%             set(gca,'FontSize',18)
%             hold off
%         end
end
end       %______________________________________________________________________________________
        %% PLOTS
        %%% Lambert W function %%%
%         figure()
%         syms x
%         fplot(lambertw(x))
%         hold on
%         fplot(lambertw(-1,x))
%         hold off
%         axis([-0.5 4 -4 2])
%         title('Lambert W function, two main branches')
%         legend('k=0','k=1','Location','best')
%         
%         figure()
%         syms x y
%         f = lambertw(x + 1i*y);
%         interval = [-100 100 -100 100];
%         fmesh(real(f),interval,'ShowContours','On')


    figure()
    axis tight
    xlabel('{\it V} [V]')
     ylabel('{\it | I-I_{exp} | / I_{sc}}');
    box on
    set(gca,'FontSize',18)
    hold on

for model = 1:3
    switch(model)
        case 1
        u = [k(8); h(8)];

        case 2
        u = [gamma(8); m(8)];    

        case 3
        u = eta(8);

    end



    switch(choose)
        case 1
    data = PSC;
  
        case 2
    % casos nuestros
  
    data = SIP;

    end
  
    
    V = data(:, 1);
    V = V';
    I_exp= data(:, 2);
    I_exp = I_exp';
    
    V_oc = max(V(:));
    I_sc = max(I_exp(:));
    
    P = I_exp.*V;
    [me, j] = max(P);
    V_mp = V(j);
    I_mp = I_exp(j);
    n = length(I_exp);

   
   if size(u,1) > 1
    error =  RECT(u,V,I_exp);
    RMSE_adim = RMSE_NonDim(u,V,I_exp);
    
    [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)] = error_PbPfun(u,V,I_exp);
   
   else 
    error =  RECT(u,V,I_exp);
    RMSE_adim = RMSE_NonDim(u,V,I_exp);
    
   [RMSE_PbP_Dim(:),RMSE_PbP_NonDim(:)] = error_PbPfun(u,V,I_exp);
   end 
   
    %Plot error en cada punto

    V_plot(:) = linspace(0,dat(4,i),size(I_exp,2));
%     plot(V_plot(:), RMSE_PbP_Dim(:),'k','LineWidth',2)

    
%     figure()
%     axis tight
%     xlabel('{\it V} [V]')
%     ylabel('{\it | I-I_{exp} | / I_{sc}}');
%     title('Error adimensional entre el modelo y las medidas experimentales')
%     box on
%     set(gca,'FontSize',18)
%     hold on
if model == 1
    plot(V_plot(:), RMSE_PbP_NonDim(:),'-sk','LineWidth',2, 'MarkerSize', 10)%,'MarkerIndices',1:20:length(V_plot))
elseif model == 2 
       plot(V_plot(:), RMSE_PbP_NonDim(:),'-ok','LineWidth',2,'MarkerSize', 10)%,'MarkerIndices',1:20:length(V_plot))
elseif model ==3
        plot(V_plot(:), RMSE_PbP_NonDim(:),'-dk','LineWidth',2,'MarkerSize', 10)%,'MarkerIndices',1:20:length(V_plot))
end 
legend('Das', 'Kalmalkar y Hannefa' , 'Pindado y Cubas')
        
  if model < 3
      stop  
    clear RMSE_PbP_Dim
    clear RMSE_PbP_NonDim
    clear V_plot
    clear u
  end

    %RMSE_PbP_NonDim(i,:) = error_PbPfun(u(:,i),V,I_exp);
  
end
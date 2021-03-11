function I_modelo = panel(u, V)
global V_oc
global I_sc
global V_mp
global I_mp

 %Ipv=u(1); I0=u(2); Rs=u(3); Rsh=u(4); a=u(5); 
 
 %I_modelo =fzero(@(I) u(1)-u(2)*(exp((V+u(3)*I)/(Vt*u(5)))-1)-(V+u(3)*I)/u(4)-I, 0);

% modelo Das
%I_modelo =  I_sc*(( 1 - ( V / V_oc )^u(1) ) / ( 1 + u(2) * ( V / V_oc )) );
    
%  Karmalkar & Hannefa's model
% gamma = u(1);
% m = u(2);
% I_modelo =  I_sc*(1 - (1- gamma)* (V/V_oc) - gamma * (V/V_oc)^m);

% Pindado & Cubas's model
eta = u(1);

if V<=V_mp
I_modelo = I_sc*(1 -(1 - I_mp/I_sc)*(V/V_mp)^(I_mp/(I_sc - I_mp)));
end

if V>=V_mp
I_modelo = I_mp*(V_mp/V)*(1 - ((V - V_mp)/(V_oc - V_mp))^eta); 
end


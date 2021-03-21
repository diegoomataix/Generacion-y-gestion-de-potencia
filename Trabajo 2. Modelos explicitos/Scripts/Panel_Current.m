function I_modelo = panel(u, V)
global V_oc
global I_sc
global V_mp
global I_mp
global model

switch(model)
case 1 
% modelo Das
I_modelo =  I_sc*(( 1 - ( V / V_oc )^u(1) ) / ( 1 + u(2) * ( V / V_oc )) );

case 2
%  Karmalkar & Hannefa's model
gamma = u(1);
m = u(2);
I_modelo =  I_sc*(1 - (1- gamma)* (V/V_oc) - gamma * (V/V_oc)^m);

case 3
% Pindado & Cubas's model
eta = u(1);

if V<=V_mp
I_modelo = I_sc*(1 -(1 - I_mp/I_sc)*(V/V_mp)^(I_mp/(I_sc - I_mp)));
end

if V>=V_mp
I_modelo = I_mp*(V_mp/V)*(1 - ((V - V_mp)/(V_oc - V_mp))^eta); 
end
end


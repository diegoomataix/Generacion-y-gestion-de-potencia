function I_modelo = panel(u, V)
global Vt
 
 %Ipv=u(1); I0=u(2); Rs=u(3); Rsh=u(4); a=u(5); 

 I_modelo = fzero(@(I) u(1)-u(2)*(exp((V+u(3)*I)/(Vt*u(5)))-1)-(V+u(3)*I)/u(4)-I, 0);



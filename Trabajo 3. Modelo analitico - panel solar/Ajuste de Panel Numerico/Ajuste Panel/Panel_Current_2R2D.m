function I_modelo = panel(u, V)
global Vt

%Ipv=u(1); I01=u(2); I02=u(3); Rs=u(4); Rsh=u(5); a1=u(6); a2=u(7);

I_modelo = fzero(@(I) u(1)-u(2)*(exp((V+u(4)*I)/(Vt*u(6)))-1)-u(3)*(exp((V+u(4)*I)/(Vt*u(7)))-1)...
 -(V+u(4)*I)/u(5) -I, 1);
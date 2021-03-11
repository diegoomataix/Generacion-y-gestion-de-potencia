function error = RECT(u, V, I_exp)
global n 
for i=1:size(V,2)
	I_modelo(i) = Panel_Current(u,V(i));
end
error = (sum((I_modelo - I_exp).^2)/n)^0.5;


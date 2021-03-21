function RMSE_adim = RMSE_NonDim(u, V, I_exp)
global n 
global I_sc
for i=1:size(V,2)
	I_modelo(i) = Panel_Current(u,V(i));
end
RMSE_adim = (1/I_sc)*((sum((I_modelo - I_exp).^2)/n)^0.5) %/sum(I_exp);

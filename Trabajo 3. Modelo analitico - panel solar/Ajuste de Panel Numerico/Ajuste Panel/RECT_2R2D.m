function error = RECT_2R2D(u, V, I_exp)
global n_dat
for i=1:size(V,2)
	I_modelo(i) = Panel_Current_2R2D(u,V(i));
end
error = (sum((I_modelo - I_exp).^2)/n_dat)^0.5;
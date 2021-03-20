function error = RECT_2R2D(u, V, I_exp)
for i=1:size(V,2)
	I_modelo(i) = Panel_Current_2R2D(u,V(i));
end
error = (sum((I_modelo - I_exp).^2))^0.5;
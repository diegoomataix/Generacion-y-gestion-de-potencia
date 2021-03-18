function error = ect_linea_2(u, x, z_exp)
for i=1:size(x,2)
	z_modelo(i) = linea(u,x(i));
end
error = (sum ((z_modelo - z_exp).^2))^0.5;


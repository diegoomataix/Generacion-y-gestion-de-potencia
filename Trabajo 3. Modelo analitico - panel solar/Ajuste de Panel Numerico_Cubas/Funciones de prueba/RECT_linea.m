function error = ect_linea(u, x, z_exp)
for i=1:size(x,2)
	z_modelo(i) = u(1) + u(2)*x(i);
end
error = (sum ((z_modelo - z_exp).^2))^0.5;


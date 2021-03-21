function [RMSE_PbP_Dim,RMSE_PbP_NonDim] = error_PbPfun(u, V, I_exp)
global I_sc
for i=1:size(V,2)
	I_modelo(i) = Panel_Current(u,V(i));
end
RMSE_PbP_Dim(:) = abs(I_modelo(:) - I_exp(:));
RMSE_PbP_NonDim(:) = (1/I_sc)*abs(I_modelo(:) - I_exp(:));
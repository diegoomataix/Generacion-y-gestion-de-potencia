function [RMSE_PbP_Dim,RMSE_PbP_NonDim,I_modelo] = error_PbPfun(u, V, I_exp)
global I_sc
for i=1:size(V,2)
	I_modelo(i) = Panel_Current(u,V(i));
    RMSE_PbP_Dim(i) = abs(I_modelo(i) - I_exp(i));
    RMSE_PbP_NonDim(i) = (1/I_sc)*abs(I_modelo(i) - I_exp(i));
end

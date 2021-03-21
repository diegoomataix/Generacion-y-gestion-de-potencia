function [RMSE_PbP_Dim,RMSE_PbP_NonDim] = error_PbPfun(I_modelo, I_exp)
I_sc = max(I_exp);
    RMSE_PbP_Dim(:) = abs(I_modelo(:) - I_exp(:));
    RMSE_PbP_NonDim(:) = (1/I_sc)*abs(I_modelo(:) - I_exp(:));
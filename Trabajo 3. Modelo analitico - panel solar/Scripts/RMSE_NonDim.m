function RMSE_adim = RMSE_NonDim(I_modelo, I_exp)
 
I_sc = max(I_exp);

RMSE_adim = (1/I_sc)*((sum((I_modelo - I_exp).^2)/size(I_exp,2))^0.5) %/sum(I_exp);

syms I_modelo I_sc V V_oc k h eta I_mp V_mp RMSE n I_exp

a = latex (RMSE == (sum((I_modelo - I_exp).^2)/n)^0.5)
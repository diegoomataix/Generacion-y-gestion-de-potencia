function error = RECT(I_modelo, I_exp)

error = ((sum((I_modelo - I_exp).^2)/size(I_exp,1))^0.5); %/sum(I_exp);


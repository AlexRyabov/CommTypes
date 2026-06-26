function [Res] = lme_results(lme)
Res.Rsquared =  lme.Rsquared.Adjusted;
Res.Res_Std =  sqrt(lme.MSE);

end
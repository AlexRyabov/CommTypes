function Posterior = ml_returnPosterior(Mdl,X)
    [~,~,~,Posterior] = predict(Mdl,X);   % Get posterior 
end
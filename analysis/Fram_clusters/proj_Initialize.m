%if ~exist('HGP', 'var')
    clearvars
    projMainVariables;
    projSetIDs;

    if ~exist('MaxFigure', 'var')
        MaxFigure = 201;  %this is the max figure number i used so far in this script
    end
    DataParams = DT_Params(DataID, root_data_folder);
    HGP = projParams(JobID, DataParams);
    [Data, DI] = dt_loadData(HGP, []);
    tSamplesNorm = Data.tSamplesNorm;
    tStPredictors = Data.tStPredictors;

    if HGP.Samples.BoxCox4NitriteAndAmmonium
        % the most General way to use the box-cox transformation (but only for not nans)
        %Box-Cox: trans_data = (data^lambda-1)/lambda;  %inverse transf (trans_data * lambda + 1).^(1/lambda) -> data
        notNaN = ~isnan(tStPredictors.nitrite_imputed);
        [tStPredictors.nitrite_imputed(notNaN),  HGP.Samples.BoxCox.lambda_nitrite] = boxcox(tStPredictors.nitrite_imputed(notNaN));
        notNaN = ~isnan(tStPredictors.ammonium_imputed);
        [tStPredictors.ammonium_imputed(notNaN), HGP.Samples.BoxCox.lambda_ammonium] = boxcox(tStPredictors.ammonium_imputed(notNaN));
        %Test variable should be transformed independently to avoid data leakage
    else
        HGP.Samples.BoxCox.lambda_nitrite = 1;
        HGP.Samples.BoxCox.lambda_ammonium = 1;
    end
%end


function indSelected = selectSpecies(Species, selectionParams, dim)
arguments
    Species 
    selectionParams 
    dim = 1  %dimension 1 or 2
end
% select species based on some criteria
% the matrix: samples (rows) x species (columns)
% selectionParams.Threshold threshold for selection 
% selectionParams.Type
% 'minMedian'  median abundance of species should be less than threshold
% 'minMedianInPresence'  median abundance when species is present should be less than threshold
% selectionParams.Type = 'nTop'  median abundance when species is present
% should be less than threshold which is provided selectionParams.nTop 
% selectionParams.Type =  'minOccup' the species should present in Threshold fraction of samples  
%
% example, select species which appear in >=5% of samples 
% selectSpecies(Species, struct('Type', 'minOccup', 'Threshold', 0.05))
% 

if dim == 2
    Species = Species';
end 
switch  selectionParams.Type
    case    'minMedian'
        %find species with median value higher than minAbund
        medianAbnd = median(Species, 1);
        indSelected  = medianAbnd >selectionParams.Threshold;
    case    'minMedianInPresence'   %minimal medium density when present
        %find species with median value higher than minAbund
        Species(Species==0) = NaN;
        medianAbnd = median(Species, 1, "omitmissing");
        indSelected  = medianAbnd >selectionParams.Threshold;
    case 'nTop'
        %find n species with top mean value
        medianAbnd = mean(Species, 1);
        B = maxk(medianAbnd,selectionParams.nTop); % returns the k largest elements of A.
        indSelected  = medianAbnd >= min(B);
    case  'minOccup'  %minimal occupancy
        Occupancy = sum(Species>0, 1)/size(Species, 1);
        indSelected  = Occupancy >selectionParams.Threshold;
end
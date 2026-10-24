function [Data, DI] = dt_loadData(HGP, Opt)
arguments
    HGP %% if necessary can later assign default value from HG_GetParams
    Opt = {'refresh_data'};  %refresh_data to reload all data from the source
end

%%list of phyla

tPhylaRoles = readtable([HGP.Data.PathMetaData, 'Phyla_roles_and_sizes.xlsx']);
tPhylaRoles = tPhylaRoles(tPhylaRoles.Active==1, :);

%% read data for samples
Samples = load(HGP.Data.fnSamplesProc);
tSamples = Samples.tSamples;
clearvars Samples

if ismember('DNA', tSamples.Properties.VariableNames)  %remove DNA column, from 2023 sample format
    tSamples.DNA = [];
end


if HGP.Samples.cleanSpeciesNames
    %remove suffixes and spaces
    taxa_list =    {'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'};
    for iT = 1:length(taxa_list)
        if ismember(taxa_list{iT}, tSamples.Properties.VariableNames)
            tSamples.(taxa_list{iT}) = string(cleanSpecNames(tSamples.(taxa_list{iT})));
        end
    end
end


%Replace Calanus_helgolandicus with Calanus spp
if HGP.Samples.renameCalanus
    ind = contains(tSamples.Species, 'Calanus_helgolandicus');
    tSamples.Species(ind) = 'Calanus_spp.';
end
%split Biecheleria
if HGP.Samples.separateBiecheleriaASV
    ind = strcmp(cleanSpecNames(tSamples.Species), "Biecheleria_sp.");
    tSamples(ind, {'ASV', 'Species'});
    tSamples.Species(ind) = strcat(tSamples.Species(ind), tSamples.ASV(ind));
end

if HGP.Samples.SplitSyndinialesDinophyceaeIsPhylum
  ind = strcmp(cleanSpecNames(tSamples.Class), "Syndiniales");
  tSamples.Phylum(ind) = "Syndiniales";
  ind = strcmp(cleanSpecNames(tSamples.Class), "Dinophyceae");
  tSamples.Phylum(ind) = "Dinophyceae";
end 

if HGP.Samples.FilterMetazoa
    %% exclude metazoa species
    indMetazoa = strcmp(tSamples.Phylum, 'p_ Metazoa') | strcmp(tSamples.Phylum, 'Metazoa');
    tSamples = tSamples(~indMetazoa, :);
    warning('%i metazoa AVS were removed from the sample list', sum(indMetazoa))
end

if HGP.Samples.removeUndefinedTaxaWhenGroupping
    %remove undefined species
    indUndef = strcmp(cleanSpecNames(tSamples.(HGP.Samples.GroupBy)), 'NA');
    tSamples = tSamples(~indUndef, :);
    if sum(indUndef)>0
        warning('%i NA AVS were removed from the sample list', sum(indUndef))
    end
end

if  HGP.Samples.RemoveNAGenus
    %remove species with undefined Genus
    indUndef = strcmp(cleanSpecNames(tSamples.('Genus')), 'NA');
    tSamples = tSamples(~indUndef, :);
    if sum(indUndef)>0
        warning('%i NA AVS without Genus were remove', sum(indUndef))
    end
end


if HGP.Samples.FillerMainPhyla
    ind = ismember(tSamples.Phylum, tPhylaRoles.Phylum);
    tSamples = tSamples(ind, :);
end

switch HGP.Samples.GroupBy
     case {'ASV'}
         tSamples.TopName = tSamples.ASV;
         tSamples.ASV = [];
         tSamples = movevars(tSamples,"TopName",'Before',"Kingdom");
          SmplStart =  find(strcmp(tSamples.Properties.VariableNames, 'Species'))+1;
         idxSamplCols = SmplStart:size(tSamples, 2);
    case {'Species'}
        if ismember('ASV', tSamples.Properties.VariableNames)  %remove ASV column
            tSamples.ASV = [];
        end
        tgrSamples = grpstats(tSamples, {'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'}, "sum");
        tgrSamples = f_TableColsCropStart(tgrSamples, 'sum_');
        tgrSamples.TopName = tgrSamples.Species;
        tgrSamples = movevars(tgrSamples,"TopName",'Before',"Kingdom");
        tgrSamples.GroupCount = [];
        tSamples = tgrSamples;
        SmplStart =  find(strcmp(tSamples.Properties.VariableNames, 'Species'))+1;
        idxSamplCols = SmplStart:size(tSamples, 2);

end
tSamples.Row = [];
% define the index of sample rows in tSamples

%remove rare species 
selectionParams = struct('Type', 'minOccup', 'Threshold', HGP.Samples.minOccurence); %select species which occur in this fraction of samples 
indSelectedSpecies1 = selectSpecies(tSamples{:, idxSamplCols}', selectionParams);

%1st normalization
dSamplesNormTmp = tSamples{:, idxSamplCols}./sum(tSamples{:, idxSamplCols});
selectionParams = struct('Type', 'minMedianInPresence', 'Threshold', HGP.Samples.minMedianInPresence); %median density when present 
indSelectedSpecies2 = selectSpecies(dSamplesNormTmp', selectionParams);

indSelectedSpecies = indSelectedSpecies1&indSelectedSpecies2;

tSamples = tSamples(indSelectedSpecies, :);
warning('%i rare taxons were removed from the sample list', sum(~indSelectedSpecies));
warning('The dataset have now %i taxons', sum(indSelectedSpecies));

%2nd normalization data (as some taxon could be removed above)
tSamplesNorm = tSamples;
tSamplesNorm{:, idxSamplCols} = tSamples{:, idxSamplCols}./sum(tSamples{:, idxSamplCols});
clearvars tSamples;  %Remove the original dataset to keep only the normalized tSamplesNorm


countTaxonsNumbers(tSamplesNorm);

%% load traits
tTraits = tSamplesNorm(:, {'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'});


%% sample metadata
%read metadata
tStPredictors = readtable([HGP.Data.fnSamplEnvImp]);
if ismember('date_time', tStPredictors.Properties.VariableNames)
    tStPredictors.date_time = [];
end

%% Sort metadata by year and month
[tStPredictors] = sortrows(tStPredictors, ["Year", "MonthNum", "doy", "Lat", "Lon", "Depth"]);

%% remove data with depth > 102 m
ind = (tStPredictors.Depth>102);
toRemove = sum(ind);
if toRemove >0
    warning("_____________________")
    warning("WARNING, DATA CLEANING")
    warning( "%i records with depth > 102 were removed", toRemove)
    warning("_____________________")
    tStPredictors = tStPredictors(~ind, :);
end

%remove old samples 
indSamples = tStPredictors.Year>2000; %include all samples  
tStPredictors = tStPredictors(indSamples, :);
tStPredictors.ID = (1:size(tStPredictors, 1))';
tStPredictors = movevars(tStPredictors,'ID','Before',1);


%resort species samples 
tSamplesNorm = [tSamplesNorm(:, 1:(SmplStart-1)), tSamplesNorm(:, tStPredictors.Sample)];
idxSamplCols = SmplStart:size(tSamplesNorm, 2);

%Add species richness and ESN
tStPredictors.divSimpsonESN = f_DiversityMetrics(tSamplesNorm{:, idxSamplCols}, [], 'SimpsonESN', [], 1)';
tStPredictors.divRichness = f_DiversityMetrics(tSamplesNorm{:, idxSamplCols}, [], 'Richness', [], 1)';


%% return
% sorted data
Data.tStPredictors = tStPredictors;  %table of predictors
Data.tSamplesNorm = tSamplesNorm;    %normalized samples (only species observed > 20 times)
Data.dSamplesNorm= tSamplesNorm{:, idxSamplCols};  %array of samples species x samples
Data.tTraits = tTraits;  %species traits (phyla, kingdom, ... )
Data.tPhylaRoles = tPhylaRoles; %phyla roles and sizes

%return data information
%sample column names
DI.SamplNames =  tSamplesNorm.Properties.VariableNames(idxSamplCols);
%sample column index
DI.idxSamplCols = idxSamplCols;
%number of species which we use for comparison
DI.nSp = size(tSamplesNorm, 1);

writenewtable(tSamplesNorm, [HGP.Data.Path2ResultsData '\tSamplesNorm.xlsx'])
writenewtable(tStPredictors, [HGP.Data.Path2ResultsData '\tStPredictors.xlsx'])


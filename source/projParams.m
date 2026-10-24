%Get a record a main Hausgarten parameters 
function HGP = projParams(JobID, DataParams)
arguments
    JobID
    DataParams
end

%Assign params of the data files
HGP = DataParams;


%set params specific fo JobID
HGP.JobID = JobID;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
HGP.Data.Path2ResultsData = fullfile(projectRoot, 'results', HGP.DataID, JobID);
if ~isfolder(HGP.Data.Path2ResultsData)
    mkdir(HGP.Data.Path2ResultsData);
end

%set default params
HGP.dm.Objects = 'Stations';
HGP.Samples.topGroupName = 'Phylum';
HGP.Samples.topGroupNameNMostDiverse = 15;  %how many most diverse taxa to select?
HGP.Samples.FilterMetazoa = true;
HGP.Samples.separateBiecheleriaASV = true;  %different species of genus 

HGP.Samples.FillerMainPhyla = true;  %include only main phyla

% Class name Syndiniales and Dinophyceae will go into Phylum name 
%to separate the diverse dinoflag group into two, when group by Phylum,
HGP.Samples.SplitSyndinialesDinophyceaeIsPhylum = true;
%Rename calanus
HGP.Samples.renameCalanus = true;
HGP.Samples.cleanSpeciesNames = false; %remove suffixes and spaces 
HGP.Samples.removeUndefinedTaxaWhenGroupping = false; %if true, NA species 
% will be removed, if group by species name, otherwise they will merged at
% higher taxonomic level
%set thresholds for filtering species abundance matrix 
if strcmp(JobID, 'ASVNoMetazoaJfewSp')
    HGP.Samples.minOccurence = 0.05; %the taxon (after grouping by HGP.Samples.GroupBy) should occur in at least minOccurence fraction of samples 
    HGP.Samples.minMedianInPresence = 1e-3; % Minimum median relative abundance of a taxon (when present) required for selection. 
    HGP.Samples.topGroupNameNMostDiverse = 10;  %how many most diverse taxa to select?
else
    HGP.Samples.minOccurence = 0.05; %the taxon (after grouping by HGP.Samples.GroupBy) should occur in at least minOccurence fraction of samples 
    HGP.Samples.minMedianInPresence = 1e-4; % Minimum median relative abundance of a taxon (when present) required for selection. 
end

HGP.Samples.RemoveNAGenus = strcmp(JobID, 'ASVNoMetazoaJRemoveNAGenus');  %remove all species where genus is NA

%Machine learning params 
%Transform columns 
HGP.Samples.BoxCox4NitriteAndAmmonium = false; % use box cox transformation for Nitrite And Ammonium as they are left biased 
HGP.ML.KFolds = 5; % use 5 fold for training models  
HGP.ML.CVBlocked = true; %if true, then cv partition will be blocked by cruises

%cluster parameters 
HGP.dmCluster.Method = 'gmm'; %Gaus clustering / other: rad -  radius clustering)
HGP.dmCluster.NumComponents = 8; %Number of clusters
HGP.dmCluster.rad_clust = 0.45; %radius for radius clustering 
HGP.dmCluster.DMCoords = 1:6; %diffusion map coordinates to find clusters

switch JobID
    case 'ASVNoMetazoaJHier'
        HGP.dmCluster.Method = 'hierarch'; %Hierarchical clustering
end


switch JobID
    case {'ASV'}
        HGP.Samples.GroupBy =  'ASV';
        HGP.Samples.FilterMetazoa = false;
    case {'ASVNoMetazoa', 'ASVNoMetazoaPears', ...
            'ASVNoMetazoaSERa', 'ASVNoMetazoaBC', ...
            'ASVNoMetazoaMC', 'ASVNoMetazoamcc', 'ASVNoMetazoaJ', 'ASVNoMetazoaJBlocked', 'ASVNoMetazoaJRemoveNAGenus', 'ASVNoMetazoaJfewSp', 'ASVNoMetazoaJHier',...
            'ASVNoMetazoaEq'}
        HGP.Samples.GroupBy =  'ASV';
        HGP.Samples.FilterMetazoa = true;
    case {'Species'}
        HGP.Samples.GroupBy =  'Species';
        HGP.Samples.FilterMetazoa = false;
    case {'SpeciesNoMetazoa'}
        HGP.Samples.GroupBy =  'Species';
        HGP.Samples.FilterMetazoa = true;
    otherwise
        error('JobID %s not found', JobID);
end


%set similarity metric
HGP.dm.BlockByYears = false;  %no blocking of similarity matrix by year

switch JobID
   case 'ASVNoMetazoaPears'
        HGP.dm.SimilarityMetric = 'pearson';
    case 'ASVNoMetazoaSERa' 
        HGP.dm.SimilarityMetric = 'SERa';
    case'ASVNoMetazoaBC'
        HGP.dm.SimilarityMetric = 'dist_BrayCurtis';
    case'ASVNoMetazoaMC'
        HGP.dm.SimilarityMetric = 'morisitaHorn';
    case 'ASVNoMetazoamcc'
        HGP.dm.SimilarityMetric = 'mcc';
    case {'ASVNoMetazoaJ', 'ASVNoMetazoaJRemoveNAGenus', 'ASVNoMetazoaJfewSp', 'ASVNoMetazoaJHier'}
        HGP.dm.SimilarityMetric = 'jaccard';
    case {'ASVNoMetazoaJBlocked'}
        HGP.dm.SimilarityMetric = 'jaccard';
        HGP.dm.BlockByYears = true; %blocking of similarity matrix by year
    case 'ASVNoMetazoaEq'
         HGP.dm.SimilarityMetric = 'euclidean';
   otherwise
    HGP.dm.SimilarityMetric = 'spearman';
end



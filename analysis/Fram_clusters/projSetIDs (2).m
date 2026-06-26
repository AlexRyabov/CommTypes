
DataID = 'HG_20092023';
%JobID = 'SpeciesNoMetazoa'; %dm of samples, remove Metazoa7
%JobID = 'ASVNoMetazoa'; %dm of samples, remove Metazoa7, Spearman for DM
%% Change the DM metric
% JobID = 'ASVNoMetazoaPears'; %dm of samples, remove Metazoa7, Pearson for dm
%JobID = 'ASVNoMetazoaSERa'; %dm of samples, remove Metazoa7, SERa for dm
% JobID = 'ASVNoMetazoaBC'; %dm of samples, remove Metazoa7, dist_BrayCurtis for dm
 %JobID = 'ASVNoMetazoaMC'; %dm of samples, remove Metazoa7, morisitaHorn for dm
% JobID = 'ASVNoMetazoamcc'; %dm of samples, remove Metazoa7, mcc for dm
 JobID = 'ASVNoMetazoaJ'; %dm of samples, remove Metazoa7, jaccard for dm
%JobID = 'ASVNoMetazoaJRemoveNAGenus'; %dm of samples, remove Metazoa, jaccard for dm and remove NA genus
% JobID = 'ASVNoMetazoaJfewSp'; %dm of samples, remove Metazoa7, jaccard for dm, use only 50 most abundant ASVs
% JobID = 'ASVNoMetazoaEq'; %dm of samples, remove Metazoa7, euclidean for dm

%% Use Jaccard metric but change clustering 
%JobID = 'ASVNoMetazoaJHier'; %dm of samples, remove Metazoa7, jaccard for dm, Hierarchical clustering

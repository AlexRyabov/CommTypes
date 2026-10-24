function HGP = DT_Params(DataID, root)
%DT_PARAMS Input paths for the manuscript analysis package.
% Paths for the bundled input data.
arguments
    DataID
    root
end

HGP.DataID = DataID;
HGP.Data.PathRoot = root;
HGP.Data.PathMetaData = [root filesep];
HGP.Data.Path2ProcData = [root filesep];
HGP.Data.fnSamplesProc = fullfile(root, 'HGSamples.mat');
HGP.Data.fnSamplEnvImp = fullfile(root, 'Samples_imp.xlsx');
HGP.Data.fnStatCoors = fullfile(root, 'Stations_coorsd.xlsx');
end

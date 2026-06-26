%Do diffusion map, clustering and community composition 

%Find diffusion map 
Fram_0100_DM;
pause
%Find clusters 
Fram_0250_Cluster
%pause
%Find communities and variability
Fram_0271_CommComposition_viol
  
%get last result
fprintf('%s\t %s \t %s\n', 'JobID', 'DivStd', 'AbndStd'); 
fprintf('%s\t%f\t%f\n', JobID, Stat.meanDivStd, Stat.meanAbndStd); 
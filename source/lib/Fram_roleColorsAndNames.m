function [roleClrs, roleShortNames] = Fram_roleColorsAndNames()
%FRAM_ROLECOLORSANDNAMES Shared ecological role colors and display names.

roleClrs = containers.Map;
roleClrs('Autotroph')   = [0.08 0.46 0.20];
roleClrs('Mixotroph')   =  [1.0000    0.7843    0.0980]; % [0.9608 0.5098 0.1882];
roleClrs('Zooplankton') = [0.8000 0.2000 0];
roleClrs('Parasite')    = [0.6667 0.4314 0.1569];
roleClrs('Recycler')    = [0.5686 0.1176 0.7059];
roleClrs('Unknown')     = [0.15 0.15 0.15];

roleShortNames = containers.Map;
roleShortNames('Autotroph')   = 'Autotroph';
roleShortNames('Mixotroph')   = 'Mixotroph';
roleShortNames('Zooplankton') = 'Zooplankton';
roleShortNames('Parasite')    = 'Parasite';
roleShortNames('Recycler')    = 'Recycler';
roleShortNames('Unknown')     = 'Unk';
end

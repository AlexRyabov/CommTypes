function f_Set_scales( Scales, varargin )
%% Alexey Ryabov 2014-2017
%example 
% f_Set_scales({'x', 'y'}) 
% f_Set_scales({'logx', 'logy'});
% f_Set_scales({'logx', 'logy', 'logz'});
% f_Set_scales({'x', 'y'}, xlabel, ylabel);
% f_Set_scales({'x', 'y'}, xlabel, ylabel, title);
      if sum(strcmp(Scales, 'logx'))
         set(gca, 'XScale','log')
      end
      if sum(strcmp(Scales, 'x'))
         set(gca, 'XScale','linear')
      end
      if sum(strcmp(Scales, 'logy'))
         set(gca, 'YScale','log')
      end
      if sum(strcmp(Scales, 'y'))
         set(gca, 'YScale','linear')
      end
      if sum(strcmp(Scales, 'logz'))
         set(gca, 'ZScale','log')
      end
      if sum(strcmp(Scales, 'z'))
         set(gca, 'ZScale','linear')
      end
      if isempty(varargin)
        f_FramePlot;
      else
        f_FramePlot(varargin);
      end

end


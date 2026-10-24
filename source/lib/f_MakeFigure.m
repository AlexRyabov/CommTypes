function [Fg] = f_MakeFigure( Fg, Position)
%create a figure, clear it and sets its position
Fg = figure(Fg); 
clf; 
ind = Position < 0;
Position(ind) = Fg.Position(ind);
set(Fg, 'Position', Position);
end


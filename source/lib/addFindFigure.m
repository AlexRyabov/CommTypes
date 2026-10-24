function [fg, iFgNrNew] = addFindFigure(PlotID, iFgNr, Position, WindowStyle, Visible, ToolBar)
% ADDFINDFIGURE  Create or retrieve a named figure and optionally adjust its properties.
%  This function is especially useful for creating docked plots in live
%  matlab scripts
%   [fg, iFgNrNew] = addFindFigure(PlotID, iFgNr) searches for an existing
%   figure with the name PlotID. If found, that figure is brought to the front.
%   Otherwise, a new figure is created. If iFgNr is positive, the figure is
%   created (or reused) with a specific numeric handle (iFgNr), and then
%   cleared before returning.
%
%   [fg, iFgNrNew] = addFindFigure(PlotID, iFgNr, Position, WindowStyle, Visible)
%   additionally allows specifying:
%       - Position    : The position of the figure window (e.g., [x, y, width, height]).
%       - WindowStyle : The window style (e.g., 'normal', 'docked', etc.).
%       - Visible     : The visibility of the figure (e.g., 'on', 'off').
%
%   Inputs:
%       PlotID      : (char or string) Unique name used to identify the figure.
%       iFgNr       : (scalar integer) Figure handle index. 
%                     - If iFgNr > 0, the function attempts to create or 
%                       retrieve a figure with handle iFgNr, and then increments 
%                       iFgNrNew by 1.
%                     - If iFgNr <= 0, the function creates (or retrieves) a 
%                       figure with an automatically assigned handle and sets 
%                       iFgNrNew to -1.
%       Position    : (1x4 numeric array, optional) [left, bottom, width, height].
%                     If omitted or empty, the figure’s current/ default 
%                     position is used.
%       WindowStyle : (char or string, optional) 'normal', 'docked', 'modal', etc.
%                     If omitted or empty, the figure’s current window style 
%                     remains unchanged.
%       Visible     : (char or string, optional) 'on' or 'off'.
%                     If omitted or empty, the figure’s current visibility 
%                     remains unchanged.
%       ToolBar     : type of toolbar 
%
%   Outputs:
%       fg      : Figure handle to the created or retrieved figure.
%       iFgNrNew: New figure index. If iFgNr was positive, this is iFgNr + 1.
%                 Otherwise, this is -1.
%
%   Example:
%       % Create a named figure 'MyPlot' with numeric handle 100, 
%       % positioned at [100, 100, 800, 600], undocked, and visible:
%       [fh, nextID] = addFindFigure('MyPlot', 100, [100, 100, 800, 600], 'normal', 'on');
%
%   See also: figure, findobj

% Define default arguments via the arguments block
arguments
    PlotID
    iFgNr
    Position = []
    WindowStyle = []  % 'docked', 'normal', ...
    Visible = []      % 'on', 'off'
    ToolBar = []  %'figure', 'auto', 'none'
end
alreadyOpen = false;

% If iFgNr > 0, create or find a figure with handle iFgNr
if iFgNr > 0
    % Check if the figure (by name) already exists
    fg = findobj('Type', 'Figure', 'Name', PlotID);
    if isempty(fg)
        fg = figure(iFgNr);
        fg.Name = PlotID;
    else
        figure(fg);  % Bring the existing figure to foreground
        alreadyOpen = true;
    end
    iFgNrNew = iFgNr + 1;
    clf;  % Clear the figure
else
    % If iFgNr <= 0, find or create a figure with an automatically assigned handle
    fg = findobj('Type', 'Figure', 'Name', PlotID);
    if isempty(fg)
        fg = figure;
        figure(fg);
        fg.Name = PlotID;
    else
        fg = fg(1);
        figure(fg)
        fg.Name = PlotID;
        alreadyOpen = true;
    end
    clf
    iFgNrNew = -1;
end

% If a position is specified, update the figure's position
if ~isempty(Position)
    if alreadyOpen
        %set only width and height
        fg.Position(3:4) = Position(3:4);
    else
        fg.Position = Position;
    end
end

% If a window style is specified, update it
if ~isempty(WindowStyle)
    fg.WindowStyle = WindowStyle;
end

% If visibility is specified, update it
if ~isempty(Visible)
    fg.Visible = Visible;
end

if ~isempty(ToolBar)
    fg.ToolBar = ToolBar;
end 
end

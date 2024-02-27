% Clear all variables and close all figures. restore all break points.

% Close all figures including those with hidden handles
close all hidden;

% Store all the currently set breakpoints in a variable
temporaryBreakpointData=dbstatus('-completenames');

% Clear functions and their persistent variables (also clears breakpoints 
% set in functions)
clear functions;

% Restore the previously set breakpoints
dbstop(temporaryBreakpointData);

% Clear global variables
clear global;

% Clear variables (including the temporary one used to store breakpoints)
clear variables;

% Clear command window
clc;
% HigherOrderClosure_v2025.m : 
% 
% Matlab version of the C version of the code I made of the Fortran version 
% of Meyers and Paw U Higher Order Closure model
% 
% 
%  Meyers TP + Paw U KT. 1987. Ag.Forest Met.
% 
% Meyers , T. P., and K. T. Paw U (1987), 
% Modeling the plant canopy micrometeorology with higher order 
% closure principles, Agricultural and Forest Meteorology, 41, 143-163.

%  
%  Dennis Baldocchi
%  Ecosystem Science Division
%  Dept of Environmental Science, Policy and Management
%  University of California, Berkeley
%  baldocchi@nature.berkeley.edu


% Feb 5, 2025

% Tilden gave us an updated version of Fortran in 2025

%  Dec. 2003 (From Tilden's 1987 Fortran version)
% 
%  Using code to compute turbulence statistics in sparse savanna canopy
% 
%  Nine equations are solved for
% 
%  <u>, <w'u'>,<u'u'>, <v'v'>, <w'w'>,<w'u'u'>, <w'v'v'>, <w'w'w'>, <w'w'u'>
% 
%  Fourth order terms are solved with the quasi-Gaussian approximation
% 
%  <w'w'w'w'>=3 <w'w'>^2
% 
%  <w'w'w'u'> = 3 <w'w'><w'u'>
% 
%  <w'u'u'u'> = 3 <u'u'> <w'u'>
% 
%  <w'w'u'u'> = <w'w'><u'u'> + 2 <w'u'><w'u'>
% 
%  <u'u'u'u'> = 3 <u'u'><u'u'>
% 
%  Computations are made for a normalized canopy with x layers

%  The grid domain of the canopy extends to three times canopy height times layers
% 
%  The turbulence profiles are also normalized in terms of friction velocity, u*
% 
%  the model assumes stationarity



HigherOrderClosure_v2025.m : 

Dennis Baldocchi
Dept of Environmental Science, Policy and Management
University of California, Berkeley
baldocchi@berkeley.edu


Feb 5, 2025

Matlab version of the C version of the code I made from the Fortran version 
of Meyers and Paw U Higher Order Closure model

Meyers , T. P., and K. T. Paw U (1987), 
Modeling the plant canopy micrometeorology with higher order 
closure principles, Agricultural and Forest Meteorology, 41, 143-163.

Tilden Meyers was my colleague at ATDD/NOAA where we worked together measuring and modeling
turbulence and turbulent exchange inside forests

Tilden gave me an updated version of Fortran in 2025, from which the most recent code has provenance too.

The main code is HigherOrderClosure_v2025.m

It calls a few function and inputs a profile of leaf area index.

Nine equations are solved for
 
<u>, <w'u'>,<u'u'>, <v'v'>, <w'w'>,<w'u'u'>, <w'v'v'>, <w'w'w'>, <w'w'u'>
 
Fourth order terms are solved with the quasi-Gaussian approximation
 
<w'w'w'w'>=3 <w'w'>^2

<w'w'w'u'> = 3 <w'w'><w'u'>
 
<w'u'u'u'> = 3 <u'u'> <w'u'>
 
<w'w'u'u'> = <w'w'><u'u'> + 2 <w'u'><w'u'>
 
<u'u'u'u'> = 3 <u'u'><u'u'>
 
Computations are made for a normalized canopy with n layers

The grid domain of the canopy extends to three times canopy height times layers
The turbulence profiles are also normalized in terms of friction velocity, u*
the model assumes stationarity

I have made a few modifications from the original code

I increased Cd to 0.2 and get a 2nd wind max!!

I refined the tests for convergence and interations. I am testing for changes in the canopy profile of turbulent kinetic
energy

Tilden's Cionco code U function was extinquishing u in the canopy too,
too fast. I needed to normalize the layers by the number of canopy layers

I am applying the Adams Bashford integration for U profiles.
It seems to work better than the approach in the FORTRAN code.
I am using trapezoid for layer 2, as it depends on 2, not 3 layers

I set the weight factor, wt, to 0.25 and I get some interesting profiles and lots of
iterations




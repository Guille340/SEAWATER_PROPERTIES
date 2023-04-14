
function z = seawaterPressure2Depth_BissetBerman1971(p,lat)

%**************************************************************************
%  z = seawaterPressure2Depth_LeroyParthiot1998(p,lat,loc)
%
%  DESCRIPTION: calculates depth from hydrostatic (gauge) pressure in 
%  seawater for all oceans and seas using the equation proposed by Bisset-
%  Berman (1971), as quoted in Leroy & Parthiot (1998).
%
%  Pressure to depth conversion is used in ocean engineering applications,
%  e.g. to calculate the depth of an autonomous underwater vehicle (AUV)
%  equipped with a pressure sensor. 
%
%  INPUT VARIABLES
%  - p: gauge pressure, above atmospheric (vector) [dbar] 
%    (1 dbar = 10 kPa corresponds to a depth increase in seawater of ~1 m)
%  - lat: latitude [deg]
%   
%   NOTE: if ¦lat¦ is not known, choose a value of ¦lat¦ = 45 deg
%        
%  OUTPUT VARIABLES
%  - z: depth below sea surface (vector) [m]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - The depth ¦z¦ is 1-5 m lower than that from Leroy & Parthiot (1998).
%  - The units used for the variables in the equations are: 
%    p [kg cm-2], z [m], lat [deg]
%
%  FUNCTION CALLS
%  1) z = seawaterPressure2Depth_BissetBerman1971(p,lat)
%
%  REFERENCES
%  - Bisset Berman (1971). ‘‘Instruction manual for salinity/temperature/
%    depth/sound velocity measuring systems models 9040’’, San Diego, CA.
%  - Leroy, C.C., and Parthiot, F. (1998). “Depth-pressure relationships 
%    in the oceans and seas”, J. Acoust. Soc. Am. 103, 1346–1352.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  29 Nov 2017
%
%**************************************************************************

p = p/9.80665; % % hydrostatic (gauge) pressure [dbar] 
z = (-2.07e-4*p + 9.7512 ./ (1 + 5.3e-3*sin(lat*pi/180).^2)).*p; % depth below sea surface [m]

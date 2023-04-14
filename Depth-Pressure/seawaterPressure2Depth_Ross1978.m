
function z = seawaterPressure2Depth_Ross1978(p)

%**************************************************************************
%  z = seawaterPressure2Depth_Ross1978(p)
%
%  DESCRIPTION: calculates depth from hydrostatic (gauge) pressure in 
%  seawater for all oceans and seas using the equation proposed by Ross
%  (1978). The formula is an inversion of p(z) equation from Leroy &  
%  Parthiot (1968) evaluated at mid-latitudes.
%
%  Pressure to depth conversion is used in ocean engineering applications, 
%  e.g. to calculate the depth of an autonomous underwater vehicle (AUV)
%  equipped with a pressure sensor. 
%
%  INPUT VARIABLES
%  - p: gauge pressure, above atmospheric (vector) [dbar] 
%    (1 dbar = 10 kPa corresponds to a depth increase in seawater of ~1 m)
%        
%  OUTPUT VARIABLES
%  - z: depth below sea surface (vector) [m]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - The equation assumes a latitude of 45 deg.
%  - Very limited accuracy.
%  - The units used for the variables in the equations are: 
%    p [kg cm-2], z [m]
%
%  FUNCTION CALLS
%  1) z = seawaterPressure2Depth_Ross1978(p)
%
%  REFERENCES
%  - Ross, D. (1978). Revised simplified formulae for calculating the 
%    speed of sound in sea water. Saclant ASW Research Centre, No. SM 107.
%  - Leroy, C.C. (1968). “Formulas for the calculation of underwater 
%    pressure in acoustics”, J. Acoust. Soc. Am. 44(2), 651-653
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  29 Nov 2017
%
%**************************************************************************

p = p/9.80665; % % hydrostatic (gauge) pressure [kg cm-2] 
z = (-2.2e-4*p +9.74).*p; % depth below sea surface [m]

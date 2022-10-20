
function p = seawaterDepth2Pressure_Leroy1988(z,lat)

%**************************************************************************
%  p = seawaterDepth2Pressure_Leroy1988(z,lat)
%
%  DESCRIPTION: calculates depth from hydrostatic (gauge) pressure in 
%  seawater for open oceans using the equation proposed by Leroy (1988),
%  as quoted in Leroy & Parthiot (1998). The equation is a compromise 
%  between those proposed by Leroy (1968) and Siess (1982).
%
%  Depth to pressure conversion is required for the calculation of the 
%  speed of sound in seawater from general empirical equations, such as 
%  the NRLII (DelGrosso, 1974), UNESCO (Chen & Millero, 1977; Fofonoff & 
%  Millard, 1983) or Wilson 1st and 2nd equations (Wilson, 1960b; Wilson, 
%  1960c). These empirical equations are fitted to precise laboratory 
%  measurements and use pressure as a parameter, instead of depth, since 
%  pressure is the measured magnitude. Profiles of salinity and temperature
%  are generally expressed as a function of depth, and a conversion from 
%  depth to pressure is then required.
%
%  INPUT VARIABLES
%  - z: depth below sea surface (vector) [m]
%  - lat: latitude [deg]
%
%   NOTE: if ¦lat¦ is not known, choose a value of ¦lat¦ = 45 deg.
%        
%  OUTPUT VARIABLES
%  - p: gauge pressure, above atmospheric (vector) [dbar] 
%    (1 dbar = 10 kPa corresponds to a depth increase in seawater of ~1 m)
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - The units used for the variables in the equations are: 
%    p [dbar], z [m], lat [deg]
%
%  FUNCTION CALLS
%  1) p = seawaterDepth2Pressure_Leroy1988(z,lat)
%
%  REFERENCES
%  - IPA Siess (1982). ‘‘Relation entre pression et immersion dans l’océan” 
%    (relationship between pressure and depth in the ocean), Annex to 
%    letter 293 EPSHOM/E/OC (in French).
%  - Leroy, C.C., and Parthiot, F. (1998). “Depth-pressure relationships 
%    in the oceans and seas”, J. Acoust. Soc. Am. 103, 1346–1352.
%  - Fofonoff, N.P., & Millard, R.C. (1983) “Algorithm for computation of 
%    fundamental properties of seawater” UNESCO Technical Papers in Marine 
%    Science, No. 44.
%  - Del Grosso, V. A. (1974). “New equation for the speed of sound in 
%    natural waters (with comparisons to other equations)” J. Acoust. Soc. 
%    Am. 56, 1084–1091.
%  - Chen, C.T., & Millero, F.J. (1977) “Sound speed in seawater at high 
%    pressures”, J. Acoust. Soc. Am. 62, 1129–1135.
%  - Wilson, W.D. (1960b). “Speed of sound in sea water as a function of
%    temperature, pressure, and salinity”, J. Acoust. Soc. Am. 32(6), 641-
%    644.
%  - Wilson, W.D. (1960c). “Equation for the speed of sound in sea water”, 
%    J. Acoust. Soc. Am. 32(10), 1357.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  29 Nov 2017
%
%**************************************************************************

p = 1.00534*(2.46e-6*z + 1 + 5.294e-3*sin(lat*pi/180)^2).*z;

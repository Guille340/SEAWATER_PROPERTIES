
function p = seawaterDepth2Pressure_Leroy1968(z,lat,loc)

%**************************************************************************
%  p = seawaterDepth2Pressure_Leroy1968(z,lat,loc)
%
%  DESCRIPTION: calculates depth from hydrostatic (gauge) pressure in 
%  seawater using the equations proposed by Ross (1978) for the main
%  oceans, Black Sea and Baltic Sea.
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
%  - loc: number specifying the geographic location (ocean or sea). 
%
%    loc    Area                   Latitude [deg]   
%  -------------------------------------------------
%    0      Common Oceans          40S - 60N        
%    1      Black Sea              41N - 47N        
%    2      Baltic Sea             53N - 63N    
%
%   NOTE: if ¦lat¦ and ¦loc¦ are not known, choose a value of ¦lat¦ = 45 
%   deg and ¦loc¦ = 0 for an approximate estimation of ¦z¦.
%        
%  OUTPUT VARIABLES
%  - p: gauge pressure, above atmospheric (vector) [dbar] 
%    (1 dbar = 10 kPa corresponds to a depth increase in seawater of ~1 m)
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - For ¦loc¦ = 1 (Black Sea) or ¦loc¦ = 2 (Baltic Sea), ¦lat¦ is ignored
%  - The units used for the variables in the equations are: 
%    p [kg cm-2], z [m], lat [deg]
%
%  FUNCTION CALLS
%  1) p = seawaterDepth2Pressure_Leroy1968(z,lat,loc)
%
%  REFERENCES
%  - Leroy, C.C. (1968). “Formulas for the calculation of underwater 
%    pressure in acoustics”, J. Acoust. Soc. Am. 44(2), 651-653.
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

switch loc
    case 0 % Common Oceans
        p = (2.524e-7*z + 0.102506*(1 + 5.28e-3*sin(lat*pi/180).^2)).*z + 1e-2;    
    case 1 % Black Sea
        p = (2.6e-7*z + 1.0168e-1).*z;
    case 2 % Baltic Sea
        p = (1.4e-6*z + 1.008e-1).*z;    
    otherwise
        error('Not a valid LOC number')  
end

p = p*9.80665; % hydrostatic (gauge) pressure [dbar]

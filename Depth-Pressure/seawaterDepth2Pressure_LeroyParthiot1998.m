
function p = seawaterDepth2Pressure_LeroyParthiot1998(z,lat,loc)

%**************************************************************************
%  p = seawaterDepth2Pressure_LeroyParthiot1998(z,lat,loc)
%
%  DESCRIPTION: calculates the hydrostatic (gauge) pressure from depth in 
%  seawater for all oceans and seas using the equation proposed by Leroy & 
%  Parthiot (1998). 
%
%  The current depth-pressure equation is based on the UNESCO algorithm 
%  (Fofonoff & Millard, 1983) and on calculations from temperature and 
%  salinity profiles. The equation leads to error in the sound speed lower
%  than +/- 0.02 m/s, compared to the +/- 0.5 to +/- 0.1 m/s of previous
%  simple equations. The improvement over simple equations is substantial, 
%  and those equations should be abandoned.
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
%    loc    Area                   Latitude [deg]   Std. Dev. [+/- dbar]
%  --------------------------------------------------------------------
%    0      Common Oceans          40S - 60N        0.8
%    1      North East Atlantic    30N - 35N        0.3
%    2      Circumpolar Antarctic  90S - 66.6S      0.1
%    3      Mediterranean Sea      30N - 46N        0.2
%    4      Red Sea                12N - 30N        0.2
%    5      Arctic Ocean           66.6N - 90N      0.1
%    6      Sea of Japan           33N - 53N        0.1 (0.8 for ¦loc¦=0)
%    7      Sulu Sea               5N - 13N         <0.1
%    8      Halmahera Basin        3S - 3N          <0.1
%    9      Celebes Basin          1N - 8N          0.2
%   10      Weber Deep             8S - 3S          0.2
%   11      Black Sea              41N - 47N        0.1
%   12      Baltic Sea             53N - 63N        0.1
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
%  - The units used for the variables in the equations are: 
%    p [MPa], z [m], lat [deg]
%
%  FUNCTION CALLS
%  1) p = seawaterDepth2Pressure_LeroyParthiot1998(z,lat,loc)
%
%  REFERENCES
%  - Leroy, C.C., and Parthiot, F. (1998). “Depth-pressure relationships 
%    in the oceans and seas”, J. Acoust. Soc. Am. 103, 1346–1352.
%  - Leroy, C.C. (1998). “Depth-pressure relationships in the oceans and 
%    seas” J. Acoust. Soc. Am. 121, 2447 Erratum.
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
%  28 Nov 2017
%
%**************************************************************************

% Error Management
switch loc
    case 0 % Common Oceans 
        latmin = -40;
        latmax = 60;
    case 1 % North East Atlantic
        latmin = 30;
        latmax = 35;
    case 2 % Circumpolar Antarctic
        latmin = -90;
        latmax = -66.6;
    case 3 % Mediterranean Sea
        latmin = 30;
        latmax = 46;
    case 4 % Red Sea
        latmin = 12;
        latmax = 30;
    case 5 % Arctic Ocean
        latmin = 66.6;
        latmax = 90;
    case 6 % Sea of Japan
        latmin = 33;
        latmax = 53;
    case 7 % Sulu Sea
        latmin = 5;
        latmax = 13;
    case 8 % Halmahera Basin
        latmin = -3;
        latmax = 3;
    case 9 % Celebes Basin
        latmin = 1;
        latmax = 8;
    case 10 % Weber Deep
        latmin = -8;
        latmax = -3;
    case 11 % Black Sea
        latmin = 41;
        latmax = 47;
    case 12 % Baltic Sea
        latmin = 53;
        latmax = 63;
    otherwise
        error('Not a valid LOC number')   
end

if ~(lat>latmin && lat<latmax)
    warning(['The latitude is outside the limits of selected sea. Check ' ...
        'that LAT and LOC are correct'])
end 

% Depth to Pressure Conversion
SIN1 = sin(lat*pi/180);
SIN2 = SIN1.*SIN1;
g = 9.780318*(1 + (-2.36e-5*SIN2 + 5.2788e-3)*SIN2); % acceleration of gravity [m s-2]
k = (g - 2e-5*z)./(9.80612 - 2e-5*z);
ps45 = (((2.8e-19*z - 1.25e-13).*z + 2.465e-8).*z + 1.00818e-2).*z;
ps = ps45.*k;

switch loc
    case 0 % Common Oceans 
        Dp = (1e-2./(z+100) + 6.2e-6).*z;
    case 1 % North East Atlantic
        Dp = (8e-3./(z+200) + 4e-6).*z;
    case 2 % Circumpolar Antarctic
        Dp = (8e-3./(z+1000) + 1.6e-6).*z;  
    case 3 % Mediterranean Sea
        Dp = (1.4e-9*z - 8.5e-6).*z;
    case {4 5} % Red Sea & Arctic Ocean
        Dp = 0;  
    case 6 % Sea of Japan
        Dp = 7.8e-6*z; 
    case 7 % Sulu Sea
        Dp = (1e-9*z + 1.6e-5 + 1e-2./(z+100)).*z ;  
    case 8 % Halmaera Basin
        Dp = (8e-3./(z+50) + 1.3e-5).*z;
    case {9 10} % Celebes Basin & Weber Deep
        Dp = (2.5e-10*z + 7e-6 + 1.2e-2./(z+100)).*z ;       
    case 11 % Black Sea
        Dp = 1.13e-4*z;
    case 12 % Baltic Sea
        Dp = 1.8e-4*z;
    otherwise
        error('Not a valid LOC number')   
end

p = (ps - Dp)*1e2; % hydrostatic (gauge) pressure [dbar] 

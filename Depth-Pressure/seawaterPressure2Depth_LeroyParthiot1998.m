
function z = seawaterPressure2Depth_LeroyParthiot1998(p,lat,loc)

%**************************************************************************
%  z = seawaterPressure2Depth_LeroyParthiot1998(p,lat,loc)
%
%  DESCRIPTION: calculates depth from hydrostatic (gauge) pressure in 
%  seawater for all oceans and seas using the equation proposed by Leroy 
%  & Parthiot (1998). 
%
%  The current pressure-depth equation is based on the UNESCO algorithm 
%  (Fofonoff & Millard, 1983) and on calculations from temperature and 
%  salinity profiles. The improvement over previous simple equations is 
%  substantial, and those equations should be abandoned.
% 
%  Pressure to depth conversion is used in ocean engineering applications,
%  e.g. to calculate the depth of an autonomous underwater vehicle (AUV)
%  equipped with a pressure sensor. 
%
%  INPUT VARIABLES
%  - p: gauge pressure, above atmospheric (vector) [dbar] 
%    (1 dbar = 10 kPa corresponds to a depth increase in seawater of ~1 m)
%  - lat: latitude [deg]
%  - loc: number specifying the geographic location (ocean or sea). 
%
%    loc    Area                   Latitude [deg]   Std. Dev. [+/- m]
%  --------------------------------------------------------------------
%    0      Common Oceans          40S - 60N        0.8
%    1      North East Atlantic    30N - 35N        0.3
%    2      Circumpolar Antarctic  90S - 66.6S      0.1
%    3      Mediterranean Sea      30N - 46N        0.2
%    4      Red Sea                12N - 30N        0.2
%    5      Arctic Ocean           66.6N - 90N      0.1
%    6      Sea of Japan           33N - 53N        0.1 (0.8 for ¦loc¦=0)
%    7      Sulu Sea               5N - 13N         0.2
%    8      Halmahera Basin        3S - 3N          0.1
%    9      Celebes Basin          1N - 8N          0.4
%   10      Weber Deep             8S - 3S          0.4
%   11      Black Sea              41N - 47N        0.1
%   12      Baltic Sea             53N - 63N        0.1
%   
%   NOTE: if ¦lat¦ and ¦loc¦ are not known, choose a value of ¦lat¦ = 45 
%   deg and ¦loc¦ = 0 for an approximate estimation of ¦z¦.
%        
%  OUTPUT VARIABLES
%  - z: depth below sea surface (vector) [m]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - The units used for the variables in the equations are: 
%    p [MPa], z [m], lat [deg]
%
%  FUNCTION CALLS
%  1) z = seawaterPressure2Depth_LeroyParthiot1998(p,lat,loc)
%
%  REFERENCES
%  - Leroy, C.C., and Parthiot, F. (1998). “Depth-pressure relationships 
%    in the oceans and seas”, J. Acoust. Soc. Am. 103, 1346–1352.
%  - Leroy, C.C. (1998). “Depth-pressure relationships in the oceans and 
%    seas” J. Acoust. Soc. Am. 121, 2447 Erratum.
%  - Fofonoff, N.P., & Millard, R.C. (1983) “Algorithm for computation of 
%    fundamental properties of seawater” UNESCO Technical Papers in Marine 
%    Science, No. 44.
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

% Pressure to Depth Conversion
p = p*1e-2; % gauge pressure [MPa]
SIN1 = sin(lat*pi/180);
SIN2 = SIN1.*SIN1;
g = 9.780318*(1 + (-2.36e-5*SIN2 + 5.2788e-3)*SIN2); % acceleration of gravity [m s-2]
zs = ((((-1.82e-7*p) + 2.279e-4).*p - 2.2512e-1).*p + 9.72659e2).*p ./ (g +  1.092e-4*p);

% Correction Factor
switch loc
    case 0 % Common Oceans 
        Dz = (1./(p + 1) + 5.7e-2).*p;
    case 1 % North East Atlantic
        Dz = (1./(p + 2) + 3e-2).*p;
    case 2 % Circumpolar Antarctic
        Dz = (-2e-4*p + 4e-2).*p;
    case 3 % Mediterranean Sea
        Dz = (2e-3*p - 7e-2).*p;
    case {4 5} % Red Sea & Arctic Ocean
        Dz = 0;
    case 6 % Sea of Japan
        Dz = 6e-2*p;
    case 7 % Sulu Sea
        Dz = (7e-4*p + 0.17 + 0.9./(p + 1)).*p;
    case 8 % Halmahera Basin
        Dz = (0.8./(p + 5) +1.25e-1).*p;
    case {9 10} % Celebes Basin & Weber Deep
        Dz = (2.2e-4*p + 6.7e-2 + 1.2./(p + 1)).*p; 
    case 11 % Black Sea
        Dz = 1.1*p;
    case 12 % Baltic Sea
        Dz = 1.8*p;
    otherwise
        error('Not valid LOC number')   
end

z = zs + Dz; % depth below sea surface [m]

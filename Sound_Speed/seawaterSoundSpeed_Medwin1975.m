
function c = seawaterSoundSpeed_Medwin1975(t,s,z)

%**************************************************************************
%  c = seawaterSoundSpeed_Medwin1975(t,s,z)
%
%  DESCRIPTION: calculates the speed of sound in seawater using the 
%  equation proposed by Medwin (1975). The expression is an approximation
%  to the NRLII equation (Del Grosso, 1974).
%
%  INPUT VARIABLES
%  - t: temperature, in ITS-90 (vector) [ºC]
%  - s: salinity (vector) [ppt]
%  - z: depth below sea surface (vector) [m]
%
%    NOTE: any input ¦t¦, ¦s¦ or ¦z¦ can either be a vector or a number.
%        
%  OUTPUT VARIABLES
%  - c: speed of sound in seawater (vector) [m s-1]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - The Medwin (1975) equation fits the NRLII equation (Del Grosso,1974) 
%    with a standard error of +/- 0.2 m/s. The range of applicability for 
%    the given standard error is: 
%
%    0 < t < 35 ºC, 0 < s < 45 ppt, 0 < z < 1000 m.
%
%  - Medwin (1975) does not mention the temperature scale considered
%    in the equation. The NRLII equation (Del Grosso, 1974), in which the 
%    Medwin's expression is based, is expressed in the International 
%    Practical Temperature Scale of 1968 (IPTS-68), and so is assumed to 
%    be the equation from Medwin. A conversion is made within the current 
%    function between the input temperature, expressed in the International 
%    Temperature Scale of 1990 (ITS-90), and the temperature in IPTS-68.
%  - The units used for the variables in the equation are: 
%    t [ºC], s [ppt], z [m]
%
%  FUNCTION CALLS
%  1) c = seawaterSoundSpeed_Medwin1975(t,s,z)
%
%  REFERENCES
%  - Medwin, H. (1975) “Sound speed in water: A simple equation for 
%    realistic parameters,” J. Acoust. Soc. Am. 58, 1318–1319.
%  - Del Grosso, V. A. (1974). “New equation for the speed of sound in 
%    natural waters (with comparisons to other equations)” J. Acoust. Soc. 
%    Am. 56, 1084–1091.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  21 Nov 2017
%
%**************************************************************************

% Convert to column vectors
t = t(:);
s = s(:);
z = z(:);
        
% Variables and Units Conversion
t = 1.00024*t; % temperature in IPTS-68 scale [ºC]

% Speed of Sound Equation
c = 1449.2 + 4.6*t - 5.5e-2*t.^2 + 2.9e-4*t.^3 + (1.34 - 1e-2*t).*(s-35) + 1.6e-2*z; % sound speed in seawater [m s-1]


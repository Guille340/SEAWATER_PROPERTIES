
function c = seawaterSoundSpeed_Ross1978(t,s,p)

%**************************************************************************
%  c = seawaterSoundSpeed_Ross1978(t,s,p)
%
%  DESCRIPTION: calculates the speed of sound in seawater using the
%  equation proposed by Ross (1978). The equation is based on a weighted 
%  average of sound speed values calculated with equations from Chen &
%  Millero (1977), DelGrosso (1974) and Wilson (1960c). 
%
%  INPUT VARIABLES
%  - t: temperature, in ITS-90 (vector) [ºC]
%  - s: salinity (vector) [ppt]
%  - p: gauge pressure, above atmospheric (vector) [dbar] 
%    (1 dbar = 10 kPa corresponds to a depth increase in seawater of ~1 m)
%
%    NOTE: any input ¦t¦, ¦s¦ or ¦p¦ can either be a vector or a number.
%        
%  OUTPUT VARIABLES
%  - c: speed of sound in seawater (vector) [m s-1]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS - General
%  - The equation is "sufficiently accurate" for realistic combinations
%    of temperature, salinity and pressure found in open oceans.
%  - The range of validity for t and s is not given, but for z > 3000 m 
%    the equation is not applicable (standard deviation > 0.6 m/s).
%  - The temeprature scale used for the equation is not specified. Chen 
%    & Millero (1977) and Del Grosso (1974) use IPTS-68 and Wilson (1960c)
%    uses IPTS-48. A conversion of the input ITS-90 temperature into
%    IPTS-68 is applied.
%  - The units used for the variables in the equations are: 
%    t [ºC], s [ppt], p [kg cm-2]
%
%  FUNCTION CALLS
%  1) c = seawaterSoundSpeed_Ross1978(t,s,p)
%
%  REFERENCES
%  - Ross, D. (1978). Revised simplified formulae for calculating the 
%    speed of sound in sea water. Saclant ASW Research Centre, No. SM 107.
%  - Chen, C.T., & Millero, F.J. (1977) “Sound speed in seawater at high 
%    pressures”, J. Acoust. Soc. Am. 62, 1129–1135.
%  - Del Grosso, V. A. (1974). “New equation for the speed of sound in 
%    natural waters (with comparisons to other equations)” J. Acoust. Soc. 
%    Am. 56, 1084–1091.
%  - Wilson, W.D. (1960c). “Equation for the speed of sound in sea water”, 
%    J. Acoust. Soc. Am. 32(10), 1357.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  24 Nov 2017
%
%**************************************************************************

% Convert to column vectors
t = t(:);
s = s(:);
p = p(:);

% Variables and Units Conversion
t = 1.00024*t; % temperature in IPTS-68 scale [ºC]
S = s - 35; % salinity referred to 35 ppt
p = p/9.80665; % pressure [kg cm2]

% Speed of Sound Equation
c0 = 1449.1;
ct = ((2.21e-4*t -5.17e-2).*t + 4.565).*t;
cs = 1.338*S;
cp = (1.25e-5*p + 1.592e-1).*p;
cstp = (-2.4e-7.*p + 2e-4).*p.*S + (-7.5e-7*p +2e-4).*t.*p ...
    + ((1e-4*t - 1.3e-2).*t + 1.338).*S;

c = c0 + ct + cs + cp + cstp;



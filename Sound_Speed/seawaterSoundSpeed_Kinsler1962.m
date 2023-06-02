
function c = seawaterSoundSpeed_Kinsler1962(t,s,z)

%**************************************************************************
%  c = seawaterSoundSpeed_Kinsler1962(t,s,z)
%
%  DESCRIPTION: calculates the speed of sound in seawater using the 
%  equation proposed by Kinsler & Frey (1962). The expression is a
%  simplification of 1st Wilson equation (Wilson, 1960b).
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
%  - The measurements used by Wilson (1960a) for his 1st equation (1960b), 
%    on which the current simple equation is based (Kinsler & Frey, 2000), 
%    have been proved inaccurate.
%  - The Kinsler & Frey (1962) equation fits 1st Wilson equation (Wilson, 
%    1962b) with a standard error of +/- 1 m/s for open oceans. The range 
%    of applicability for the given standard error is: 
%
%    -4 < t < 30 ºC, 33 < s < 37 ppt, 0 < z < 1000 m.
%
%  - Kinsler & Frey (1962) do not mention the temperature scale considered
%    in the equation. The equations from Wilson are expressed in the 
%    International Practical Temperature Scale of 1948 (IPTS-48), and so
%    is assumed to be the equation from Kinsler & Frey. A conversion is
%    made within the current function between the input temperature, 
%    expressed in the International Temperature Scale of 1990 (ITS-90),
%    and the temperature in IPTS-48. In any case, errors due to temperature
%    scale are expected to be small (order of 1/100 m/s).
%  - The units used for the variables in the equation are: 
%    t [ºC], s [ppt], z [m]
%
%  FUNCTION CALLS
%  1) c = seawaterSoundSpeed_Kinsler1962(t,s,z)
%
%  REFERENCES
%  - Kinsler, L.E., & Frey, A.R. (1962). Fundamentals of Acoustics. 2nd 
%    Ed., John Wiley & Sons: New York, London
%  - Wilson, W.D. (1960b). “Speed of sound in sea water as a function of 
%    temperature, pressure, and salinity”, J. Acoust. Soc. Am. 32(6), 641-
%    644.
%  - Wilson, W.D. (1960a) “Ultrasonic measurement of the velocity of sound 
%    in distilled and sea water,” Naval Ordnance Report 6746, US Naval 
%    Ordnance Laboratory, White Oak, Maryland.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  20 Nov 2017
%
%**************************************************************************
       
% Variables and Units Conversion
t = (-0.99956 + sqrt(0.9991202 + 1.76e-5*1.00024*t))/(8.8e-6); % temperature in IPTS-48 scale [ºC]

% Speed of Sound Equation
c = 1449 + 4.6*t - 5.5e-2*t.^2 + 3e-4*t.^3 + (1.39-1.2e-2*t).*(s-35) + 1.7e-2*z; % sound speed in seawater [m s-1]


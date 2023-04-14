
function c = seawaterSoundSpeed_Kinsler2000(t,s,p)

%**************************************************************************
%  c = seawaterSoundSpeed_Kinsler2000(t,s,p)
%
%  DESCRIPTION: calculates the speed of sound in seawater using the 
%  equation proposed by Kinsler et al. (2000). The expression is an
%  approximation to the NRLII equation (Del Grosso, 1974).
%
%  INPUT VARIABLES
%  - t: temperature, in ITS-90 (vector) [ºC]
%  - s: salinity (vector) [ppt]
%  - p: pressure above atmospheric pressure (vector) [dbar] (1 dbar = 10
%    kPa corresponds to a depth increase in seawater of ~1 m).
%
%    NOTE: any input ¦t¦, ¦s¦ or ¦p¦ can either be a vector or a number.
%        
%  OUTPUT VARIABLES
%  - c: speed of sound in seawater (vector) [m s-1]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - The Kinsler et al. (2000) equation fits Del Grosso's (1974) with a 
%    standard error of +/- 0.1 m/s for open oceans (+/- 0.2 m/s for 
%    atypical regions such as Sulu Sea, Halmaera Basin, Caribbean Basin, 
%    and East Indian Basins), with reduced accuracy (sigma = +/- 0.6 m/s) 
%    for Black and Baltic seas. The range of applicability for the given 
%    standard errors is: 
%
%    -2 < t < 30 ºC, 25 < s < 40 ppt, 0 < z < 6000 m.
%
%  - Kinsler et al. (2000) do not mention the temperature scale considered
%    in the equation. The equation from Del Grosso (1974) is expressed in 
%    the International Practical Temperature Scale of 1968 (IPTS-48), and 
%    so is assumed to be the equation from Kinsler et al. A conversion is
%    made within the current function between the input temperature, 
%    expressed in the International Temperature Scale of 1990 (ITS-90),
%    and the temperature in IPTS-68. In any case, errors due to temperature
%    scale are expected to be small (order of 1/100 m/s).
%  - The units used for the variables in the equation are: 
%    t [ºC], s [ppt], p [atm] (1 atm = 10.132501 dbar)
%
%  FUNCTION CALLS
%  1) c = seawaterSoundSpeed_Kinsler2000(t,s,p)
%
%  REFERENCES
%  - Kinsler, L.E., Frey,A.,Coppens, A.B., & Sanders, J.V. (2000). 
%    Fundamentals of Acoustics. 4th Ed., John Wiley & Sons: New York.
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

% Convert to column vectors
t = t(:);
s = s(:);
p = p(:);
        
% Variables and Units Conversion
t = 1.00024*t; % temperature in IPTS-68 scale [ºC]
p = p/10.132501; % pressure [atm]

% Speed of Sound Equation
c = 1449.08 + 4.57*t.*exp(-t/86.9 -(t/360).^2) + 1.33*(s-35).*exp(-t/120) ...
    + 0.1522*p.*exp(t/1200 + (s-35)/400) + 1.46e-5*p.^2.*exp(-t/20 + (s-35)/10); % sound speed in seawater [m s-1]




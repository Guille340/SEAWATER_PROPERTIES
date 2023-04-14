
function c = seawaterSoundSpeed_Coppens1980(t,s,z,lat,varargin)

%**************************************************************************
%  c = seawaterSoundSpeed_Coppens(t,s,z,varargin)
%
%  DESCRIPTION: calculates the speed of sound in seawater using Coppens 
%  formulation (Coppens, 1980). Coppens provides two simple equations to
%  fit the more general Lovett's 3rd equation (Lovett, 1978). The two
%  equations provide different degrees of fitting accuracy and only 
%  differ in the depth-dependent correction term. The most accurate of the
%  two equations is valid for all Neptunian waters (see definition below).
% 
%  Neptunian waters: all open ocean waters and seas connected by natural 
%  waterways navigable by medium-sized ships (Leroy, 1969).
%   
%  INPUT VARIABLES
%  - t: temperature, in ITS-90 (vector) [ºC]
%  - s: salinity (vector) [ppt]
%  - z: depth below sea surface (vector) [m]
%  - lat: latitude [deg]. The latitude affects the depth values used to 
%    calculate ¦c¦. If the latitude is not known, use ¦lat¦ = 45 is an
%    approximate estimation of the sound speed (for ¦lat¦ = 45 the exact 
%    depth input ¦z¦ is used in the calculations).
%  - eq (varargin{1}): string specifying the sound speed equation. Options
%    ¬ 'com': uses the complete depth-dependent correction term (DEFAULT).
%       This option provides a fitting accuracy of 0.03 m/s to 3rd Lovett
%       equation for s < 45 ppt, t < 35 ºC and z < 4 km. Valid for all
%       Neptunian waters.
%    ¬ 'bas': uses the basic depth-dependent correction term. 
%       This option provides a fitting accuracy of 0.1 m/s to 3rd Lovett
%       equation for s < 40 ppt, t < 30 ºC and z < 4 km. Red Sea,
%       Mediterranean Sea and Persian Gulf are not covered.
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
%  - Although the equation from Coppens (1980) is a generally accepted
%    simple expression for the calculation of speed of sound in seawater, 
%    the Lovett equations (Lovett, 1980) it is based on are fitted to
%    measurements made by Wilson (1960a), which have been proved 
%    inaccurate. 
%  - The expression from Coppens (1980) fits 3rd Lovett equation (Lovett, 
%    1978) with a standard error of +/- 0.1 m/s for the basic equation
%    and +/- 0.03 m/s for the complete equation. The corresponding ranges 
%    of applicability for both approximations are: 
%
%    s < 45 ppt, t < 35 ºC and z < 4 km - complete eq, all Neptunian waters
%
%    s < 40 ppt, t < 30 ºC and z < 4 km - basic eq, all Neptunian waters
%                                         except Red Sea, Mediterranean Sea
%                                         and Persian Gulf
%    
%  - Coppens (1980) does not mention the temperature scale considered
%    in the equation. The 3rd equation from Lovett, in which Coppens' 
%    expression is based, is expressed in the International Practical 
%    Temperature Scale of 1948 (IPTS-48), and so is assumed to be the 
%    equation Coppens. A conversion is made within the current function 
%    between the input temperature, expressed in the International 
%    Temperature Scale of 1990 (ITS-90), and the temperature in IPTS-48. 
%    Lovett (1978) confirms that errors due to temperature scale are 
%    expected to be small (0.03 m/s for input ¦t¦ in IPTS-68).
%  - The units used for the variables in the equation are: 
%    t [ºC], s [ppt], z [km], lat [deg]
%
%  FUNCTION CALLS
%  1) c = seawaterSoundSpeed_Coppens(t,s,z,lat)
%     ¬ eq = 'com'
%  2) c = seawaterSoundSpeed_Coppens(t,s,z,lat,eq)
%
%  REFERENCES
%  - Coppens, A.B. (1981) “Simple equations for the speed of sound in 
%    Neptunian waters” J. Acoust. Soc. Am. 69, 862–863.
%  - Lovett, J.R. (1978) “Merged seawater sound-speed equations” J. Acoust.
%    Soc. Am. 63(6), 1713-1718.
%  - Leroy, C.C. (1969) “Development of simple equations for accurate and 
%    more realistic calculation of the speed of sound in seawater” J. 
%    Acoust. Soc. Am. 46, 216–226.
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
switch nargin
    case {0 1 2 3}
        error('Not enough input arguments')
    case 4
        eq = 'com'; % complete expression
    case 5
        eq = varargin{1};
    otherwise
        error('Too many input arguments')
end

% Convert to column vectors
t = t(:);
s = s(:);
z = z(:);
        
% Variables and Units Conversion
z = z*1e-3; % depth [km]
lat = lat*pi/180; % latitude [rad]
z = z*(1-2.6e-3*cos(2*lat)); % corrected depth (z@45deg = z*(1))[km]
t = (-0.99956 + sqrt(0.9991202 + 1.76e-5*1.00024*t))/(8.8e-6); % temperature in IPTS-48 scale [ºC]
T = t/10; % normalised temperature [0.1 ºC]
c0 = 1449.05 + 45.7*T - 5.21*T.^2 + 0.23*T.^3 + (1.333 - 0.126*T + 0.009*T.^2).*(s-35); % seawater sound speed at atmospheric pressure (z = 0)

switch eq
    case 'com'
        Dc = (16.23 +  0.253*T).*z + (0.213 - 0.1*T).*z.^2 + (0.016 + 0.0002*(s-35)).*(s-35).*T.*z; % complete depth-dependent correction term [m s-1]
    case 'bas'
        Dc = 16.3*z + 0.18*z.^2; % basic depth-dependent correction term [m s-1]
end

c = c0 + Dc; % total sound speed in seawater [m s-1]

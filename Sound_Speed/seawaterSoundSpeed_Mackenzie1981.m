
function c = seawaterSoundSpeed_Mackenzie1981(t,s,z)

%**************************************************************************
%  c = seawaterSoundSpeed_Mackenzie1981(t,s,z)
%
%  DESCRIPTION: calculates the speed of sound in seawater using the 
%  equation proposed by Mackenzie (1981). The nine-term equation is based
%  on good quality oceanographic data gathered at 15 representative
%  worldwide stations. The fitting was performed over 14135 data points for
%  25 combinations of t,s,z.
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
%  - The fitting to measurements of Mackenzie's (1981) equation is achieved
%    with a standard error of +/- 0.07 m/s. The range of applicability is: 
%
%    -2 < t < 30 ºC, 25 < s < 40 ppt, 0 < z < 8000 m.
%
%  - In Mackenzie (1981) the temperature scale is not taken into account,
%    due to the small impact it will have on the sound speed calculations.
%    The temperature measurements used to derive the equations could be in
%    either IPTS-48 or IPTS-68. To increase precision, input temperature, 
%    expressed in ITS-90, is converted to IPTS-68.
%  - The units used for the variables in the equation are: 
%    t [ºC], s [ppt], z [m]
%
%  FUNCTION CALLS
%  1) c = seawaterSoundSpeed_Mackenzie1981(t,s,z)
%
%  REFERENCES
%  - Mackenzie, K.V. (1981) “Nine term equation for the sound speed in the 
%    oceans” J. Acoust. Soc. Am. 70, 807–812.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  21 Nov 2017
%
%**************************************************************************
        
% Variables and Units Conversion
t = 1.00024*t; % temperature in IPTS-68 scale [ºC]

% Speed of Sound Equation
c = 1448.96 + 4.591*t - 5.304e-2*t.^2 + 2.374e-4*t.^3 + 1.340*(s-35) + 1.630e-2*z ...
    + 1.675e-7*z.^2 - 1.025e-2*t.*(s-35) - 7.139e-13*t.*z.^3; % sound speed in seawater [m s-1]

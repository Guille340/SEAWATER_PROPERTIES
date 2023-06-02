
function c = seawaterSoundSpeed_Wilson1960(t,s,p,varargin)

%**************************************************************************
%  c = seawaterSoundSpeed_Wilson1960(t,s,p,varargin)
%
%  DESCRIPTION: calculates the speed of sound in seawater using the
%  equations proposed by Wilson (1960). The first equation (Wilson,1960b)
%  follows the general form of Mackenzie's equation (1960). The second 
%  equation (Wilson, 1960c) was developed to extend the validity of
%  the first equation to 0 < s < 37 ppt, after finding it was inaccurate 
%  for salinities outside the range of 33-37 ppt. Both equations are fitted
%  to laboratory measurements of seawater under pressure (Wilson, 1960a).
%
%  INPUT VARIABLES
%  - t: temperature, in ITS-90 (vector) [ºC]
%  - s: salinity (vector) [ppt]
%  - p: absolute pressure, including atmospheric (vector) [dbar] 
%    (1 dbar = 10 kPa corresponds to a depth increase in seawater of ~1 m)
%  - eq (varargin{1}): version of Wilson's equation. Options
%    ¬ 1: first equation (Wilson,1960b)
%    ¬ 2: second equation (Wilson,1960c) [DEFAULT]
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
%  - The temperature is expressed in the International Practical Salinity
%    Scale of 1948 (IPTS-48). A conversion of the input temperature, 
%    expressed in the International Temperature Scale of 1990 (ITS-90),
%    into IPTS-48 is made within the current function.
%  - The units used for the variables in the equations are: 
%    t [ºC], s [ppt], p [kg cm-2] (1 kg/cm2 = 9.80665 dbar)
%
%  CONSIDERATIONS & LIMITATIONS - 1st Equation (Wilson, 1960b)
%  - The measurements used to fit the equation consisted of 581 sound
%    speed values obtained from a combination of 15 temperatures, 5
%    salinities and 8 pressures (Wilson, 1960a).
%  - The equation fits the measurements with an standard error of 
%    +/- 0.22 m/s for the following tsp range, within which measurements 
%    where confined to:
%    
%    -4 < t < 30 ºC, 33 < s < 37 ppt, 1.033 < p < 1000 kg/cm2
%
%  CONSIDERATIONS & LIMITATIONS - 2nd Equation (Wilson, 1960c)
%  - The equation fits the measurements with an standard error of +/- 0.3
%    m/s for the following tsp range:
%
%    -4 < t < 30 ºC, 0 < s < 37 ppt, 1.033 < p < 1000 kg/cm2
%
%  FUNCTION CALLS
%  1) c = seawaterSoundSpeed_Wilson1960(t,s,p)
%     ¬ eq = 2
%  2) c = seawaterSoundSpeed_Wilson1960(t,s,p,eq)
%
%  REFERENCES
%  - Wilson, W.D. (1960a) “Ultrasonic measurement of the velocity of sound 
%    in distilled and sea water,” Naval Ordnance Report 6746, US Naval 
%    Ordnance Laboratory, White Oak, Maryland.
%  - Wilson, W.D. (1960b). “Speed of sound in sea water as a function of 
%    temperature, pressure, and salinity”, J. Acoust. Soc. Am. 32(6), 641-
%    644.
%  - Wilson, W.D. (1960c). “Equation for the speed of sound in sea water”, 
%    J. Acoust. Soc. Am. 32(10), 1357.
%  - Mackenzie, K.V. (1960). “Formulas for the computation of sound speed 
%    in sea water”, J. Acoust. Soc. Am. 32(100), 100-104.
%  - Mackenzie, K.V. (1960E). Errata: [“Formulas for the computation of 
%    sound speed in sea water”, J. Acoust. Soc. Am. 32(100), 100-104].
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  22 Nov 2017
%
%**************************************************************************

switch nargin
    case {0 1 2}
        error('Not enough input arguments')
    case 3
        if ~isnumeric([t s p]), error('t, s and p must be numbers'), end
        eq = 2;
    case 4
        if ~isnumeric([t s p]), error('t, s and p must be numbers'), end
        eq = varargin{1};
    otherwise
        error('Too many input arguments')
end

% Variables and Units Conversion
t = (-0.99956 + sqrt(0.9991202 + 1.76e-5*1.00024*t))/(8.8e-6); % temperature in IPTS-48 scale [ºC]
p = p/9.80665; % pressure [kg cm-2]
S = s - 35; % salinity referred to 35 ppt 

switch eq
    case 1        
        % Coefficients
        c0 = 1449.22; % sound speed in seawater at atmospheric pressure, 35 ppt salinity and 0ºC
        A1 =   4.62330e+00;
        A2 =  -5.45850e-02;
        A3 =   2.82200e-04;
        A4 =  -5.07000e-07;
        B1 =   1.60518e-01;
        B2 =   1.02790e-05;
        B3 =   3.45100e-09;
        B4 =  -3.50300e-12;
        C1 =   1.39100e+00;
        C2 =  -7.80000e-02;
        E01 = -1.19700e-02;
        E02 =  2.61000e-04;
        E03 = -1.96000e-07;
        E04 = -2.09000e-06;
        E11 = -2.79600e-04;
        E12 =  1.33020e-05;
        E13 = -6.64400e-08;
        E21 = -2.39100e-07;
        E22 =  9.28600e-10;
        E31 = -1.74500e-10;
        
        % Equation
        Dct = (((A4*t + A3).*t + A2).*t + A1).*t;
        Dcp = (((B4*p + B3).*p + B2).*p + B1).*p;
        Dcs = (C2*S + C1).*S;
        Dctps = ((E31*t.*p + (E22*t + E21).*t).*p + ((E13*t + E12).*t + E11).*t).*p ...
            + S.*((E04*t + E03*p + E02).*p + E01*t);     
        
    case 2
        c0 = 1449.14; % sound speed in seawater at atmospheric pressure, 35 ppt salinity and 0ºC
        A1 =   4.57210e+00;
        A2 =  -4.45320e-02;
        A3 =  -2.60450e-04;
        A4 =   7.98510e-06;
        B1 =   1.60272e-01;
        B2 =   1.02680e-05;
        B3 =   3.52160e-09;
        B4 =  -3.36030e-12;
        C1 =   1.39799e+00;
        C2 =   1.69202e-03;
        E01 = -1.12440e-02;
        E02 =  7.77110e-07;
        E03 =  7.70160e-05;
        E04 = -1.29430e-07;
        E05 =  3.15800e-08;
        E06 =  1.57900e-09;
        E11 = -1.86070e-04;
        E12 =  7.48120e-06;
        E13 =  4.52830e-08;
        E21 = -2.52940e-07;
        E22 =  1.85630e-09;
        E31 = -1.96460e-10;
        
        Dct = (((A4*t + A3).*t + A2).*t + A1).*t;
        Dcp = (((B4*p + B3).*p + B2).*p + B1).*p;
        Dcs = (C2*S + C1).*S;
        Dctps = ((E31*t.*p + (E22*t + E21).*t).*p + ((E13*t + E12).*t + E11).*t).*p ...
            + S.*(((E06*t + E05).*t + E04*p + E03).*p + (E02*t + E01).*t);
        
    otherwise
        error('Worng value for input EQ')
end
       
% Speed of Sound Equation
c = c0 + Dct + Dcp + Dcs + Dctps; % sound speed in seawater [m s-1]

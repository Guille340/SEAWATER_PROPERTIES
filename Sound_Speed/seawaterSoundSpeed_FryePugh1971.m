
function c = seawaterSoundSpeed_FryePugh1971(t,s,p)

%**************************************************************************
%  c = seawaterSoundSpeed_FryePugh1971(t,s,p)
%
%  DESCRIPTION: calculates the speed of sound in seawater using the
%  equation proposed by Frye & Pugh (1971). The equation fits the
%  measurements used by Wilson for his 1st equation (Wilson,1960a) for 
%  realistic combinations of temperature, salinity and pressure in 
%  Neptunian waters (the extended dataset used by Wilson for his 2nd 
%  equation was rejected by Frye & Pugh because it was not from the same 
%  statistical population as the original data). The formula is similar 
%  to Leroy's equation (Leroy, 1969) without selectable accuracy (single 
%  expression), and provides a better fit to Wilson measurements than 
%  Wilson 1st and 2nd equations (Wilson, 1960b; Wilson, 1960c). 
%
%  INPUT VARIABLES
%  - t: temperature, in ITS-90 (vector) [ºC]
%  - s: salinity (vector) [ppt]
%  - p: absolute pressure, including atmospheric (vector) [dbar] 
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
%  - The standard error of the equation with respect to Wilson data is 
%    0.1 m/s, lower than the 0.22 m/s and 0.3 m/s of 1st and 2nd Wilson
%    equations, respectively.
%  - The equation is applicable to 99.5% of the entire ocean volume for
%    the given standard deviation. The exact range of validity is the
%    same as that from Wilson 1st equation (Wilson, 1960b):
%
%    -3 < t < 30 ºC, 33 < s < 37 ppt, 1.033 < p < 1000 kg/cm2
%
%  - Frye & Pugh (1971) do not mention the temperature scale considered
%    in their equation. The Wilson (1960a), in which the current equation
%    is based, is expressed in the International Practical Temperature 
%    Scale of 1948 (IPTS-48). Accordingly, Frye & Pugh's equation is
%    assumed to be expressed in IPTS-48. A conversion is made within the 
%    current function between the input temperature, expressed in the 
%    International Temperature Scale of 1990 (ITS-90), and the temperature 
%    in IPTS-48. Errors caused by the use of the wrong temperature scale 
%    are generally < 0.05 m/s.
%  - The units used for the variables in the equations are: 
%    t [ºC], s [ppt], p [kg cm-2]
%
%  FUNCTION CALLS
%  1) c = seawaterSoundSpeed_FryePugh1971(t,s,p)
%
%  REFERENCES
%  - Frye, H.W., & Pugh, J.D. (1971) “A new equation for the speed of sound 
%    in seawater”, J. Acoust. Soc. Am. 50(1), 384-386.
%  - Wilson, W.D. (1960a) “Ultrasonic measurement of the velocity of sound 
%    in distilled and sea water,” Naval Ordnance Report 6746, US Naval 
%    Ordnance Laboratory, White Oak, Maryland.
%  - Leroy, C.C. (1969) “Development of simple equations for accurate and 
%    more realistic calculation of the speed of sound in seawater” J. 
%    Acoust. Soc. Am. 46, 216–226.
%  - Wilson, W.D. (1960b). “Speed of sound in sea water as a function of 
%    temperature, pressure, and salinity”, J. Acoust. Soc. Am. 32(6), 641-
%    644.
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
t = (-0.99956 + sqrt(0.9991202 + 1.76e-5*1.00024*t))/(8.8e-6); % temperature in IPTS-48 scale [ºC]
S = s - 35; % salinity referred to 35 ppt
p = p/9.80665; % pressure [kg cm2]

% Coefficients
c0   =   1.4493e+03;
cp1  =   1.5848e-01;
cp2  =   1.5720e-05;
cp4  =  -3.4600e-12;
ct1  =   4.5870e+00;
ct2  =  -5.3560e-02;
ct3  =   2.6040e-04;
cs1  =   1.1900e+00;
cs3  =   9.6000e-02;
ct2p =   1.3540e-05;
ctp2 =  -7.1900e-07;
cst  =  -1.2000e-02;

% Terms
cp = (((cp4*p).*p + cp2).*p + cp1).*p;
ct = ((ct3*t + ct2).*t + ct1).*t;
cs = ((cs3*S).*S + cs1).*S;
ctsp = ct2p*t.^2.*p + ctp2*t.*p.^2 + cst*S.*t;

% Equation
c =  c0 + ct + cs + cp + ctsp; % sound speed in seawater [m s-1]

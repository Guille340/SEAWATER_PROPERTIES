
function c = seawaterSoundSpeed_Lovett1978(t,s,p,varargin)

%**************************************************************************
%  c = seawaterSoundSpeed_Lovett1978(t,s,p,varargin)
%
%  DESCRIPTION: calculates the speed of sound in seawater using the
%  equations proposed by Lovett (1978). The three proposed equations apply
%  Anderson's (1971) pressure dependence, based on Wilson's (1960a) data,
%  to either DelGrosso (1974) equation or to DelGrosso-Mader (1972)
%  data at atmospheric pressure. In words of Lovett, data from DelGrosso-
%  Mader (1972) provides the best relationship of sound speed with 
%  temperature and salinity, whilst data from Wilson (1960a) provides the 
%  best relationship with pressure.
%
%  INPUT VARIABLES
%  - t: temperature, in ITS-90 (vector) [ºC]
%  - s: salinity (vector) [ppt]
%  - p: pressure above atmospheric pressure (vector) [dbar] (1 dbar = 10
%    kPa corresponds to a depth increase in seawater of ~1 m).
%  - eq (varargin{1}): version of Wilson's equation. Options
%    ¬ 1: first equation
%    ¬ 2: second equation [DEFAULT]
%    ¬ 3: third equation
%
%    * NOTE: for details about the characteristics and range of validity 
%            of each equation see "Considerations & Limitations".
%        
%  OUTPUT VARIABLES
%  - c: speed of sound in seawater (vector) [m s-1]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS - General
%  - The range of validity for the equations is not specified, but is
%    expected to be the same as that DelGrosso (1974) equation:
%
%    0 < t  < 30 ºC, 30 < s < 40 ppt, 0 < p < 10000 dbar
%
%  - The temperature in the 3 equations is expressed in the International
%    Practical Salinity Scale of 1948 (IPTS-48). A conversion of the input 
%    temperature, expressed in the International Temperature Scale of 1990
%    (ITS-90), into IPTS-48 is made within the current function.
%  - The units used for the variables in the equations are: 
%    t [ºC], s [ppt], p [dbar]

%
%  CONSIDERATIONS & LIMITATIONS - 1st Equation
%  - Anderson's pressure dependence is incorporated to the original
%    DelGrosso (1974) equation just by modifying 3 coefficients, in 
%    terms p, p^2 and p^3.
%  - The 1st Lovett equation fits DelGrosso (1974) equation with a 
%    standard error of 0.049 m/s.
%  - Term ct2p2 appears to be wrong (ct2p2 = 2.760566e-02 in Lovett, 1978).
%    An error in the exponential is the most likely cause; changing the
%    exponential from e-02 to e-10 seems to fix the problem.
%  - Check value: t =  2ºC, s =  34.7 ppt, p = 6000 dbar, 1559.462 m/s
%    * NOTE: perform the check without the temperature scale conversion, 
%      i.e. CTRL+R on line 111)
%
%  CONSIDERATIONS & LIMITATIONS - 2nd Equation
%  - Anderson's pressure dependence is incorporated to the original
%    DelGrosso-Mader (1972) data, resulting in a 24-variables equation.
%  - The 2nd Lovett equation fits DelGrosso (1974) data with a standard 
%    error of 0.035 m/s.
%  - Check value: t =  2ºC, s =  34.7 ppt, p = 6000 dbar, 1559.393 m/s
%    * NOTE: perform the check without the temperature scale conversion, 
%      i.e. CTRL+R on line 111)
%
%  CONSIDERATIONS & LIMITATIONS - 3rd Equation
%  - Simplification of the 24-variable equation, reduced to 12 variables.
%  - The 3rd Lovett equation fits DelGrosso (1974) data with a standard 
%    error of 0.063 m/s (0.048 m/s for combinations of tsp representing
%    real ocean conditions)
%  - Check value: t =  2ºC, s =  34.7 ppt, p = 6000 dbar, 1559.499 m/s
%    * NOTE: perform the check without the temperature scale conversion, 
%      i.e. CTRL+R on line 111)
%
%  FUNCTION CALLS
%  1) c = seawaterSoundSpeed_Lovett1978(t,s,p)
%     ¬ eq = 2
%  2) c = seawaterSoundSpeed_Lovett1978(t,s,p,eq)
%
%  REFERENCES
%  - Lovett, J.R. (1978) “Merged seawater sound-speed equations” J. Acoust.
%    Soc. Am. 63(6), 1713-1718.
%  - Wilson, W.D. (1960a) “Ultrasonic measurement of the velocity of sound 
%    in distilled and sea water,” Naval Ordnance Report 6746, US Naval 
%    Ordnance Laboratory, White Oak, Maryland.
%  - Anderson, E.R. (1971). Sound Speed in Seawater as a Function of 
%    Realistic Temperature, Salinity, Pressure. Naval Undersea Research 
%    and Development Centre (NUC), rep. NUC TP 243, 51 pp. 
%  - Del Grosso, V. A. (1974). “New equation for the speed of sound in 
%    natural waters (with comparisons to other equations)” J. Acoust. Soc. 
%    Am. 56, 1084–1091.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  23 Nov 2017
%
%**************************************************************************

switch nargin
    case {0 1 2}
        error('Not enough input arguments')
    case 3
        eq = 2;
    case 4
        eq = varargin{1};
    otherwise
        error('Too many input arguments')
end

t = (-0.99956 + sqrt(0.9991202 + 1.76e-5*1.00024*t))/(8.8e-6); % temperature in IPTS-48 scale [ºC]

switch eq
    case 1        
        % Coefficients
        c0    =  1402.392;
        ct1   =  5.011094e+00;
        ct2   = -5.509468e-02;
        ct3   =  2.215360e-04;
        cs1   =  1.329523e+00;
        cs2   =  1.289558e-04;
        cp1   =  1.598938e-02;
        cp2   =  2.478901e-07;
        cp3   = -8.485727e-12;
        cts   = -1.275628e-02;
        ctp   =  6.477152e-04;
        ct2p2 =  2.760566e-10; % error in the original paper (ct2p2 = 2.760566e-02 in Lovett, 1978, p. 1715)
        ctp2  = -1.656950e-08;
        ctp3  =  5.536118e-13;
        ct3p  = -4.466674e-08;
        cs2p2 = -1.681126e-11;
        ct2s  =  9.684032e-05;
        cts2p =  4.952146e-07;
        ctsp  = -3.473123e-05;
        
        % Equation Terms
        ct = ((ct3*t + ct2).*t + ct1).*t;
        cs = (cs2*s + cs1).*s;
        cp = ((cp3*p + cp2).*p + cp1).*p;       
        ctsp = cts*t.*s + ctp*t.*p + ct2p2.*t.^2.*p.^2 + ctp2*t.*p.^2 + ctp3*t.*p.^3 ...
            + ct3p*t.^3.*p + cs2p2*s.^2.*p.^2 + ct2s*t.^2.*s + cts2p*t.*s.^2.*p + ctsp*t.*s.*p;       
    case 2
        % Coefficients
        c0   = 1402.394;
        ct1  = 5.028849e+00;
        ct2  = -5.723758e-02;
        ct3  = 2.858485e-04;
        ct5  = -1.404216e-08;
        cs1  = 1.280746e+00;
        cs2  = 2.830167e-03;
        cs3  = -3.787896e-05;
        cp1  = 1.594777e-02;
        cp2  = 2.778778e-07;
        cp5  = 7.069489e-21;
        cts  = -1.280698e-02;
        ct2s = 1.040167e-04;
        ct3s3 = -9.301259e-11;
        ctp = 9.466535e-05;
        ctp2 = -1.23743e-08;
        ct2p = -7.100174e-06;
        ct2p3 = 8.592724e-14;
        ct3p = -9.02519e-08;
        ct3p2 = -2.70148e-11;
        csp3 = -7.816551e-13;
        cs2p3 = 1.303142e-14;
        cs3p2 = -6.265617e-13;
        ctsp = -2.238383e-06;
        ct2sp = 2.85346e-07;
        
        % Equation Terms
        ct = ((((ct5*t).*t + ct3).*t + ct2).*t + ct1).*t;
        cs = ((cs3*s + cs2).*s + cs1).*s;
        cp = ((((cp5*p).*p).*p +  cp2).*p + cp1).*p;
        ctsp = cts*t.*s + ct2s*t.^2.*s + ct3s3*t.^3.*s.^3 + ctp*t.*p + ctp2*t.*p.^2 ...
            + ct2p*t.^2.*p + ct2p3*t.^2.*p.^3 + ct3p*t.^3.*p + ct3p2*t.^3.*p.^2 ...
            + csp3*s.*p.^3 + cs2p3*s.^2.*p.^3 + cs3p2*s.^3.*p.^2 + ctsp*t.*s.*p + ct2sp*t.^2.*s.*p;
    case 3
        % Coefficients  
        c0 = 1402.394;
        ct1 = 5.01132e+00;
        ct2 = -5.513036e-02;
        ct3 = 2.221008e-04;
        cs1 = 1.332947e+00;
        cp1 = 1.605336e-02;
        cp2 = 2.12448e-07;
        cts = -1.266383e-02;
        ct2s = 9.543664e-05;
        ctp2 = -1.052396e-08;
        ctp3 = 2.183988e-13;
        csp3 = -2.253828e-13;
        cts2p = 2.062107e-08;
        
        % Equation Terms
        ct = ((ct3.*t + ct2).*t + ct1).*t;
        cs = cs1*s;
        cp = (cp2.*p + cp1).*p;
        ctsp = cts*t.*s + ct2s*t.^2.*s + ctp2*t.*p.^2 + ctp3*t.*p.^3 ...
            + csp3*s.*p.^3 + cts2p*t.*s.^2.*p; 
        
    otherwise
        error('Wrong value for input EQ')
end

c = c0 + ct + cs + cp + ctsp; % sound speed in seawater [m s-1]


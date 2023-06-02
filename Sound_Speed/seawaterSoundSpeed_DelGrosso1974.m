
function c = seawaterSoundSpeed_DelGrosso1974(t,s,p,varargin)

%**************************************************************************
%  c = seawaterSoundSpeed_DelGrosso1974(t,s,p,varargin)
%
%  DESCRIPTION: calculates the speed of sound in seawater using the NRLII
%  equation (DelGrosso, 1974) or its revised version (Wong & Zhu, 1995).
%
%  INPUT VARIABLES
%  - t: temperature, in ITS-90 (vector or matrix) [ºC]
%  - s: salinity (vector or matrix) [ppt]
%  - p: gauge pressure, above atmospheric pressure (vector or matrix) [dbar] 
%    (1 dbar = 10 kPa corresponds to a depth increase in seawater of ~1 m).
%  - eq (varargin{1}): string specifying the sound speed equation. Options
%    ¬ 'gro': uses the original formulation from Del Grosso, 1974 (IPTS-68)
%    ¬ 'won': uses the Wong & Zhu (1995) formulation of Del Grosso equation
%      (revised coefficients, ITS-90) [DEFAULT]
%
%    NOTE: any input ¦t¦, ¦s¦ or ¦p¦ can either be a vector or a number.
%        
%  OUTPUT VARIABLES
%  - c: speed of sound in seawater (vector or matrix) [m s-1]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS - Original DelGrosso eq. (Del Grosso, 1974)
%  - The Del Grosso equation, also known as NRLII, is based on absolute, 
%    extremely precise measurements of sound speed in seawater.
%  - The measurements used to derive the equation cover 110 combinations of
%    tsp, with 627 data points in total. Samples were limited to:
%
%    a) z < 2000 m --------------- 30 < s < 41 ppt
%    b) z > 2000 m, t < 15 ºC ---- 33 < s < 38 ppt
%
%  - No measurements were made at intermediate salinities under pressure 
%    (i.e. 0 < s < 30 ppt, p > 0 dbar). The NRLII is inaccurate for low
%    salinity waters at depth, due to the extrapolation of sound speed 
%    calculations outside the domain of measurements.
%  - The equation has a standard error of +/-0.05 m/s with respect to 
%    measurements. The general domain of applicability for the equation is:
%
%    0 < t  < 30 ºC, 30 < s < 40 ppt, 0 < p < 1000 kg/cm2
%
%  - The temperature is expressed in the International Practical Salinity
%    Scale of 1968 (IPTS-68). A conversion of the input temperature, 
%    expressed in the International Temperature Scale of 1990 (ITS-90),
%    into IPTS-68 is made within the current function.
%  - The units used for the variables in the equation are: 
%    t [ºC], s [ppt], p [kg cm-2] (1 kg/cm2 = 9.80665 dbar)
%
%  CONSIDERATIONS & LIMITATIONS - Revised DelGrosso eq. (Wong & Zhu, 1995)
%  - The equation is fitted on 1331 sound speed values evaluated with 
%    DelGrosso (1974) equation, and obtained from all possible tsp
%    combinations, each consisting of 11 values equally spaced within the
%    domain of validity.
%  - The domain of validity of the equation is exactly the same as that 
%    from DelGrosso's:
%
%    0 < t  < 30 ºC, 30 < s < 40 psu, 0 < p < 1000 kg/cm2
%
%  - Wong & Zhu (1995) revision provides updated coefficients for the
%    application of the ITS-90, but does not alter the original DelGrosso
%    equation. No temperature scale conversion is required.
%  - The units used for the variables in the equation are: 
%    t [ºC], s [psu], p [kg cm-2] (1 kg/cm2 = 9.80665 dbar)
%  - NOTE: ppt (parts per thousand), psu (practical salinity unit), and
%    g/kg are all equivalent salinity units.
%
%  FUNCTION CALLS
%  1) c = seawaterSoundSpeed_DelGrosso1974(t,s,p)
%     ¬ eq = 'won'
%  2) c = seawaterSoundSpeed_DelGrosso1974(t,s,p,eq)
%
%  REFERENCES
%  - Del Grosso, V. A. (1974). “New equation for the speed of sound in 
%    natural waters (with comparisons to other equations)” J. Acoust. Soc.
%    Am. 56, 1084–1091.
%  - Wong, G.S.K., Zhu, S. (1995) “Speed of sound in seawater as a function
%    of salinity, temperature and pressure” J. Acoust. Soc. Am. 97, 2235–
%    2237.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  21 Nov 2017
%
%**************************************************************************

switch nargin
    case {0 1 2}
        error('Not enough input arguments')
    case 3
        if ~isnumeric([t s p]), error('t, s and p must be numbers'), end
        eq = 'won'; % revised DelGrosso equation (Wong & Zhu, 1995)
    case 4
        if ~isnumeric([t s p]), error('t, s and p must be numbers'), end
        eq = varargin{1};
    otherwise
        error('Too many input arguments')
end

switch eq
    case 'won'       
        c0 =     1402.392; % speed of sound in pure water (s=0) at atmospheric pressure (p=0) and 0 ºC (t=0) [m s-1]
        ct1 =    5.012285e+00; 
        ct2 =   -5.511840e-02;
        ct3 =    2.216490e-04;
        cs1 =    1.329530e+00;
        cs2 =    1.288598e-04;
        cp1 =    1.560592e-01;
        cp2 =    2.449993e-05;
        cp3 =   -8.833959e-09;
        cts =   -1.275936e-02;
        ctp =    6.353509e-03;
        ct2p2 =  2.656174e-08;
        ctp2 =  -1.593895e-06;
        ctp3 =   5.222483e-10;
        ct3p =  -4.383615e-07;
        cs2p2 = -1.616745e-09;
        ct2s =   9.688441e-05;
        cts2p =  4.857614e-06;
        ctsp =  -3.406824e-04;
          
    case 'gro'
        t = 1.00024*t; % temperature in IPTS-68 scale [ºC]
        
        c0 =     1402.392; % speed of sound in pure water (s=0) at atmospheric pressure (p=0) and 0 ºC (t=0) [m s-1]
        ct1 =    5.01109398873e+00; 
        ct2 =   -5.50946843172e-02;
        ct3 =    2.21535969240e-04;
        cs1 =    1.32952290781e+00;
        cs2 =    1.98955756844e-04;
        cp1 =    1.56059257041e-01;
        cp2 =    2.44998688441e-05;
        cp3 =   -8.83392332513e-09;
        cts =   -1.27562783426e-02;
        ctp =    6.35191613389e-03;
        ct2p2 =  2.65484716608e-08;
        ctp2 =  -1.59349479045e-06;
        ctp3 =   5.22116437235e-10;
        ct3p =  -4.38031096213e-07;
        cs2p2 = -1.61674495909e-09;
        ct2s =   9.68403156410e-05;
        cts2p =  4.85639620015e-06;
        ctsp =  -3.40597039004e-04; 
    otherwise
        error('Invalid string for input parameter ¦eq¦')
        
end
       
% Variables and Units Conversion
p = p/9.80665; % pressure [kg cm-2]    

% Speed of Sound Equation
Dct = ((ct3*t + ct2).*t + ct1).*t; % t-dependent correction term [m s-1]
Dcs = (cs2*s + cs1).*s; % s-dependent correction term [m s-1]
Dcp = ((cp3*p + cp2).*p + cp1).*p; % p-dependent correction term [m s-1]
Dcstp = (ctp + ct3p*t.^2 +ctp2*p + ct2p2*t.*p + ctp3*p.^2).*t.*p ... 
+ (cts + ct2s*t + ctsp*p + cts2p*s.*p).*s.*t + (cs2p2.*s.*p).*s.*p; % tsp-dependent correction term [m s-1]
c = c0 + Dct + Dcs + Dcp + Dcstp; % total sound speed in seawater [m s-1]


function c = seawaterSoundSpeed_ChenMillero1977(t,s,p,varargin)

%**************************************************************************
%  c = seawaterSoundSpeed_ChenMillero1977(t,s,p,varargin)
%
%  DESCRIPTION: calculates the speed of sound in seawater using the 
%  equation proposed by Chen & Millero (1977) and included in UNESCO 
%  Techical Paper No. 44 (Fofonoff & Millard, 1983), or the revised 
%  version proposed by Wung & Zhu (1995).
%
%  INPUT VARIABLES
%  - t: temperature, in ITS-90 (vector or matrix) [ºC]
%  - s: salinity (vector or matrix) [ppt]
%  - p: pressure above atmospheric pressure (vector or matrix) [dbar] 
%    (1 dbar = 10 kPa corresponds to a depth increase in seawater of ~1 m).
%  - eq (varargin{1}): string specifying the sound speed equation. Options
%    ¬ 'une': uses the original formulation from Chen & Millero, as it
%      appears in the UNESCO Technical Paper No. 44 (IPTS-68).
%    ¬ 'won': uses the Wong & Zhu (1995) formulation of the UNESCO equation
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
%  CONSIDERATIONS & LIMITATIONS - Original UNESCO eq. (Fofonoff & Millard, 
%                                                      1983)
%
%  - The Chen & Millero equation, also knows as the UNESCO eq., is based 
%    on less accurate measurements than those from DelGrosso, 1974
%    (commercial velocimeter, non-repeated measurements).
%  - The Chen & Millero equation provides the difference in sound speeds
%    between saline and pure water. To calculate the sound speed in 
%    seawater c(s,t,p) an equation for the speed of sound in pure water 
%    c(0,t,p) is required. Two equations for c(0,t,p) were available at
%    the time: 1) Barlow & Yazgan (1967), with most temperatures out of
%    the domain of interest; 2) Wilson (1959), known to contain unreliable
%    data. Chen & Millero incorporated a corrected version of Wilson data
%    into the UNESCO equation, but the correction still produced inaccurate
%    results. A more accurate, recent equation for the speed of sound in
%    pure water under pressure (Belogol'skii et al., 1999) was used by 
%    Leroy et al. (2008) to improve UNESCO equation at depth (see
%    #seawaterSoundSpeed_Leroy2008.m).
%  - The domain of validity for the UNESCO equation is:
%
%    0 < t < 40 ºC, 0 < s < 40 ppt, 0 < p < 1000 bar
%
%  - The temperature is expressed in the International Practical Salinity
%    Scale of 1968 (IPTS-68). A conversion of the input temperature, 
%    expressed in the International Temperature Scale of 1990 (ITS-90),
%    into IPTS-68 is made within the current function.
%  - The units used for the variables in the equation are: 
%    t [ºC], s [ppt], p [bar] (1 bar = 10 dbar)
%
%  CONSIDERATIONS & LIMITATIONS - Revised UNESCO eq. (Wong & Zhu, 1995)
%  - The equation is fitted on 1331 sound speed values evaluated with 
%    UNESCO (1983) equation, and obtained from all possible tsp
%    combinations, each consisting of 11 values equally spaced within the
%    domain of validity.
%  - The domain of validity of the equation is exactly the same as that 
%    from UNESCO's:
%
%    0 < t < 40 ºC, 0 < s < 40 psu, 0 < p < 1000 bar
%
%  - Wong & Zhu (1995) revision provides updated coefficients for the
%    application of the ITS-90, but does not alter the original UNESCO
%    equation. No temperature scale conversion is required.
%  - The units used for the variables in the equation are: 
%    t [ºC], s [psu], p [bar] (1 bar = 10 dbar)
%  - NOTE: ppt (parts per thousand), psu (practical salinity unit), and
%    g/kg are all equivalent salinity units.
%
%  FUNCTION CALLS
%  1) c = seawaterSoundSpeed_ChenMillero1977(t,s,p)
%     ¬ eq = 'won'
%  2) c = seawaterSoundSpeed_ChenMillero1977(t,s,p,eq)
%
%  REFERENCES
%  - Chen, C.T., & Millero, F.J. (1977) “Sound speed in seawater at high 
%    pressures”, J. Acoust. Soc. Am. 62, 1129–1135.
%  - Fofonoff, N.P., & Millard, R.C. (1983) “Algorithm for computation of
%    fundamental properties of seawater” UNESCO Technical Papers in Marine 
%    Science, No. 44.
%  - Wong, G.S.K., Zhu, S. (1995) “Speed of sound in seawater as a function
%    of salinity, temperature and pressure” J. Acoust. Soc. Am. 97, 2235–
%    2237.
%  - Belogol’skii, V.A., Sekoyan, S.S., Samorukova, L.M., Stefanov, S.R., 
%    & Levtsov, V.I. (1999). “Pressure dependence of the sound velocity in 
%    distilled water,” Meas. Tech. 42, 406–413.
%  - Leroy, C.C., Robinson, S.P., & Goldsmith, M.J. (2008) “A new equation 
%    for the accurate calculation of sound speed in all oceans”, J. Acoust.
%    Soc. Am. 124(5), 2774-2782.
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
        C00 = 1402.388;
        C01 =  5.03830e+00;
        C02 = -5.81090e-02;
        C03 =  3.34320e-04;
        C04 = -1.47797e-06;
        C05 =  3.14190e-09;
        C10 =  1.53563e-01;
        C11 =  6.89990e-04; 
        C12 = -8.18290e-06;
        C13 =  1.36320e-07;
        C14 = -6.12600e-10;
        C20 =  3.12600e-05;
        C21 = -1.71110e-06;
        C22 =  2.59860e-08;
        C23 = -2.53530e-10;
        C24 =  1.04150e-12;
        C30 = -9.77290e-09;
        C31 =  3.85130e-10;
        C32 = -2.36540e-12;
        A00 =  1.38900e+00;
        A01 = -1.26200e-02;
        A02 =  7.16600e-05;
        A03 =  2.00800e-06;
        A04 = -3.21000e-08;
        A10 =  9.47420e-05;
        A11 = -1.25830e-05;
        A12 = -6.49280e-08;
        A13 =  1.05150e-08;
        A14 = -2.01420e-10;
        A20 = -3.90640e-07;
        A21 =  9.10610e-09;
        A22 = -1.60090e-10;
        A23 =  7.99400e-12;
        A30 =  1.10000e-10;
        A31 =  6.65100e-12;
        A32 = -3.39100e-13;
        B00 = -1.92200e-02;
        B01 = -4.42000e-05;
        B10 =  7.36370e-05;
        B11 =  1.79500e-07;
        D00 =  1.72700e-03;
        D10 = -7.98360e-06;
    case 'une'
        t = 1.00024*t; % temperature in IPTS-68 scale [ºC]
        
        C00 = 1402.388;
        C01 =  5.03711e+00;
        C02 = -5.80852e-02;
        C03 =  3.34200e-04;
        C04 = -1.47800e-06;
        C05 =  3.14640e-09;
        C10 =  1.53563e-01;
        C11 =  6.89820e-04; 
        C12 = -8.17880e-06;
        C13 =  1.36210e-07;
        C14 = -6.11850e-10;
        C20 =  3.12600e-05;
        C21 = -1.71070e-06;
        C22 =  2.59740e-08;
        C23 = -2.53350e-10;
        C24 =  1.04050e-12;
        C30 = -9.77290e-09;
        C31 =  3.85040e-10;
        C32 = -2.36430e-12;
        A00 =  1.38900e+00;
        A01 = -1.26200e-02;
        A02 =  7.16400e-05;
        A03 =  2.00600e-06;
        A04 = -3.21000e-08;
        A10 =  9.47420e-05;
        A11 = -1.25800e-05;
        A12 = -6.48850e-08;
        A13 =  1.05070e-08;
        A14 = -2.01220e-10;
        A20 = -3.90640e-07;
        A21 =  9.10410e-09;
        A22 = -1.60020e-10;
        A23 =  7.98800e-12;
        A30 =  1.10000e-10;
        A31 =  6.64900e-12;
        A32 = -3.38900e-13;
        B00 = -1.92200e-02;
        B01 = -4.42000e-05;
        B10 =  7.36370e-05;
        B11 =  1.79450e-07;
        D00 =  1.72700e-03;
        D10 = -7.98360e-06;
        
    otherwise
        error('Invalid string for input parameter EQ')
end
       
% Variables and Units Conversion
p = p/10; % pressure [bar]  

% Equation for the Speed of Sound
% Term s^2
D = D00 + D10*p;

% Term s^(3/2)
B1 = B10 + B11*t;
B0 = B00 + B01*t;
B  = B0 + B1.*p;

% Term s^1
A3 =   (A32*t + A31).*t + A30;
A2 =  ((A23*t + A22).*t + A21).*t + A20;
A1 = (((A14*t + A13).*t + A12).*t + A11).*t + A10;
A0 = (((A04*t + A03).*t + A02).*t + A01).*t + A00;
A  = ((A3.*p + A2).*p + A1).*p + A0;

% Term s^0
C3 =    (C32.*t + C31).*t + C30;
C2 =  (((C24.*t + C23).*t + C22).*t + C21).*t + C20;
C1 =  (((C14.*t + C13).*t + C12).*t + C11).*t + C10;
C0 = ((((C05.*t + C04).*t + C03).*t + C02).*t + C01).*t + C00;
C  = ((C3.*p + C2).*p + C1).*p + C0;

% Sound Speed 
c = C + (A + B.*sqrt(s) + D.*s).*s; % total sound speed in seawater [m s-1]

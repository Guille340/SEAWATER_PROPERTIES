
function rho = seawaterDensity_FofonoffMillard1983(t,s,p)

%**************************************************************************
%  rho = seawaterDensity_FofonoffMillard1983(t,s,p)
%
%  DESCRIPTION: calculates the density of seawater using the UNESCO
%  international equation of state of seawater EOS1980 (Fofonoff & Millard,
%  1983).
%
%  INPUT VARIABLES
%  - t: temperature, in ITS-90 (vector or matrix) [ºC]
%  - s: salinity (vector or matrix) [ppt]
%  - p: gauge pressure, above atmospheric pressure (vector or matrix) [dbar] 
%    (1 dbar = 10 kPa corresponds to a depth increase in seawater of ~1 m)
%
%    NOTE: any input ¦t¦, ¦s¦ or ¦p¦ can either be a vector or a number.
%        
%  OUTPUT VARIABLES
%  - rho: density of seawater (vector or matrix) [kg m-3]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - The EOS80 algorithm is valid for: 
%    0 < s < 42 ppt, -2 < t < 40 °C, 0 < p < 10000 dbar
%  - The units used for the variables in the equations are: 
%    t [°C, IPTS-68], s [ppt, PSS-78], p [bar], rho [kg m-3]
%  - Checking values (*): 
%        t [°C]       s [ppt]      p [dbar]      rho [kg m-3]
%       ----------------------------------------------------
%          5            0               0         999.96675
%          5           35               0        1027.67547
%         25           35           10000        1062.53817
%
%    * NOTE: the checking must be done directly with temperatures in 
%      IPTS-68 (i.e. do CTRL+R on line 52).
%
%  FUNCTION CALLS
%  1) rho = seawaterDensity_FofonoffMillard1983(t,s,p)
%
%  REFERENCES
%  - Fofonoff, N.P., & Millard, R.C. (1983) “Algorithm for computation of 
%    fundamental properties of seawater” UNESCO Technical Papers in Marine 
%    Science, No. 44.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  29 Nov 2017
%
%**************************************************************************

% Variables and Units Conversion
t = 1.00024*t; % temperature conversion from ITS-90 to IPTS-68 [°C]
p = p*1e-1; % gauge pressure [bar]

% Calculation of the Secant Bulk Modulus k - Coefficients
E = (((-5.155288e-5*t + 1.360477e-2).*t - 2.327105).*t + 148.4206).*t + 19652.21;
F = ((-6.1670e-5*t + 1.09987e-2).*t - 0.603459).*t + 54.6746;
G = (-5.3009e-4*t + 1.6483e-2).*t + 7.944e-2;
H = ((-5.77905e-7*t + 1.16092e-4).*t + 1.43713e-3).*t + 3.239908;
I = (-1.6078e-6*t - 1.0981e-5).*t + 2.2838e-3;
J =  1.91075e-4;
K = (5.2787e-8*t - 6.12293e-6).*t + 8.50935e-5;
M = (9.1697e-10*t + 2.0816e-8).*t - 9.9348e-7;
Ak = H + (I + J.*sqrt(s)).*s;
Bk = K + M.*s;

% Calculation of the Secant Bulk Modulus k - Equation
k0 = E + (F + G.*sqrt(s)).*s; % secant bulk modulus at surface (p = 0) [bar]
k = k0 + (Ak + Bk.*p).*p; % secant bulk modulus at depth z [bar]

% Calculation of Sea Water Density - Coefficients
A = ((((6.536332e-9*t - 1.120083e-6).*t + 1.001685e-4).*t - 9.09529e-3).*t + 6.793952e-2).*t + 999.842594; % density of reference pure water (IUPAC, 1976) [kg m-3]
B = (((5.3875e-9*t - 8.2467e-7).*t + 7.6438e-5).*t - 4.0899e-3).*t + 8.24493e-1;
C = (-1.6546e-6*t + 1.0227e-4).*t - 5.72466e-3;
D = 4.8314e-4;

% Calculation of Sea Water Density - Equation
rho0 = A + (B + C.*sqrt(s) + D.*s).*s; % sea water density at surface (p = 0) [kg m-3]
rho = rho0./(1 - p./k); % density of seawater at depth [kg m-3]

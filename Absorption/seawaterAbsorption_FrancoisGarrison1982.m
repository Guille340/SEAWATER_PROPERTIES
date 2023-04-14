
function alpha = seawaterAbsorption_FrancoisGarrison1982(t,s,pH,z,f,out)

%**************************************************************************
%  alpha = seawaterAbsorption_FrancoisGarrison1982(t,s,pH,z,f,out)
%
%  DESCRIPTION: calculates the sound absorption coefficient in seawater 
%  for a given frequency using the equation proposed by Francois & 
%  Garrison (1982a, 1982b).
%
%  The equation follows the nomenclature of Fisher & Simmons (1977). The
%  total absorption is the result of contributions from pure water (H2O)
%  and small amounts of salts in solution, namely magnesium sulphate 
%  (MgSO4) and boric acid (H3BO3). The absorption from pure water is 
%  caused by viscous attenuation; absorption from salts is related to 
%  chemical relaxation processes.Contribution from other components are 
%  considered negligible; if they exist, their effect is included as a 
%  residual error in the boric acid term.
%
%  Francois & Garrison (1982) equation represents an improvement of 10-20%
%  compared to the absorption predicted by equations from Fisher & Simmons 
%  (1977) and Schulkin & Marsh (1978). 
%
%  INPUT VARIABLES
%  - t: temperature, in ITS-90 (vector) [ºC]
%  - s: salinity (vector) [ppt]
%  - pH: coefficient of acidity (vector) [-] (pH inc., acidity dec.)
%  - z: depth below water surface (vector) [m]
%  - f: frequency (vector) [kHz]
%  - out: string specifying the type of output. Two options:
%    ¬ 'tot': total absorption from all contributors (M x 1 column vector, 
%       [H3BO3 + MgS04 + H20])
%    ¬ 'par': partial absorption from each contributor (M x 3 matrix,
%       [H3BO3 MgSO4 H20]).
%
%    NOTE: any input ¦t¦, ¦s¦, ¦pH¦, ¦z¦ or ¦f¦ can either be a vector
%          or a single number.
%        
%  OUTPUT VARIABLES
%  - alpha: absorption coefficient of seawater at frequency ¦f¦ [dB km-1].
%    Divide ¦alpha¦ by 1000 and multiply by the wavelength (lambda [m]) 
%    to express the absorption in [dB/lambda]. Lambda = c/f, with c the 
%    speed of sound in [m/s] and f the frequency in [Hz].
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - The equation is valid for: 
%    200 Hz < f < 1 MHz, z < 5000 m
%  - For an accuracy to within 5% the domain of applicability is:
%
%     MgSO4 Domain (10 < f < 500 kHz)    H2O Domain ( f > 500 kHz)
%     -----------------------------------------------------------
%      -2 < t < 22 °C                     0 < t < 30 °C
%      30 < s < 35 ppt                    0 < s < 40 ppt
%       0 < z < 3.5 km                    0 < z < 10 km
%
%  - The units used for the variables in the equations are: 
%    t [°C, IPTS-68], s [ppt], z [m], pH [-], f [kHz]
%  - Check values (*):
%    
%      t [°C]   s [ppt]   f [kHz]    pH    z [m]   alpha [dB km-1]
%    --------------------------------------------------------------------
%      0        30        1          8     0          6.10e-02 
%      0        35        1          8     0          6.40e-02
%     10        35        1          8     0          6.00e-02
%      0        30       10          8     0          1.12e+00
%      0        35       10          8     0          1.29e+00
%     10        35       10          8     0          9.63e-01
%      0        30      100          8     0          2.10e+01
%      0        35      100          8     0          2.34e+01
%     10        35      100          8     0          3.36e+01
%     
%   (*) NOTE: to do the check, the temperature scale conversion must be
%       disabled (temperature directly specified in IPTS-68; replace
%       t = 1.00024*t by t = 1*t).
%
%  FUNCTION CALLS
%  1) alpha = seawaterAbsorption_FrancoisGarrison1982(t,s,pH,z,f,out)
%
%     Examples
%     --------------------
%     a) Multiple Oceanographic Values, Fixed Frequency:
%        z  = [   0  100  500 1000 5000];
%        t  = [  20   15    5    3    2];
%        s  = [  35   35   35   35   35]; --> can be replaced by s = 35
%        pH = [   8    8    8    8    8];
%        f  = 10;
%
%     b) Multiple Frequencies, Fixed Oceanographic Values
%        z  = 0;
%        t  = 20;
%        s  = 35;
%        pH = 8;
%        f  = 0.1:0.1:1000; --> from 100 Hz to 1 MHz
%
%  NOTES ON CONTRIBUTION FROM PURE WATER, MAGNESIUM SULPHATE AND BORIC ACID
%  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
%   Contributor     Dominant         Absorption    Relaxation 
%                   Frequencies      Increases     Frequency
%                   [kHz]            as...         Increases as...
%  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%   Pure water      f > 200          t  increases   N/A
%   H2O                              p  increases
%                                    f  increases
%  
%   Magnesium       10 < f < 200     t  increases   t increases
%   Sulphate                         s  increases
%   MgSO4                            p  decreases
%                                    f  increases 
%                                       (f < f2)
% 
%   Boric          0.2 < f < 10      pH increases   t increases
%   Acid                             f  increases   s increases
%   H3BO3                                    (f < f1)
%                                           
%   Scattering     f < 0.2
%  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%
%  REFERENCES
%  - Francois, R.E., & Garrison, G.R. (1982a). “Sound absorption based on 
%    ocean measurements. Part I: Pure water and magnesium sulfate 
%    contributions”, J. Acoust. Soc. Am., 72(3), 896-907.
%  - Francois, R.E., & Garrison, G.R. (1982b). “Sound absorption based on 
%    ocean measurements. Part II: Boric acid contribution and equation for 
%    total absorption”, J. Acoust. Soc. Am., 72(6), 1879-1890.
%  - Fisher, F.H., & Simmons, V.P. (1977). “Sound absorption in sea water”, 
%    J. Acoust. Soc. Am. 62(3), 558-564
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  29 Nov 2017
%  
%  REVISION 1.1 (30 Nov 2017)
%  - Equation modified to accept vectors of ¦t¦, ¦s¦, ¦pH¦ and ¦p¦.
%  - Error control for input vectors of different size.
%
%  REVISION 1.2 (30 Nov 2017)
%  - Depth ¦z¦ replaces ¦p¦ as input parameter, using the approximation
%    described by Francois & Garrison (1982a) (¦p¦ = ¦z¦/10, with ¦p¦ in
%    [atm] and ¦z¦ in [m]). To obtain a better result both pressure ¦p¦ 
%    and speed of sound c should be calculated using more accurate
%    equations (e.g. Leroy & Parthiot, 1998; Chen & Millero, 1974). For
%    this function it has been opted for the simplest, original formulas.
%  - Equation modified to accept vectors of ¦f¦, along with ¦t¦, ¦s¦, ¦pH¦
%    and ¦z¦. Error control on input vectors of different size removed
%    to accept single values.
%
%**************************************************************************

% Convert to column vectors
t = t(:);
s = s(:);
z = z(:);
pH = pH(:);
f = f(:);

% Variables and Units Conversion
t = 1.00024*t; % temperature conversion from ITS-90 to IPTS-68 [°C]
T = t + 273; % temperature in IPTS-68 [K]
c = 1412 + 3.21*t + 1.19*s + 1.67e-2*z; % speed of sound in seawater [m s-1]

% Boric Acid (H3BO3)
P1 = 1; % depth-dependent term [-]
A1 = 8.86*10.^(0.78*pH - 5)./c; % general absorption term at sea surface [dB km-1 kHz-1]
f1 = 2.8*sqrt(s/35).*10.^(4 - 1245./T); % relaxation frequency [kHz]
alpha1 = A1.*P1.*f1.*f.^2./(f.^2 + f1.^2); % partial absorption coefficient [dB km-1]

% Magnesium Sulphate (MgSO4)
P2 = (6.2e-9*z - 1.37e-4).*z + 1; % depth-dependent term [-]
A2 = 21.44*s.*(1 + 2.5e-2*t)./c; % general absorption term at sea surface [dB km-1 kHz-1]
f2 = (8.17*10.^(8 - 1990./T))./(1 + 1.8e-3*(s-35)); % relaxation frequency [kHz]
alpha2 = A2.*P2.*f2.*f.^2./(f.^2 + f2.^2); % partial absorption coefficient [dB km-1]

% Pure Water (H2O)
P3 = (4.9e-10*z - 3.83e-5).*z + 1; % depth-dependent term [-]
ind_01 = t > 20; % logical vector, 1 for t > 20 °C
ind_02 = ~ind_01; % logical vector, 1 for t <= 20 °C
t_01 = t(ind_01); % vector of temperatures, t > 20 °C
t_02 = t(ind_02); % vector of temperatures, t <= 20 °C
A3_01 = ((-6.5e-10*t_01 + 1.45e-7).*t_01 - 1.146e-5).*t_01 + 3.964e-4; % general absorption term at sea surface for t > 20 ºC [dB km-1 kHz-1]
A3_02 = ((-1.5e-8*t_02 + 9.11e-7).*t_02 - 2.59e-5).*t_02 + 4.937e-4; % general absorption term at sea surface for t <= 20 ºC [dB km-1 kHz-1]
A3(ind_01,1) = A3_01; % general absorption term at sea surface for t > 20 ºC (original positions) [dB km-1 kHz-1]
A3(ind_02,1) = A3_02; % general absorption term at sea surface for t <= 20 ºC (original positions) [dB km-1 kHz-1]
alpha3 = A3.*P3.*f.^2; % partial absorption coefficient [dB km-1]
    
% Absorption Coefficient
switch out
    case 'tot'
        alpha = alpha1 + alpha2 + alpha3;
    case 'par'
        alpha = [alpha1 alpha2 alpha3];
    otherwise
        error('Wrong string for OUT input')
end

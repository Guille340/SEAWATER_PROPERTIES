
function alpha = seawaterAbsorption_AinslieMcColm1998(t,s,pH,z,f,out)

%**************************************************************************
%  alpha = seawaterAbsorption_AinslieMcColm1998(t,s,pH,z,f,out)
%
%  DESCRIPTION: calculates the sound absorption coefficient in seawater 
%  for a given frequency using the equation proposed by Ainslie & McColm
%  (1998).
%
%  The equation is a simplified version of that from Francois & Garrison 
%  (1982a, 1982b). As the latter, it includes the contributions from pure 
%  water H2O (viscous absorption, f > 100 kHz), magnesium sulphate MgSO4 
%  (chemical relaxation, f > 10 kHz) and boric acid H3BO3 (chemical 
%  relaxation, f < 10 kHz). The absorption caused by other components 
%  (e.g. magnesium carbonate MgCO3) are not explicitly included in 
%  Francois & Garrison equation, but their effects are implicit in the 
%  boric acid term as a residual error.
% 
%  The equation fits Francois & Garrison's (1982b) to within 10% for:
%  100 Hz < f < 1 MHz, -6 < t < 35 ºC, 7.7 < pH < 8.3, 5 < s < 50 ppt,
%  0 < z < 7 km (calculated at t = 10 ºC, s = 35 ppt, z = 0 km, pH = 8).
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
%  - The equation is valid for (10% error): 
%    0.1 < f < 1000 kHz, -6 < t < 35 ºC, 7.7 < pH < 8.3, 5 < s < 50 ppt,
%    0 < z < 7 km
%  - ¦t¦ assumed to be in IPTS-68, same as Francois & Garrison (1982b)
%  - The units used for the variables in the equations are: 
%    t [°C, IPTS-68], s [ppt], z [km], pH [-], f [kHz]
%
%  FUNCTION CALLS
%  1) alpha = seawaterAbsorption_AinslieMcColm1998(t,s,pH,z,f,out)
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
%  NOTES ON EFFECTS OF WATER PARAMETERS ON THE ABSORPTION
%  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
%  Parameter  As Parameter     Absorption        Frequency
%             ...              ...               Range
%  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%   t         Increases        Increases         f1 and f2
%                       
%                              Decreases         All f except
%                                                f1 and f2
%   
%   s         Increases        Increases         f > f1
%                              Decreases         f << f1
%
%   pH        Increases        Increases         f < f1
%                              Constant          f > f1
%
%   z         Increases        Decreases         f > f1
%                              Constant          f < f1
%  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%
%  REFERENCES
%  - Ainslie, M.A., & McColm, J.G. (1998). “A simplified formula for 
%    viscous and chemical absorption in sea water”, J. Acoust. Soc. Am. 
%    103(3), 1671-1672.
%  - Francois, R.E., & Garrison, G.R. (1982a). “Sound absorption based on 
%    ocean measurements. Part I: Pure water and magnesium sulfate 
%    contributions”, J. Acoust. Soc. Am., 72(3), 896-907.
%  - Francois, R.E., & Garrison, G.R. (1982b). “Sound absorption based on 
%    ocean measurements. Part II: Boric acid contribution and equation for 
%    total absorption”, J. Acoust. Soc. Am., 72(6), 1879-1890.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  30 Nov 2017
%
%**************************************************************************

% Variables and Units Conversion
t = 1.00024*t; % temperature conversion from ITS-90 to IPTS-68 [°C]
z = z*1e-3; % depth below sea surface [km]

% Boric Acid (H3BO3)
f1 = 0.78*sqrt(s/35).*exp(t/26); % relaxation frequency [kHz]
alpha1 = 0.106*f1.*f.^2./(f.^2 + f1.^2).*exp((pH - 8)/0.56); % partial absorption coefficient [dB km-1]

% Magnesium Sulphate (MgSO4)
f2 = 42*exp(t/17); % relaxation frequency [kHz]
alpha2 = 0.52*(1 + t/43).*s/35.*f2.*f.^2./(f.^2 + f2.^2).*exp(-z/6); % partial absorption coefficient [dB km-1]

% Pure Water (H2O)
alpha3 = 4.9e-4*f.^2.*exp(-t/27 - z/17); % partial absorption coefficient [dB km-1]

% Absorption Coefficient
switch out
    case 'tot'
        alpha = alpha1 + alpha2 + alpha3;
    case 'par'
        alpha = [alpha1 alpha2 alpha3];
    otherwise
        error('Wrong string for OUT input')
end


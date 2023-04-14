
function alpha = seawaterAbsorption_FisherSimmons1977(t,z,f,out)

%**************************************************************************
%  alpha = seawaterAbsorption_FisherSimmons1977(t,z,f,out)
%
%  DESCRIPTION: calculates the sound absorption coefficient in seawater 
%  for a given frequency using the equation proposed by Fisher & Simmons
%  (1977).
%
%  The equation is based on: 1) laboratory measurements from Simmons (PhD
%  thesis) on the sound absorption in seawater at atmospheric pressure
%  due to magnesium sulphate and boric acid; 2)reexamination of the
%  pressure effect on the absorption in magnesium sulphate and pure water.
%  The equation is a function of frequency, temperature and pressure, and
%  unlike the equations from Francois & Garrison (1982a, 1982b), Ainslie &
%  McColm (1998) or Kinsler et al. (2000), is only applicable for a fixed
%  salinity and acidity coefficient (s = 35 ppt, pH = 8).
%
%  The general expression shows that sound absorption in seawater is the 
%  result of contributions from pure water (H2O), magnesium sulphate
%  (MgSO4) and boric acid (H3BO3).
%
%  INPUT VARIABLES
%  - t: temperature, in ITS-90 (vector) [ºC]
%  - z: depth below water surface (vector) [m]
%  - f: frequency (vector) [kHz]
%  - out: string specifying the type of output. Two options:
%    ¬ 'tot': total absorption from all contributors (M x 1 column vector, 
%       [H3BO3 + MgS04 + H20])
%    ¬ 'par': partial absorption from each contributor (M x 3 matrix,
%       [H3BO3 MgSO4 H20]).
%
%    NOTE: any input ¦t¦, ¦z¦ or ¦f¦ can either be a vector or a number.
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
%  - The equation is only valid for s = 35 ppt and pH = 8.
%  - The units used for the variables in the equations are: 
%    t [°C, IPTS-68], p [atm], f [Hz]
%  - The pressure is approximated to p = ¦z¦/10, as in Francois & Garrison
%    (1982a, 1982b) and other absorption equations based on it (e.g.
%    Ainslie & McColm, 1998; Kinsler et al., 2000).
%
%  FUNCTION CALLS
%  1) alpha = seawaterAbsorption_FisherSimmons1977(t,z,f,out)
%
%     Examples
%     --------------------
%     a) Multiple Oceanographic Values, Fixed Frequency:
%        z  = [   0  100  500 1000 5000];
%        t  = [  20   15    5    3    2];
%        f  = 10;
%
%     b) Multiple Frequencies, Fixed Oceanographic Values
%        z  = 0;
%        t  = 20;
%        f  = 0.1:0.1:1000; --> from 100 Hz to 1 MHz
%
%  FREQUENCY RANGE OF CONTRIBUTIONS FROM H2O, MgSO4 and H3BO3
%  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%   Contributor            Dominant Frequencies [kHz]
%  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%   Pure water                  f > 200         
%   H2O                              
%  
%   Magnesium Sulphate     10 < f < 200                            
%   MgSO4                           
%                                    
%   Boric Acid                  f < 10                                
%   H3BO3
%  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%
%  REFERENCES
%  - Fisher, F.H., & Simmons, V.P. (1977). “Sound absorption in sea water”,
%    J. Acoust. Soc. Am. 62(3), 558-564
%  - Francois, R.E., & Garrison, G.R. (1982a). “Sound absorption based on 
%    ocean measurements. Part I: Pure water and magnesium sulfate 
%    contributions”, J. Acoust. Soc. Am., 72(3), 896-907.
%  - Francois, R.E., & Garrison, G.R. (1982b). “Sound absorption based on 
%    ocean measurements. Part II: Boric acid contribution and equation for 
%    total absorption”, J. Acoust. Soc. Am., 72(6), 1879-1890.
%  - Ainslie, M.A., & McColm, J.G. (1998). “A simplified formula for 
%    viscous and chemical absorption in sea water”, J. Acoust. Soc. Am. 
%    103(3), 1671-1672.
%  - Kinsler, L.E., Frey,A.,Coppens, A.B., & Sanders, J.V. (2000). 
%    Fundamentals of Acoustics. 4th Ed., John Wiley & Sons: New York.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  30 Nov 2017
%
%**************************************************************************

% Convert to column vectors
t = t(:);
z = z(:);
f = f(:);

% Variables and Units Conversion
p = z/10; % gauge pressure [atm]
f = f*1e3; % frequency [Hz]
t = 1.00024*t; % temperature conversion from ITS-90 to IPTS-68 [°C]
T = t + 273.1; % temperature [K]

% Boric Acid (H3BO3)
P1 = 1; % depth-dependent term [-]
A1 = (-5.22e-12*t + 2.36e-10).*t + 1.03e-8; % general absorption term at sea surface [Np m-1 Hz-1] (std. dev. +/- 2%)
f1 = 1315*T.*exp(-1700./T); % relaxation frequency [Hz] (std. dev. +/- 1%)
alpha1 = A1.*P1.*f1.*f.^2./(f1.^2 + f.^2); % partial absorption coefficient [Np m-1]

% Magnesium Sulphate (MgSO4)
P2 = (3.7e-7*p - 1.03e-3).*p + 1; % depth-dependent term [-] (std. dev. +/- 1%)
A2 = 7.52e-10*t + 5.62e-8; % general absorption term at sea surface [Np m-1 Hz-1] (std. dev. +/- 5%)
f2 = 1.55e7*T.*exp(-3052./T); % relaxation frequency [Hz] (std. dev. +/- 4%)
alpha2 =  A2.*P2.*f2.*f.^2./(f2.^2 + f.^2); % partial absorption coefficient [Np m-1]

% Pure Water (H2O)
P3 = (7.57e-8*p - 3.84e-4).*p + 1; % depth-dependent term [-] (std. dev. +/- 1%)
A3 = (((-3.48e-4*t + 4.77e-2).*t - 2.37).*t + 55.9)*1e-15;  % general absorption term at sea surface [Np m-1 Hz-1] (std. dev. +/- 4%)
alpha3 = A3.*P3.*f.^2; % partial absorption coefficient [Np m-1]

% Absorption Coefficient
switch out
    case 'tot'
        alpha = (alpha1 + alpha2 + alpha3)*8686; % total absorption coefficient [dB km-1]
    case 'par'
        alpha = [alpha1 alpha2 alpha3]*8686; % partial absorption coefficients [dB km-1]
    otherwise
        error('Wrong string for OUT input')
end

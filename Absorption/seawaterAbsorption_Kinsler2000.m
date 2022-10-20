
function alpha = seawaterAbsorption_Kinsler2000(t,s,pH,z,f,out)

%**************************************************************************
%  alpha = seawaterAbsorption_Kinsler2000(t,s,pH,z,f,out)
%
%  DESCRIPTION: calculates the sound absorption coefficient in seawater 
%  for a given frequency using the equation proposed by Kinsler et al.
%  (2000).
%
%  The equation includes the contributions from pure water H2O (structural
%  relaxation, f > 100 kHz), magnesium sulphate MgSO4 (chemical relaxation,
%  f > 10 kHz) and boric acid H3BO3 (chemical relaxation, f < 10 kHz). 
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
%    100 Hz < f < 1 MHz, z < 6000 m
%  - ¦t¦ assumed to be in IPTS-90 (no scale conversion required)
%  - The units used for the variables in the equations are: 
%    t [°C, IPTS-90], s [ppt], z [km], pH [-], f [kHz]
%
%  FUNCTION CALLS
%  1) alpha = seawaterAbsorption_Kinsler2000(t,s,pH,z,f,out)
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
%  NOTES ON CAUSES OF ABSORPTION IN SEAWATER
%  a) Structural Relaxation
%     - Structural relaxation processes are the main cause of sound
%       absorption in pure water.
%     - A relaxational dissipation occurs when a passing compressional
%       wave transfers the water molecules from the equilibrium state to
%       a high energy state.
%     - Structural relaxation in pure water can be taken into account by
%       a coefficient of bulk viscosity.
%
%  b) Chemical Relaxation
%     - Chemical relaxation processes are the main cause of sound
%       absorption from dissolved salts in seawater. Magnesium sulphate
%       (MgSO4) and boric acid (H3BO3) are the major cause of absorption 
%       associated to chemical relaxation. Magnesium carbonate MgCO3 and
%       other components are expected to have a much lower contribution.
%     - The relaxational dissipation occurs when the acoustic process
%       changes the concentration of associated and dissociated ions of
%       MgSO4 or H3BO3. The relaxation time of this process results in
%       sound energy absorption.
%
%  FREQUENCY RANGE OF CONTRIBUTIONS FROM H2O, MgSO4 and H3BO3
%  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%   Contributor            Dominant Frequencies [kHz]
%  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%   Pure water                  f > 100         
%   H2O                              
%  
%   Magnesium Sulphate     10 < f < 100                            
%   MgSO4                           
%                                    
%   Boric Acid              1 < f < 10                                
%   H3BO3
%
%   Scattering                  f < 1
%  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%
%  REFERENCES
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
s = s(:);
z = z(:);
pH = pH(:);
f = f(:);

% Variables and Units Conversion
z = z*1e-3; % depth below sea surface [km]
f = f*1e3; % frequency [Hz]

% Boric Acid (H3BO3)
f1 = 780.*exp(t/29); % relaxation frequency [Hz]
alpha1 = 0.083*(s/35).*exp(t/31 - z/91 + 1.8*(pH - 8)).*f.^2./(f.^2 + f1.^2); % partial absorption coefficient [dB km-1]

% Magnesium Sulphate (MgSO4)
f2 = 42e3*exp(t/18); % relaxation frequency [Hz]
alpha2 = 22*(s/35).*exp(t/14 - z/6).*f.^2./(f.^2 + f2.^2); % partial absorption coefficient [dB km-1]

% Pure Water (H2O)
alpha3 = 4.9e-10*f.^2.*exp(-t/26 + z/25); % partial absorption coefficient [dB km-1]

% Absorption Coefficient
switch out
    case 'tot'
        alpha = alpha1 + alpha2 + alpha3;
    case 'par'
        alpha = [alpha1 alpha2 alpha3];
    otherwise
        error('Wrong string for OUT input')
end

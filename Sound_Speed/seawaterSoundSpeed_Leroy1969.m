
function c = seawaterSoundSpeed_Leroy1969(t,s,z,lat,varargin)

%**************************************************************************
%  c = seawaterSoundSpeed_Leroy1969(t,s,z,lat,varargin)
%
%  DESCRIPTION: calculates the speed of sound in seawater using the
%  equations proposed by Leroy (1969). The 1st equation is an approximation
%  to Wilson's 2nd equation (Wilson, 1960c). The 2nd equation fits the 
%  measurements from Wilson (1960a), and provides a more accurate and 
%  simpler expression than Wilson's 2nd equation (Wilson, 1960c). The two
%  equations proposed by Leroy provide different levels of accuracy
%  depending on the number of correction terms considered.
%
%  INPUT VARIABLES
%  - t: temperature, in ITS-90 (vector) [ºC]
%  - s: salinity (vector) [ppt]
%  - z: depth below sea surface (vector) [m]
%  - lat: latitude [deg]
%  - eq (varargin{1}): version of Leroy's equation. Options
%    ¬ 1: first equation, approximation to Wilson's 2nd eq. (Wilson,1960c)
%    ¬ 2: second equation, fitting to data from Wilson (1960a) [DEFAULT]
%  - ac (varargin{2}): string specifying the accuracy of selected equation.
%    ¬ 'sim': simplified, 1/5 terms (c = c0)
%    ¬ 'bas': basic, 3/5 terms (c = c0 + ca + cb)
%    ¬ 'com': complete, 5/5 terms (c = c0 + ca + cb + cc + cd) [DEFAULT]
%
%    NOTE: any input ¦t¦, ¦s¦ or ¦z¦ can either be a vector or a number.
%        
%  OUTPUT VARIABLES
%  - c: speed of sound in seawater (vector) [m s-1]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS - General
%  - Leroy (1969) does not mention the temperature scale considered
%    in his equations. The 2nd Wilson equation  (Wilson, 1960c) and his
%    data on sound speed in seawater under pressure (Wilson, 1960a), in 
%    which Leroy's equations are based, are expressed in the International
%    Practical Temperature Scale of 1948 (IPTS-48). Accordingly, Leroy's
%    equations are assumed to be the expressed in IPTS-48. A conversion 
%    is made within the current function between the input temperature, 
%    expressed in the International Temperature Scale of 1990 (ITS-90), 
%    and the temperature in IPTS-48. Errors caused by the use of the wrong
%    temperature scale are generally < 0.05 m/s.
%  - The units used for the variables in the equations are: 
%    t [ºC], s [ppt], z [m], lat [deg].
%
%  CONSIDERATIONS & LIMITATIONS - 1st Equation
%  - Approximation to Wilson's 2nd equation (Wilson, 1960c).
%  - Maximum deviation of 0.2 m/s from 2nd Wilson eq. for combinations of
%    tsp found in Neptunian waters.
%  - For the standard deviation of 0.2 m/s, the accuracy term ¦ac¦ limits
%    the domain of validity to:
%    
%    ¦ac¦        t [ºC]          s [ppt]       z [m]
%   ---------------------------------------------------------------
%    'sim'       -2 to 24.5      30 to 42      0 to 1000
%    'bas'       -2 to 34        25 to 42      0 to 8000
%    'com'       -2 to 34        20 to 42      all depths
%
%  CONSIDERATIONS & LIMITATIONS - 2nd Equation (Wilson, 1960c)
%  - Fitting to measurements from Wilson on sound speed in seawater under
%    pressure (Wilson, 1960a). The resulting equation is more accurate and
%    simpler than Wilson's 2nd equation (Wilson, 1960c).
%  - The standard fitting error of Leroy's eq. to Wilson's measurements
%    is not specified, but is expected to be small.
%  - The domain of validity depends on ¦ac¦ and is assumed to be the same
%    as for 1st Leroy's equation:
%    
%    ¦ac¦        t [ºC]          s [ppt]       z [m]
%   ---------------------------------------------------------------
%    'sim'       -2 to 24.5      30 to 42      0 to 1000
%    'bas'       -2 to 34        25 to 42      0 to 8000
%    'com'       -2 to 34        20 to 42      all depths
%
%  FUNCTION CALLS
%  1) c = seawaterSoundSpeed_Leroy1969(t,s,z,lat)
%     ¬ eq = 2, ac = 'com'
%  2) c = seawaterSoundSpeed_Leroy1969(t,s,z,lat,eq)
%     ¬ ac = 'com'
%  2) c = seawaterSoundSpeed_Leroy1969(t,s,z,lat,eq,ac)
%
%  REFERENCES
%  - Leroy, C.C. (1969) “Development of simple equations for accurate and 
%    more realistic calculation of the speed of sound in seawater” J. 
%    Acoust. Soc. Am. 46, 216–226.
%  - Wilson, W.D. (1960a) “Ultrasonic measurement of the velocity of sound 
%    in distilled and sea water,” Naval Ordnance Report 6746, US Naval 
%    Ordnance Laboratory, White Oak, Maryland.
%  - Wilson, W.D. (1960c). “Equation for the speed of sound in sea water”, 
%    J. Acoust. Soc. Am. 32(10), 1357.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  22 Nov 2017
%
%**************************************************************************

switch nargin
    case {0 1 2 3}
        error('Not enough input arguments')
    case 4
        eq = 2;
        ac = 'com';
    case 5
        eq = varargin{1};
        ac = 'com';
    case 6
        eq = varargin{1};
        ac = varargin{2};       
    otherwise
        error('Too many input arguments')
end

% Variables and Units Conversion
t = (-0.99956 + sqrt(0.9991202 + 1.76e-5*1.00024*t))/(8.8e-6); % temperature in IPTS-48 scale [ºC]
Z = z/1000; % depth [km]
S = s - 35; % salinity referred to 35 ppt

% Speed of Sound Equation
switch eq
    case 1
        c0 = 1492.9 + 3*(t-10) - 6e-3*(t-10).^2 - 4e-2*(t-18).^2 ...
            + 1.2*S - 1e-2*(t-18).*S + z/61;
        switch ac
            case 'sim'
                ca = 0;
                cb = 0;
                cc = 0;
                cd = 0;
            case 'bas'
                ca = 1e-1*Z.^2 +  2e-4*Z.^2.*(t-18).^2 + 1e-1*Z.*lat/90;
                cb = 2e-7*t.*(t-10).^4;
                cc = 0;
                cd = 0;
            case 'com'
                ca = 1e-1*Z.^2 +  2e-4*Z.^2.*(t-18).^2 + 1e-1*Z.*lat/90;
                cb = 2e-7*t.*(t-10).^4;
                cc = -5e-4*Z.^2.*(Z-6).^2;
                cd = 1.5e-3*S.^2.*(1-Z);
            otherwise
                error('Wrong string for input ¦ac¦')
        end 
        
    case 2
        c0 = 1493 + 3*(t-10) - 6e-3*(t-10).^2 - 4e-2*(t-18).^2 ...
            + 1.2*S - 1e-2*(t-18).*S + z/61;
        switch ac
            case 'sim'
                ca = 0;
                cb = 0;
                cc = 0;
                cd = 0;
            case 'bas'
                ca = 1e-1*Z.^2 + 2e-4*Z.^2.*(t-18).^2 + 1e-1*Z.*lat/90;
                cb = 2.6e-4*t.*(t-5).*(t-25);
                cc = 0;
                cd = 0;
            case 'com'
                ca = 1e-1*Z.^2 + 2e-4*Z.^2.*(t-18).^2 + 1e-1*Z.*lat/90;
                cb = 2.6e-4*t.*(t-5).*(t-25);
                cc = -1e-3*Z.^2.*(Z-4).*(Z-8);
                cd = 1.5e-3*S.^2.*(1-Z) + 3e-6*t.^2.*(t-30).*S;
            otherwise
                error('Wrong string for input AC')
        end
    otherwise
        error('Wrong value for input EQ')
end

c = c0 + ca + cb + cc + cd; % sound speed in seawater [m s-1]

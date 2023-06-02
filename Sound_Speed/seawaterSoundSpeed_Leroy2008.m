
function c = seawaterSoundSpeed_Leroy2008(t,s,z,lat)

%**************************************************************************
%  c = seawaterSoundSpeed_Leroy2008(t,s,z,varargin)
%
%  DESCRIPTION: calculates the speed of sound in seawater using the
%  equation proposed by Leroy et al. (2008). The expression fits the NRLII
%  equation (DelGrosso, 1974) as revised by Wong & Zhu (1995) for medium-
%  high salinities (s > 30 ppt), and the so-called "merged equation" for 
%  medium-los salinities (s < 30 ppt), using realistic combinations of 
%  temperature, salinity and depth found in the oceans and seas.
%
%  The "merged equation" is an improved version of the UNESCO equation 
%  (Fofonoff & Millard, 1983; Chen & Millero, 1977), where the incorrect 
%  expression from Wilson (1959)for pure water under pressure is replaced 
%  by a more recent, accurate equation (Belogol'skii et al., 1999).
%
%  The NRLII equation is inaccurate for low saline water at depth (s < 30
%  ppt, z > 0 m), and the accuracy of the UNESCO equation is compromised 
%  by the use of an inaccurate expression for the calculation of the speed 
%  of sound in pure water under pressure (Wilson, 1959). 
% 
%  The 14-term equation of Leroy et al (2008) is obtained from the
%  combination of the NRLII and "merged equation" for their respective
%  domains of applicability (NRLII for s > 30 ppt, "merged equation" for 
%  s < 30 ppt). The equation is a function of temperature t [ºC], salinity
%  s [ppt], depth z [m] and latitude [deg].
%
%  INPUT VARIABLES
%  - t: temperature, in ITS-90 (vector) [ºC]
%  - s: salinity (vector) [ppt]
%  - z: depth below sea surface (vector) [m]
%  - lat: latitude [deg]. If ¦lat¦ is not known, choose a value of 45 deg
%    for an approximate estimation of sound speed.
%
%    NOTE: any input ¦t¦, ¦s¦ or ¦z¦ can either be a vector or a number.
%        
%  OUTPUT VARIABLES
%  - c: speed of sound in seawater (vector) [m s-1]
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - The equation (Leroy et al., 2008) agrees to +/- 0.2 m/s with the two
%    reference complex equations (NRLII for s > 30 ppt, "merged equation"
%    for s < 30 ppt) anywhere in the oceans and seas (including Baltic Sea,
%    Black Sea and Red Sea) and to the greatest depths. It is valid for
%    salinities as high as 42 ppt. The only exceptions to its validity are 
%    abnormal waters with extremely high salinities (e.g. some inland seas 
%    and hot brine spots at the bottom of some seas).
%
%  - The temperature is expressed in the International Salinity Scale of 
%    1990 (IPTS-90). Since the input temperature is expressed in ITS-90, 
%    no conversion is required.
%
%  - The units used for the variables in the equation are: 
%    t [ºC], s [ppt], z [m], lat [deg].
%
%  FUNCTION CALLS
%  1) c = seawaterSoundSpeed_Leroy2008(t,s,z)
%     ¬ lat = 45;
%  2) c = seawaterSoundSpeed_Leroy2008(t,s,z,lat)
%
%  REFERENCES
%  - Leroy, C.C., Robinson, S.P., & Goldsmith, M.J. (2008) “A new equation 
%    for the accurate calculation of sound speed in all oceans”, J. Acoust. 
%    Soc. Am. 124(5), 2774-2782.
%  - Del Grosso, V. A. (1974). “New equation for the speed of sound in 
%    natural waters (with comparisons to other equations)” J. Acoust. Soc.
%    Am. 56, 1084–1091.
%  - Wong, G.S.K., Zhu, S. (1995) “Speed of sound in seawater as a function
%    of salinity, temperature and pressure” J. Acoust. Soc. Am. 97, 2235–
%    2237.
%  - Fofonoff, N.P., & Millard, R.C. (1983) “Algorithm for computation of 
%    fundamental properties of seawater” UNESCO Technical Papers in Marine 
%    Science, No. 44.
%  - Chen, C.T., & Millero, F.J. (1977) “Sound speed in seawater at high 
%    pressures”, J. Acoust. Soc. Am. 62, 1129–1135.
%  - Wilson, W.D. (1959) “Speed of sound in distilled water as a function 
%    of temperature and pressure” J. Acoust. Soc. Am. 31, 1067–1072.
%  - Belogol’skii, V.A., Sekoyan, S.S., Samorukova, L.M., Stefanov, S.R., 
%    & Levtsov, V.I. (1999). “Pressure dependence of the sound velocity in 
%    distilled water,” Meas. Tech. 42, 406–413.
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  22 Nov 2017
%
%**************************************************************************

% Speed of Sound Equation
c = 1402.5 + 5*t - 5.44e-2*t.^2 + 2.1e-4*t.^3 + 1.33*s - 1.23e-2*s.*t ...
    + 8.7e-5*s.*t.^2 + 1.56e-2*z + 2.55e-7*z.^2 - 7.3e-12*z.^3 + 1.2e-6*z.*(lat - 45)...
    - 9.5e-13*t.*z.^3 + 3e-7*t.^2.*z + 1.43e-5*s.*z; % sound speed in seawater [m s-1]

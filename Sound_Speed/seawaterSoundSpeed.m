

function seawaterSoundSpeed(t,s,p)

%    Recommended alternatives are the full equations NRLII (Del Grosso, 
%    1974; Wong & Zhu, 1995) and UNESCO (Chen & Millero, 1977; Fofonoff & 
%    Millard, 1983l Wong & Zhu, 1995), and simple equations from Mackenzie 
%    (1981) and Leroy et al. (2008). 

%  REFERENCES
%  - Del Grosso, V. A. (1974). “New equation for the speed of sound in 
%    natural waters (with comparisons to other equations)” J. Acoust. Soc. 
%    Am. 56, 1084–1091.
%  - Chen, C.T., & Millero, F.J. (1977) “Sound speed in seawater at high 
%    pressures”, J. Acoust. Soc. Am. 62, 1129–1135.
%  - Fofonoff, N.P., & Millard, R.C. (1983) “Algorithm for computation of 
%    fundamental properties of seawater” UNESCO Technical Papers in Marine 
%    Science, No. 44.
%  - Wong, G.S.K., Zhu, S. (1995) “Speed of sound in seawater as a function 
%    of salinity, temperature and pressure” J. Acoust. Soc. Am. 97, 2235–
%    2237.
%  - Mackenzie, K.V. (1981) “Nine term equation for the sound speed in the
%    oceans” J. Acoust. Soc. Am. 70, 807–812.
%  - Leroy, C.C., Robinson, S.P., & Goldsmith, M.J. (2008) “A new equation 
%    for the accurate calculation of sound speed in all oceans”, J. Acoust. 
%    Soc. Am. 124(5), 2774-2782.
%  - Leroy, C.C. (1969) “Development of simple equations for accurate and 
%    more realistic calculation of the speed of sound in seawater” J. 
%    Acoust. Soc. Am. 46, 216–226.
%  - Wilson, W.D. (1960) “Ultrasonic measurement of the velocity of sound 
%    in distilled and sea water,” Naval Ordnance Report 6746, US Naval 
%    Ordnance Laboratory, White Oak, Maryland.
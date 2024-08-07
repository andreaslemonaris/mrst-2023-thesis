-- This reservoir simulation deck is made available under the Open Database
-- License: http://opendatacommons.org/licenses/odbl/1.0/. Any rights in
-- individual contents of the database are licensed under the Database Contents
-- License: http://opendatacommons.org/licenses/dbcl/1.0/

-- Copyright (C) 2015 Statoil
-- Copyright (C) 2015 SINTEF ICT, Applied Mathematics

NOECHO
PLYVISC
0      1.0
0.5    4.0
1.0    8.0
1.5   13.0
2.0   26.0
3.0   52.0 /


PLYROCK
--IPV  RRF  dens  AI  max ads
 0.05  1.3  2600  2  0.000025 /


PLYADS
-- conc  adsorb-conc
 0.000   0.000000
 0.250   0.000014
 0.500   0.000018
 0.750   0.000021
 1.000   0.000022
 1.250   0.000023
 1.500   0.000023
 1.750   0.000027
 3.000   0.000027 /

PLMIXPAR
 1.0  /

PLYMAX
 3.00  0  /

-- PLYSHEAR: defined by a viscosity reduction factor (VRF) as a function of water phase velocity (Vw).
-- Omit shear as CPU-time is significantly reduced and shear thinning is only significant for injectivity which
--  will have to be studied by other tools than Eclipse
PLYSHEAR
-- ---- Vw.            VRF
0           1
0.00002   0.9
0.00008   0.8
0.0004    0.74
0.0008    0.69
0.0017    0.67
0.0043    0.63
0.011     0.60
0.015     0.58
0.019     0.58
0.22      0.48
2.2       0.41   /

ECHO

makeRiver_cr.f90 (ver 1.0)

1. purpose: output river files based on creager-type regression of cross section

2. input file : makeRiverCr.txt

3. format of input file

   L1 : [in] acc file
   L2 : [out] width file
   L3 : [out] depth file
   L4 : [out] height file
   L5 : [in] s_c
   L6 : [in] wd_a, wd_b
   L7 : [in] thresh
   L8 : [in] utm (1: utm, 0: latlon)

4. Formula

1) s = s_c * area ** (area ** (-0.05))
2) wd = wd_a * area ** wd_b

-> width = sqrt( wd * s )
-> depth = sqrt( s / wd )

where

area: catchment area [km2]
s: cross_section_area [m2]
wd: width / depth

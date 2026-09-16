evalPeak.f90 (ver 1.0)

1. purpose: calculate inundation evaluation index

2. input file : evalPeak.txt

3. format of input file

   L1 : [in] simulation data file (RRI output hs)
   L2 : [in] observation data file (ESRI/ASCII Format)
   L3 : [in] evaluation range file

4. format of evaluation range file

   (example.)
   353 410   : imin, imax
   479 571   : jmin, jmax
   0   0     : ishift, jshift

   i : y-direction from top
   j : x-direction from left
   ishift : shift the simulation result on y-direction [defalut : 0]
   jshift : shift the simulation result on x-direction [defalut : 0]

5. note that the threshold to identify the extent of
   simulation result is defined within the program
   (ex. more than 1.0 m in hs is considered to be inundated)

6. index is calculated as

   eff = (S (A) O) / (S (U) O)

   where (A) : product set
   where (U) : union set

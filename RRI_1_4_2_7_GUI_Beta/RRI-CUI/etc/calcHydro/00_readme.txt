calcHydro.f90 (ver 1.0)

1. purpose: calculate hydrograph from RRI output files

2. input file : calcHydro.txt

3. format of input file

   L1 : [in] location file
   L2 : [in] RRI output file (ex. ./out/qr_)
   L3 : [out] hydrograph file

4. format of location

   name(1) loc_i(1) loc_j(1)
   name(2) loc_i(2) loc_j(2)
   ...
   (up to any number of location)

5. Note(1)

   loc_i : row (from top)
   loc_j : column (from left)

   loc_i, loc_j can be calculated from 
   coordinate with "coordinate.xlsx"

rainBasin.f90 (ver 1.0)

1. purpose: 1) calcluate basin avarage rainfall (hyetograph)
            2) calculate total rainfall distribution

2. input file : rainBasin.txt

3. format of input file

   L1 : [in] rainfall file (RRI format) [mm/h]
   L2 : [in] dem file
   L3 : [in] rainfall xll corner
   L4 : [in] rainfall yll corner
   L5 : [in] rainfall cellsize (x, y)
   L6 : [out] hyetograph [mm/h]
   L7 : [out] total rainfall distribution map [mm]
   L8 : [out] cumulative rainfall [mm]
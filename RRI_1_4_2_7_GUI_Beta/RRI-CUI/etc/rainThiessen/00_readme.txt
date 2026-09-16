rainThiessen.f90 (ver 1.1)

1. purpose: Interpolate gauged rainfall data 
            by Thiessen Polygon to prepare 2D rainfall

2. input file : rainThiessen.txt

3. format of input file

   L1 : [in] gauged rainfall data file [mm/h] or [mm/d]
   L2 : [in] divide parameter (set 1 if input is [mm/h], 24 if input is [mm/d])
   L3 : [out] outfile
   L4 : [out] outfile_map
   L5 : [in] ncols ***
   L6 : [in] nrows ***
   L7 : [in] xll ***
   L8 : [in] yll ***
   L9 : [in] cellsize ***

where *** must be replaced by values

4. format of gauged rainfall data file specified in L1

   See Section 4 in RRI_Manual.pdf for the detail

5. notice

When the input data is daily data, set "divide = 24", 
so that the rainThiessen program divides the output data 
by 24 to obtain rainfall intensity in (mm/h) for RRI Model.


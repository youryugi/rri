scaleUp.f90 (ver 1.1)

1. purpose: scale up DEM, DIR, ACC data

2. input file : scaleUp.txt

3. format of input file

   L1 : [in] original DEM file
   L2 : [in] original DIR file
   L3 : [in] original ACC file
   L4 : [in] ups : a number of grid-cell to be integrated (ex. 2, 3, ...)
   L5 : [out] upscaled DEM file
   L6 : [out] upscaled DIR file
   L7 : [out] upscaled ACC file

5. notice

The program upscale a original flow network to coaser spatial resolution. The upscaled resolution become "ups" times larger than the original resolution. The algorithm follows the method proposed by Masutani et al. (2006).
.

Note that if a target river basin has upstream catchments, extra "acc" is added automatically (after ver 1.1).
[Reference]

K. Masutani, K. Akai, J. Magome: A New Scaling Algorithm of Gridded River Networks, J. Japan Soc. Hydrol. & Water Resour., Vol. 19, No. 2, 2006.

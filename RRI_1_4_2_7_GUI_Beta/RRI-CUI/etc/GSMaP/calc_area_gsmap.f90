! Compute suitable range of GSMaP data for a target catchment

integer temporal_resolution

integer ncols, nrows
real xllcorner, yllcorner, cellsize
integer itop, ibottom, jleft, jright
integer xll_negative ! added on Dec 7, 2022 (1: negatie, 0: positive)

! IMPORTANT: Spatial Resolution of Rainfall Data 
real cellsize_rain
!parameter( cellsize_rain = 0.1 ) ! Hourly
!parameter( cellsize_rain = 0.25 ) ! Daily

real xul_data, yul_data ! Origin of the Data
!parameter( xul_data = 0.05 )   ! Hourly-0.1deg
!parameter( yul_data = 59.95 )  ! Hourly-0.1deg
!parameter( xul_data = 0.125 )   ! Daily-0.25deg
!parameter( yul_data = 59.875 )  ! Daily-0.25deg

real xur, yur
real xll_rain, yll_rain
real xur_rain, yur_rain
real xllcorner_rain, yllcorner_rain

character*256 infile

! read
open(10, file = "calc_area_gsmap.txt", status = 'old')
read(10, *) cellsize_rain
read(10, *) temporal_resolution
read(10, *) ncols
read(10, *) nrows
read(10, *) xll
read(10, *) yll
read(10, *) cellsize

if( xll .lt. 0 )then
 xll_negative = 1
 xll = xll + 360.0 ! the range of topographic data: xll = -180 ~ 180
                   ! the range of GSMaP rainfall data: xll = 0 ~ 360
else
 xll_negative = 0
endif

xul_data = cellsize_rain / 2.
yul_data = 60. - cellsize_rain / 2.

! calc extent

!jleft = (xll - xul_data) / cellsize_rain
! modified on Aug 20, 2021
jleft = (xll - xul_data) / cellsize_rain + 1
xll_rain = xul_data + jleft * cellsize_rain

!ibottom = (yll - yul_data) / cellsize_rain
! modified on Aug 20, 2021
ibottom = (yll - yul_data) / cellsize_rain - 1
if(ibottom.lt.0) ibottom = ibottom - 1
yll_rain = yul_data + ibottom * cellsize_rain

xur = xll + ncols * cellsize
jright = (xur - xul_data) / cellsize_rain + 1
xur_rain = xul_data + jright * cellsize_rain

yur = yll + nrows * cellsize
itop = (yur - yul_data) / cellsize_rain
if(ibottom.gt.0) ibottom = ibottom + 1
yur_rain = yul_data + itop * cellsize_rain

! output
write(*,*) "horizontal_resolution [d] : ", cellsize_rain
write(*,*) "temporal_resolution [h]   : ", temporal_resolution
write(*,*) 

write(*,*) "xll : ", xll
write(*,*) "yll : ", yll
write(*,*) "xur : ", xur
write(*,*) "yur : ", yur
write(*,*) 

write(*,*) "xll_rain : ", xll_rain
write(*,*) "yll_rain : ", yll_rain
write(*,*) "xur_rain : ", xur_rain
write(*,*) "yur_rain : ", yur_rain
write(*,*) 

write(*,*) "jleft   : ", jleft
write(*,*) "ibottom : ", -ibottom
write(*,*) "jright  : ", jright
write(*,*) "itop    : ", -itop
write(*,*) 

if( xll_negative .eq. 0 )then
 write(*,*) "xllcorner_rain (raster): ", xll_rain - cellsize_rain * 0.5
else
 write(*,*) "xllcorner_rain (raster): ", xll_rain - cellsize_rain * 0.5 - 360.0 ! modified on Dec. 7, 2022
endif
write(*,*) "yllcorner_rain (raster): ", yll_rain - cellsize_rain * 0.5
write(*,*) "cellsize_rain : ", cellsize_rain
write(*,*) 

open(20, file = "out_by_calc_area_gsmap.txt")
write(20,*) "horizontal_resolution [d] : ", cellsize_rain
write(20,*) "temporal_resolution [h]   : ", temporal_resolution
write(20,*) 

write(20,*) "xll : ", xll
write(20,*) "yll : ", yll
write(20,*) "xur : ", xur
write(20,*) "yur : ", yur
write(20,*) 

write(20,*) "xll_rain : ", xll_rain
write(20,*) "yll_rain : ", yll_rain
write(20,*) "xur_rain : ", xur_rain
write(20,*) "yur_rain : ", yur_rain
write(20,*) 

write(20,*) "jleft   : ", jleft
write(20,*) "ibottom : ", -ibottom
write(20,*) "jright  : ", jright
write(20,*) "itop    : ", -itop
write(20,*) 

if( xll_negative .eq. 0 )then
 write(20,*) "xllcorner_rain (raster): ", xll_rain - cellsize_rain * 0.5
else
 write(20,*) "xllcorner_rain (raster): ", xll_rain - cellsize_rain * 0.5 - 360.0 ! modified on Dec. 7, 2022
endif
write(20,*) "yllcorner_rain (raster): ", yll_rain - cellsize_rain * 0.5
write(20,*) "cellsize_rain : ", cellsize_rain
write(20,*) 

end

! Calculate extraction area of 3B42RT data

integer ncols, nrows
real xllcorner, yllcorner, cellsize
integer itop, ibottom, jleft, jright

real cellsize_rain ! spatial resolution of rainfall data
!parameter( cellsize_rain = 0.25 )

real xul_data, yul_data ! origin of rainfall data
!parameter( xul_data = 0.125 )
!parameter( yul_data = 59.875 ) 

real xur, yur
real xll_rain, yll_rain
real xur_rain, yur_rain
real xllcorner_rain, yllcorner_rain

character*256 infile

! read
open(10, file = "calc_area.txt", status = 'old')
read(10, *) cellsize_rain
read(10, *) xul_data
read(10, *) yul_data
read(10, *) ncols
read(10, *) nrows
read(10, *) xll
read(10, *) yll
read(10, *) cellsize

! calc extent

jleft = (xll - xul_data) / cellsize_rain
xll_rain = xul_data + jleft * cellsize_rain

ibottom = (yll - yul_data) / cellsize_rain
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

write(*,*) "xllcorner_rain (raster): ", xll_rain - cellsize_rain * 0.5
write(*,*) "yllcorner_rain (raster): ", yll_rain - cellsize_rain * 0.5
write(*,*) "cellsize_rain : ", cellsize_rain
write(*,*) 

open(20, file = "out_by_calc_area.txt")
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

write(20,*) "xllcorner_rain (raster): ", xll_rain - cellsize_rain * 0.5
write(20,*) "yllcorner_rain (raster): ", yll_rain - cellsize_rain * 0.5
write(20,*) "cellsize_rain : ", cellsize_rain
write(20,*) 

end

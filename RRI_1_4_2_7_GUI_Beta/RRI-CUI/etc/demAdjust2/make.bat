del *.obj *.mod
ifort /O3 demAdjust2.f90
editbin /stack:40000000 demAdjust2.exe
del *.obj *.mod

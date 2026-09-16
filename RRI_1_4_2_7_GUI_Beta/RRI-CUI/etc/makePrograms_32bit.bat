ifort /o ../bin32/calcHydro.exe ./calcHydro/calcHydro.f90 
ifort /o ../bin32/calcPeak.exe ./calcPeak/calcPeak.f90
ifort /o ../bin32/calcTecplot.exe ./calcTecplot/calcTecplot.f90
ifort /o ../bin32/demAdjust2.exe ./demAdjust2/demAdjust2.f90
ifort /o ../bin32/evalHydro.exe ./evalHydro/evalHydro.f90
ifort /o ../bin32/evalPeak.exe ./evalPeak/evalPeak.f90
ifort /o ../bin32/rainBasin.exe ./rainBasin/rainBasin.f90
ifort /o ../bin32/rainThiessen.exe ./rainThiessen/rainThiessen.f90
ifort /o ../bin32/scaleUp.exe ./scaleUp/scaleUp.f90
ifort /o ../bin32/calc_area_gsmap.exe ./GSMaP/calc_area_gsmap.f90
ifort /o ../bin32/read_gsmap.exe ./GSMaP/read_gsmap.f90
ifort /o ../bin32/makeRiver3.exe ./makeRiver3/makeRiver3.f90
ifort /o ../bin32/flowDirection.exe ./flowDirection/flowDirection.f90

del *.obj *.mod

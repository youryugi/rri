ifort /o ../bin/calcHydro.exe ./calcHydro/calcHydro.f90 
ifort /o ../bin/calcPeak.exe ./calcPeak/calcPeak.f90
ifort /o ../bin/calcTecplot.exe ./calcTecplot/calcTecplot.f90
ifort /o ../bin/demAdjust2.exe ./demAdjust2/demAdjust2.f90
ifort /o ../bin/evalHydro.exe ./evalHydro/evalHydro.f90
ifort /o ../bin/evalPeak.exe ./evalPeak/evalPeak.f90
ifort /o ../bin/rainBasin.exe ./rainBasin/rainBasin.f90
ifort /o ../bin/rainThiessen.exe ./rainThiessen/rainThiessen.f90
ifort /o ../bin/scaleUp.exe ./scaleUp/scaleUp.f90
ifort /o ../bin/calc_area_gsmap.exe ./GSMaP/calc_area_gsmap.f90
ifort /o ../bin/read_gsmap.exe ./GSMaP/read_gsmap.f90
ifort /o ../bin/makeRiver3.exe ./makeRiver3/makeRiver3.f90
ifort /o ../bin/flowDirection.exe ./flowDirection/flowDirection.f90

del *.obj *.mod


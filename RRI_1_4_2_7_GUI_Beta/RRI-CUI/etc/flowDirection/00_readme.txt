flowDirection.f90 (ver 1.0)

1. purpose: calculate flow direction from dem

2. input file : flowDirection.txt

3. format of input file

   L1 : [in] original dem file
   L2 : [out] dir file
   L3 : [out] acc file
   L4 : [in] [1] read riv file [0] not read riv file
   L5 : [in] riv file basin

4. Note

If riv file is read with L4 = 1, 
flow directiosns on rivers are decided first, then 
flod directions on other grid-cells are decided.

In riv file,
100 : river grid-cells, whose directions will be calculated by this program.
0   : river outlet cell.
1, 2, 4, 8, 16, 32, 64, 128 : flow directions are pre-defined (no change)



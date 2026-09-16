demAdjust2.f90

1. purpose: 1) DEM adjustment by carving and filling, and
            2) zero setting of flow direction at outlet cells

2. input file : demAdjust2.txt

3. format of input file

   L1 : [in] dem file
   L2 : [in] dir file
   L3 : [in] acc file
   L4 : [out] adjusted dem file (adem)
   L5 : [out] modified dir file (adir)

4. About DEM adjustment

   The program adjusts DEM by carving and filling to remove pits along flow line.
   This DEM adjustment is necessary to avoid unrealistic discontinuity of flow.
   The "adem2" should be used as the input of RRI model.

5. About Flow Direction Zero Setting

   RRI requires to specify outlet cells.
   By setting zero values at outlet cells, RRI recognizes the locations.
   demAdjustment2 program sets flow direction zero when it finds a grid that
   is directed to the outlet of the simulation domain.
   The "adem2" must be used as the input of RRI model.

6. Algorithm of DEM adjustment

   Based on the flow direction, demAdjust2 finds upstream cells (i.e. cells with no inflow).

   Among the detected upstream cells, searching order is determined from the total length of the flow paths from each upstream cell to its most downstream cell.

   Following the above decided order, demAdjust2 adjusts elevations based on the following procedures.
  1) The negative elevation is set to be zero.

  2) Lifting: If a single cell is extremely low (likely as a noise error) compared to its upstream and downstream cells, the cellÅfs elevation will be replaced by the same elevation as the upstream cell. The parameter ÅgliftÅh is used as the threshold to detect sudden drop and its default value is set to be 500 m.

  3) Carving: If the elevation suddenly increases along the flow direction, the cellÅfs elevation will be replaced by the same elevation as the upstream cell. The parameter ÅgcarveÅh is used as the threshold to detect the sudden increase and its default value is 5 m.

  4) Lifting and Carving: By searching from the most upstream, it finds a cell whose downstream elevation is higher than that cell (point L). By searching from point L toward downstream, it finds a cell whose downstream is lower than that cell (point H). 
The point L is lifted and point H is carved by the parameter ÅgincrementÅh, whose default is 0.01 m.

  The demAdjust2 program conducts each of the above procedure repeatedly for each flow path ways from all the detected upstream cells until all negative slopes are removed. Note that the above procedure does not change flow direction.

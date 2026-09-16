        !COMPILER-GENERATED INTERFACE MODULE: Sat Apr 02 21:16:31 2016
        MODULE DOWN__genmod
          INTERFACE 
            SUBROUTINE DOWN(DIR,DEM,NX,NY,I,J,LENGTH,II,JJ,DIS,SLO,WID)
              INTEGER(KIND=4) :: NY
              INTEGER(KIND=4) :: NX
              INTEGER(KIND=4) :: DIR(NY,NX)
              REAL(KIND=8) :: DEM(NY,NX)
              INTEGER(KIND=4) :: I
              INTEGER(KIND=4) :: J
              REAL(KIND=8) :: LENGTH
              INTEGER(KIND=4) :: II
              INTEGER(KIND=4) :: JJ
              REAL(KIND=8) :: DIS
              REAL(KIND=8) :: SLO
              REAL(KIND=8) :: WID
            END SUBROUTINE DOWN
          END INTERFACE 
        END MODULE DOWN__genmod

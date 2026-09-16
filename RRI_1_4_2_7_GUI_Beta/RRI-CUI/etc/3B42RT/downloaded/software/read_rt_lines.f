        PROGRAM READ_RT_LINES
C*******************************************************************
C* READ_RT_LINES  Read by lines for real-time TRMM and Other Data  *
C*                                                                 *
C* This program provides an example of line-by-line reading for    *
C* each of the three data files produced by the real-time TRMM and *
C* Other Data precipitation estimation system.  It is set up to    *
C* prompt the user for a file name and data type, then the header  *
C* and values of all fields are output for the four boxes          *
C* surrounding the intersection of the Date Line and Equator.      *
C* NOTES:                                                          *
C* 1. The units for specifying FORTRAN unformatted direct access   *
C*    record size are system-dependent.  On SGI it's in 4-byte     *
C*    words, but other machines use the number of bytes.  These	   *
C*    are accommodated by setting RECLSZ in the PARAMETER state-   *
C*    ment to 'WORD' or 'BYTE' respectively.		           *
C* 2. The data are written big-endian, so little-endian machines   *
C*    might require byte swapping, unless gzip has successfully    *
C*    handled it.                                                  *
C**                                                                *
C* Log:                                                            *
C* G.Huffman/SSAI       12/01                                      *
C* G.Huffman/SSAI	12/05	Revise setting RECL		   *
C* G.Huffman/SSAI	12/08	Corrected use of RECLSZ; add       *
C*				commented-out byte swapping        *
C* G.Huffman/SSAI	1/09	Upgrade to Feb. 2009 TMPA-RT       *
C*******************************************************************
        IMPLICIT        NONE
        INTEGER*4       LN, LT90, LT60
	CHARACTER*4	RECLSZ
        PARAMETER       ( LN     = 1440,
     +                    LT90   = 720,
     +                    LT60   = 480,
     +			  RECLSZ = 'WORD' )
C
C       The arrays are dimensioned for fully global fields to accomodate
C       3B40RT.  3B41RT and 3B42RT only cover the latitude band 60N-S,
C       which is accounted for by reading only 480 rows.
C       Arrays: p       - estimated precipitation (mm/hr)
C               e       - estimated precipitation random error (mm/hr)
C               npixels - number of original data pixels contributing to
C                         gridbox estimates
C               nambig  - number of original data pixels that are
C                         classified as "ambiguous"
C               nrain   - number of original data pixels that are
C                         classified as raining
C               nsource - source of estimate (<0, 0, 100 = none, HQ, IR)
C               p_uncal - estimated precipitation without climatological
C                         calibration to 3B42 V.6 (mm/hr)
C
        REAL*4          p       (LN, LT90), e       (LN, LT90),
     +                  p_uncal (LN, LT90), r4dum
        INTEGER*1       npixel  (LN, LT90), nambig  (LN, LT90),
     +                  nrain   (LN, LT90), nsource (LN, LT90),
     +                  i1dum
        INTEGER*4       iret, i, j, i1, j1
        CHARACTER*256   type, name
        CHARACTER*2880  head
C
C       Prompt for file name.
C
        WRITE (*, *)  'Enter file name:'
        READ (*, 10) name
   10   FORMAT ( A256 )
C
C       Prompt for data type.
C
   15   WRITE (*, *) 'Enter data type (3B40RT, 3B41RT, 3B42RT):'
        READ (*, 10) type
        IF  ( type (1:6) .NE. '3B40RT' .AND. 
     +        type (1:6) .NE. '3B41RT' .AND.
     +        type (1:6) .NE. '3B42RT' )  THEN
            WRITE (*, *) type (1:6), 'isn''t a choice'
            GO TO 15
        END IF
C
C       Call the example reading routine and write out the values of
C       all five fields for the four boxes surrounding the intersection
C       of the Date Line and the Equator.  The undimensioned dummy
C       arguments i1dum and r4dum work because the subroutine is
C       guaranteed to never touch them.
C
C       Case 1: 3B40RT - merged microwave, or high quality (HQ)
C
        IF  ( type (1:6) .EQ. '3B40RT' )  THEN
            CALL GET_ARRS_LINES ( name,    type,   LN,    LT90, 
     +				  RECLSZ,  head,   p, e,
     +                            npixel,  nambig, nrain, 
     +                            nsource, r4dum,  iret )
            IF  ( iret .LT. 0 )  STOP
C
C         Case 2: 3B41RT - variable rainrate infrared (VAR).
C
          ELSE IF  ( type (1:6) .EQ. '3B41RT' )  THEN
            CALL GET_ARRS_LINES ( name,   type,  LN,    LT60,  
     +				  RECLSZ, head,  p, e,
     +                            npixel, i1dum, i1dum, 
     +                            i1dum,  r4dum, iret )
            IF  ( iret .LT. 0 )  STOP
C
C         Case 3: 3B42RT - VAR merged with HQ (VARHQ).
C
          ELSE IF  ( type (1:6) .EQ. '3B42RT' )  THEN
            CALL GET_ARRS_LINES ( name,    type,    LN,    LT60,  
     +				  RECLSZ,  head,    p, e,
     +                            nsource, i1dum,   i1dum, 
     +                            i1dum,   p_uncal, iret )
            IF  ( iret .LT. 0 )  STOP
C
C         Case 4: Any other input is not recognized.
C
          ELSE
            WRITE (*, *) 'Error - don''t recognize data type ',
     +                   type (1:6)
            STOP
        END IF
C
C       Write out the header and values at the points around the
C       intersection of the Date Line and Equator.
C
        WRITE (*, *) 'Header:'
        WRITE (*, *) head
C
        IF  ( type (1:6) .EQ. '3B40RT' )  THEN
            WRITE (*, *) 'i, j, p, e, npixel, nambig, nrain, nsource:'
            i1 = LN   / 2
            j1 = LT90 / 2
            DO  j = j1, j1+1
            DO  i = i1, i1+1
                WRITE (*, 100) i, j, p (i, j), e (i, j), npixel (i, j),
     +                         nambig (i, j), nrain (i, j), 
     +                         nsource (i, j)
  100           FORMAT ( 2I5, 2F10.2, 4I4 )
            END DO
            END DO
          ELSE IF  ( type (1:6) .EQ. '3B41RT' )  THEN
            WRITE (*, *) 'i, j, p, e, npixel:'
            i1 = LN   / 2
            j1 = LT60 / 2
            DO  j = j1, j1+1
            DO  i = i1, i1+1
                WRITE (*, 110) i, j, p (i, j), e (i, j), npixel (i, j)
  110           FORMAT ( 2I5, 2F10.2, I4 )
            END DO
            END DO
          ELSE
            WRITE (*, *) 'i, j, p, e, nsource, p_uncal:'
            i1 = LN   / 2
            j1 = LT60 / 2
            DO  j = j1, j1+1
            DO  i = i1, i1+1
                WRITE (*, 120) i, j, p (i, j), e (i, j), 
     +                         nsource (i, j), p_uncal (i, j)
  120           FORMAT ( 2I5, 2F10.2, I4, F10.2 )
            END DO
            END DO
        END IF
C
        STOP
        END
C
        SUBROUTINE GET_ARRS_LINES  ( name, type,  LN, LT,  RECLSZ,
     +				     head, p, e,  np1,     np2,
     +                               np3,  np4,   p_uncal, iret )
C********************************************************************
C* GET_VAR_LINES  Do line-by-line read of a real-time file          *
C*                                                                  *
C* This subroutine opens the specified TRMM and Other Data          *
C* precipitation estimate file and reads the specified arrays a     *
C* line at a time, converting the first two from scaled I*2 to R*4. *
C*                                                                  *
C* GET_ARRS_LINES  ( name, type, LN,  LT,  RECLSZ,  head, p, e,     *
C*                   np1,  np2,  np3, np4, p_uncal, iret )          *
C*                                                                  *
C* Input parameters:                                                *
C*   name       CHAR*   File name                                   *
C*   type       CHAR*   Product I.D. ('3B40RT', '3B41RT', '3B42RT') *
C*   LN         INT     X size                                      *
C*   LT         INT     Y size                                      *
C*   RECLSZ	CHAR*4	Size of a RECL unit; 4-byte word, 1 byte    *
C*			are denoted by 'WORD', 'BYTE' respectively  *
C*                                                                  *
C* Output parameters:                                               *
C*   head       CHAR*   File header                                 *
C*   p          REAL*   Precip field in mm/hr (size LN,LT)          *
C*   e          REAL*   Precip error field in mm/hr (size LN,LT)    *
C*   np1        BYTE*   Precip companion field 1 (number or source) *
C*                      (size LN,LT)                                *
C*   np2        BYTE*   Precip companion field 2 (number or source) *
C*                      (size LN,LT)                                *
C*   np3        BYTE*   Precip companion field 3 (number or source  *
C*                      (size LN,LT)                                *
C*   np4        BYTE*   Precip companion field 4 (number or source  *
C*                      (size LN,LT)                                *
C*   p_uncal    REAL*   Precip field without climatological         *
C*                      calibration to 3B42V.6 in mm/hr (size LN,LT)*
C*   iret       INT     Return code: 0 = normal                     *
C*                                  -1 = local array too small      *
C*                                  -2 = I/O error                  *
C*                                  -3 = header array too short     *
C**                                                                 *
C* Log:                                                             *
C* G.Huffman/SSAI       12/01   Adapt GET_ARRS.F                    *
C* G.Huffman/SSAI	1/09	Upgrade to Feb. 2009 TMPA-RT        *
C********************************************************************
        IMPLICIT        NONE
        REAL*4          RMISS
        INTEGER*4       LNL, LTL
        INTEGER*2       I2MISS
        PARAMETER       ( LNL    = 1440,
     +                    LTL    = 720,
     +                    RMISS  = -99999.,
     +                    I2MISS = -31999 )
C
        INTEGER*4       LN, LT, iret
        REAL*4          p   (LN, LT), e   (LN, LT), p_uncal (LN, LT)
        BYTE            np1 (LN, LT), np2 (LN, LT), np3 (LN, LT),
     +                  np4 (LN, LT)
        CHARACTER*(*)   name, type, RECLSZ, head
C
        INTEGER*4       iunit, l_n, i, j, l_head, lrec, iostat, joff
        INTEGER*2       i2  (LNL, LTL)
        CHARACTER*1     dum
        LOGICAL         in_use
cc
cc	If byte swapping is required, uncomment this section.
cc
cc	INTEGER*2	ivarin, ivar
cc	CHARACTER*1	cvarin (2), cvar (2)
cc	EQUIVALENCE (cvarin, ivarin)
cc	EQUIVALENCE (cvar,   ivar)
cc
C       Check that the local array is big enough.
C
        IF  ( LNL .LT. LN .OR. LTL .LT. LT )  THEN
            WRITE (*, *) 'Error: local work array too small, with',
     +                   LNL, ',', LTL, ' .LT.', LN, ',', LT
            iret = -1
            RETURN
        END IF
C
C       Open the file.  Start by finding a free FORTRAN unit number.
C       The logical record is one row of I*2, or LN * 2 bytes long.
C       ==> The units for specifying FORTRAN unformatted direct   <==
C       ==> access record size are system-dependent.  On SGI it's <==
C       ==> in 4-byte words, but other machines use the number of <==
C       ==> bytes.  These choices are denoted by RECLSZ values    <==
C	==> 'WORD','BYTE' respectively.                           <==
C
        in_use = .TRUE.
        iunit  = 2
        DO  WHILE ( in_use )
            iunit = iunit + 1
            INQUIRE ( iunit, OPENED=in_use )
        END DO
C
        lrec = LN / 2
	IF  ( RECLSZ .EQ. 'BYTE' )  lrec = lrec * 4
        OPEN  ( UNIT = iunit, ACCESS = 'DIRECT',
     +          FILE = name,  FORM   = 'UNFORMATTED',
     +          RECL = lrec,  IOSTAT = iostat,
     +          ERR  = 15 )
        GO TO 16
   15       WRITE (*, *) 'Error', iostat, ' opening ', name
            iret = -2
            RETURN
   16   CONTINUE
C
C       Get a 1-row header. It's an error if the receiving array is
C       short.
C
        l_head = LEN ( head )
        IF  ( l_head .GE. LN * 2 )  THEN
            READ ( iunit, REC=1, ERR=100, IOSTAT=iostat ) head
          ELSE
            WRITE (*, *) 'Error: header size (', l_head,
     +                   ' bytes) .LT. one row (', LN*2, ' bytes)'
            iret = -3
            RETURN
        END IF
C
C       Read and convert the precip from scaled (by 100) I*2.
C
        joff = 1
        DO  j = 1, LT
            READ ( iunit, REC=j+joff, ERR=100, IOSTAT=iostat )
     +          ( i2 (i, j), i = 1, LN )
        END DO
cc
cc	If byte swapping is required, uncomment this section.
cc
cc	DO  j = 1, LT
cc	DO  i = 1, LN
cc	    ivarin    = i2 (i, j)
cc	    cvar (1)  = cvarin (2)
cc	    cvar (2)  = cvarin (1)
cc	    i2 (i, j) = ivar
cc	END DO
cc	END DO
cc
        DO  j = 1, LT
        DO  i = 1, LN
            IF  ( i2 ( i, j ) .EQ. I2MISS )  THEN
                p (i, j) = RMISS
              ELSE
                p (i, j) = i2 ( i, j ) * .01
            END IF
        END DO
        END DO
C
C       Read and convert the error from scaled (by 100) I*2.
C
        joff = 1 + LT
        DO  j = 1, LT
            READ ( iunit, REC=j+joff, ERR=100, IOSTAT=iostat )
     +          ( i2 (i, j), i = 1, LN )
        END DO
cc
cc	If byte swapping is required, uncomment this section.
cc
cc	DO  j = 1, LT
cc	DO  i = 1, LN
cc	    ivarin    = i2 (i, j)
cc	    cvar (1)  = cvarin (2)
cc	    cvar (2)  = cvarin (1)
cc	    i2 (i, j) = ivar
cc	END DO
cc	END DO
cc
        DO  j = 1, LT
        DO  i = 1, LN
            IF  ( i2 ( i, j ) .EQ. I2MISS )  THEN
                e (i, j) = RMISS
              ELSE
                e (i, j) = i2 ( i, j ) * .01
            END IF
        END DO
        END DO
C
C       Proceed to read the companion I*1 fields, as requested.  Now
C       each READ pulls in two rows.  All types have the first, but
C       the next 3 are for 3B40RT only.
C
        joff = 1 + 2 * LT
        DO  j = 1, LT / 2
            READ ( iunit, REC=j+joff, ERR=100, IOSTAT=iostat )
     +          ( np1 (i, j*2-1), i = 1, LN ),
     +          ( np1 (i, j*2  ), i = 1, LN )
        END DO
C
        IF  ( type (1:6) .EQ. '3B40RT' )  THEN
            joff = 1 + 2.5 * LT
            DO  j = 1, LT / 2
                READ ( iunit, REC=j+joff, ERR=100, IOSTAT=iostat )
     +              ( np2 (i, j*2-1), i = 1, LN ),
     +              ( np2 (i, j*2  ), i = 1, LN )
            END DO
C
            joff = 1 + 3 * LT
            DO  j = 1, LT / 2
                READ ( iunit, REC=j+joff, ERR=100, IOSTAT=iostat )
     +              ( np3 (i, j*2-1), i = 1, LN ),
     +              ( np3 (i, j*2  ), i = 1, LN )
            END DO
C
            joff = 1 + 3.5 * LT
            DO  j = 1, LT / 2
                READ ( iunit, REC=j+joff, ERR=100, IOSTAT=iostat )
     +              ( np4 (i, j*2-1), i = 1, LN ),
     +              ( np4 (i, j*2  ), i = 1, LN )
            END DO
        END IF
C
C       Type 3B42RT has p_uncal following one companion field, so
C       the offset looks like that for companion field 2, above.
C       Read and convert p_uncal from scaled (by 100) I*2.
C
        IF  ( type (1:6) .EQ. '3B42RT' )  THEN
            joff = 1 + 2.5 * LT
            DO  j = 1, LT
                READ ( iunit, REC=j+joff, ERR=100, IOSTAT=iostat )
     +              ( i2 (i, j), i = 1, LN )
            END DO
cc
cc	If byte swapping is required, uncomment this section.
cc
cc	DO  j = 1, LT
cc	DO  i = 1, LN
cc	    ivarin    = i2 (i, j)
cc	    cvar (1)  = cvarin (2)
cc	    cvar (2)  = cvarin (1)
cc	    i2 (i, j) = ivar
cc	END DO
cc	END DO
cc
            DO  j = 1, LT
            DO  i = 1, LN
                IF  ( i2 ( i, j ) .EQ. I2MISS )  THEN
                    p_uncal (i, j) = RMISS
                  ELSE
                    p_uncal (i, j) = i2 ( i, j ) * .01
                END IF
            END DO
            END DO
        END IF
C
        CLOSE ( iunit )
        iret = 0
        RETURN
C
C       I/O error return.
C
  100   WRITE (*, *)  'Read error', iostat, ' on', name
        CLOSE ( iunit )
        iret = -2
        RETURN
        END

        PROGRAM READ_HEADER
C********************************************************************
C* READ_HEADER  Read header for real-time TRMM Multi-Satellite Data *
C*                                                                  *
C* This program reads the header record and file arrays of KEYWORD  *
C* and VALUE.                                                       *
C                                                                   *
C* The header is written in a "PARAMETER=VALUE" format, where       *
C* PARAMETER is a string without embedded blanks that gives the     *
C* parameter name, VALUE is a string that gives the value of the    *
C* parameter, and blanks separate each "PARAMETER=VALUE" set.  To   *
C* prevent ambiguity, no spaces or "=" are permitted as characters  *
C* in either PARAMETER or VALUE.                                    *
C*                                                                  *
C* The data arrays are dimensioned large enough that we don't have  *
C* to be careful about overflows; they could be reduced if space    *
C* is short.                                                        *
C**                                                                 *
C* Log:                                                             *
C* G.Huffman/SSAI       12/01                                       *
C* G.Huffman/SSAI        2/02   Documentation                       *
C* G.Huffman/SSAI	12/05	Revise setting RECL		    *
C********************************************************************
        IMPLICIT        NONE
	CHARACTER*4	RECLSZ
	PARAMETER	( RECLSZ = 'WORD' )
C
        CHARACTER*2880  header
        CHARACTER*80    keywd (50), value (50), file
        INTEGER         neq   (50), kstrt (50), nvend (50)
        INTEGER         iret, i, l_header, ipt, in, numkey, j,
     +			lrec
C
C       Open the data file (using January 1997 as an example) with a
C       RECL of 1 data row.
C	          ==>>    WARNING WARNING WARNING    <<==
C       The RECL is defined differently on different machines; it isn't
C       specified in the FORTRAN77 standard.  On SGI it's in 4-B words.
C       If you find that you only get some good values and then garbage
C       (either all zeros or random values), your machine probably
C       wants RECL in bytes, and you should say RECLSZ='BYTE' in the
C       preceding PARAMETER statement.  [The header is 2880 bytes.]
C
        file = '3B40RT.2002020500.bin'
	lrec = 720
	IF  ( LRECSZ .EQ. 'BYTE' )  lrec = lrec * 4
        OPEN  ( UNIT=10, FILE=file, ACCESS='DIRECT',
     +          FORM='UNFORMATTED', STATUS='OLD', RECL=lrec,
     +          IOSTAT=iret )
        IF  ( iret .NE. 0 )  THEN
            WRITE (*, *) 'Error: open error', iret,
     +                   ' on file ', file
            STOP
        END IF
C
C       Read the header (the first record) and close the file.
C
        READ ( UNIT=10, REC=1, IOSTAT=iret )  header
        IF  ( iret .NE. 0 )  THEN
            WRITE (*, *) 'Error: read error', iret,
     +                   ' on file ', file
            STOP
        END IF
        CLOSE ( UNIT=10 )
C
C       Find the actual length of the header (as opposed to the
C       declared FORTRAN size) by parsing back from the end for the
C       first non-blank character (it was written blank-filled).
C
        DO  10 i = 1, 2880
            IF  ( header (2881-i:2881-i) .NE. ' ' )  GO TO 20
   10   CONTINUE
        WRITE (*, *) 'Error: found no non-blanks in the header'
        STOP
   20   l_header = 2881 - i
C
C       Parse for "=".
C
        ipt = 1
        DO  30 i = 1, l_header
            in = INDEX ( header (ipt:l_header), '=' )
            IF  ( in .EQ. 0 )  THEN
                GO TO 40
              ELSE
                neq (i) = ipt + in - 1
                ipt     = ipt + in
            END IF
   30   CONTINUE
        WRITE (*, *) 'Error: ran through header without ending parsing'
        STOP
   40   CONTINUE
        numkey = i - 1
C
C       Now find corresponding beginning of each keyword by parsing
C       backwards for " ".  The first automatically starts at 1.  We
C       assume that there are at least 2 keywords!
C
        kstrt (1) = 1
        DO  60 i = 2, numkey
            DO  50 j = 1, neq (i) - 1
                IF  ( header (neq(i)-j:neq(i)-j) .EQ. ' ' )  GO TO 55
   50       CONTINUE
   55       kstrt (i) = neq (i) - j + 1
   60   CONTINUE
C
C       The end of the value string is the 2nd character before the
C       start of the next keyword, except the last is at l_header.
C
        DO  70 i = 1, numkey - 1
            nvend (i) = kstrt (i+1) - 2
   70   CONTINUE
        nvend (numkey) = l_header
C
C       Now use these indices to load the arrays.  We assume that null
C       strings will not be encountered.
C
        DO  80 i = 1, numkey
            keywd (i) = header (kstrt(i):neq(i)-1)
            value (i) = header (neq(i)+1:nvend(i))
   80   CONTINUE
C
C       Now there are "numkey" keywords with corresponding values ready
C       to be manipulated, printed, etc.  For example, print them:
C
        DO  85 i = 1, numkey
            WRITE (*, *) '"', keywd (i) (1:neq(i)-kstrt(i)), '" = "',
     +                   value (i) (1:nvend(i)-neq(i)), '"'
   85   CONTINUE
        STOP
        END

! calcPeak.f90
!
! coded by T.Sayama on May 2010
!
! calculate peak water level from 2D-RRI output
!
program calcPeak
implicit none

! cndfile
character*256 cndfile
parameter( cndfile = "calcPeak.txt" )

! input files
character*256 infile_dem, infile1, infile
!parameter( infile_dem = "../../Model/infile/indus/adem_indus_lg.txt" )
!parameter( infile1 = "../../Model/out_s/h_" )

! output file
character*256 outfile
!parameter( outfile = "hpeak.txt" )
!parameter( outfile = "../../summary/miyagawa/case13/hpeak.txt" )

! parameter
integer maxt
!parameter( maxt = 96 )

! variable definitions
integer nx, ny
real xllcorner, yllcorner, cellsize, nodata
real, allocatable :: h(:,:), hp(:,:)
character*6 t_char

integer i, j, t
character*256 ctemp

! STEP 0 : open files
open(1, file = cndfile, status = "old")
read(1, '(a)') infile_dem
read(1, '(a)') infile1
read(1, '(i72)') maxt
read(1, '(a)') outfile

open(10, file = infile_dem, status = "old")
open(30, file = outfile)

! STEP 1 : reading files
read(10, *) ctemp, nx
read(10, *) ctemp, ny
read(10, *) ctemp, xllcorner
read(10, *) ctemp, yllcorner
read(10, *) ctemp, cellsize
read(10, *) ctemp, nodata

rewind(10)

read(10, '(a30)') ctemp
write(30, '(a30)') ctemp
read(10, '(a30)') ctemp
write(30, '(a30)') ctemp
read(10, '(a30)') ctemp
write(30, '(a30)') ctemp
read(10, '(a30)') ctemp
write(30, '(a30)') ctemp
read(10, '(a30)') ctemp
write(30, '(a30)') ctemp
read(10, '(a30)') ctemp
write(30, '(a30)') ctemp
close(10)

! STEP 2 : calc peak
allocate( h(ny, nx), hp(ny, nx) )

hp = -1000000.
do t = 1, maxt
 write(*,'(i5, "/", i5)') t, maxt
 call int2char(t, t_char)
 infile = trim(infile1) // trim(t_char) // ".out"
 open(20, file = infile )

 do i = 1, ny
  read(20, *) (h(i, j), j = 1, nx)
 enddo

 do i = 1, ny
  do j = 1, nx
   if( h(i,j) .gt. hp(i,j) ) hp(i,j) = h(i,j)
  enddo
 enddo

 close(20)
enddo

where( hp.lt.0.d0 ) hp = nodata

!where( hp.ge.4.5d0 ) hp = 4.5d0

! STEP 3 : output
do i = 1, ny
 !write(30, '(<nx>f14.5)') (hp(i, j), j = 1, nx)
 write(30, '(<nx>f14.2)') (hp(i, j), j = 1, nx)
enddo
close(30)

end program calcPeak



! numbers to characters
subroutine int2char( num, cwrk )
implicit none

integer num, j
character*6 cwrk
cwrk = ' '
write( cwrk, '( I6 )' ) num
do j = 1, 6
if( cwrk( j:j ) .eq. ' ') cwrk(j:j) = '0'
enddo
end subroutine int2char

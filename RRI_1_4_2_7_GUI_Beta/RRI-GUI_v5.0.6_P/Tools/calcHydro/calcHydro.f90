! calcHydro.f90
!
! obtaining hydrograph information from RRI model output

program calcHydro
implicit none

! condition file
character*256 condfile
parameter( condfile = "calcHydro.txt" )

! location list
character*256 location
!parameter( location = "location.txt" )

! input folder
character*256 infile1, infile
!parameter( infile1 = "./out_q/q_" )

! outfile
character*256 outfile1, outfile
!parameter( outfile1 = "hydro_" )

! target point
character*256 name(100)
integer target_i(100), target_j(100), maxloc
integer col, row

! variable definitions
real(8), allocatable :: disc(:,:), hydro(:,:)
real(8) rtemp
integer maxt
integer t, i, j, ios
character*6 ct

! read condition file
open(1, file = condfile, status = "old")
read(1, '(a)') location
read(1, '(a)') infile1
read(1, '(a)') outfile1

! read location file
open( 5, file = location, status = "old")
i = 1
do
 read(5,*,iostat=ios) name(i), target_i(i), target_j(i)
 if(ios.ne.0) exit
 i = i + 1
enddo
maxloc = i-1
close(5)

! row, col
col = maxval( target_j(:) )
row = maxval( target_i(:) )
allocate( disc(row, col) )
allocate( hydro(maxloc, 100000) )

! read rainfall file
t = 0
do

 t = t + 1
 call int2char( t, ct )
 infile = trim(infile1) // trim(ct) // ".out"

 ! read only the necessary part
 open(20, file = infile, status = "old", err = 1000)

 write(*,'(a)') trim(infile)

 do i = 1, row
  read(20, *) (disc(i, j), j = 1, col)
 enddo
 close(20)

 do i = 1, maxloc
  hydro(i, t) = disc( target_i(i), target_j(i) )
 enddo

enddo
1000 continue
maxt = t - 1

do i = 1, maxloc

 outfile = trim(outfile1) // trim(name(i)) // ".txt"
 open(30, file = outfile)
 do t = 1, maxt
  !write(30, '(i5, f13.5)') t, hydro(i, t)
  write(30, '(i5, e17.8)') t, hydro(i, t)
 enddo

enddo

end program calcHydro



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

! rainThiessen.f90
!
! coded by T.Sayama
! v1.1 Jan 11, 2013
!
! rainThiessen Program
!
implicit none

!!!!!!!!!! NOTICE !!!!!!!!!!!!!!!!!!!!!!!!!!
! Divide !!! When the input data is daily data,
!        !!! set "divide = 24", so that
!        !!! it divides the output data by 24
!        !!! to obtain rainfall intensity (mm/h).

character*256 infile_cond
parameter( infile_cond = "rainThiessen.txt" )

character*256 infile, outfile, outfile_map
real divide, xll, yll, cellsize
integer ncols, nrows

! Variables
integer dt, maxt, num
real, allocatable :: rain_2d(:,:,:), rain_gauge(:,:)
real, allocatable :: lon(:), lat(:)
integer, allocatable :: loc_x(:), loc_y(:), map(:,:)

real dis_min
integer kmin

integer i, j, k, t, ios, itemp, jtemp
real rtemp
character*256 ctemp

! STEP 0: Open files
open(1, file = infile_cond, status = "old" )
read(1, '(a)') infile
read(1, *) divide
read(1, '(a)') outfile
read(1, '(a)') outfile_map
read(1, *) ctemp, ncols
read(1, *) ctemp, nrows
read(1, *) ctemp, xll
read(1, *) ctemp, yll
read(1, *) ctemp, cellsize

! STEP 1: Reading file

open(10, file = infile, status = 'old')
read(10, *) num
read(10, *) ctemp
read(10, *) ctemp

t = 0
do
 read(10, *, iostat = ios) itemp
 if(ios.lt.0) exit
 if(t .eq. 1) dt = itemp
 t = t + 1
enddo
maxt = t - 1

allocate( rain_gauge(0:maxt, num) )
allocate( lon(num), lat(num) )

rewind(10)

read(10, *) num
read(10, *) ctemp, (lat(k), k = 1, num)
read(10, *) ctemp, (lon(k), k = 1, num)
do t = 0, maxt
 read(10, *) itemp, (rain_gauge(t, k), k = 1, num)
 if( itemp .ne. t * dt ) then
  stop "error"
 endif
enddo

close(10)

write(*, *) "dt: ", dt
write(*, *) "maxt: ", maxt


! STEP 2: Calc position

allocate( loc_x(num), loc_y(num) ) 

do k = 1, num
 loc_x(k) = int((lon(k) - xll) / real(cellsize)) + 1
 loc_y(k) = nrows - int((lat(k) - yll) / real(cellsize))
 if( loc_x(k) .le. 0 .or. loc_x(k) .gt. ncols .or. loc_y(k) .le. 0 .or. loc_y(k) .gt. nrows ) then
  loc_x(k) = 0
  loc_y(k) = 0
  write(*,"(i5, ': (', 2f10.3, ')', '  are out of the range.')") k, lon(k), lat(k)
 else
  write(*,"(i5, ': (', 2f10.3, ')', '  (', 2i10, ')')") k, lon(k), lat(k), loc_x(k), loc_y(k)
 endif
enddo

! STEP 3: find closest station

allocate( map(nrows, ncols) )

do i = 1, nrows
 do j = 1, ncols
  dis_min = 1.d15
  kmin = 10000
  do k = 1, num
   if( loc_x(k) .eq. 0 ) cycle
   rtemp = (abs( loc_y(k) - (i - 0.5) )) ** 2. + (abs( loc_x(k) - (j - 0.5) )) ** 2.
   if( rtemp .lt. dis_min ) then
    dis_min = rtemp
    kmin = k
   endif
  enddo
  map(i,j) = kmin
 enddo
enddo

! STEP 4: Interpolation

allocate( rain_2d(0:maxt, nrows, ncols) )

do t = 0, maxt
 do i = 1, nrows
  do j = 1, ncols

   k = map(i, j)

   ! missing data operation (gauge)
   if( rain_gauge(t, k) .lt. 0 ) then
    dis_min = 1.d15
    kmin = 10000
    do k = 1, num
     if( loc_x(k) .eq. 0 ) cycle
     if( rain_gauge(t, k) .lt. 0 ) cycle
     rtemp = (abs( loc_y(k) - (i - 0.5) )) ** 2. + ( abs( loc_x(k) - (j - 0.5) )) ** 2.
     if( rtemp .lt. dis_min ) then
      dis_min = rtemp
      kmin = k
     endif
    enddo
    k = kmin
   endif

   rain_2d(t, i, j) = rain_gauge(t, k)

  enddo
 enddo
enddo

! STEP 5: output

open(20, file = outfile)
do t = 0, maxt
write(20,'(i15, 2i6)') t * dt, ncols, nrows
 do i = 1, nrows
  write(20, '(<ncols>f12.3)') (rain_2d(t, i, j) / divide, j = 1, ncols)  
 enddo
enddo
close(20)

open(30, file = outfile_map)
write(30, '("ncols", i10)') ncols
write(30, '("nrows", i10)') nrows
write(30, '("xllcorner", f20.5)') xll
write(30, '("yllcorner", f20.5)') yll
write(30, '("cellsize", f20.5)') cellsize
write(30, '("NODATA_value", i10)') -9999
do i = 1, nrows
 write(30, '(<ncols>i10)') (map(i, j), j = 1, ncols)
enddo
close(30)

end

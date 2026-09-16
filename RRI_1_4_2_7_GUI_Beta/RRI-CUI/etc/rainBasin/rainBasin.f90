! rainBasin.f90
!
! coded by T.Sayama
! ver 1.0 : Jan. 3, 2011
!
program rainBasin
implicit none

! variable definition
character*256 infile
parameter( infile = "rainBasin.txt" )

character*256 rainfile
character*256 demfile
real(8) xllcorner_rain
real(8) yllcorner_rain
real(8) cellsize_rain_x, cellsize_rain_y
character*256 outfile
character*256 distfile
character*256 cumfile

! topographic variable
integer nx, ny
real(8) xllcorner, yllcorner, cellsize, nodata
real(8), allocatable :: dem(:,:)
integer, allocatable :: domain(:,:)

! rainfall variable
integer, allocatable :: rain_i(:), rain_j(:)
integer maxt
integer, allocatable :: time(:)
integer nx_rain, ny_rain
real(8), allocatable :: qp_org(:,:), qp(:,:,:)
real(8), allocatable :: qp_series(:), qp_dist(:,:), qp_cum(:)

integer num_cell, total_time

! other variable
integer i, j, t, ios, idummy
real(8) rdummy
character*256 ctemp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! STEP 0: Input File Reading
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(10, file = infile, status = 'old')
read(10, '(a)') rainfile
read(10, '(a)') demfile
read(10, *) xllcorner_rain
read(10, *) yllcorner_rain
read(10, *) cellsize_rain_x, cellsize_rain_y
read(10, '(a)') outfile
read(10, '(a)') distfile
read(10, '(a)') cumfile
close(10)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! STEP 1: Dem File Reading
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(20, file = demfile, status = "old")
read(20,*) ctemp, nx
read(20,*) ctemp, ny
read(20,*) ctemp, xllcorner
read(20,*) ctemp, yllcorner
read(20,*) ctemp, cellsize
read(20,*) ctemp, nodata
allocate(dem(ny, nx), domain(ny, nx))
do i = 1, ny
 read(20, *) (dem(i, j), j = 1, nx)
enddo

! domain setting
domain(:,:) = 0
where( dem(:,:) .gt. -10.0 ) domain(:,:) = 1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! STEP 2: Reading Rainfall Data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(30, file = rainfile, status = 'old')

allocate (rain_i(ny), rain_j(nx))
rain_i(:) = 0
rain_j(:) = 0

t = 0
do
 read(30, *, iostat = ios) idummy, nx_rain, ny_rain
 do i = 1, ny_rain
  read(30, *, iostat = ios) (rdummy, j = 1, nx_rain)
 enddo
 if( ios.lt.0 ) exit
 t = t + 1
enddo
maxt = t - 1

allocate( time(0:maxt), qp_org(ny_rain, nx_rain), qp(0:maxt, ny, nx))
rewind(30)
qp = 0.d0
qp_org = 0.d0

do j = 1, nx
 rain_j(j) = int( (xllcorner + (dble(j) - 0.5d0) * cellsize - xllcorner_rain) / cellsize_rain_x ) + 1
enddo
do i = 1, ny
 rain_i(i) = ny_rain - int( (yllcorner + (dble(ny) - dble(i) + 0.5d0) * cellsize - yllcorner_rain) / cellsize_rain_y )
enddo

do t = 0, maxt
 read(30, *) time(t), nx_rain, ny_rain
 do i = 1, ny_rain
  read(30, *) (qp_org(i, j), j = 1, nx_rain)
 enddo

 do i = 1, ny
  if(rain_i(i) .lt. 1 .or. rain_i(i) .gt. ny_rain ) cycle
  do j = 1, nx
   if(rain_j(j) .lt. 1 .or. rain_j(j) .gt. nx_rain ) cycle
   qp(t, i, j) = qp_org(rain_i(i), rain_j(j))
  enddo
 enddo
enddo

close(30)
write(*,*) "done: reading rain file"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! STEP 3: Rain Sum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(qp_series(0:maxt), qp_dist(ny, nx), qp_cum(0:maxt))

qp_series(:) = 0.d0
qp_dist(:,:) = 0.d0
qp_cum(:) = 0.d0

num_cell = 0
do i = 1, ny
 do j = 1, nx
  if(domain(i, j) .eq. 1) then
   num_cell = num_cell + 1
  endif
 enddo
enddo

write(*,*) "grid : ", num_cell

do t = 0, maxt
 do i = 1, ny
  do j = 1, nx
   if(domain(i, j) .eq. 1) then
    qp_series(t) = qp_series(t) + qp(t, i, j) ! qp_series : (mm/h)
    if(t.ne.0) qp_dist(i, j) = qp_dist(i, j) + qp(t, i, j) / 3600.d0 * (time(t) - time(t-1)) ! qp_dist : (mm)
   endif
  enddo
 enddo
 qp_series(t) = qp_series(t) / real(num_cell)
enddo

do t = 1, maxt
 qp_cum(t) = qp_cum(t-1) + qp_series(t) / 3600.d0 * (time(t) - time(t-1)) ! qp_cum : (mm)
enddo
    
where(domain .ne. 1) qp_dist(:,:) = nodata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! STEP 4: Output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! output
open(100, file = outfile)
open(110, file = distfile)
open(120, file = cumfile)

do t = 0, maxt
 write(100,*) time(t), qp_series(t)
enddo
close(100)

do t = 0, maxt
 write(120,*) time(t), qp_cum(t)
 !if(t.ne.0 .and. mod(time(t), 86400) .eq. 0) write(120,*) time(t) / 86400, qp_cum(t)
enddo
close(120)

rewind(20)
read(20, '(a30)') ctemp
write(110, '(a30)') ctemp
read(20, '(a30)') ctemp
write(110, '(a30)') ctemp
read(20, '(a30)') ctemp
write(110, '(a30)') ctemp
read(20, '(a30)') ctemp
write(110, '(a30)') ctemp
read(20, '(a30)') ctemp
write(110, '(a30)') ctemp
read(20, '(a30)') ctemp
write(110, '(a30)') ctemp
close(20)

do i = 1, ny
 write(110,'(10000f14.5)') (qp_dist(i, j), j = 1, nx)
enddo
close(110)

end program rainBasin


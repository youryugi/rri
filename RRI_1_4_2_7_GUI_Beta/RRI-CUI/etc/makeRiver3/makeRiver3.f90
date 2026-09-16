! makeRiver3.f90
!
! coded by T.Sayama on Nov. 6, 2014
!
program makeRiver3
implicit none

character*256 rri_file, width_file, depth_file, height_file
character*256 dem_file, dir_file, acc_file, riv_file

integer thresh, height_limit, nx, ny, utm
real(8) w_c, w_s, d_c, d_s, height_param
real(8) xllcorner, yllcorner, cellsize, nodata
real(8) x1, x2, y1, y2, d1, d2, d3, d4, dx, dy, unit_area

real(8), dimension(:,:), allocatable :: dem, area, width, depth, height
integer, dimension(:,:), allocatable :: dir, acc, riv

integer i, j
character*256 ctemp

open(1, file = "makeRiver3.txt", status = "old")
read(1, '(a30)') rri_file
read(1, '(a30)') width_file
read(1, '(a30)') depth_file
read(1, '(a30)') height_file
read(1, '(a30)') riv_file
close(1)

open(2, file = rri_file, status = "old")
do i = 1, 3 ! L1-L3
 read(2, '(a)') ctemp
enddo
read(2, '(a)') dem_file ! L4
read(2, '(a)') acc_file ! L5
read(2, '(a)') dir_file ! L6
read(2, '(a)') ctemp    ! L7
read(2, *) utm          ! L8
do i = 9, 37
 read(2, '(a)') ctemp   ! L9-L37
enddo
read(2, *) thresh       ! L38
read(2, *) w_c          ! L39
read(2, *) w_s          ! L40
read(2, *) d_c          ! L41
read(2, *) d_s          ! L42
read(2, *) height_param ! L43
read(2, *) height_limit ! L44
close(2)

write(*,'("dem: ", a)') trim(dem_file)
write(*,'("acc: ", a)') trim(acc_file)
write(*,'("dir: ", a)') trim(dir_file)
write(*,*)
write(*,'("utm: ", i10)') utm
write(*,'("thresh: ", i7)') thresh
write(*,'("w_c: ", f10.3)') w_c
write(*,'("w_s: ", f10.3)') w_s
write(*,'("d_c: ", f10.3)') d_c
write(*,'("d_s: ", f10.3)') d_s
write(*,'("height_param: ", f10.3)') height_param
write(*,'("height_limit: ", i10)') height_limit
write(*,*)

open(3, file = dem_file, status = "old")
read(3, *) ctemp, nx
read(3, *) ctemp, ny
read(3, *) ctemp, xllcorner
read(3, *) ctemp, yllcorner
read(3, *) ctemp, cellsize
read(3, *) ctemp, nodata

allocate( dem(ny, nx), dir(ny, nx), acc(ny, nx) )
allocate( width(ny, nx), depth(ny, nx), height(ny, nx), area(ny, nx), riv(ny, nx) )

rewind(3)

open(4, file = dir_file)
open(5, file = acc_file)
open(6, file = width_file)
open(7, file = depth_file)
open(8, file = height_file)
open(9, file = riv_file)

do i = 1, 6
 read(3, '(a30)') ctemp
 read(4, '(a30)') ctemp
 read(5, '(a30)') ctemp
 write(6, '(a30)') ctemp
 write(7, '(a30)') ctemp
 write(8, '(a30)') ctemp
 write(9, '(a30)') ctemp
enddo
do i = 1, ny
 read(3, *) (dem(i, j), j = 1, nx)
 read(4, *) (dir(i, j), j = 1, nx)
 read(5, *) (acc(i, j), j = 1, nx)
enddo


! dx, dy calc
! d1: south side length
x1 = xllcorner
y1 = yllcorner
x2 = xllcorner + nx * cellsize
y2 = yllcorner
if( utm.eq.0 ) call hubeny_sub( x1, y1, x2, y2, d1 )

! d2: north side length
x1 = xllcorner
y1 = yllcorner + ny * cellsize
x2 = xllcorner + nx * cellsize
y2 = yllcorner + ny * cellsize
if( utm.eq.0 ) call hubeny_sub( x1, y1, x2, y2, d2 )

! d3: west side length
x1 = xllcorner
y1 = yllcorner
x2 = xllcorner
y2 = yllcorner + ny * cellsize
if( utm.eq.0 ) call hubeny_sub( x1, y1, x2, y2, d3 )

! d4: east side length
x1 = xllcorner + nx * cellsize
y1 = yllcorner
x2 = xllcorner + nx * cellsize
y2 = yllcorner + ny * cellsize
if( utm.eq.0 ) call hubeny_sub( x1, y1, x2, y2, d4 )

if( utm.eq.1 ) then
 dx = cellsize
 dy = cellsize
else
 dx = (d1 + d2) / 2.d0 / real(nx)
 dy = (d3 + d4) / 2.d0 / real(ny)
endif
write(*,*) "dx [m] : ", dx, "dy [m] : ", dy

! length and area of each cell
unit_area = (dx / 1000.d0) * (dy / 1000.d0) ! [km2]

! calculation
width(:,:) = -9999
depth(:,:) = -9999
height(:,:) = -9999
riv(:,:) = -9999
area(:,:) = -9999
where( dem(:,:) .ge. -100.d0 ) riv(:,:) = -1

do i = 1, ny
 do j = 1, nx
  if( dem(i, j) .le. -100.d0 ) cycle
  area(i, j) = acc(i, j) * unit_area
  if( acc(i, j) .le. thresh ) cycle
  riv(i, j) = 100
  if( dir(i, j) .eq. 0 ) riv(i, j) = 0
  width(i, j) = w_c * area(i, j) ** w_s
  depth(i, j) = d_c * area(i, j) ** d_s
  if( acc(i, j) .gt. height_limit ) height(i, j) = height_param
 enddo
enddo

do i = 1, ny
 write(6, '(10000f9.2)') ( width(i,j), j = 1, nx )
 write(7, '(10000f9.2)') ( depth(i,j), j = 1, nx )
 write(8, '(10000f9.2)') ( height(i,j), j = 1, nx )
 write(9, '(10000i6)') (riv(i, j), j = 1, nx)
enddo

end program makeRiver3



! Hubeny_sub.f90
subroutine hubeny_sub( x1_deg, y1_deg, x2_deg, y2_deg, d )
implicit none

real(8) x1_deg, y1_deg, x2_deg, y2_deg
real(8) x1, y1, x2, y2
real(8) pi, dx, dy, mu, a, b, e, W, N, M, d

pi = 3.1415926535897d0

x1 = x1_deg * pi / 180.d0
y1 = y1_deg * pi / 180.d0
x2 = x2_deg * pi / 180.d0
y2 = y2_deg * pi / 180.d0

dy = y1 - y2
dx = x1 - x2
mu = (y1 + y2) / 2.

a = 6378137.000d0 ! Semi-Major Axis
b = 6356752.314d0 ! Semi-Minor Axis

e = sqrt((a**2.d0 - b**2.d0) / (a**2.d0))

W = sqrt(1. - e**2.d0 * (sin(mu))**2.d0)

N = a / W

M = a * (1. - e ** 2.d0) / W**(3.d0)

d = sqrt((dy * M) ** 2.d0 + (dx * N * cos(mu)) ** 2.d0)

end subroutine

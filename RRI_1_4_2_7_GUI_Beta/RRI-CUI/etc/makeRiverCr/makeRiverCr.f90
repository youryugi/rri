! makeRiver4.f90
!
! coded by T.Sayama on Dec 8, 2021

program makeRiverCr
implicit none

character*256 acc_file, width_file, depth_file, height_file

integer thresh, nx, ny, utm
real(8) s_c, wd_a, wd_b
real(8) xllcorner, yllcorner, cellsize, nodata
real(8) x1, x2, y1, y2, d1, d2, d3, d4, dx, dy, unit_area

real(8), dimension(:,:), allocatable :: area, width, depth, height, s, wd
integer, dimension(:,:), allocatable :: acc

integer i, j
character*256 ctemp

open(1, file = "makeRiverCr.txt", status = "old")
read(1, '(a256)') acc_file
read(1, '(a256)') width_file
read(1, '(a256)') depth_file
read(1, '(a256)') height_file
read(1, *) s_c
read(1, *) wd_a, wd_b
read(1, *) thresh
read(1, *) utm
close(1)

write(*,'("acc: ", a)') trim(acc_file)

open(1, file = acc_file, status = "old")
read(1, *) ctemp, nx
read(1, *) ctemp, ny
read(1, *) ctemp, xllcorner
read(1, *) ctemp, yllcorner
read(1, *) ctemp, cellsize
read(1, *) ctemp, nodata

allocate( acc(ny, nx) )
allocate( width(ny, nx), depth(ny, nx), height(ny, nx), area(ny, nx) )
allocate( s(ny, nx), wd(ny, nx) )

rewind(1)

open(1, file = acc_file)
open(2, file = width_file)
open(3, file = depth_file)
open(4, file = height_file)

do i = 1, 6
 read(1, '(a30)') ctemp
 write(2, '(a30)') ctemp
 write(3, '(a30)') ctemp
 write(4, '(a30)') ctemp
enddo
do i = 1, ny
 read(1, *) (acc(i, j), j = 1, nx)
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
write(*,*) "unit_area [km2]: ", unit_area

! calculation
width(:,:) = -9999
depth(:,:) = -9999
height(:,:) = -9999
area(:,:) = -9999

do i = 1, ny
 do j = 1, nx
  if( acc(i, j) .lt. thresh ) cycle
  area(i, j) = acc(i, j) * unit_area
  s(i, j) = s_c * area(i,j) ** (area(i,j) ** (-0.05))
  wd(i, j) = wd_a * area(i, j) ** wd_b
  width(i, j) = sqrt( wd(i, j) * s(i, j) )
  depth(i, j) = sqrt( s(i, j) / wd(i, j) )
 enddo
enddo

do i = 1, ny
 write(2, '(10000f9.2)') ( width(i,j), j = 1, nx )
 write(3, '(10000f9.2)') ( depth(i,j), j = 1, nx )
 write(4, '(10000f9.2)') ( 0.0, j = 1, nx )
enddo

end program makeRiverCr



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

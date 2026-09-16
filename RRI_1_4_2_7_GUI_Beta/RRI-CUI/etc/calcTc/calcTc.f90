! calcTc.f90
!
! coded by T.Sayama on Mar 30, 2016
!
implicit none

character*256 infile, infile_dem, infile_dir, infile_acc, infile_width
character*256 outfile_tc

parameter( infile = "calcTc.txt" )

integer nx, ny, nxy
real(8) xllcorner, yllcorner, cellsize, nodata

real(8) d1, d2, d3, d4
real(8) x1, x2, y1, y2, dx, dy
integer utm, manning, acc_thresh
real(8) rain, nr, ns, ka
real(8) alph, m

real(8) length, dis, tl, maxlen, riv_ratio
real(8) slo_temp, wid_temp
integer i, j, k, l, ii, jj, iii, jjj, ios

integer, dimension(:,:), allocatable :: dir, acc
real(8), dimension(:,:), allocatable :: dem, width
real(8), dimension(:,:), allocatable :: slo, wid
real(8), dimension(:,:), allocatable :: tc, len
real(8), dimension(:,:), allocatable :: q, h, c, t


character*256 ctemp
character*20 ctemp2

! STEP 0 : Open Files
open(1, file = infile, status = "old")
read(1, '(a)') infile_dem
read(1, '(a)') infile_dir
read(1, '(a)') infile_acc
read(1, '(a)') infile_width
read(1, *) utm
read(1, *) manning
read(1, *) acc_thresh
read(1, *) rain
read(1, *) nr
read(1, *) ns
read(1, *) ka
read(1, '(a)') outfile_tc

open(10, file = infile_dem, status = "old")
open(20, file = infile_dir, status = "old")
open(30, file = infile_acc, status = "old")
open(40, file = infile_width, status = "old")
open(50, file = outfile_tc)

write(*,*) "Done STEP 0"

! STEP 1 : Reading File
read(10, *) ctemp, nx
read(10, *) ctemp, ny
read(10, *) ctemp, xllcorner
read(10, *) ctemp, yllcorner
read(10, *) ctemp, cellsize
read(10, *) ctemp, nodata

allocate( dem(ny, nx), dir(ny, nx), acc(ny, nx), width(ny, nx) )
allocate( slo(ny, nx), wid(ny, nx) )
allocate( tc(ny, nx), len(ny, nx) )
allocate( q(ny, nx), h(ny, nx), c(ny, nx), t(ny, nx))

rewind(10)

read(10, '(a30)') ctemp
write(50, '(a30)') ctemp

read(10, '(a30)') ctemp
write(50, '(a30)') ctemp

read(10, '(a30)') ctemp
write(50, '(a30)') ctemp

read(10, '(a30)') ctemp
write(50, '(a30)') ctemp

read(10, '(a30)') ctemp
write(50, '(a30)') ctemp

read(10, '(a30)') ctemp
write(50, '(a30)') ctemp

do i = 1, ny
 read(10, *) (dem(i, j), j = 1, nx)
enddo

read(20, *) ctemp, nx
read(20, *) ctemp, ny
read(20, *) ctemp, xllcorner
read(20, *) ctemp, yllcorner
read(20, *) ctemp, cellsize
read(20, *) ctemp, nodata

do i = 1, ny
 read(20, *) (dir(i, j), j = 1, nx)
enddo

read(30, *) ctemp, nx
read(30, *) ctemp, ny
read(30, *) ctemp, xllcorner
read(30, *) ctemp, yllcorner
read(30, *) ctemp, cellsize
read(30, *) ctemp, nodata

do i = 1, ny
 read(30, *) (acc(i, j), j = 1, nx)
enddo

read(40, *) ctemp, nx
read(40, *) ctemp, ny
read(40, *) ctemp, xllcorner
read(40, *) ctemp, yllcorner
read(40, *) ctemp, cellsize
read(40, *) ctemp, nodata

do i = 1, ny
 read(40, *) (width(i, j), j = 1, nx)
enddo

write(*,*) "Done STEP 1"

! STEP 2 : dx, dy

! d1 : South
x1 = xllcorner
y1 = yllcorner
x2 = xllcorner + nx * cellsize
y2 = yllcorner
if( utm.eq.0 ) call hubeny_sub( x1, y1, x2, y2, d1 )

! d2 : North
x1 = xllcorner
y1 = yllcorner + ny * cellsize
x2 = xllcorner + nx * cellsize
y2 = yllcorner + ny * cellsize
if( utm.eq.0 ) call hubeny_sub( x1, y1, x2, y2, d2 )

! d3 : West
x1 = xllcorner
y1 = yllcorner
x2 = xllcorner
y2 = yllcorner + ny * cellsize
if( utm.eq.0 ) call hubeny_sub( x1, y1, x2, y2, d3 )

! d4 : East
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

length = (dx + dy) / 2.d0

write(*,*) "Done STEP 2"

! STEP 3

q(:,:) = nodata
h(:,:) = nodata
c(:,:) = nodata
t(:,:) = nodata
slo(:,:) = nodata
wid(:,:) = nodata

do i = 1, ny
 do j = 1, nx

  if( dem(i, j) .lt. -100.d0 .or. dir(i, j) .eq. 0 ) cycle
  call down(dir, dem, nx, ny, i, j, length, ii, jj, dis, slo_temp, wid_temp)

  slo(i, j) = slo_temp
  wid(i, j) = wid_temp

  if( manning .eq. 1 ) then
   alph = sqrt( slo(i, j) ) / ns
   m = 5.d0 / 3.d0
  else
   alph = ka * slo(i, j)
   m = 1.d0
  endif

  if( width(i, j) .gt. 0.d0 ) then
   alph = sqrt( slo(i, j) ) / nr
   m = 5.d0 / 3.d0
   wid(i, j) = width(i, j)
  endif

  q(i, j) = acc(i, j) * length * length * rain / 1000.d0 / 3600.d0
  h(i, j) = (q(i, j) / wid(i, j) / alph) ** (1.d0 / m)
  c(i, j) = alph * m * h(i, j) ** (m - 1.d0)
  t(i, j) = 0.d0
  if( c(i, j) .gt. 0.d0 ) t(i, j) = length / c(i, j) / 3600.d0 ! sec -> hr
 enddo
enddo

! STEP 4
tc(:,:) = nodata
do iii = 1, ny
 write(*, *) iii, "/", ny
 do jjj = 1, nx
  if( dem(iii, jjj) .lt. -100.d0 .or. dir(iii, jjj) .eq. 0 ) cycle
  i = iii
  j = jjj
  tl = 0.d0
  do
   call down(dir, dem, nx, ny, i, j, length, ii, jj, dis, slo_temp, wid_temp)
   if( dir(ii, jj) .eq. 0 .or. dem(ii, jj) .lt. -100.d0 .or. ii.eq.0 .or. jj.eq.0 .or. ii.gt.ny .or. jj.gt.nx ) then
    tc(iii, jjj) = tl
    exit
   else
    tl = tl + t(i, j)
    i = ii
    j = jj
   endif
  enddo
 enddo
enddo

! STEP 5
do i = 1, ny
 write(50, '(10000f10.3)') ( tc(i, j), j = 1, nx )
 !write(50, '(10000f10.3)') ( width(i, j), j = 1, nx )
enddo

close(50)

end


! Search Downstream
subroutine down(dir, dem, nx, ny, i, j, length, ii, jj, dis, slo, wid)

integer nx, ny, i, j, ii, jj
real(8) dem(ny, nx), length, dis, slo, wid
integer dir(ny, nx)

! right
if( dir(i,j).eq.1 ) then
 ii = i
 jj = j + 1
 dis = length
 wid = length
! right down
elseif( dir(i,j).eq.2 ) then
 ii = i + 1
 jj = j + 1
 dis = length * sqrt(2.d0)
 wid = length / sqrt(2.d0)
! down
elseif( dir(i,j).eq.4 ) then
 ii = i + 1
 jj = j
 dis = length
 wid = length
! left down
elseif( dir(i,j).eq.8 ) then
 ii = i + 1
 jj = j - 1
 dis = length * sqrt(2.d0)
 wid = length / sqrt(2.d0)
! left
elseif( dir(i,j).eq.16 ) then
 ii = i
 jj = j - 1
 dis = length
 wid = length
! left up
elseif( dir(i,j).eq.32 ) then
 ii = i - 1
 jj = j - 1
 dis = length * sqrt(2.d0)
 wid = length / sqrt(2.d0)
! up
elseif( dir(i,j).eq.64 ) then
 ii = i - 1
 jj = j
 dis = length
 wid = length
! right up
elseif( dir(i,j).eq.128 ) then
 ii = i - 1
 jj = j + 1
 dis = length * sqrt(2.d0)
 wid = length / sqrt(2.d0)
! zero
else
 ii = i
 jj = j
 dis = length
 wid = length
endif

dif = abs( dem(i, j) - dem(ii, jj) )
slo = dif / sqrt((dif ** 2.d0 + dis ** 2.d0))
if( slo .lt. 0.01d0 ) slo = 0.01d0

end subroutine


! Hubeny_sub.f90
subroutine hubeny_sub( x1_deg, y1_deg, x2_deg, y2_deg, d )
implicit none

real(8) x1_deg, y1_deg, x2_deg, y2_deg
real(8) x1, y1, x2, y2
real(8) pi, dx, dy, mu, a, b, e, W, N, M, d

!read(*,*) x1_deg, y1_deg, x2_deg, y2_deg ! in degree

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

!write(*,*) d ! in m

end subroutine

! section.f90
!
! coded by T.Sayama on March 18, 2012
! updated on Nov. 8, 2017

implicit none

character*256 infile, outfile
real, allocatable :: x(:), y(:)
real, allocatable :: sec_depth(:), sec_area(:), sec_peri(:), sec_b(:)
real datum, dy, len, h1, h2, b, ratio, rdummy, height, depth, ns_river
integer div, startx, endx, maxi, i, j, ios

! STEP 0: File Prep
open(1, file = "section.txt", status = "old")
read(1, '(a)') infile
read(1, '(a)') outfile
read(1, *) ns_river
read(1, *) div
read(1, *) datum
read(1, *) startx
read(1, *) endx

! STEP 1: File Reading
open(2, file = infile, status = "old")

do i = 1, startx - 1
 read(2, *) rdummy, rdummy
enddo

maxi = endx - startx
allocate( x(0:maxi), y(0:maxi) )

do i = 0, maxi
 read(2, *) x(i), y(i)
enddo
close(2)

! STEP 2: Calculation
allocate( sec_depth(0:div), sec_peri(div), sec_area(div), sec_b(div) )
sec_depth(:) = 0.
sec_peri(:) = 0.
sec_area(:) = 0.
sec_b(:) = 0.

depth = datum - minval(y(:))
height = min( y(0), y(maxi) )
height = height - datum

y(:) = y(:) - minval(y(:))
if( height .le. 0 ) height = 0.

dy = maxval(y(:)) / real(div)
do j = 0, div
 sec_depth(j) = j * dy
enddo

do j = 1, div

 do i = 1, maxi

  h1 = sec_depth(j) - y(i-1)
  h2 = sec_depth(j) - y(i)
  b = x(i) - x(i-1)
  len = sqrt( (abs(h1 - h2)) ** 2. + b ** 2. )

  if( y(i-1) .le. sec_depth(j) .and. y(i) .le. sec_depth(j) ) then

   sec_peri(j) = sec_peri(j) + len
   !sec_area(j) = sec_area(j) + (h1 + h2) * b / 2.
   sec_b(j) = sec_b(j) + b

  elseif( y(i-1) .le. sec_depth(j)  ) then

   ratio = h1 / (h1 - h2)
   sec_peri(j) = sec_peri(j) + len * ratio
   !sec_area(j) = sec_area(j) + h1 * (b * ratio) / 2.
   sec_b(j) = sec_b(j) + b * ratio

  elseif( y(i) .le. sec_depth(j) ) then

   ratio = h2 / (h2 - h1)
   sec_peri(j) = sec_peri(j) + len * ratio
   !sec_area(j) = sec_area(j) + h2 * (b * ratio) / 2.
   sec_b(j) = sec_b(j) + b * ratio

  endif

 enddo

enddo

! STEP 3: Output
open(10, file = outfile)
write(10, '(i13, 1000f13.5)') div, depth, height
do j = 1, div
 write(10, '(1000f13.5)') sec_depth(j), sec_peri(j), sec_b(j), ns_river
enddo
close(10)

end



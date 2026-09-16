! evalPeak.f90
!
! coded by T.Sayama on May 2011
!
! evaluation of flood inundation by RRI
!
program evalPeak
implicit none

! infile
character*256 infile
parameter( infile = "evalPeak.txt" )

! infile_sim, infile_obs
character*256 infile_sim, infile_obs, infile_range
!parameter( infile_sim = "infile_sim.txt" )
!parameter( infile_obs = "infile_obs.txt" )
!parameter( infile_range = "infile_range.txt" )

! variable definitions
integer nx, ny
real xllcorner, yllcorner, cellsize, nodata
real, allocatable :: sim(:,:), obs(:,:)
character*6 t_char
integer imin, imax, jmin, jmax, ishift, jshift
integer sim_and_obs, sim_or_obs

integer i, j, count, ii, jj
character*256 ctemp

! STEP 0 : open files
open(1, file = infile, status = "old")
read(1, '(a)') infile_sim
read(1, '(a)') infile_obs
read(1, '(a)') infile_range

open(10, file = infile_sim, status = "old")
open(20, file = infile_obs, status = "old")
open(30, file = infile_range, status = "old")

! STEP 1 : reading files

read(10, *) ctemp, nx
read(10, *) ctemp, ny
read(10, *) ctemp, xllcorner
read(10, *) ctemp, yllcorner
read(10, *) ctemp, cellsize
read(10, *) ctemp, nodata

read(20, *) ctemp, nx
read(20, *) ctemp, ny
read(20, *) ctemp, xllcorner
read(20, *) ctemp, yllcorner
read(20, *) ctemp, cellsize
read(20, *) ctemp, nodata

allocate( sim(ny, nx), obs(ny, nx) )

do i = 1, ny
 read(10, *) (sim(i, j), j = 1, nx)
enddo

do i = 1, ny
 read(20, *) (obs(i, j), j = 1, nx)
enddo

read(30, *) imin, imax
read(30, *) jmin, jmax
read(30, *) ishift, jshift

! STEP 2 : Evaluate only within the range

if( imin .lt. 1 ) imin = 1
if( imax .gt. ny ) imax = ny
if( jmin .lt. 1 ) jmin = 1
if( jmax .gt. nx ) jmax = nx

write(*,*) "i : ", imin, imax
write(*,*) "j : ", jmin, jmax

sim_or_obs = 0
sim_and_obs = 0

do i = imin, imax
 do j = jmin, jmax
  if( sim(i, j) .lt. 0. ) cycle
  if( obs(i, j) .lt. 0. ) cycle
  if( sim(i, j) .ge. 1.0 .or. obs(i - ishift, j - jshift) .ge. 0.0 ) sim_or_obs = sim_or_obs + 1
  if( sim(i, j) .ge. 1.0 .and. obs(i - ishift, j - jshift) .ge. 0.0 ) sim_and_obs = sim_and_obs + 1
 enddo
enddo

write(*,*) sim_and_obs, sim_or_obs, real(sim_and_obs) / real(sim_or_obs)

close(10)
close(20)
close(30)

pause

end program evalPeak

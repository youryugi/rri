! Tecplot_animation .f90
! obtaining a data file for Tecplot animation for the results of RRI model 

program Tec_anime
implicit none

! condition file
character*256 infile_dem, condfile, outfile
character*256 infile_rain
character*256 outfile_hs, outfile_hr, outfile_hg, outfile_qr, outfile_qu, outfile_qv
character*256 outfile_gu, outfile_gv, outfile_gampt_ff
character*256 outfile_hs2, outfile_hr2, outfile_hg2, outfile_qr2, outfile_qu2, outfile_qv2
character*256 outfile_gu2, outfile_gv2, outfile_gampt_ff2
character*256 a,b,c,d,e,f,g,h,p,u,v,w,x,y,yy
character*256 g1,g2,g3,g4,g5,g6,aa
character*6 ctemp,ct
parameter( condfile = "calcTecplot.txt" )

! variable definitions
real(8), allocatable :: dem (:,:),slope(:,:), river(:,:), qr(:,:), rain(:,:,:), rain_(:,:)
real(8), allocatable :: hg(:,:),qu(:,:), qv(:,:), gu(:,:), gv(:,:), gampt(:,:), all(:,:,:), qp(:,:,:)
real(8) rdummy
integer, allocatable :: rain_i(:), rain_j(:)
integer tt_max_rain
integer, allocatable :: t_rain(:)
real(8) xllcorner, yllcorner, cellsize, time, itemp, jtemp, tt
real(8) xllcorner_rain, yllcorner_rain, cellsize_rain_x, cellsize_rain_y
integer  t,i,j,lasth, maxt, nx,ny, recno, ii, fn,nx_rain,ny_rain, hen, k, trim_i, it, BB
integer  outswitch_hs, outswitch_hr, outswitch_hg, outswitch_qr, outswitch_qu, outswitch_qv, outswitch_gu, outswitch_gv, outswitch_gampt_ff
integer year, month , day, hour
real(8)sirial
fn=50

!----read maketec_input file ---------------------------------- 
open(1, file = condfile, status = "old")
read(1, *) year, month, day , hour
read(1, *) lasth
read(1, *) maxt
read(1, *)
read(1, '(a)') infile_rain
read(1, '(a)') infile_dem
read(1, *)
read(1, *) xllcorner_rain
read(1, *) yllcorner_rain
read(1, *) cellsize_rain_x, cellsize_rain_y
read(1, *)
read(1, *) outswitch_hs, outswitch_hr, outswitch_hg, outswitch_qr, outswitch_qu, outswitch_qv, outswitch_gu, outswitch_gv, outswitch_gampt_ff
write(*, '(20I3)') outswitch_hs, outswitch_hr, outswitch_hg, outswitch_qr, outswitch_qu, outswitch_qv, outswitch_gu, outswitch_gv, outswitch_gampt_ff
read(1, '(a)') outfile_hs
read(1, '(a)') outfile_hr
read(1, '(a)') outfile_hg
read(1, '(a)') outfile_qr
read(1, '(a)') outfile_qu
read(1, '(a)') outfile_qv
read(1, '(a)') outfile_gu
read(1, '(a)') outfile_gv
read(1, '(a)') outfile_gampt_ff
read(1, *)
read(1, '(a)') outfile
read(1, *)
read(1, '(a)') a
read(1, '(a)') b
read(1, '(a)') c
read(1, '(a)') d
read(1, '(a)') y   !rainfall
read(1, '(a)') e
read(1, '(a)') f
read(1, '(a)') g
read(1, '(a)') g1
read(1, '(a)') g2
read(1, '(a)') g3
read(1, '(a)') g4
read(1, '(a)') g5
read(1, '(a)') g6
read(1, *)
read(1, '(a)') h
read(1, '(a)') p
read(1, '(a)') u
read(1, '(a)') v
read(1, '(a)') w
close(1)
call calcsirial(year, month, day, hour, lasth, maxt,sirial)
!open(unit=fn,file='Tec_anime.bin',status='unknown', access='direct', recl=24)

!------------------------------------------------------
!write output file-------------------------------------
outfile = trim(outfile)
open(unit=5,file=outfile)
write(5,'(a)') trim(a)
write(5,'(a)') trim(b)
write(5,'(a)') trim(c)
write(5,'(a)') trim(d)
write(5,'(a)') trim(y)
if(outswitch_hs .ne. 0) write(5,'(a)') trim(e)
if(outswitch_hr .ne. 0) write(5,'(a)') trim(f)
if(outswitch_hg .ne. 0) write(5,'(a)') trim(g)
if(outswitch_qr .ne. 0) write(5,'(a)') trim(g1)
if(outswitch_qu .ne. 0) write(5,'(a)') trim(g2)
if(outswitch_qv .ne. 0) write(5,'(a)') trim(g3)
if(outswitch_gu .ne. 0) write(5,'(a)') trim(g4)
if(outswitch_gv .ne. 0) write(5,'(a)') trim(g5)
if(outswitch_gampt_ff .ne. 0)write(5,'(a)') trim(g6)
recno=0

call write_char256(fn, recno, a)!
call write_char256(fn, recno, b)
call write_char256(fn, recno, c)
call write_char256(fn, recno, d)
call write_char256(fn, recno, e)
call write_char256(fn, recno, f)
call write_char256(fn, recno, g)
call write_char256(fn, recno, g1)
call write_char256(fn, recno, g2)
call write_char256(fn, recno, g3)
call write_char256(fn, recno, g4)
call write_char256(fn, recno, g5)
call write_char256(fn, recno, y)


!---reading dem file -------------------------------------
infile_dem = trim(infile_dem)
open(10, file = infile_dem, status = "old")
read(10, *) ctemp, nx
read(10, *) ctemp, ny
read(10, *) ctemp, xllcorner
read(10, *) ctemp, yllcorner
read(10, *) ctemp, cellsize
read(10, *) ctemp

    allocate( dem(ny, nx) ) 
    do i = 1, ny
        read(10, *) (dem(i, j), j = 1, nx)
      do j = 1,nx
        if(dem(i,j).le.-9000) dem(i,j) = -10
      end do
    enddo
!--------------------------------------------------------
  recno = 4 + outswitch_hs+outswitch_hr+outswitch_hg+outswitch_qr+outswitch_qu+outswitch_qv+outswitch_gu+outswitch_gv+outswitch_gampt_ff
  hen = recno-4

  allocate( slope(ny,nx),river(ny,nx), qr(ny,nx), rain(ny,nx,100), rain_(ny,nx), rain_i(ny), rain_j(nx))
  allocate( hg(ny, nx),qu(ny, nx), qv(ny, nx), gu(ny, nx), gv(ny, nx), gampt(ny, nx))
  allocate( all(ny, nx, hen))

! reading rainfall data --------------------------------------------
infile_rain = trim(infile_rain)
open(41, file = infile_rain, status = "old")   ! Rainfall

tt = 0

do
 read(41, *, end = 99) t, nx_rain, ny_rain
 do i = 1, ny_rain
  read(41, *, end = 99) (rdummy, j = 1, nx_rain)
 enddo
 tt = tt + 1
enddo
99 continue
tt_max_rain = tt - 1
write(*,*)tt_max_rain, t

allocate( t_rain(0:tt_max_rain), qp(0:tt_max_rain, ny_rain, nx_rain))
rewind(41)

qp = 0.d0
do tt = 0, tt_max_rain
 read(41, *) t_rain(tt), nx_rain, ny_rain
 do i = 1, ny_rain
  read(41, *) (qp(tt, i, j), j = 1, nx_rain)
 enddo
enddo

  do j = 1, nx
   rain_j(j) = int( (xllcorner + (dble(j) - 0.5d0) * cellsize - xllcorner_rain) / cellsize_rain_x ) + 1
  enddo
  do i = 1, ny
  rain_i(i) = ny_rain - int( (yllcorner + (dble(ny) - dble(i) + 0.5d0) * cellsize - yllcorner_rain) / cellsize_rain_y )
  enddo
  close(41)

write(*,*)'read OK'
!--------------------------------------------------------
! read output data---------------------------------------

  t = 0
 do  

 t = t + 1
 time = real(t) * real(lasth) / real(maxt) * 3600
! write(*,'(a10,f10.3,a5)') 'time =', time/3600, '(h)' 
  if(t.gt.maxt) exit 
    call int2char( t, ct )
   
     write (*, '(i5,a3,i5)')t, ' / ', maxt
        outfile_hs2 = trim(outfile_hs) // trim(ct) // ".out"
        outfile_hr2 = trim(outfile_hr) // trim(ct) // ".out"
        outfile_hg2 = trim(outfile_hg) // trim(ct) // ".out"
        outfile_qr2 = trim(outfile_qr) // trim(ct) // ".out"
        outfile_qu2 = trim(outfile_qu) // trim(ct) // ".out"
        outfile_qv2 = trim(outfile_qv) // trim(ct) // ".out"
        outfile_gu2 = trim(outfile_gu) // trim(ct) // ".out"
        outfile_gv2 = trim(outfile_gv) // trim(ct) // ".out"
        outfile_gampt_ff2 = trim(outfile_gampt_ff) // trim(ct) // ".out"

if(outswitch_hs .ne. 0) open(20, file = outfile_hs2, status = "old") 
if(outswitch_hr .ne. 0) open(30, file = outfile_hr2, status = "old") 
if(outswitch_hg .ne. 0) open(35, file = outfile_hg2, status = "old") 
if(outswitch_qr .ne. 0) open(40, file = outfile_qr2, status = "old") 
if(outswitch_qu .ne. 0) open(70, file = outfile_qu2, status = "old") 
if(outswitch_qv .ne. 0) open(71, file = outfile_qv2, status = "old") 
if(outswitch_gu .ne. 0) open(72, file = outfile_gu2, status = "old") 
if(outswitch_gv .ne. 0) open(73, file = outfile_gv2, status = "old") 
if(outswitch_gampt_ff .ne. 0) open(74, file = outfile_gampt_ff2, status = "old") 
 
 it = 1
if(outswitch_hs .ne. 0) then
      read(20, *) ((slope(i, j), j = 1, nx),i=1,ny)
      all(:,:,it) = slope(:, :)
      it = it + 1
end if
if(outswitch_hr .ne. 0) then
      read(30, *) ((river(i, j), j = 1, nx),i=1,ny)
      all(:,:,it) = river(:, :)
      it = it + 1
end if
if(outswitch_hg .ne. 0) then
      read(35, *) ((hg(i, j), j = 1, nx),i=1,ny)
      all(:,:,it) = hg(:, :)
      it = it + 1
end if
if(outswitch_qr .ne. 0) then
      read(40, *) ((qr(i, j), j = 1, nx),i=1,ny)
      all(:,:,it) = qr(:, :)
      it = it + 1
end if
if(outswitch_qu .ne. 0) then
      read(70, *) ((qu(i, j), j = 1, nx),i=1,ny)
      all(:,:,it) = qu(:, :)
      it = it + 1
end if
if(outswitch_qv .ne. 0) then
      read(71, *) ((qv(i, j), j = 1, nx),i=1,ny)
      all(:,:,it) = qv(:, :)
      it = it + 1
end if
if(outswitch_gu .ne. 0)then
      read(72, *) ((gu(i, j), j = 1, nx),i=1,ny)
      all(:,:,it) = gu(:, :)
      it = it + 1
end if
if(outswitch_gv .ne. 0) then
      read(73, *) ((gv(i, j), j = 1, nx),i=1,ny)
      all(:,:,it) = gv(:, :)
      it = it + 1
end if
if(outswitch_gampt_ff .ne. 0) then
      read(74, *) ((gampt(i, j), j = 1, nx),i=1,ny)
      all(:,:,it) = gampt(:, :)
      it = it + 1
end if


! remake the time scale of rainfall---------------------------------
  itemp = -1
  do jtemp = 1, tt_max_rain
   if( t_rain(jtemp-1) .lt. (time ) .and. (time ) .le. t_rain(jtemp) ) itemp = jtemp
  enddo
  do i = 1, ny
   do j = 1, nx
    rain_(i, j) = qp(itemp, rain_i(i), rain_j(j))
   enddo
  enddo

!   write(55,*) 'itemp =' ,itemp , t_rain(itemp)
!  do i = 1, ny
!   write(55,'(<nx>f6.2)') (rain_(i, j) , j= 1, nx)
!  enddo
!------------------------------------------------------
       write (5,*) trim(h)
       write (5,*) trim(p),sirial + real(lasth)/real(maxt)*(t-1)/24
!       write (5,*) trim(u)
       write (5,'(1x,a3,i4,a5,i4,a25)') 'I =', nx, ', J =', ny,', K=1, ZONETYPE=Ordered'
       write (5,*) trim(v)
       AA = repeat('SINGLE ',hen+4)
       BB = len_trim(AA)
       write (5,'(A5,A<BB>,A3)') 'DT=(', trim(AA), ' )'
!c       write (5,'(a)') trim(w)
       call write_char256(fn, recno, h)
       call write_char256(fn, recno, p)
       recno=recno+1
!       write (fn,rec=recno) t
       call write_char256(fn, recno, u)
       call write_char256(fn, recno, v)
       call write_char256(fn, recno, w)
      
     i = 1
      
      
  do i = 1, ny
     do j = 1, nx 
    
    write(5, '(i4,x,i4,x,f8.2,x,f8.4,x,f8.4,x,10f13.5)') i,j,dem(i,j),rain_(i, j), (all(i,j,k), k = 1,hen)
   
    recno=recno+1 
    ! write(fn, rec=recno) i,j,real(dem(i,j)),real(slope(i,j)), real(river(i,j)), real(qr(i,j) )
    
    enddo
  enddo

! write(*, *)'1'
102 continue
 end do
  close(fn)
  pause
end program Tec_anime


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

subroutine write_char256(fn, recno, a)
    character*256 :: a
    integer, intent(in) :: fn
    integer, intent(inout) :: recno
    integer ii
    do ii = 1, 10
    recno=recno+1
!    write(fn,rec=recno) a(1+24*(ii-1):24*ii)
    enddo
    recno=recno+1
!    write(fn,rec=recno) a(241:256)
end subroutine 

subroutine calcsirial(year, month, day, hour, lasth, maxt,sirial)
    implicit none
    integer year, month , day, hour
    integer lasth, maxt
    integer ii,jj 
    real(8) sirial
    
    sirial = 1
    do ii = 1900, year-1
    if(mod(ii,4).eq.0) then
      sirial = sirial + 366
    else
      sirial = sirial + 365
    end if
    end do

if (month.eq.1)go to 111    
    do ii = 1, month-1
     if(( ii .eq. 4 ).or.( ii .eq. 6 ).or.( ii .eq. 9 ).or.(ii .eq. 11 )) sirial = sirial + 30
     if(( ii .eq. 1 ).or.( ii .eq. 3 ).or.( ii .eq. 5 ).or.(ii .eq. 7 ).or.(ii .eq. 8).or.(ii .eq. 10)) sirial = sirial + 31      
     if(( ii .eq. 2 ).and.(mod(year,4).eq.0)) sirial = sirial + 29      
     if(( ii .eq. 2 ).and.(mod(year,4).ne.0)) sirial = sirial + 28
    end do
111 continue

if (day.eq.1)go to 112    
    sirial = sirial + real(day-1)
112 continue

    sirial = sirial + real(hour) /24.


end subroutine 





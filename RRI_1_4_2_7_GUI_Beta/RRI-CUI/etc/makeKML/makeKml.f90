! Program to create kml for animation
program makekml
!implicit none

character*4 csy, cey
character*17 startt, endt
character*2 csm, csd, csh, cem, ced, ceh, csmin, cemin
integer sy, sm, sd, sh, em, ed, eh, ey, smin, emin
character*6 ct
character*100 demfile, outfile, ctemp
integer t, tmax, dtm, dth
real(8) dt, dtm0
real(8) xllcorner, yll_corner, xllcorner2, yll_corner2, cellsize
!Read input file

open(90, file = "Kml_input.txt")
read(90, *) sy, sm, sd, sh, smin
read(90, *) tmax
read(90, *) dt
read(90,'(a)')demfile
write(*,'("demfile : ", a)') trim(adjustl(demfile))
read(90, *)
read(90,'(a)')outfile
write(*,'("outfile : ", a)') trim(adjustl(outfile))
open(95, file = demfile)
read(95,*) ctemp, ncols
read(95,*) ctemp, nrows
read(95,*) ctemp, xllcorner
read(95,*) ctemp, yllcorner
read(95,*) ctemp, cellsize
  xllcorner2 = xllcorner + real(ncols) * cellsize
  yllcorner2 = yllcorner + real(nrows) * cellsize
  dth = int(aint(dt))
  dtm0 = mod(dt, 1.0)
  dtm = int(aint(dtm0 * 60))

open(100, file = outfile)

write(100,'("<Folder>")')

do t = 1, tmax

 call tplus(sy, sm, sd, sh, smin, ey, em, ed, eh, emin, dth, dtm)
 
 call int4char4(sy, csy)
 call int2char2(sm, csm)
 call int2char2(sd, csd)
 call int2char2(sh, csh)
 call int2char2(smin, csmin)

 call int4char4(ey, cey)
 call int2char2(em, cem)
 call int2char2(ed, ced)
 call int2char2(eh, ceh)
 call int2char2(emin, cemin)

 call int2char(t, ct)
 
 startt = csy // "-" // csm // "-" // csd // "T" // csh // ":" // csmin // "Z"
 endt   = cey // "-" // cem // "-" // ced // "T" // ceh // ":" // cemin // "Z"

 write(*,*) t, startt, endt
 write(100,'(" <GroundOverlay>")')
 write(100,'("  <TimeSpan>")')
 write(100,'("   <begin>", a17, "</begin>")') startt
 write(100,'("   <end>", a17, "</end>")') endt
 write(100,'("  </TimeSpan>")')
 write(100,'("  <Icon>")')
 write(100,'("    <href>hs_kml/hs_", a6, ".gif</href>")') ct
 write(100,'("  </Icon>")')
 write(100,'("  <LatLonBox>")')
 write(100,'("   <north>", f10.5 , "</north>")') yllcorner2
 write(100,'("   <south>", f10.5 , "</south>")') yllcorner
 write(100,'("   <east>", f10.5 , "</east>")') xllcorner2
 write(100,'("   <west>", f10.5 , "</west>")') xllcorner
 write(100,'("  </LatLonBox>")')
 write(100,'(" </GroundOverlay>")')
 sy = ey
 sm = em
 sd = ed
 sh = eh
 smin = emin
enddo
write(100,'("</Folder>")')
end


! Plus Timestep
subroutine tplus(sy, sm, sd, sh, smin, ey, em, ed, eh, emin, dth, dtm)
integer sy, sm, sd, sh, smin, ey, em, ed, eh, emin, dth, dtm

ey = sy
em = sm
ed = sd
eh = sh + dth

if(dtm.ge.1)then
  emin = smin + dtm
  if(emin .ge. 60) then
    emin = emin - 60
    eh = eh +1
  end if
end if

if( eh .ge. 24 ) then
99 continue 
 eh = eh - 24
 ed = ed + 1
 if(eh.ge.24) go to 99
endif


if( ed .ge. 29 ) then
 if(( em .eq. 2 ).and.(mod(ey,4).ne.0)) then
    em = em + 1
    ed = ed - 29 + 1
 end if
end if
if( ed .ge. 30 )then
 if(( em .eq. 2 ).and.(mod(ey,4).eq.0)) then
    em = em + 1
    ed = ed - 30 + 1
 end if
end if

if( ed .ge. 31 ) then
 if(( em .eq. 4 ).or.( em .eq. 6 ).or.( em .eq. 9 ).or.( em .eq. 11 ))then 
   em = em + 1
   ed = ed - 31 + 1
 endif
endif

if( ed .ge. 32 ) then
 if(( em .eq. 1 ).or.( em .eq. 3 ).or.( em .eq. 5 ).or.( em .eq. 7 ).or.( em .eq. 8 ).or.( em .eq. 10 ))then 
   em = em + 1
   ed = ed - 32 + 1
 elseif( em .eq. 12 ) then
   em = 1
   ed = ed - 32 + 1
   ey = ey + 1
 endif
end if
end subroutine tplus


! numbers to characters
subroutine int4char4( num, cwrk )
implicit none

integer num, j
character*4 cwrk
cwrk = ' '
write( cwrk, '( I4 )' ) num
do j = 1, 4
if( cwrk( j:j ) .eq. ' ') cwrk(j:j) = '0'
enddo
end subroutine int4char4

! numbers to characters
subroutine int2char2( num, cwrk )
implicit none

integer num, j
character*2 cwrk
cwrk = ' '
write( cwrk, '( I2 )' ) num
do j = 1, 2
if( cwrk( j:j ) .eq. ' ') cwrk(j:j) = '0'
enddo
end subroutine int2char2

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


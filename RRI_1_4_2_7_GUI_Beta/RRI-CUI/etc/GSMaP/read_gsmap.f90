implicit none
! Read GSMaP binary data and convert it to ASCII format,
! which can be read by RRI Model
!
! (Note) i and j are switched as compared to the original
! read_GSMaP_MVK_0.1deg program
!
! list.txt file must be prepared in advance by using make.bat
!

character*256 listfile
parameter( listfile = "list.txt" )

! input folder
character*256 infile1, infile
parameter( infile1 = "./infile/" )

character*256 cutfile1, cutfile
parameter( cutfile1 = "./cutfile/" )

character*256 rainfile
parameter( rainfile = "rain.data" )

character*256 sumfile
parameter( sumfile = "sumRain.data" )

character*256 areafile
parameter( areafile = "out_by_calc_area_gsmap.txt" )

! Target Area
integer jleft, ibottom, jright, itop

! The range to delineate
! The values can be calculated by calc_area.f90 program
! for a target catchment
!parameter( jleft   = 1101 )
!parameter( ibottom = 683 )
!parameter( jright  = 1130 )
!parameter( itop    = 665 )

! Time Step (sec) <86400>
integer tstep
!parameter( tstep = 86400 )

! Other Variables
integer i,j,k,n,ios,t
integer idim,jdim
!parameter(idim=1200, jdim=3600) ! 0.1 deg
!parameter(idim=480, jdim=1440) ! 0.25 deg

character*256 list_day, list_time, ctemp

!real*4 rain(idim, jdim), sumRain(idim, jdim)
real*4 rain(1200, 3600), sumRain(1200, 3600)

real horizontal_resolution
integer temporal_resolution

open(0, file = areafile, status = 'old')
read(0, '(a28, f20.5)') ctemp, horizontal_resolution
read(0, '(a28, i20)') ctemp, temporal_resolution
write(*, '("holirontal_resolution [d] : ", f20.5)') horizontal_resolution
write(*, '("temporal_resolution [h]    : ", f20.5)') temporal_resolution
do i = 1, 11
 read(0, '(a)') ctemp
enddo
read(0, '(a12, i15)') ctemp, jleft
read(0, '(a12, i15)') ctemp, ibottom
read(0, '(a12, i15)') ctemp, jright
read(0, '(a12, i15)') ctemp, itop
write(*,'("jleft   : ", i15)') jleft
write(*,'("ibottom : ", i15)') ibottom
write(*,'("jright  : ", i15)') jright
write(*,'("itop    : ", i15)') itop

if( horizontal_resolution .eq. 0.1 ) then
 idim = 1200
 jdim = 3600
elseif( horizontal_resolution .eq. 0.25 ) then
 idim = 480
 jdim = 1440
else
 write(*,*) "Error : horizontal_resolution must be 0.1 or 0.25"
 stop
endif

tstep = 3600 * temporal_resolution

open(1, file = listfile, status = 'old')
open(1000, file = rainfile)
open(1100, file = sumfile)

t = 0
sumRain = 0.d0
do

 read(1,*,iostat=ios) list_day
 if(ios.ne.0) exit

 infile = trim(infile1) // list_day
 cutfile = trim(cutfile1) // list_day

 open(10,file= infile, &
  form='unformatted',access='direct',recl=4*idim*jdim, &
  status='old') 

 read(10,rec=1)((rain(i,j),j=1,jdim),i=1,idim) ! Modified

 close(10)

 do i = itop, ibottom
  do j = jleft, jright
   if(rain(i,j).lt.0) rain(i,j) = 0.0
  enddo
 enddo

 open(100, file = cutfile)
 do i = itop, ibottom
  write(100,'(100000f7.2)') (rain(i,j), j = jleft, jright)
 enddo
 close(100)

 write(1000,*) t * tstep, (jright - jleft) + 1, (ibottom - itop) + 1
 do i = itop, ibottom
  write(1000,'(100000f7.2)') (rain(i,j), j = jleft, jright)
 enddo
 t = t + 1

 sumRain = sumRain + rain * temporal_resolution

enddo
close(1000)

do i = itop, ibottom
 write(1100,'(100000f8.2)') (sumRain(i,j), j = jleft, jright)
enddo

end

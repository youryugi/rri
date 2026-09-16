! makeRiver2.f90
!
! coded by T.Sayama on April 21, 2010
!
implicit none

character*256 infile_acc, outfile_width, outfile_depth, outfile_height

!parameter( infile_acc = "../../Model/infile/indus/acc_indus_60s_lg.txt" )
!parameter( outfile_width = "../../Model/infile/indus/width_indus_60s_lg.txt" )
!parameter( outfile_depth = "../../Model/infile/indus/depth_indus_60s_lg_6b.txt" )
!parameter( outfile_height = "../../Model/infile/indus/height_indus_60s_lg.txt" )

!parameter( infile_acc = "../../Model/infile/kabul/acc_pak.txt" )
!parameter( outfile_width = "../../Model/infile/kabul/width_pak.txt" )
!parameter( outfile_depth = "../../Model/infile/kabul/depth_pak.txt" )
!parameter( outfile_height = "../../Model/infile/kabul/height_pak.txt" )

!parameter( infile_acc = "../../Model/infile/thai/acc_thai_60s.txt" )
!parameter( outfile_width = "../../Model/infile/thai/width_thai_60s.txt" )
!parameter( outfile_depth = "../../Model/infile/thai/depth_thai_60s.txt" )
!parameter( outfile_height = "../../Model/infile/thai/height_thai_60s.txt" )

parameter( infile_acc = "../../Model/infile/chao/acc_chao_30s.txt" )
parameter( outfile_width = "../../Model/infile/chao/width_chao_30s.txt" )
parameter( outfile_depth = "../../Model/infile/chao/depth_chao_30s.txt" )
parameter( outfile_height = "../../Model/infile/chao/height_chao_30s.txt" )

! 変数定義
integer ncols, nrows
real xllcorner, yllcorner, cellsize, nodata

integer i, j, k

integer, dimension(:,:), allocatable :: acc
real, dimension(:,:), allocatable :: width, depth, height, area

character*256 ctemp

! STEP 0 : Open Files
open(10, file = infile_acc, status = "old")
open(30, file = outfile_width)
open(40, file = outfile_depth)
open(50, file = outfile_height)

! STEP 1 : Reading File
read(10, *) ctemp, ncols
read(10, *) ctemp, nrows
read(10, *) ctemp, xllcorner
read(10, *) ctemp, yllcorner
read(10, *) ctemp, cellsize
read(10, *) ctemp, nodata

allocate( acc(nrows, ncols), width(nrows, ncols), depth(nrows, ncols), height(nrows, ncols), area(nrows, ncols) )

acc = 0.
width = -9999.
depth = -9999.
height = -9999.
area = -9999.

rewind(10)

read(10, '(a30)') ctemp
write(30, '(a30)') ctemp
write(40, '(a30)') ctemp
write(50, '(a30)') ctemp

read(10, '(a30)') ctemp
write(30, '(a30)') ctemp
write(40, '(a30)') ctemp
write(50, '(a30)') ctemp

read(10, '(a30)') ctemp
write(30, '(a30)') ctemp
write(40, '(a30)') ctemp
write(50, '(a30)') ctemp

read(10, '(a30)') ctemp
write(30, '(a30)') ctemp
write(40, '(a30)') ctemp
write(50, '(a30)') ctemp

read(10, '(a30)') ctemp
write(30, '(a30)') ctemp
write(40, '(a30)') ctemp
write(50, '(a30)') ctemp

read(10, '(a30)') ctemp
write(30, '(a30)') ctemp
write(40, '(a30)') ctemp
write(50, '(a30)') ctemp

do i = 1, nrows
 read(10, *) (acc(i, j), j = 1, ncols)
enddo

! STEP 2 : Widthの計算

do i = 1, nrows
 do j = 1, ncols

  !if( acc(i, j) .ge. 10 ) then ! Kabul
  !if( acc(i, j) .ge. 100 ) then ! Indus, thai
  if( acc(i, j) .ge. 1000 ) then ! Chao

   ! (Scaleupを実行しても、acc は元の値になっているので、Scaleup前のグリッドサイズで面積を計算)
   !area(i, j) = acc(i, j) * 0.761 * 0.924  ! [km2] (Kabul)
   !area(i, j) = acc(i, j) * 0.46 * 0.46  ! [km2] (Solo)
   !area(i, j) = acc(i, j) * 0.820 * 0.923  ! [km2] (Indus / Guddu / Thai)
   area(i, j) = acc(i, j) * 0.887 * 0.922  ! [km2] (Chao)
 
   ! Kabul
   !width(i, j) = 2.5 * area(i, j) ** 0.4  ! [m]
   !depth(i, j) = 0.1 * area(i, j) ** 0.4  ! [m]
   !height(i, j) = 0.0
   !if(depth(i, j) .gt. 5.0) depth(i, j) = 5.0

   ! Solo
   !width(i, j) = 5. * area(i, j) ** 0.35  ! [m]
   !depth(i, j) = 0.95 * area(i, j) ** 0.2  ! [m]
   !depth(i, j) = 3.0 * area(i, j) ** 0.2  ! [m]
   !height(i, j) = 0.0
   !if( acc(i, j) .ge. 47150 ) height(i, j) = 4.0
   !if( acc(i, j) .ge. 6600 .and. acc(i, j) .lt. 47150 ) height(i, j) = 4.0

   ! Indus / Guddu
   !width(i, j) = 3.5 * area(i, j) ** 0.4   ! [m]
   !depth(i, j) = 0.1 * area(i, j) ** 0.4   ! [m]
   !if(depth(i, j) .gt. 5.0) depth(i, j) = 5.0
   !height(i, j) = 0.0

   ! Thai
   width(i, j) = 16.931 * area(i, j) ** 0.186   ! [m]
   !if(width(i, j) .gt. 100.0) width(i, j) = 100.0
   depth(i, j) = 2.481 * area(i, j) ** 0.1197   ! [m]
   !if(depth(i, j) .gt. 6.0) depth(i, j) = 6.0
   height(i, j) = 0.0

   !height(i, j) = 10.0
   !if(acc(i, j) .ge. 1013301 .and. acc(i, j) .le. 1018029 ) height (i, j) = 0.0

   ! Main river
   !if(acc(i, j) .ge. 490320 ) depth(i, j) = 2.0

  endif

 enddo
enddo

! STEP 3 : output
do i = 1, nrows
 write(30, '(<nrows>f9.2)') ( width(i,j), j = 1, ncols )
 write(40, '(<nrows>f9.2)') ( depth(i,j), j = 1, ncols )
 write(50, '(<nrows>f9.2)') ( height(i,j), j = 1, ncols )
enddo

end

! flowDirection
!
! originally developed by T.Kojima on Apr. 10, 2003 (mesh1v3)
! revised by T.Sayama on Oct. 31, 2014 as flowDirection for RRI
!
implicit integer(f)
character*256 dem_file,acc_file,dir_file,riv_file
integer gostrt,check, mmm, nnn, dir_zero, riv_check, thresh
integer x1, y1
doubleprecision e1, n1, ee0, nn0
integer,dimension(:,:), allocatable :: muse, ksuse
real,dimension(:,:), allocatable :: u
integer,dimension(:), allocatable :: nsuc, nda, nsuctmp
integer,dimension(:), allocatable :: catch_temp, catch, riv

character*256 ctemp
dimension v(9)

! mmm : number of column 
! nnn : number of row 

open(1, file = "flowDirection.txt", status = 'old')
read(1, '(a)') dem_file
read(1, '(a)') dir_file
read(1, '(a)') acc_file
read(1, *) riv_switch
if(riv_switch .eq. 1) then
 read(1, '(a)') riv_file
endif
close(1)
g = 1000000

open(10, file = dem_file, status = 'old')
read(10, *) ctemp, mds
read(10, *) ctemp, nds
read(10, *) ctemp, e1
read(10, *) ctemp, n1
read(10, *) ctemp, xscale
read(10, *) ctemp, nodata

rewind(10)

open(11, file = dir_file)
open(12, file = acc_file)

do i = 1, 6
 read(10, '(a40)') ctemp
 write(11, '(a40)') ctemp
 write(12, '(a40)') ctemp
enddo

yscale = xscale

ee0 = e1
nn0 = n1
mmm=mds
nnn=nds

allocate(u(-2:mmm+3,-2:nnn+3))
allocate(nsuc(1:mmm*nnn))
allocate(nda(1:mmm*nnn))
allocate(muse(-2:mmm+3,-2:nnn+3))
allocate(ksuse(-2:mmm+3,-2:nnn+3))
allocate(nsuctmp(1:mmm))
allocate(catch_temp(1:mmm*nnn))
allocate(catch(1:mmm*nnn))
allocate(riv(1:mmm*nnn))

yschi=xscale/yscale
mn=mds*nds
s2=sqrt(1+yschi*yschi)

u(:,:) = -9999
do n = 1, nds
 read (10,*) (u(m,n), m = 1, mds)
enddo
close(10)

do j = 1, mn
 nda(j) = 1
 nsuc(j) = -9999
 riv(j) = 0
enddo

riv_check = 0
if(riv_switch .eq. 1) then
 open(13, file = riv_file, status = 'old')
 read(13, *) ctemp, mds
 read(13, *) ctemp, nds
 read(13, *) ctemp, e1
 read(13, *) ctemp, n1
 read(13, *) ctemp, xscale
 read(13, *) ctemp, nodata

 noi=1
 do j=1,nds
  read (13,*) (nsuctmp(m), m=1,mds)
  do i=1,mds
   if(nsuctmp(i)==0) then
    x1 = i
    y1 = j
    nsuc(noi) = 0
   else if(nsuctmp(i)==100) then ! river cell
    x1 = i
    y1 = j
    nsuc(noi) = -9999
    riv(noi) = 1
    riv_check = 1
   else if(nsuctmp(i)==1) then
    x1 = i + 1
    y1 = j
    nsuc(noi) = fnno(x1,y1,mds)
   else if(nsuctmp(i)==2) then
    x1 = i + 1
    y1 = j + 1
    nsuc(noi) = fnno(x1,y1,mds)
   else if(nsuctmp(i)==4) then
    x1 = i
    y1 = j + 1
    nsuc(noi) = fnno(x1,y1,mds)
   else if(nsuctmp(i)==8) then
    x1 = i - 1
    y1 = j + 1
    nsuc(noi) = fnno(x1,y1,mds)
   else if(nsuctmp(i)==16) then
    x1 = i - 1
    y1 = j
    nsuc(noi) = fnno(x1,y1,mds)
   else if(nsuctmp(i)==32) then
    x1 = i - 1
    y1 = j - 1
    nsuc(noi) = fnno(x1,y1,mds)
   else if(nsuctmp(i)==64) then
    x1 = i
    y1 = j - 1
    nsuc(noi) = fnno(x1,y1,mds)
   else if(nsuctmp(i)==128) then
    x1 = i + 1
    y1 = j - 1
    nsuc(noi) = fnno(x1,y1,mds)
   else
    !write(*,*) "error in input error: ", nsuctmp(i)
   endif
   noi = noi + 1
  enddo
 enddo
close(13)
endif

mm1 = 0
mp1 = 0
nm1 = 0
np1 = 0

do j = 1, nds
 do i = 1, mds
  ! エッジのグリッドセルは nsuc = -1 にする
  if(j.eq.1 .or. j.eq.nds .or. i.eq.1 .or. i.eq.mds) then
   now=fnno(i,j,mds)
   nsuc(now) = -1
  endif
 enddo
enddo

1001 continue

! main loop
do j=1, nds

 if(riv_check.eq.1) then
  write(*,*) 'River Processing j=',j,'/',nds
 else
  write(*,*) 'Processing j=',j,'/',nds
 endif

 do i=1,mds
  m=i
  n=j
  now=fnno(m,n,mds)

  if( riv_check .eq. 1 .and. riv(now) .ne. 1 ) cycle ! -> 河川セルでない

  ! 領域外もしくは nsuc = 0 の点
  if( nsuc(now) .eq. -1 ) cycle ! エッジグリッドセル
  if( nsuc(now) .eq. 0 .or. u(i,j).lt.-100 ) cycle

  if(nsuc(now) .gt. 0) then
   ! nsuc(now)が既に決まっている
   mm=fnlocx(nsuc(now),mds)
   nn=fnlocy(nsuc(now),mds)
   goto 1003
  end if

  muse(:,:) = 1
  noi = 1
  do jj = 1, nds
   do ii = 1, mds
    if( riv_check .eq. 0 ) then
     if( u(ii,jj) .ge. -100. ) muse(ii,jj) = 0 ! 領域内だけ muse = 0 にする
    else
     if( u(ii,jj) .ge. -100. .and. riv(noi) .eq. 1 )  muse(ii,jj) = 0 ! 河川セルだけ muse = 0 にする
     if( u(ii,jj) .ge. -100. .and. nsuc(noi) .eq. 0 )  muse(ii,jj) = 0 ! 河口セルも muse = 0 にする
     if( u(ii,jj) .ge. -100. .and. nsuc(noi) .ge. 1 )  muse(ii,jj) = 0 ! 既に決まっているセルも muse = 0 にする
    endif
    noi = noi + 1
   enddo
  enddo

  ! loop
  1002 continue

  muse(m,n)=1 ! 通ってきた道に戻らないように muse(m, n) = 1 にする

  25 gostrt=0
  ! 通ってきた道と斜めでクロスしないように ksuse = 1 にする
  call pit(m,n,nsuc,muse,mp1,mm1,np1,nm1,mmm,nnn,mds,gostrt,ksuse)
  if(gostrt.eq.1) then
   ! (m, n)からは8方向どこにも行き場がない
   ksuse(mm1,nm1)=0
   ksuse(mp1,nm1)=0
   ksuse(mm1,np1)=0
   ksuse(mp1,np1)=0

   call back(m,n,mds,nda,nsuc,mmm,nnn) ! nsuc(m, n) = 0 にして、(m, n)を一つ戻す

   now=fnno(m,n,mds)
   goto 25
  end if

  ! (m, n) から 最急勾配方向(minn) と 勾配(sminn) を検索
  call calgrd(m,n,u,muse,mm1,ksuse,mp1,nm1,np1,s2,v,yschi,minn,smin,mmm,nnn,g,mds)

  if(smin.gt.0.) then ! 全て逆勾配
   ! 周囲２セル分みて最急勾配のセルに移動 (m, n) -> (mm, nn)
   call hosei(m,n,u,muse,g,s2,yschi,mm,nn,mmm,nnn,mds,ksuse,mp1,mm1,np1,nm1)
   nsuc(now)=fnno(mm,nn,mds)
   goto 1003
  endif

  io=0
  do k=1,9
   if(v(k).eq.smin) io = io + 1
  enddo

  if(io.eq.1) then
   ! 最急勾配のセルに移動 (m, n) -> (mm, nn)
   call onmin(minn,m,n,mm,nn)
   nsuc(now)=fnno(mm,nn,mds)
   goto 1003
  else
   ! 同一勾配のセルが複数ある
   call samegr(m,n,muse,u,mm,nn,mmm,nnn,g,smin,v,ksuse,mds,mp1,mm1,np1,nm1)
   nsuc(now)=fnno(mm,nn,mds)
  end if

  1003 continue

  ksuse(mm1,nm1)=0
  ksuse(mp1,nm1)=0
  ksuse(mm1,np1)=0
  ksuse(mp1,np1)=0

  next=fnno(mm,nn,mds)
  check=0

  1004 nda(next)=nda(next)+1

  if( nsuc(next) .eq. 0) cycle ! 出口グリッドセルに到達
  if( nsuc(next) .eq. -1 ) cycle ! エッジグリッドセルに到達

  if(nsuc(next) .gt. 0) then
   ! nsuc(next)が0以外で、既に決まっている場合は、ndaを一つ足すために(mm, nn)を進める
   next=nsuc(next)
   mm=fnlocx(next,mds)
   nn=fnlocy(next,mds)
   check=check+1
   if(check.gt.10000) write(6,*) mm,nn
   if(check.gt.10010) then
    write(*,*) "Error: incomplete"
    stop
   endif
   goto 1004
  end if

  m=mm
  n=nn
  now=next

  goto 1002

 enddo
enddo

if( riv_check .eq. 1 ) then
 riv_check = 0
 go to 1001
endif

dir_zero = 0
noi = 1
do j = 1, nds
 do i = 1, mds
  if( nsuc(noi) .eq. 0 ) dir_zero = 1
  noi = noi + 1
 enddo
enddo

if( dir_zero .eq. 0 ) then
 catch(:) = 1
else
 ! dir = 0 の上流域を検索
 catch(:) = 0
 do j = 1, nds
  write(*,*) 'delineate j=',j,'/',nds
  do i = 1, mds
   catch_temp(:) = 0
   m = i
   n = j
   do
    now=fnno(m,n,mds)
    if( nsuc(now) .lt. 0 ) exit ! エッジグリッドセル
    catch_temp(now) = 1
    if( nsuc(now) .eq. 0 ) then
     catch(:) = catch(:) + catch_temp(:)
     exit
    elseif( nsuc(now) .eq. -1 ) then
     catch_temp(now) = 0
     exit
    endif
    now=nsuc(now)
    m=fnlocx(now,mds)
    n=fnlocy(now,mds)
   enddo
  enddo
 enddo
 where(catch(:).ge.1) catch(:) = 1
endif

noi = 1
do j=1, nds
 do i=1, mds
  x1 = fnlocx(nsuc(noi),mds)
  y1 = fnlocy(nsuc(noi),mds)
  if( u(i,j) .lt. -100.  .or. catch(noi) .eq. 0 ) then
   nsuctmp(i) = -9999
  elseif( nsuc(noi) .eq. 0 .or. nsuc(noi) .eq. -1 )then
   nsuctmp(i) = 0
  else
   if(x1.lt.i) then
    if(y1.lt.j) nsuctmp(i) = 32
    if(y1.gt.j) nsuctmp(i) = 8
    if(y1.eq.j) nsuctmp(i) = 16
   else if(x1.gt.i) then
    if(y1.lt.j) nsuctmp(i) = 128
    if(y1.gt.j) nsuctmp(i) = 2
    if(y1.eq.j) nsuctmp(i) = 1
   else if(x1.eq.i) then
    if(y1.lt.j) nsuctmp(i) = 64
    if(y1.gt.j) nsuctmp(i) = 4
    if(y1.eq.j) nsuctmp(i) = 0
   endif
  endif
  noi = noi+1
 enddo
 write(11, '(10000i8)') (nsuctmp(i),i=1,mds)
enddo

noi = 1
do j=1, nds
 do i=1, mds
  nsuctmp(i) = nda(noi)
  if( u(i,j) .lt. -100. .or. catch(noi) .eq. 0 ) nsuctmp(i) = -9999 
  noi = noi+1
 enddo
 write(12, '(10000i8)') (nsuctmp(i),i=1,mds)
enddo

close(11)
close(12)

stop
end


subroutine onmin(min,m,n,mm,nn)

! ArcGISの流れ方向の定義とは違うので要注意
! 1 2 3
! 4 5 6
! 7 8 9

if(min.eq.1) then
 mm=m-1
 nn=n-1
else if(min.eq.2) then
 mm=m
 nn=n-1
else if(min.eq.3) then
 mm=m+1
 nn=n-1
else if(min.eq.4) then
 mm=m-1
 nn=n
else if(min.eq.5) then
 mm=m
 nn=n
else if(min.eq.6) then
 mm=m+1
 nn=n
else if(min.eq.7) then
 mm=m-1
 nn=n+1
else if(min.eq.8) then
 mm=m
 nn=n+1
else if(min.eq.9) then
 mm=m+1
 nn=n+1
end if
return
end

subroutine hosei(m,n,u,muse,g,s2,yschi,mm,nn,mmm,nnn,mds,ksuse,mp1,mm1,np1,nm1)
implicit integer(f)
real u,uu
dimension u(-2:mmm+3,-2:nnn+3),muse(-2:mmm+3,-2:nnn+3),ksuse(-2:mmm+3,-2:nnn+3)
dimension w(25)
!fnno(ml,nl,mds)=(nl-1)*mds+ml
!fnlocx(no,mds)=mod((no-1),mds)+1
!fnlocy(no,mds)=(no-1)/mds+1
!
!  1  2  3  4  5
!  6  7  8  9 10
! 11 12 13 14 15
! 16 17 18 19 20
! 21 22 23 24 25
!
uu=u(m,n)
mp2=m+2
mm2=m-2
np2=n+2
nm2=n-2
m1=muse(mm2,nm2)
m2=muse(mm1,nm2)
m3=muse(m,nm2)
m4=muse(mp1,nm2)
m5=muse(mp2,nm2)
m6=muse(mm2,nm1)
m7=muse(mm1,nm1)
m8=muse(m,nm1)
m9=muse(mp1,nm1)
m10=muse(mp2,nm1)
m11=muse(mm2,n)
m12=muse(mm1,n)
m14=muse(mp1,n)
m15=muse(mp2,n)
m16=muse(mm2,np1)
m17=muse(mm1,np1)
m18=muse(m,np1)
m19=muse(mp1,np1)
m20=muse(mp2,np1)
m21=muse(mm2,np2)
m22=muse(mm1,np2)
m23=muse(m,np2)
m24=muse(mp1,np2)
m25=muse(mp2,np2)
!*---------   cross check   ----------
if(ksuse(mm1,nm1).eq.1) m7=1
if(ksuse(mp1,nm1).eq.1) m9=1
if(ksuse(mm1,np1).eq.1) m17=1
if(ksuse(mp1,np1).eq.1) m19=1
!*------------------------------------
if(m1.eq.0 .and. m7.eq.0) then
 w(1)=real(u(mm2,nm2))
else
 w(1)=g
end if
if(m2.eq.0 .and. (m7.eq.0 .or. m8.eq.0)) then
 w(2)=real(u(mm1,nm2))
else
 w(2)=g
end if
if(m3.eq.0 .and. (m7.eq.0 .or. m8.eq.0 .or. m9.eq.0))then
 w(3)=real(u(m,nm2))
else
 w(3)=g
end if
if(m4.eq.0 .and. (m8.eq.0 .or. m9.eq.0)) then
 w(4)=real(u(mp1,nm2))
else
 w(4)=g
end if
if(m5.eq.0 .and. m9.eq.0) then
 w(5)=real(u(mp2,nm2))
else
 w(5)=g
end if
if(m6.eq.0 .and. (m7.eq.0 .or. m12.eq.0)) then
 w(6)=real(u(mm2,nm1))
else
 w(6)=g
end if
if(m7.eq.0) then
 w(7)=real(u(mm1,nm1))
 y7=(w(7)-uu)/s2
else
 w(7)=g
 y7=g
end if
if(m8.eq.0) then
 w(8)=real(u(m,nm1))
 y8=w(8)-uu
else
 w(8)=g
 y8=g
end if
if(m9.eq.0) then
 w(9)=real(u(mp1,nm1))
 y9=(w(9)-uu)/s2
else
 w(9)=g
 y9=g
end if
if(m10.eq.0 .and. (m9.eq.0 .or. m14.eq.0)) then
 w(10)=real(u(mp2,nm1))
else
 w(10)=g
end if
if(m11.eq.0 .and. (m7.eq.0 .or. m12.eq.0 .or. m17.eq.0))then
 w(11)=real(u(mm2,n))
else
 w(11)=g
end if
if(m12.eq.0) then
 w(12)=real(u(mm1,n))
 y12=(w(12)-uu)/yschi
else
 w(12)=g
 y12=g
end if
w(13)=1.e8
y13=w(13)
if(m14.eq.0) then
 w(14)=real(u(mp1,n))
 y14=(w(14)-uu)/yschi
else
 w(14)=g
 y14=g
end if
if(m15.eq.0 .and. (m9.eq.0 .or. m14.eq.0 .or. m19.eq.0)) then
 w(15)=real(u(mp2,n))
else
 w(15)=g
end if
if(m16.eq.0 .and. (m12.eq.0 .or. m17.eq.0)) then
 w(16)=real(u(mm2,np1))
else
 w(16)=g
end if
if(m17.eq.0) then
 w(17)=real(u(mm1,np1))
 y17=(w(17)-uu)/s2
else
 w(17)=g
 y17=g
end if
if(m18.eq.0) then
 w(18)=real(u(m,np1))
 y18=w(18)-uu
else
 w(18)=g
 y18=g
end if
if(m19.eq.0) then
 w(19)=real(u(mp1,np1))
 y19=(w(19)-uu)/s2
else
 w(19)=g
 y19=g
end if
if(m20.eq.0 .and. (m14.eq.0 .or. m19.eq.0)) then
 w(20)=real(u(mp2,np1))
else
 w(20)=g
end if
if(m21.eq.0 .and. m17.eq.0) then
 w(21)=real(u(mm2,np2))
else
 w(21)=g
end if
if(m22.eq.0 .and. (m17.eq.0 .or. m18.eq.0)) then
  w(22)=real(u(mm1,np2))
else
  w(22)=g
end if
if(m23.eq.0 .and. (m17.eq.0 .or. m18.eq.0 .or. m19.eq.0)) then
  w(23)=real(u(m,np2))
else
  w(23)=g
end if
if(m24.eq.0 .and. (m18.eq.0 .or. m19.eq.0)) then
  w(24)=real(u(mp1,np2))
else
  w(24)=g
end if
if(m25.eq.0 .and. m19.eq.0) then
  w(25)=real(u(mp2,np2))
else
  w(25)=g
end if
wmin=1.e10
do k=1,25
  if(wmin.lt.w(k)) cycle
  wmin=w(k)
  k1=k
enddo
if(k1.eq.1 .or. k1.eq.7) goto 500
if(k1.eq.5 .or. k1.eq.9) goto 502
if(k1.eq.17 .or. k1.eq.21) goto 506
if(k1.eq.19 .or. k1.eq.25) goto 508

if(k1.eq.2) goto 400
if(k1.eq.3) goto 401
if(k1.eq.4) goto 402
if(k1.eq.6) goto 403
if(k1.eq.8) goto 501
if(k1.eq.10) goto 404
if(k1.eq.11) goto 405
if(k1.eq.12) goto 503
if(k1.eq.13) goto 504
if(k1.eq.14) goto 505
if(k1.eq.15) goto 407
if(k1.eq.16) goto 406
if(k1.eq.18) goto 507
if(k1.eq.20) goto 408
if(k1.eq.22) goto 409
if(k1.eq.23) goto 410
if(k1.eq.24) goto 411

400 if((y7.lt.y8 .and. m7.eq.0) .or. m8.ne.0) then
 goto 500
else
 goto 501
end if
401 if((y7.lt.y8 .and. m7.eq.0 .and. y7.lt.y9).or.(m8.ne.0.and. m9.ne.0)) then
 goto 500
else
 goto 402
end if
402 if((y8.lt.y9 .and. m8.eq.0).or. m9.ne.0) then
 goto 501
else
 goto 502
end if
403 if((y7.lt.y12 .and. m7.eq.0).or. m12.ne.0) then
 goto 500
else
 goto 503
end if
404 if((y9.lt.y14 .and. m9.eq.0).or. m14.ne.0) then
 goto 502
else
 goto 505
end if
405 if((y7.lt.y12 .and. m7.eq.0 .and. y7.lt.y17) .or.(m12.ne.0 .and. m17.ne.0)) then
 goto 500
else
 goto 406
end if
406 if((y12.lt.y17 .and. m12.eq.0).or. m17.ne.0) then
 goto 503
else
 goto 506
end if
407 if((y9.lt.y14 .and. y9.lt.y19 .and. m9.eq.0) .or.(m14.ne.0 .and. m19.ne.0)) then
 goto 502
else
 goto 408
end if
408 if((y14.lt.y19 .and. m14.eq.0).or. m19.ne.0) then
 goto 505
else
 goto 508
end if
409 if((y17.lt.y18 .and. m17.eq.0).or. m18.ne.0) then
 goto 506
else
 goto 507
end if
410 if((y17.lt.y18 .and. y17.lt.y19 .and. m17.eq.0) .or.(m18.ne.0 .and. m19.ne.0)) then
 goto 506
else
 goto 411
end if
411 if((y18.lt.y19 .and. m18.eq.0).or. m19.ne.0) then
 goto 507
else
 goto 508
end if
500 mm=m-1
nn=n-1
if(muse(mm,nn).eq.1) write(6,*) '1',m,n,mm,nn
goto 600
501 mm=m
nn=n-1
if(muse(mm,nn).eq.1) write(6,*) '2',m,n,mm,nn
goto 600
502 mm=m+1
nn=n-1
if(muse(mm,nn).eq.1) write(6,*) '3',m,n,mm,nn
goto 600
503 mm=m-1
nn=n
if(muse(mm,nn).eq.1) write(6,*) '4',m,n,mm,nn
goto 600
504 mm=m
nn=n
if(muse(mm,nn).eq.1) write(6,*) '5',m,n,mm,nn
goto 600
505 mm=m+1
nn=n
if(muse(mm,nn).eq.1) write(6,*) '6',m,n,mm,nn
goto 600
506 mm=m-1
nn=n+1
if(muse(mm,nn).eq.1) write(6,*) '7',m,n,mm,nn
goto 600
507 mm=m
nn=n+1
if(muse(mm,nn).eq.1) write(6,*) '8',m,n,mm,nn
goto 600
508 mm=m+1
nn=n+1
if(muse(mm,nn).eq.1) write(6,*) '9',m,n,mm,nn
600 return
end

subroutine samegr(m,n,muse,u,mm,nn,mmm,nnn,g,smin,v,ksuse,mds,mp1,mm1,np1,nm1)
implicit integer(f)
real u
dimension u(-2:mmm+3,-2:nnn+3),muse(-2:mmm+3,-2:nnn+3),ksuse(-2:mmm+3,-2:nnn+3)
dimension x(12),t(9),v(9)
!fnno(ml,nl,mds)=(nl-1)*mds+ml
!fnlocx(no,mds)=mod((no-1),mds)+1
!fnlocy(no,mds)=(no-1)/mds+1
mm3=m-3
mm2=m-2
mp2=m+2
mp3=m+3
nm3=n-3
nm2=n-2
np2=n+2
np3=n+3
x(1)=real(u(mm3,nm2)+u(mm2,nm2)+u(mm1,nm2)+u(mm3,nm1)+u(mm2,nm1)+u(mm1,nm1))
x(2)=real(u(mm2,nm3)+u(mm2,nm2)+u(mm2,nm1)+u(mm1,nm3)+u(mm1,nm2)+u(mm1,nm1))
x(3)=real(u(mm1,nm2)+u(mm1,nm1)+u(m,nm2)+u(m,nm1)+u(mp1,nm2)+u(mp1,nm1))
x(4)=real(u(mp1,nm3)+u(mp1,nm2)+u(mp1,nm1)+u(mp2,nm3)+u(mp2,nm2)+u(mp2,nm1))
x(5)=real(u(mp1,nm2)+u(mp1,nm1)+u(mp2,nm2)+u(mp2,nm1)+u(mp3,nm2)+u(mp3,nm1))
x(6)=real(u(mp1,nm1)+u(mp1,n)+u(mp1,np1)+u(mp2,nm1)+u(mp2,n)+u(mp2,np1))
x(7)=real(u(mp1,np1)+u(mp1,np2)+u(mp2,np1)+u(mp2,np2)+u(mp3,np1)+u(mp3,np2))
x(8)=real(u(mp1,np1)+u(mp1,np2)+u(mp1,np3)+u(mp2,np1)+u(mp2,np2)+u(mp2,np3))
x(9)=real(u(mm1,np1)+u(mm1,np2)+u(m,np1)+u(m,np2)+u(mp1,np1)+u(mp1,np2))
x(10)=real(u(mm2,np1)+u(mm2,np2)+u(mm2,np3)+u(mm1,np1)+u(mm1,np2)+u(mm1,np3))
x(11)=real(u(mm3,np1)+u(mm3,np2)+u(mm2,np1)+u(mm2,np2)+u(mm1,np1)+u(mm1,np2))
x(12)=real(u(mm2,nm1)+u(mm2,n)+u(mm2,np1)+u(mm1,nm1)+u(mm1,n)+u(mm1,np1))
if(muse(mm3,nm2).ne.0 .or. muse(mm2,nm2).ne.0.or. muse(mm1,nm2).ne.0 .or. muse(mm3,nm1).ne.0.or. muse(mm2,nm1).ne.0 .or. muse(mm1,nm1).ne.0)then
 x(1)=x(1)+10000.
end if
if(muse(mm2,nm3).ne.0 .or. muse(mm2,nm2).ne.0.or. muse(mm2,nm1).ne.0 .or. muse(mm1,nm3).ne.0.or. muse(mm1,nm2).ne.0 .or. muse(mm1,nm1).ne.0)then
 x(2)=x(2)+10000.
end if
if(muse(mp1,nm1).ne.0 .or. muse(mm1,nm2).ne.0.or. muse(m,nm2).ne.0 .or. muse(m,nm1).ne.0.or. muse(mp1,nm2).ne.0 .or. muse(mm1,nm1).ne.0)then
 x(3)=x(3)+10000.
end if
if(muse(mp1,nm3).ne.0 .or. muse(mp1,nm2).ne.0.or. muse(mp1,nm1).ne.0 .or. muse(mp2,nm3).ne.0.or. muse(mp2,nm2).ne.0 .or. muse(mp2,nm1).ne.0)then
 x(4)=x(4)+10000.
end if
if(muse(mp1,nm2).ne.0 .or. muse(mp1,nm1).ne.0.or. muse(mp2,nm2).ne.0 .or. muse(mp2,nm1).ne.0.or. muse(mp3,nm2).ne.0 .or. muse(mp3,nm1).ne.0)then
 x(5)=x(5)+10000.
end if
if(muse(mp1,nm1).ne.0 .or. muse(mp1,n).ne.0.or. muse(mp1,np1).ne.0 .or. muse(mp2,nm1).ne.0.or. muse(mp2,n).ne.0 .or. muse(mp2,np1).ne.0)then
 x(6)=x(6)+10000.
end if
if(muse(mp1,np1).ne.0 .or. muse(mp1,np2).ne.0.or. muse(mp2,np1).ne.0 .or. muse(mp2,np2).ne.0.or. muse(mp3,np1).ne.0 .or. muse(mp3,np2).ne.0)then
 x(7)=x(7)+10000.
end if
if(muse(mp1,np1).ne.0 .or. muse(mp1,np2).ne.0.or. muse(mp1,np3).ne.0 .or. muse(mp2,np1).ne.0.or. muse(mp2,np2).ne.0 .or. muse(mp2,np3).ne.0)then
 x(8)=x(8)+10000.
end if
if(muse(mm1,np1).ne.0 .or. muse(mm1,np2).ne.0.or. muse(m,np1).ne.0 .or. muse(m,np2).ne.0.or. muse(mp1,np1).ne.0 .or. muse(mp1,np2).ne.0)then
 x(9)=x(9)+10000.
end if
if(muse(mm2,np1).ne.0 .or. muse(mm2,np2).ne.0.or. muse(mm2,np3).ne.0 .or. muse(mm1,np1).ne.0.or. muse(mm1,np2).ne.0 .or. muse(mm1,np3).ne.0)then
 x(10)=x(10)+10000.
end if
if(muse(mm3,np1).ne.0 .or. muse(mm3,np2).ne.0.or. muse(mm2,np1).ne.0 .or. muse(mm2,np2).ne.0.or. muse(mm1,np1).ne.0 .or. muse(mm1,np2).ne.0)then
 x(11)=x(11)+10000.
end if
if(muse(mm2,nm1).ne.0 .or. muse(mm2,n).ne.0.or. muse(mm2,np1).ne.0 .or. muse(mm1,nm1).ne.0.or. muse(mm1,n).ne.0 .or. muse(mm1,np1).ne.0)then
 x(12)=x(12)+10000.
end if
if(x(1).gt.x(2)) then
 t(1)=x(2)
else
 t(1)=x(1)
end if
if(muse(mm1,nm1).eq.1 .or. ksuse(mm1,nm1).eq.1) t(1)=g
t(2)=x(3)
if(muse(mm1,nm1).eq.1 .and. muse(m,nm1).eq.1 .and. muse(mp1,nm1).eq.1) t(2)=g
if(x(4).gt.x(5)) then
 t(3)=x(5)
else
 t(3)=x(4)
end if
if(muse(mp1,nm1).eq.1 .or. ksuse(mp1,nm1).eq.1) t(3)=g
t(5)=g
t(6)=x(6)
if(muse(mp1,nm1).eq.1 .and. muse(mp1,n).eq.1 .and. muse(mp1,np1).eq.1) t(6)=g
if(x(7).gt.x(8)) then
 t(9)=x(8)
else
 t(9)=x(7)
end if
if(muse(mp1,np1).eq.1 .or. ksuse(mp1,np1).eq.1) t(9)=g
t(8)=x(9)
if(muse(mm1,np1).eq.1 .and. muse(m,np1).eq.1 .and. muse(mp1,np1).eq.1) t(8)=g
if(x(10).gt.x(11)) then
 t(7)=x(11)
else
 t(7)=x(10)
end if
if(muse(mm1,np1).eq.1 .or. ksuse(mm1,np1).eq.1) t(7)=g
 t(4)=x(12)
if(muse(mm1,nm1).eq.1 .and. muse(mm1,n).eq.1 .and. muse(mm1,np1).eq.1) t(4)=g
tmin=1.e10
do l=1,9
 if(v(l).eq.smin .and. tmin.ge.t(l)) then
  tmin=t(l)
  minl=l
 end if
enddo
call onmin(minl,m,n,mm,nn)
return
end

! 最急勾配を探す( muse, ksuse がゼロでなければ検索対象から外す )
! 現在の問題点 行き先の標高が-9999の場合でも検索対象として残っている
subroutine calgrd(m,n,u,muse,mm1,ksuse, mp1,nm1,np1,s2,v,yschi,minn,smin,mmm,nnn,g,mds)
implicit integer(f)
real u,uu
dimension u(-2:mmm+3,-2:nnn+3),muse(-2:mmm+3,-2:nnn+3),ksuse(-2:mmm+3,-2:nnn+3)
dimension v(9)
!fnno(ml,nl,mds)=(nl-1)*mds+ml
!fnlocx(no,mds)=mod((no-1),mds)+1
!fnlocy(no,mds)=(no-1)/mds+1
uu=u(m,n)
if(muse(mm1,nm1).eq.0 .and. ksuse(mm1,nm1).eq.0) then
 v(1)=(u(mm1,nm1)-uu)/s2
else
 v(1)=g
end if
if(muse(m,nm1).eq.0) then
 v(2)=real(u(m,n-1)-uu)
else
 v(2)=g
end if
if(muse(mp1,nm1).eq.0 .and. ksuse(mp1,nm1).eq.0) then
 v(3)=(u(m+1,n-1)-uu)/s2
else
 v(3)=g
end if
if(muse(mm1,n).eq.0) then
 v(4)=(u(m-1,n)-uu)/yschi
else
 v(4)=g
end if
v(5)=g
if(muse(mp1,n).eq.0) then
 v(6)=(u(m+1,n)-uu)/yschi
else
 v(6)=g
end if
if(muse(mm1,np1).eq.0 .and. ksuse(mm1,np1).eq.0) then
 v(7)=(u(m-1,n+1)-uu)/s2
else
 v(7)=g
end if
if(muse(m,np1).eq.0) then
 v(8)=real(u(m,n+1)-uu)
else
 v(8)=g
end if
if(muse(mp1,np1).eq.0 .and. ksuse(mp1,np1).eq.0) then
 v(9)=(u(m+1,n+1)-uu)/s2
else
 v(9)=g
end if
smin=g
do k = 1, 9
 if(v(k).eq.g) cycle
 if(smin.lt.v(k)) cycle
 smin=v(k)
 minn=k
enddo
return
end

subroutine pit(m,n,nsuc,muse,mp1,mm1,np1,nm1,mmm,nnn,mds,gostrt,ksuse)
implicit integer(f)
integer gostrt
dimension nsuc(1:mmm*nnn),muse(-2:mmm+3,-2:nnn+3),ksuse(-2:mmm+3,-2:nnn+3)
!fnno(lm,ln,mds)=(ln-1)*mds+lm
!fnlocx(no,mds)=mod((no-1),mds)+1
!fnlocy(no,mds)=(no-1)/mds+1
!*  +++++++  pit no kisei  +++++++++l
!    no1
!no2  *  no4
!    no3
mp1=m+1
mm1=m-1
np1=n+1
nm1=n-1
no1=fnno(m,nm1,mds)
no2=fnno(mm1,n,mds)
no3=fnno(m,np1,mds)
no4=fnno(mp1,n,mds)
if(nsuc(no1).eq.no2 .or. nsuc(no2).eq.no1) ksuse(mm1,nm1)=1
if(nsuc(no1).eq.no4 .or. nsuc(no4).eq.no1) ksuse(mp1,nm1)=1
if(nsuc(no2).eq.no3 .or. nsuc(no3).eq.no2) ksuse(mm1,np1)=1
if(nsuc(no3).eq.no4 .or. nsuc(no4).eq.no3) ksuse(mp1,np1)=1
if((muse(mm1,nm1).eq.1 .or. ksuse(mm1,nm1).eq.1) .and. (muse(mp1,nm1).eq.1 .or. ksuse(mp1,nm1).eq.1) .and.(muse(mm1,np1).eq.1 .or. ksuse(mm1,np1).eq.1) &
    .and. (muse(mp1,np1).eq.1 .or. ksuse(mp1,np1).eq.1) .and. muse(m,nm1).eq.1 .and. muse(mm1,n).eq.1 .and. muse(mp1,n).eq.1 .and. muse(m,np1).eq.1) then
 ! (m, n)は周囲どこにも行き場がない
 write(*,*)'** 12    all muse = 1 **'
 write(*,*)'m=',m,'   n=',n
 gostrt=1
end if
return
end

subroutine back(m,n,mds,nda,nsuc,mmm,nnn)
implicit integer(f)
dimension nda(1:mmm*nnn),nsuc(1:mmm*nnn)
!fnno(mm,nn,mds)=(nn-1)*mds+mm
!fnlocx(no,mds)=mod((no-1),mds)+1
!fnlocy(no,mds)=(no-1)/mds+1
nda(fnno(m,n,mds))=nda(fnno(m,n,mds))-1
do i=1, mmm*nnn
 if(nsuc(i).eq.fnno(m,n,mds)) then
  nsuc(i)=-9999
  m=fnlocx(i,mds)
  n=fnlocy(i,mds)
  return
 end if
enddo
end

function fnno(m,n,mds)
integer fnno
fnno=(n-1)*mds+m
end function fnno

function fnlocx(no,mds)
integer fnlocx
fnlocx=mod((no-1),mds)+1
end function fnlocx

function fnlocy(no,mds)
integer fnlocy
fnlocy=(no-1)/mds+1
end function fnlocy

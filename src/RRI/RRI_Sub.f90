! RRI_Sub

! river index setting
subroutine riv_idx_setting
    use globals
    use dam_mod    ! only:dam_switch, dam_num, damflg,dam_ix, dam_iy,dam_loc,dam_qin !This line is added for RSR model 20240724
    implicit none
    integer i, j, ii, jj, k, kk, kkk, n, nn, f, l, ll, n_count ! add v1.4(k, kk) This line is added for RSR model 20240724
    real(8) distance

    riv_count = 0
    do i = 1, ny
        do j = 1, nx
            if (domain(i, j) .ne. 0 .and. riv(i, j) .eq. 1) riv_count = riv_count + 1
        end do
    end do

    allocate (riv_idx2i(riv_count), riv_idx2j(riv_count))   ! This line is added for RSR model 20240724
    allocate (down_riv_idx(riv_count), domain_riv_idx(riv_count))
    allocate (width_idx(riv_count), depth_idx(riv_count))
    allocate (height_idx(riv_count), area_ratio_idx(riv_count))
    allocate (zb_riv_idx(riv_count), dis_riv_idx(riv_count))
    allocate (dif_riv_idx(riv_count))
    allocate (sec_map_idx(riv_count), len_riv_idx(riv_count)) ! add v1.4
    allocate (damflg(riv_count))  ! This line is added for RSR model 20240724

    riv_count = 0
    riv_ij2idx(:, :) = 0
    do i = 1, ny
        do j = 1, nx

            if (domain(i, j) .eq. 0 .or. riv(i, j) .ne. 1) cycle
            ! domain(i, j) = 1 or 2 and riv(i, j) = 1
            riv_count = riv_count + 1

            riv_idx2i(riv_count) = i
            riv_idx2j(riv_count) = j
            domain_riv_idx(riv_count) = domain(i, j)
            ! added to avoid extreme small river channel 20250716
            if(chan_capa_switch==1)then
            if(width(i,j).le.min_wid) width(i,j) = min_wid   
            if(depth(i,j).le.min_dep) depth(i,j) = min_dep 
            endif
            width_idx(riv_count) = width(i, j)
            depth_idx(riv_count) = depth(i, j)
            height_idx(riv_count) = height(i, j)
            area_ratio_idx(riv_count) = area_ratio(i, j)
            zb_riv_idx(riv_count) = zb_riv(i, j)
            riv_ij2idx(i, j) = riv_count
            dif_riv_idx(riv_count) = dif(land(i, j))
            sec_map_idx(riv_count) = sec_map(i, j)
            len_riv_idx(riv_count) = len_riv(i, j) ! add v1.4

        end do
    end do

!---added for RSR model 20240724
    damflg(:) = 0
    do f = 1, dam_num
        !for dam
        if(riv(dam_iy(f),dam_ix(f))==0)then
             write(*,*) "Error: Please set dam location on a river cell" 
             stop
        endif
        dam_loc(f) = riv_ij2idx(dam_iy(f), dam_ix(f))
        damflg(dam_loc(f)) = f
    end do
!---until here

! search for downstream gridcell (down_idx)
    riv_count = 0
    do i = 1, ny
        do j = 1, nx

            if (domain(i, j) .eq. 0 .or. riv(i, j) .ne. 1) cycle
            ! domain(i, j) = 1 or 2 and riv(i, j) = 1
            riv_count = riv_count + 1

            ! right
            if (dir(i, j) .eq. 1) then
                ii = i
                jj = j + 1
                distance = dx
                ! right down
            elseif (dir(i, j) .eq. 2) then
                ii = i + 1
                jj = j + 1
                distance = sqrt(dx*dx + dy*dy)
                ! down
            elseif (dir(i, j) .eq. 4) then
                ii = i + 1
                jj = j
                distance = dy
                ! left down
            elseif (dir(i, j) .eq. 8) then
                ii = i + 1
                jj = j - 1
                distance = sqrt(dx*dx + dy*dy)
                ! left
            elseif (dir(i, j) .eq. 16) then
                ii = i
                jj = j - 1
                distance = dx
                ! left up
            elseif (dir(i, j) .eq. 32) then
                ii = i - 1
                jj = j - 1
                distance = sqrt(dx*dx + dy*dy)
                ! up
            elseif (dir(i, j) .eq. 64) then
                ii = i - 1
                jj = j
                distance = dy
                ! right up
            elseif (dir(i, j) .eq. 128) then
                ii = i - 1
                jj = j + 1
                distance = sqrt(dx*dx + dy*dy)
            elseif (dir(i, j) .eq. 0) then
                ii = i
                jj = j
            else
             !   write (*, *) "dir(i, j) is error (", i, j, ")", dir(i, j)
             write (*, *) "dir(i, j) is error (", j, ny+1-i, ")", dir(i, j) !convert i,j in IRIC GUI 20250505
                stop
            end if

            ! If the downstream cell is outside the domain, set domain(i, j) = 2
            if (ii .lt. 1 .or. ii .gt. ny .or. jj .lt. 1 .or. jj .gt. nx) then
                domain(i, j) = 2
                dir(i, j) = 0
                ii = i
                jj = j
            end if
            if (domain(ii, jj) .eq. 0) then
                domain(i, j) = 2
                dir(i, j) = 0
                ii = i
                jj = j
            end if

            if (riv(ii, jj) .eq. 0) then
 !               write (*, *) "riv(ii, jj) should be 1 (", i, j, ")", "(", ii, jj, ")"
                write (*, *) "riv(ii, jj) should be 1 (", j, ny+1-i, ")", "(", jj, ny+1-ii, ")" !convert i,j in IRIC GUI 20250505
                stop
            end if

            dis_riv_idx(riv_count) = distance
            down_riv_idx(riv_count) = riv_ij2idx(ii, jj)

        end do
    end do

! add v1.4
    if (sec_length_switch .eq. 1) then
        do k = 1, riv_count
            kk = down_riv_idx(k)
            dis_riv_idx(k) = (len_riv_idx(k) + len_riv_idx(kk))/2.d0
        end do
    end if

end subroutine riv_idx_setting

! 2D -> 1D (ij2idx)
subroutine sub_riv_ij2idx(a, a_idx)
    use globals
    implicit none

    real(8) a(ny, nx), a_idx(riv_count)
    integer k

    do k = 1, riv_count
        a_idx(k) = a(riv_idx2i(k), riv_idx2j(k))
    end do

end subroutine sub_riv_ij2idx

! 1D -> 2D (idx2ij)
subroutine sub_riv_idx2ij(a_idx, a)
    use globals
    implicit none

    real(8) a_idx(riv_count), a(ny, nx)
    integer k

    a(:, :) = 0.d0
    do k = 1, riv_count
        a(riv_idx2i(k), riv_idx2j(k)) = a_idx(k)
    end do

end subroutine sub_riv_idx2ij

! slope index setting
subroutine slo_idx_setting
    use globals
    use sediment_mod !this line is added for RSR 20240724
    implicit none
    integer i, j, ii, jj, k, kk, l, n_count !this line is modified for RSR 20240724
    real(8) distance, len, l1, l2, l3
    real(8) l1_kin, l2_kin, l3_kin

    slo_count = 0
    do i = 1, ny
        do j = 1, nx
            if (domain(i, j) .ne. 0) slo_count = slo_count + 1
        end do
    end do

    allocate (slo_idx2i(slo_count), slo_idx2j(slo_count), slo_ij2idx(ny, nx))
    allocate (down_slo_idx(i4, slo_count), domain_slo_idx(slo_count))
    allocate (zb_slo_idx(slo_count), dis_slo_idx(i4, slo_count), len_slo_idx(i4, slo_count), acc_slo_idx(slo_count))
    allocate (down_slo_1d_idx(slo_count), dis_slo_1d_idx(slo_count), len_slo_1d_idx(slo_count))
    allocate (up_slo_gather_count(slo_count), up_slo_gather_src(up_slo_gather_max, slo_count), &
              up_slo_gather_dir(up_slo_gather_max, slo_count))
    allocate (land_idx(slo_count))

    allocate (dif_slo_idx(slo_count))
    allocate (ns_slo_idx(slo_count))
    allocate (soildepth_idx(slo_count))
    allocate (gammaa_idx(slo_count))

    allocate (ksv_idx(slo_count), faif_idx(slo_count), infilt_limit_idx(slo_count),min_wc4latflow_idx(slo_count)) !added 20250314
    allocate (ka_idx(slo_count), gammam_idx(slo_count), beta_idx(slo_count), da_idx(slo_count), dm_idx(slo_count))
    allocate (ksg_idx(slo_count), kgv_idx(slo_count), gammag_idx(slo_count), tg_idx(slo_count), fpg_idx(slo_count), init_cond_gw_idx(slo_count) )!this line is modified for RSR 20240724
    allocate ( slo_riv_idx(slo_count), up_slo_idx(slo_count,8) )  !this line is modified for RSR 20240724

    slo_count = 0
    slo_ij2idx(:, :) = 0
    do i = 1, ny
        do j = 1, nx
            if (domain(i, j) .eq. 0) cycle
            slo_count = slo_count + 1
            slo_idx2i(slo_count) = i
            slo_idx2j(slo_count) = j
            domain_slo_idx(slo_count) = domain(i, j)
            zb_slo_idx(slo_count) = zb(i, j)
            acc_slo_idx(slo_count) = acc(i, j)
            slo_ij2idx(i, j) = slo_count
            land_idx(slo_count) = land(i, j)

            dif_slo_idx(slo_count) = dif(land(i, j))
            ns_slo_idx(slo_count) = ns_slope(land(i, j))
            soildepth_idx(slo_count) = soildepth(land(i, j))
            gammaa_idx(slo_count) = gammaa(land(i, j))
 
            ksv_idx(slo_count) = ksv(land(i, j))
            faif_idx(slo_count) = faif(land(i, j))
            infilt_limit_idx(slo_count) = infilt_limit(land(i, j))!added 20250314
            ka_idx(slo_count) = ka(land(i, j))
            min_wc4latflow_idx(slo_count) = min_wc4latflow(land(i,j))
            gammam_idx(slo_count) = gammam(land(i, j))
            beta_idx(slo_count) = beta(land(i, j))
            da_idx(slo_count) = da(land(i, j))
            dm_idx(slo_count) = dm(land(i, j))
            ksg_idx(slo_count) = ksg(land(i, j))
            kgv_idx(slo_count) = kgv(land(i, j)) !this line is modified for RSR 20240724
            gammag_idx(slo_count) = gammag(land(i, j))
            tg_idx(slo_count) = tg(land(i, j))   !this line is modified for RSR 20240724
!            kg0_idx(slo_count) = kg0(land(i, j))  !20240724
!            fpg_idx(slo_count) = fpg(land(i, j))  !20240724
!            rgl_idx(slo_count) = rgl(land(i, j))  !20240724
            if(riv(i,j) == 0)then    !this if is modified for RSR 20240724
                slo_riv_idx(slo_count) = 0
            else
                slo_riv_idx(slo_count) = riv_ij2idx(i,j)  !added 2021/11/5 connect slope number to river number
            end if

        end do
    end do
    
    if (eight_dir .eq. 1) then
        ! Hromadka etal, JAIH2006 (USE THIS AS A DEFAULT)
        ! 8-direction
        lmax = 4
        l1 = dy/2.d0
        l2 = dx/2.d0
        l3 = sqrt(dx**2.d0 + dy**2.d0)/4.d0
    else if (eight_dir .eq. 0) then
        ! 4-direction
        lmax = 2
        l1 = dy
        l2 = dx
        l3 = 0.d0
    else
        stop "error: eight_dir should be 0 or 1."
    end if

! search for downstream gridcell (down_slo_idx)
    slo_count = 0
    down_slo_idx(:, :) = -1

    do i = 1, ny
        do j = 1, nx

            if (domain(i, j) .eq. 0) cycle
            ! domain(i, j) = 1 or 2
            slo_count = slo_count + 1

            ! 8-direction: lmax = 4, 4-direction: lmax = 2
            do l = 1, lmax ! (1: rightC2: down, 3: right down, 4: left down)

                if (l .eq. 1) then
                    ii = i
                    jj = j + 1
                    distance = dx
                    len = l1
                elseif (l .eq. 2) then
                    ii = i + 1
                    jj = j
                    distance = dy
                    len = l2
                elseif (l .eq. 3) then
                    ii = i + 1
                    jj = j + 1
                    distance = sqrt(dx*dx + dy*dy)
                    len = l3
                else
                    ii = i + 1
                    jj = j - 1
                    distance = sqrt(dx*dx + dy*dy)
                    len = l3
                end if

                if (ii .gt. ny) cycle
                if (jj .gt. nx) cycle
                if (ii .lt. 1) cycle
                if (jj .lt. 1) cycle
                if (domain(ii, jj) .eq. 0) cycle

                down_slo_idx(l, slo_count) = slo_ij2idx(ii, jj)
                dis_slo_idx(l, slo_count) = distance
                len_slo_idx(l, slo_count) = len

            end do
        end do
    end do

! search for downstream gridcell (down_slo_1d_idx) (used only for kinematic with 1-direction)
    slo_count = 0

    down_slo_1d_idx(:) = -1
    dis_slo_1d_idx(:) = l1
    len_slo_1d_idx(:) = dx

    l1_kin = dy
    l2_kin = dx
    l3_kin = dx*dy/sqrt(dx**2.d0 + dy**2.d0)

    do i = 1, ny
        do j = 1, nx

            if (domain(i, j) .eq. 0) cycle
            ! domain(i, j) = 1 or 2
            slo_count = slo_count + 1

            ! right
            if (dir(i, j) .eq. 1) then
                ii = i
                jj = j + 1
                distance = dx
                len = l1_kin
                ! right down
            elseif (dir(i, j) .eq. 2) then
                ii = i + 1
                jj = j + 1
                distance = sqrt(dx*dx + dy*dy)
                len = l3_kin
                ! down
            elseif (dir(i, j) .eq. 4) then
                ii = i + 1
                jj = j
                distance = dy
                len = l2_kin
                ! left down
            elseif (dir(i, j) .eq. 8) then
                ii = i + 1
                jj = j - 1
                distance = sqrt(dx*dx + dy*dy)
                len = l3_kin
                ! left
            elseif (dir(i, j) .eq. 16) then
                ii = i
                jj = j - 1
                distance = dx
                len = l1_kin
                ! left up
            elseif (dir(i, j) .eq. 32) then
                ii = i - 1
                jj = j - 1
                distance = sqrt(dx*dx + dy*dy)
                len = l3_kin
                ! up
            elseif (dir(i, j) .eq. 64) then
                ii = i - 1
                jj = j
                distance = dy
                len = l2_kin
                ! right up
            elseif (dir(i, j) .eq. 128) then
                ii = i - 1
                jj = j + 1
                distance = sqrt(dx*dx + dy*dy)
                len = l3_kin
            elseif (dir(i, j) .eq. 0) then
                ii = i
                jj = j
                distance = dx
                len = l1_kin
            else
              !  write (*, *) "dir(i, j) is error (", i, j, ")", dir(i, j)
             write (*, *) "dir(i, j) is error (", j, ny+1-i, ")", dir(i, j) !convert i,j in IRIC GUI 20250505                
                stop
            end if

            if (ii .gt. ny) cycle
            if (jj .gt. nx) cycle
            if (ii .lt. 1) cycle
            if (jj .lt. 1) cycle
            if (domain(ii, jj) .eq. 0) cycle

            down_slo_1d_idx(slo_count) = slo_ij2idx(ii, jj)
            dis_slo_1d_idx(slo_count) = distance
            len_slo_1d_idx(slo_count) = len

        end do
    end do

! Build an upstream list for gather-style slope conservation updates.
! Kinematic cells use down_slo_1d_idx and store their flux in direction 1.
    up_slo_gather_count(:) = 0
    up_slo_gather_src(:, :) = 0
    up_slo_gather_dir(:, :) = 0

    do k = 1, slo_count
        if (dif_slo_idx(k) .eq. 0) then
            kk = down_slo_1d_idx(k)
            if (kk .eq. -1) cycle
            n_count = up_slo_gather_count(kk) + 1
            if (n_count .gt. up_slo_gather_max) stop "error: up_slo_gather_max is too small."
            up_slo_gather_count(kk) = n_count
            up_slo_gather_src(n_count, kk) = k
            up_slo_gather_dir(n_count, kk) = 1
        else
            do l = 1, lmax
                kk = down_slo_idx(l, k)
                if (kk .eq. -1) cycle
                n_count = up_slo_gather_count(kk) + 1
                if (n_count .gt. up_slo_gather_max) stop "error: up_slo_gather_max is too small."
                up_slo_gather_count(kk) = n_count
                up_slo_gather_src(n_count, kk) = k
                up_slo_gather_dir(n_count, kk) = l
            end do
        end if
    end do

!---------modified for RSR model : from here to the end of this subroutine
! searching for upstream slope cell; added by Qin 20220803
if(slope_ero_switch==1)then    !this loop sometimes not work correctly
slo_count = 0
do i = 1, ny
 do j =1, nx
   n_count = 0
  if(domain(i,j).eq.0) cycle
  slo_count = slo_count + 1

  ! right
   ii = i
   jj = j + 1

!pause'right'

  if( dir(ii,jj).eq.16 ) then
  if(domain(ii,jj).eq.0) goto 10
   n_count = n_count + 1 
   up_slo_idx(slo_count,n_count) = slo_ij2idx(ii, jj)
!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'right'
endif

10 continue

  ! right down
   ii = i + 1
   jj = j + 1

!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'right down'

  if( dir(ii,jj).eq.32 ) then
if(domain(ii,jj).eq.0 ) goto 20
   n_count = n_count + 1 
   up_slo_idx(slo_count,n_count) =slo_ij2idx(ii, jj)
!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'right down'
  endif

20 continue

  ! down
   ii = i + 1
   jj = j

!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'down'

  if( dir(ii,jj).eq.64 ) then
if(domain(ii,jj).eq.0 .or. riv(ii,jj).ne.1) goto 30
   n_count = n_count + 1 
   up_slo_idx(slo_count,n_count) = slo_ij2idx(ii, jj)
!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'down'
  endif

30 continue

  ! left down
   ii = i + 1
   jj = j - 1

!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'left down'

  if( dir(ii,jj).eq.128 ) then
if(domain(ii,jj).eq.0 .or. riv(ii,jj).ne.1) goto 40
   n_count = n_count + 1 
   up_slo_idx(slo_count,n_count) = slo_ij2idx(ii, jj)
!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'left down'
  endif

40 continue

  ! left
   ii = i
   jj = j - 1

!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'left'

  if( dir(ii,jj).eq.1 ) then
if(domain(ii,jj).eq.0 .or. riv(ii,jj).ne.1) goto 50
   n_count = n_count + 1 
   up_slo_idx(slo_count,n_count) = slo_ij2idx(ii, jj)
!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'left'
  endif

50 continue

  ! left up
   ii = i - 1
   jj = j - 1

!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'left up'


  if( dir(ii,jj).eq.2 ) then
if(domain(ii,jj).eq.0 .or. riv(ii,jj).ne.1) goto 60
   n_count = n_count + 1 
   up_slo_idx(slo_count,n_count) = slo_ij2idx(ii, jj)
!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'left up'
  endif

60 continue

  ! up
   ii = i - 1
   jj = j

!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'up'

  if( dir(ii,jj).eq.4 ) then
if(domain(ii,jj).eq.0 .or. riv(ii,jj).ne.1) goto 70
   n_count = n_count + 1 
   up_slo_idx(slo_count,n_count) = slo_ij2idx(ii, jj)
!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'up'
  endif

70 continue

  ! right up
   ii = i - 1
   jj = j + 1

!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'right up'

  if( dir(ii,jj).eq.8 ) then
if(domain(ii,jj).eq.0 .or. riv(ii,jj).ne.1) goto 80
   n_count = n_count + 1 
   up_slo_idx(slo_count,n_count) = slo_ij2idx(ii, jj)
!write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
!pause'right up'
  endif

80 continue

enddo
enddo


end if
write(*,*) "Finished slope idx setting"
end subroutine slo_idx_setting

! 2D -> 1D (ij2idx)
subroutine sub_slo_ij2idx(a, a_idx)
    use globals
    implicit none

    real(8) a(ny, nx), a_idx(slo_count)
    integer k

    do k = 1, slo_count
        a_idx(k) = a(slo_idx2i(k), slo_idx2j(k))
    end do

end subroutine sub_slo_ij2idx

! 1D -> 2D (idx2ij)
subroutine sub_slo_idx2ij(a_idx, a)
    use globals
    implicit none

    real(8) a_idx(slo_count), a(ny, nx)
    integer k

    a(:, :) = 0.d0
    do k = 1, slo_count
        a(slo_idx2i(k), slo_idx2j(k)) = a_idx(k)
    end do

end subroutine sub_slo_idx2ij

! 2D -> 1D (ij2idx)
subroutine sub_slo_idx2ij4(a_idx, a)
    use globals   
    implicit none

    real(8) a_idx(i4, slo_count), a(i4, ny, nx)
    integer i, k

    a(:, :, :) = 0.d0
    do i = 1, i4
        do k = 1, slo_count
            a(i, slo_idx2i(k), slo_idx2j(k)) = a_idx(i, k)
        end do
    end do

end subroutine sub_slo_idx2ij4

! storage calculation
subroutine storage_calc(hs, hr, hg, ss, sr, si, sg)
    use globals
    implicit none

    real(8) hs(ny, nx), hr(ny, nx), hg(ny, nx)
    real(8) ss, sr, si, sg
    real(8) vr_temp ! add v1.4
    integer i, j, k

    ss = 0.d0
    sr = 0.d0
    si = 0.d0
    sg = 0.d0
    do i = 1, ny
        do j = 1, nx
            if (domain(i, j) .eq. 0) cycle
            ss = ss + hs(i, j)*area
            !if(riv_thresh.ge.0 .and. riv(i,j).eq.1) sr = sr + hr(i,j) * area * area_ratio(i,j)
            ! modified v1.4
            if (riv_thresh .ge. 0 .and. riv(i, j) .eq. 1) then
                call hr2vr(hr(i, j), riv_ij2idx(i, j), vr_temp)
                sr = sr + vr_temp
            end if
            si = si + gampt_ff(i, j)*area
            sg = sg - hg(i, j)*gammag_idx(slo_ij2idx(i, j))*area ! storage deficit
        end do
    end do

end subroutine storage_calc

! numbers to characters
subroutine int2char(num, cwrk)
    implicit none

    integer num, j
    character*6 cwrk
    cwrk = ' '
    write (cwrk, '( I6 )') num
    do j = 1, 6
        if (cwrk(j:j) .eq. ' ') cwrk(j:j) = '0'
    end do
end subroutine int2char

! Hubeny_sub.f90
subroutine hubeny_sub(x1_deg, y1_deg, x2_deg, y2_deg, d)
    implicit none

    real(8) x1_deg, y1_deg, x2_deg, y2_deg
    real(8) x1, y1, x2, y2
    real(8) pi, dx, dy, mu, a, b, e, W, N, M, d

    pi = 3.1415926535897d0

    x1 = x1_deg*pi/180.d0
    y1 = y1_deg*pi/180.d0
    x2 = x2_deg*pi/180.d0
    y2 = y2_deg*pi/180.d0

    dy = y1 - y2
    dx = x1 - x2
    mu = (y1 + y2)/2.

    a = 6378137.000d0 ! Semi-Major Axis
    b = 6356752.314d0 ! Semi-Minor Axis

    e = sqrt((a**2.d0 - b**2.d0)/(a**2.d0))

    W = sqrt(1.-e**2.d0*(sin(mu))**2.d0)

    N = a/W

    M = a*(1.-e**2.d0)/W**(3.d0)

    d = sqrt((dy*M)**2.d0 + (dx*N*cos(mu))**2.d0)

end subroutine

! reading gis data (integer)
subroutine read_gis_int(fi, gis_data)
    use globals
    implicit none

    integer gis_data(ny, nx)
    character*256 fi

    integer itemp, i, j
    real(8) rtemp
    character*256 ctemp

    open (10, file=fi, status="old")

    read (10, *) ctemp, itemp
    if (nx .ne. itemp) stop "error in gis input data (int)"
    read (10, *) ctemp, itemp
    if (ny .ne. itemp) stop "error in gis input data (int)"
    read (10, *) ctemp, rtemp
    if (abs(xllcorner - rtemp) .gt. 0.01) stop "error in gis input data (int)"
    read (10, *) ctemp, rtemp
    if (abs(yllcorner - rtemp) .gt. 0.01) stop "error in gis input data (int)"
    read (10, *) ctemp, rtemp
    if (abs(cellsize - rtemp) .gt. 0.01) stop "error in gis input data (int)"
    read (10, *) ctemp, rtemp ! nodata

    do i = 1, ny
        read (10, *) (gis_data(i, j), j=1, nx)
    end do
    close (10)
end subroutine read_gis_int

! reading gis data (real)
subroutine read_gis_real(fi, gis_data)
    use globals
    implicit none

    real(8) gis_data(ny, nx)
    character*256 fi

    integer itemp, i, j
    real(8) rtemp
    character*256 ctemp

    open (10, file=fi, status="old")

    read (10, *) ctemp, itemp
    if (nx .ne. itemp) stop "error in gis input data (real)"
    read (10, *) ctemp, itemp
    if (ny .ne. itemp) stop "error in gis input data (real)"
    read (10, *) ctemp, rtemp
    if (abs(xllcorner - rtemp) .gt. 0.01) stop "error in gis input data (real)"
    read (10, *) ctemp, rtemp
    if (abs(yllcorner - rtemp) .gt. 0.01) stop "error in gis input data (real)"
    read (10, *) ctemp, rtemp
    if (abs(cellsize - rtemp) .gt. 0.01) stop "error in gis input data (real)"
    read (10, *) ctemp, rtemp ! nodata

    do i = 1, ny
        read (10, *) (gis_data(i, j), j=1, nx)
    end do
    close (10)
end subroutine read_gis_real

! RRI_Slope

! variable definition (slope)

subroutine funcs(hs_idx, qp_t_idx, fs_idx, qs_idx)
    use globals
    implicit none

    real(8) hs_idx(slo_count), qp_t_idx(slo_count), fs_idx(slo_count)
    real(8) qs_idx(i4, slo_count)

    integer k, n, ku, lu, itemp, jtemp

    fs_idx(:) = 0.d0
    qs_idx(:, :) = 0.d0

! Keep one OpenMP team across the slope flux calculation and the local
! right-hand-side update.  The conservation update gathers upstream fluxes,
! so each thread writes only its own fs_idx(k).
!$omp parallel private(k, itemp, jtemp)
    call qs_calc_do(hs_idx, qs_idx)

! boundary condition for slope (discharge boundary)
!$omp single
    if (bound_slo_disc_switch .ge. 1) then
        !itemp = time / dt_bound_slo + 1
        itemp = -1
        do jtemp = 1, tt_max_bound_slo_disc
            if (t_bound_slo_disc(jtemp - 1) .lt. (time + ddt) .and. (time + ddt) .le. t_bound_slo_disc(jtemp)) itemp = jtemp
        end do
        do k = 1, slo_count
            if (bound_slo_disc_idx(itemp, k) .le. -100.0) cycle ! not boundary
            ! right
            if (dir(slo_idx2i(k), slo_idx2j(k)) .eq. 1) then
                qs_idx(1, k) = bound_slo_disc_idx(itemp, k)/area
                ! right down
            elseif (dir(slo_idx2i(k), slo_idx2j(k)) .eq. 2) then
                qs_idx(3, k) = bound_slo_disc_idx(itemp, k)/area
                ! down
            elseif (dir(slo_idx2i(k), slo_idx2j(k)) .eq. 4) then
                qs_idx(2, k) = bound_slo_disc_idx(itemp, k)/area
                ! left down
            elseif (dir(slo_idx2i(k), slo_idx2j(k)) .eq. 8) then
                qs_idx(4, k) = bound_slo_disc_idx(itemp, k)/area
                ! left
            elseif (dir(slo_idx2i(k), slo_idx2j(k)) .eq. 16) then
                qs_idx(1, k) = -bound_slo_disc_idx(itemp, k)/area
                ! left up
            elseif (dir(slo_idx2i(k), slo_idx2j(k)) .eq. 32) then
                qs_idx(3, k) = -bound_slo_disc_idx(itemp, k)/area
                ! up
            elseif (dir(slo_idx2i(k), slo_idx2j(k)) .eq. 64) then
                qs_idx(2, k) = -bound_slo_disc_idx(itemp, k)/area
                ! right up
            elseif (dir(slo_idx2i(k), slo_idx2j(k)) .eq. 128) then
                qs_idx(4, k) = -bound_slo_disc_idx(itemp, k)/area
            end if
        end do
    end if
!$omp end single

! qs_idx > 0 --> discharge flowing out from a cell

!$omp do private(n, ku, lu)
    do k = 1, slo_count
        fs_idx(k) = qp_t_idx(k) - (qs_idx(1, k) + qs_idx(2, k) + qs_idx(3, k) + qs_idx(4, k))
        do n = 1, up_slo_gather_count(k)
            ku = up_slo_gather_src(n, k)
            lu = up_slo_gather_dir(n, k)
            fs_idx(k) = fs_idx(k) + qs_idx(lu, ku)
        end do
    end do
!$omp end do
!$omp end parallel

end subroutine funcs

! lateral discharge (slope)
      subroutine qs_calc(hs_idx, qs_idx) !20231226modified 
         use globals
         implicit none

         real(8) hs_idx(slo_count)
         real(8) qs_idx(i4, slo_count)

         qs_idx = 0.d0

!$omp parallel
         call qs_calc_do(hs_idx, qs_idx)
!$omp end parallel

      end subroutine qs_calc

! lateral discharge workshare (slope)
      subroutine qs_calc_do(hs_idx, qs_idx)
         use globals
         implicit none

         real(8) hs_idx(slo_count)
         real(8) qs_idx(i4, slo_count), q
         integer k, kk, l
         real(8) zb_p, hs_p, ns_p, ka_p, da_p, dm_p, b_p
         real(8) zb_n, hs_n, ns_n, ka_n, da_n, dm_n, b_n
         real(8) dh, distance
         real(8) lev_p, lev_n
         real(8) len, hw
         real(8) dhdx, dhdy
         integer dif_p, dif_n
         !real(8) emb

!$omp do private(kk,zb_p,hs_p,ns_p,ka_p,da_p,dm_p,b_p,dif_p,l,distance,len, &
!$omp                     zb_n,hs_n,ns_n,ka_n,da_n,dm_n,b_n,dif_n,lev_p,lev_n,dh,hw,dhdx,dhdy,q) !modified 20231226
!added q into private 20231026
         do k = 1, slo_count

             zb_p = zb_slo_idx(k)
             hs_p = hs_idx(k)
             ns_p = ns_slo_idx(k)
             ka_p = ka_idx(k)
             da_p = da_idx(k)
             dm_p = dm_idx(k)
             b_p  = beta_idx(k)
             dif_p = dif_slo_idx(k)

             dhdx = 0.d0
             dhdy = 0.d0

 ! this loop was added by Harada to compute resultant
             if( eight_dir .eq. 0)then
                do l = 1, lmax ! (1: right�C2: down, 3: right down, 4: left down)
                    if( dif_p .eq. 0 .and. l .eq. 2 ) exit ! kinematic -> 1-direction
                    kk = down_slo_idx(l, k)
                    if( dif_p .eq. 0 ) kk = down_slo_1d_idx(k)
                    if( kk .eq. -1 ) cycle
   
                    distance = dis_slo_idx(l, k)
                    len = len_slo_idx(l, k)
                    if( dif_p .eq. 0 ) distance = dis_slo_1d_idx(k)
                    if( dif_p .eq. 0 ) len = len_slo_1d_idx(k)
   
     ! information of the destination cell
                    zb_n = zb_slo_idx(kk)
                    hs_n = hs_idx(kk)
                    ns_n = ns_slo_idx(kk)
                    ka_n = ka_idx(kk)
                    da_n = da_idx(kk)
                    dm_n = dm_idx(kk)
                    b_n = beta_idx(kk)
                    dif_n = dif_slo_idx(kk)
   
                    call h2lev(hs_p, k, lev_p)
                    call h2lev(hs_n, kk, lev_n)
   
     ! diffusion wave
                    dh = ((zb_p + lev_p) - (zb_n + lev_n)) / distance
     ! 1-direction : kinematic wave
                    if( dif_p .eq. 0 ) dh = max( (zb_p - zb_n) / distance, 0.001 )
     !-----added to compute resultant
                    !to avoid very small value
                    if(dh >= 0.) dh = max(dh , 0.00001)
                    if(dh < 0.) dh = min(dh , -0.00001)

                    if( l.eq.1 ) then
                        dhdx = dh
                    else     ! l=2
                        dhdy = dh
                    end if
                    
                end do
             end if

 ! 8-direction: lmax = 4, 4-direction: lmax = 2
             do l = 1, lmax ! (1: right�C2: down, 3: right down, 4: left down)
                 if( dif_p .eq. 0 .and. l .eq. 2 ) exit ! kinematic -> 1-direction
                 kk = down_slo_idx(l, k)
                 if( dif_p .eq. 0 ) kk = down_slo_1d_idx(k)
                 if( kk .eq. -1 ) cycle

                 distance = dis_slo_idx(l, k)
                 len = len_slo_idx(l, k)
                 if( dif_p .eq. 0 ) distance = dis_slo_1d_idx(k)
                 if( dif_p .eq. 0 ) len = len_slo_1d_idx(k)

  ! information of the destination cell
                 zb_n = zb_slo_idx(kk)
                 hs_n = hs_idx(kk)
                 ns_n = ns_slo_idx(kk)
                 ka_n = ka_idx(kk)
                 da_n = da_idx(kk)
                 dm_n = dm_idx(kk)
                 b_n = beta_idx(kk)
                 dif_n = dif_slo_idx(kk)

                 if( eight_dir .eq. 1)then

                    call h2lev(hs_p, k, lev_p)
                    call h2lev(hs_n, kk, lev_n)

  ! diffusion wave
                    dh = ((zb_p + lev_p) - (zb_n + lev_n)) / distance

  ! 1-direction : kinematic wave
                    if( dif_p .eq. 0 ) dh = max( (zb_p - zb_n) / distance, 0.001 )

  !-----modified by Harada 2021/10/22
                 else  !eight_dir = 0
                    if( l .eq. 1 ) then
                        dh = dhdx*abs(dhdx) / sqrt ( dhdx**2. + dhdy**2. )
                    else  ! l=2 
                        dh = dhdy*abs(dhdy) / sqrt ( dhdx**2. + dhdy**2. ) 
                    end if
                 end if
  !------
                    
  ! embankment
  !if(emb_switch.eq.1) then
  ! if(l.eq.1) emb = emb_r_idx(k)
  ! if(l.eq.2) emb = emb_b_idx(k)
  ! if(l.eq.3) emb = max( emb_r_idx(k), emb_b_idx(k) )
  ! if(l.eq.4) emb = max( emb_r_idx(kk), emb_b_idx(k) )
  !endif

  ! water coming in or going out?
                 if( dh .ge. 0.d0 ) then
   ! going out
                     hw = hs_p
   !if(emb .gt. 0.d0) hw = max(hs_p - emb, 0.d0)
                     if( zb_p .lt. zb_n ) hw = max(0.d0, zb_p + hs_p - zb_n)
                     call hq(ns_p, ka_p, da_p, dm_p, b_p, hw, dh, len, q,k) !modified 20250314 
                     qs_idx(l,k) = q
                 else
   ! coming in
                     hw = hs_n
   !if(emb .gt. 0.d0) hw = max(hs_n - emb, 0.d0)
                    ! dh = abs(dh) !modified 20250317
                     if( zb_n .lt. zb_p ) hw = max(0.d0, zb_n + hs_n - zb_p)
                     call hq(ns_n, ka_n, da_n, dm_n, b_n, hw, dh, len, q,k) !modified 20250314 
                     qs_idx(l,k) = -q
                 endif

             enddo
         enddo
!$omp end do


      end subroutine qs_calc_do

! water depth and discharge relationship
subroutine hq(ns_p, ka_p, da_p, dm_p, b_p, h, dh, len, q,k) !modified 20250314
    use globals
    implicit none

    real(8) ns_p, da_p, dm_p, ka_p, b_p, h, dh, len, q
    real(8) km, vm, va, al, m, min_hs
    integer k
    !real(8) min_hs_idx(slo_count)
    min_hs =0.d0
    if(ka_p>0.d0) min_hs = soildepth_idx(k)*gammaa_idx(k)*min_wc4latflow_idx(k) !modified 20250314
    !modified 20250317
    if (dh>0.d0 .and. h.le.min_hs)then 
        q=0.d0
    else   
        dh = dabs(dh)
        if (b_p .gt. 0.d0) then
            km = ka_p/b_p
        else
            km = 0.d0
        end if
        vm = km*dh

        if (da_p .gt. 0.d0) then
            va = ka_p*dh
        else
            va = 0.d0
        end if

        !if (dh .lt. 0) dh = 0.d0
        al = sqrt(dh)/ns_p
        m = 5.d0/3.d0

        if (h .lt. dm_p) then
            q = vm*dm_p*(h/dm_p)**b_p
        elseif (h .lt. da_p) then
            q = vm*dm_p + va*(h - dm_p)
        else
            q = vm*dm_p + va*(h - dm_p) + al*(h - da_p)**m
        end if
    endif

! discharge per unit area
! (q multiply by width and divide by area)
    q = q*len/area

! water depth limitter (1 mm)
! note: it can be set to zero
!if( h.le.0.001 ) q = 0.d0

end subroutine hq

! water depth (h) to actual water level (lev)
subroutine h2lev(h, k, lev)
    use globals
    implicit none

    integer k
    real(8) h, lev
    real(8) rho
    real(8) da_temp

    da_temp = soildepth_idx(k)*gammaa_idx(k)

    if (soildepth_idx(k) .eq. 0.d0) then
        lev = h
    elseif (h .ge. da_temp) then ! including da = 0
        lev = soildepth_idx(k) + (h - da_temp) ! surface water
    else
        if (soildepth_idx(k) .gt. 0.d0) rho = da_temp/soildepth_idx(k)
        lev = h/rho
    end if

end subroutine h2lev

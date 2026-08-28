!  RRI_Riv.f90

! variable definition (river)

      subroutine funcr( vr_idx, fr_idx, qr_idx )
         use globals
         use sediment_mod
         use dam_mod
         implicit none

         real(8) hr_idx(riv_count), vr_idx(riv_count), fr_idx(riv_count), qr_idx(riv_count)

         real(8) qr_sum_idx(riv_count), qr_div_idx(riv_count)
         integer i, k, kk, itemp, jtemp , l
         integer kk_div
	     integer t

         fr_idx(:) = 0.d0
         qr_idx(:) = 0.d0
         qr_sum_idx(:) = 0.d0
         qr_div_idx(:) = 0.d0

! add v1.4
         do k = 1, riv_count
             call vr2hr(vr_idx(k), k, hr_idx(k)) 
         enddo

! boundary condition for river (water depth boundary)
         if( bound_riv_wlev_switch .ge. 1 ) then
             itemp = -1
             do jtemp = 1, tt_max_bound_riv_wlev
                 if( t_bound_riv_wlev(jtemp-1) .lt. (time + ddt) .and. (time + ddt) .le. t_bound_riv_wlev(jtemp) ) itemp = jtemp
             enddo
             do k = 1, riv_count
                 if( bound_riv_wlev_idx(itemp, k) .le. -100.0 ) cycle ! not boundary
                 hr_idx(k) = bound_riv_wlev_idx(itemp, k)
                 call hr2vr(hr_idx(k), k, vr_idx(k)) ! add v1.4
             enddo
         endif

         call qr_calc(hr_idx, qr_idx)

! dam control
         if( dam_switch .eq. 1 ) then
             call dam_prepare(qr_idx) ! calculate inflow to dam
             do i = 1, dam_num
                 if( dam_volmax(i) .ge. 0.d0 ) then
   ! dam
                     call dam_operation(i)
                     qr_idx( dam_loc(i) ) = dam_qout(i)
                     !hr_idx( dam_loc(i) ) = 0.d0
                 else
   ! barrage
                     if( hr_idx( dam_loc(i) ) .le. dam_floodq(i) ) then
                         qr_idx( dam_loc(i) ) = 0.d0
                     else
             continue
                     endif
                 endif
             enddo
            call qr_calc(hr_idx, qr_idx) !modified 20240407
         endif
         
         if( div_switch > 0 ) call RRI_Div(qr_idx, hr_idx,qr_div_idx)

! boundary condition for river (discharge boundary)
         if( bound_riv_disc_switch .ge. 1 ) then
             do jtemp = 1, tt_max_bound_riv_disc
                 if( t_bound_riv_disc(jtemp-1) .lt. (time + ddt) .and. (time + ddt) .le. t_bound_riv_disc(jtemp) ) itemp = jtemp
             enddo
             do k = 1, riv_count
                 if( bound_riv_disc_idx(itemp, k) .le. -100.0 ) cycle ! not boundary
                 !qr_idx(k) = bound_riv_disc_idx(itemp, k) / area  ! qr_idx: discharge per unit area
                 qr_idx(k) = bound_riv_disc_idx(itemp, k) ! modified v1.4
  ! linear interpolation of the boundary condition
  !qr_idx(k) = bound_riv_disc_idx(itemp-1, k) * (t_bound_riv_disc(itemp) - (time + ddt)) &
  !           + bound_riv_disc_idx(itemp, k) * ((time + ddt) - t_bound_riv_disc(itemp-1))
  !qr_idx(k) = qr_idx(k) / (t_bound_riv_disc(itemp) - t_bound_riv_disc(itemp-1))
                 hr_idx(k) = 0.d0
                 vr_idx(k) = 0.d0 ! add v1.4
             enddo
         endif
         
! qr_sum > 0 --> discharge flowing out from a cell
         do k = 1, riv_count
 ! outflow from (k)
             qr_sum_idx(k) = qr_sum_idx(k) + qr_idx(k)
             kk = down_riv_idx(k)
             if(domain_riv_idx(kk).eq.0) cycle
 ! qr_sum minus (flowing into) discharge at the destination cell
             qr_sum_idx(kk) = qr_sum_idx(kk) - qr_idx(k)
         enddo

! diversion
         if( div_switch > 0 ) then
             do l = 1, div_id_max
 ! outflow from (k)
                 k = div_org_idx(l)
        
                 qr_sum_idx(k) = qr_sum_idx(k) + qr_div_idx(k)
                 kk = div_dest_idx(l)
                 if(domain_riv_idx(kk).eq.0) cycle
  ! qr_sum minus (flowing into) discharge at the destination cell
                 qr_sum_idx(kk) = qr_sum_idx(kk) - qr_div_idx(k) 
             enddo
         endif

! qr_sum divide by area_ratio
! (area_ratio = ratio of river area against total cell area) 
!fr_idx = -qr_sum_idx / area_ratio_idx
         fr_idx = - qr_sum_idx ! modified v1.4

      end subroutine funcr


! lateral discharge (river)
!---------------------------------------------------------------------------------------------
      subroutine qr_calc(hr_idx, qr_idx)
         use globals
         use sediment_mod
         use dam_mod!, only: damflg
         implicit none

         real(8) hr_idx(riv_count), qr_idx(riv_count), qr_div_idx(riv_count)

         integer k, kk
         real(8) zb_p, hr_p, hr_p_eff
         real(8) zb_n, hr_n, hr_n_eff
         real(8) dh, distance
         real(8) qr_temp, hw
         real(8) qr_raw, qr_capacity, hw_eff, hcap_k, hcap_kk, hcap_link
         real(8) block_fac, block_fac_k, block_fac_kk
         real(8) depth_ratio_k, depth_ratio_kk
         integer dif_p, dif_n

         qr_idx(:) = 0.d0
         qr_div_idx(:) = 0.d0
         if(allocated(qr_raw_idx)) qr_raw_idx(:) = 0.d0
         if(allocated(qr_capacity_limited_idx)) qr_capacity_limited_idx(:) = 0.d0
         if(allocated(qr_effective_depth_idx)) qr_effective_depth_idx(:) = 0.d0
         if(allocated(qr_capacity_depth_idx)) qr_capacity_depth_idx(:) = 0.d0
         if(allocated(qr_blockage_factor_idx)) qr_blockage_factor_idx(:) = 1.d0

!$omp parallel do private(kk,zb_p,hr_p,hr_p_eff,distance,zb_n,hr_n,hr_n_eff,dh,hw,qr_temp,qr_raw, &
!$omp& qr_capacity,hw_eff,hcap_k,hcap_kk,hcap_link,block_fac,block_fac_k,block_fac_kk, &
!$omp& depth_ratio_k,depth_ratio_kk,dif_p,dif_n)
! aaded dif_p,dif_n to private 20231026
         do k = 1, riv_count

             if(domain_riv_idx(k) .eq. 2) cycle

             zb_p = zb_riv_idx(k)
             hr_p = hr_idx(k)
             dif_p = dif_riv_idx(k)

             distance = dis_riv_idx(k)

 ! information of the destination cell
             kk = down_riv_idx(k)
             zb_n = zb_riv_idx(kk)
             hr_n = hr_idx(kk)
             dif_n = dif_riv_idx(kk)

             hcap_k = max(depth_idx(k) + height_idx(k), 0.d0)
             hcap_kk = max(depth_idx(kk) + height_idx(kk), 0.d0)
             hcap_link = min(hcap_k, hcap_kk)
             hr_p_eff = hr_p
             hr_n_eff = hr_n
             if(channel_capacity_qr_switch > 0)then
                 hr_p_eff = min(max(hr_p, 0.d0), hcap_k)
                 hr_n_eff = min(max(hr_n, 0.d0), hcap_kk)
             endif

! diffusion wave
             dh = ((zb_p + hr_p_eff) - (zb_n + hr_n_eff)) / distance ! diffussion

! kinematic wave
             if( dif_p .eq. 0 ) dh = max( (zb_p - zb_n) / distance, 0.001 )

! the destination cell is outlet (domain = 2)
             !if( domain_riv_idx(kk) .eq. 2 ) dh = 0.0   !check tentative 20210418
!             if( domain_riv_idx(kk) .eq. 2 ) dh = (zb_p + hr_p - zb_n - hr_n) / distance ! yorozuya added

             if( domain_riv_idx(kk) .eq. 2 )then
                if(DBC_switch==0) then
                dh = (zb_p + hr_p_eff - zb_n) / distance ! kinematic wave (+hr_p)
                elseif(DBC_switch==1)then
                dh = (zb_p - zb_n) / distance ! free flow, modified 20250312
                endif 
             endif
!             if( domain_riv_idx(kk) .eq. 2 ) then
!               dh = max( (zb_p - zb_n) / distance, 0.0001) ! same as kinematic wave
!               hr_n = hr_p + zb_p - zb_n
!              endif

 ! the destination cell is outlet (domain = 2 ) with water depth boundary water table
             if( domain_riv_idx(kk) .eq. 2 .and. bound_riv_wlev_switch .ge. 1 ) then
                ! ver 1.4.2 mod by T.Sayama on June 24, 2015
                !if( bound_riv_wlev_idx(1, kk) .le. -100.0 ) cycle ! not boundary
                !dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance ! diffussion
                if( bound_riv_wlev_idx(1, kk) .gt. -100.0 ) dh = ((zb_p + hr_p_eff) - (zb_n + hr_n_eff)) / distance ! diffussion 
               endif

 ! the cell or the destination cell is dam (damflg(k) or damflg(kk) > 0)
             !if ( damflg(k) .gt. 0 .or. damflg(kk) .gt. 0 ) dh = max((zb_p - zb_n) / distance, 0.001) ! kinematic wave
            !off 20240406
            if ( damflg(k) .gt. 0 ) then
            dh = max((zb_p - zb_n) / distance, 0.001)
            endif
            !modified 20240406
            if ( damflg(kk) .gt. 0 ) then
            dh = max(dh, 0.001)
            hr_p = max(hr_p, 0.d0)
            !dh = min(dh, 0.001)
            endif
 ! from a tributary (levee height =< 0) to a main river (levee height > 0) : no reverse flow
 !if( height_idx(k) .le. 0.d0 .and. height_idx(kk) .gt. 0.d0 ) dh = max( (zb_p - zb_n) / distance, 0.001 )

             if( dh .ge. 0.d0 ) then
                 hw = hr_p
                 if( zb_p .lt. zb_n ) hw = max(0.d0, zb_p + hr_p - zb_n)
  !call hq_riv(hw, dh, width_idx(k), qr_temp)
                 call hq_riv(hw, dh, k, width_idx(k), qr_temp)
                 qr_raw = qr_temp
                 qr_capacity = qr_raw
                 hw_eff = hw
             else
  ! reverse flow
                 hw = hr_n
                 if( zb_n .lt. zb_p ) hw = max(0.d0, zb_n + hr_n - zb_p)
                 dh = abs(dh)
  !call hq_riv(hw, dh, width_idx(k), qr_temp)
                 call hq_riv(hw, dh, kk, width_idx(k), qr_temp)
                 qr_raw = -qr_temp
                 qr_capacity = qr_raw
                 hw_eff = hw
             endif

             if(channel_capacity_qr_switch > 0)then
                 if(qr_raw .ge. 0.d0)then
                     hw_eff = hr_p_eff
                     if( zb_p .lt. zb_n ) hw_eff = max(0.d0, zb_p + hr_p_eff - zb_n)
                     hw_eff = min(hw_eff, hcap_link)
                     call hq_riv(hw_eff, dh, k, width_idx(k), qr_temp)
                     qr_capacity = qr_temp
                 else
                     hw_eff = hr_n_eff
                     if( zb_n .lt. zb_p ) hw_eff = max(0.d0, zb_n + hr_n_eff - zb_p)
                     hw_eff = min(hw_eff, hcap_link)
                     call hq_riv(hw_eff, dh, kk, width_idx(k), qr_temp)
                     qr_capacity = -qr_temp
                 endif
             endif

             block_fac = 1.d0
             if(channel_blockage_switch > 0)then
                 block_fac_k = 1.d0
                 if(depth_idx_ini(k) > 0.d0)then
                     depth_ratio_k = max(0.d0, min(1.d0, depth_idx(k)/depth_idx_ini(k)))
                     block_fac_k = channel_blockage_leakage + &
                         (1.d0 - channel_blockage_leakage) * depth_ratio_k**channel_blockage_exponent
                 endif
                 block_fac_kk = 1.d0
                 if(depth_idx_ini(kk) > 0.d0)then
                     depth_ratio_kk = max(0.d0, min(1.d0, depth_idx(kk)/depth_idx_ini(kk)))
                     block_fac_kk = channel_blockage_leakage + &
                         (1.d0 - channel_blockage_leakage) * depth_ratio_kk**channel_blockage_exponent
                 endif
                 block_fac = min(block_fac_k, block_fac_kk)
             endif

             qr_idx(k) = qr_capacity * block_fac
             if(allocated(qr_raw_idx)) qr_raw_idx(k) = qr_raw
             if(allocated(qr_capacity_limited_idx)) qr_capacity_limited_idx(k) = qr_capacity
             if(allocated(qr_effective_depth_idx)) qr_effective_depth_idx(k) = hw_eff
             if(allocated(qr_capacity_depth_idx)) qr_capacity_depth_idx(k) = hcap_link
             if(allocated(qr_blockage_factor_idx)) qr_blockage_factor_idx(k) = block_fac

         enddo           

      end subroutine qr_calc
! water depth and discharge relationship
!subroutine hq_riv(h, dh, w, q)
subroutine hq_riv(h, dh, k, w, q)
    use globals
    implicit none

    real(8) h, dh, w, q
    real(8) a, m
    integer k ! add v1.4

    a = sqrt(abs(dh))/ns_river
    m = 5.d0/3.d0
! modified v1.4
    q = a*h**m*w

    if (sec_map_idx(k) .gt. 0) then
        call sec_hq_riv(h, dh, k, q)
    end if

! discharge per unit area
! (q multiply by width and divide by area)
!q = q * w / area

end subroutine hq_riv

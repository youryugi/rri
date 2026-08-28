! RRI_RivSlo
! river and slope interaction

subroutine funcrs(hr, hs)
    use globals
    implicit none

    real(8) hr(ny, nx), hs(ny, nx)

    call funcrs_dt(hr, hs, dble(dt))

end subroutine funcrs

subroutine funcrs_dt(hr, hs, rs_dt)
    use globals
    use sediment_mod
    implicit none

    real(8) hr(ny, nx), hs(ny, nx)
    real(8) rs_dt

    real(8) hrs ! discharge amount from slope to river [m/s]
    real(8) hrs_self
    real(8) hr_top, hs_top, h1, h2
    real(8) height_top
    real(8) mu1, mu2, mu3
    real(8) hr_new, len, b, ar ! add v1.4
    integer i, j, k, count ! add v1.4 (k, count)

    qrs(:, :) = 0.d0
    if (rs_dt .le. 0.d0) return

    mu1 = (2.d0/3.d0)**(3.d0/2.d0)
    mu2 = 0.35d0
    mu3 = 0.91d0

    do i = 1, ny
        do j = 1, nx

            if (domain(i, j) .eq. 0) cycle
            if (riv(i, j) .eq. 0) cycle

            hs_top = hs(i, j)
            hr_top = hr(i, j) - depth(i, j)
            height_top = max(height(i, j), 0.d0)

            k = riv_ij2idx(i, j)
            len = len_riv_idx(k)

            ! (Case a) : (height = 0 and hr_top < 0) or (height > 0 and hr_top < 0 and hs_top <= height)
            ! -> From slope to river : step fall (hrs : positive)

            if ((height_top .eq. 0.d0 .and. hr_top .lt. 0.d0) .or. &
                (height_top .gt. 0.d0 .and. hr_top .lt. 0.d0 .and. hs_top .le. height_top)) then

                !hrs = mu1 * hs_top * sqrt( 9.81d0 * hs_top ) * dt / length
                hrs = mu1*hs_top*sqrt(9.81d0*hs_top)*rs_dt*len/area
                if (hrs .gt. hs(i, j)) hrs = hs(i, j)
                hs(i, j) = hs(i, j) - hrs
                !hr(i,j) = hr(i,j) + hrs / area_ratio(i,j)
                call hr_update(hr(i, j), hrs*area, k, hr_new) ! add v1.4
                hr(i, j) = hr_new ! add v1.4
                qrs(i, j) = hrs

                ! avoid the situation of hr_top > hs_top
                hs_top = hs(i, j)
                hr_top = hr(i, j) - depth(i, j)
                if (hr_top .ge. -0.00001d0 .and. hr_top .gt. hs_top) then
                    do count = 1, 10
                        call sec_h2b(hr(i, j), k, b)
                        ar = len*b/area
                        hrs = (hs_top - hr_top)/(1.d0 + 1.d0/ar)
                        hs(i, j) = hs(i, j) - hrs
                        !hr(i,j) = hr(i,j) + hrs / area_ratio(i,j)
                        call hr_update(hr(i, j), hrs*area, k, hr_new) ! add v1.4
                        hr(i, j) = hr_new ! add v1.4
                        qrs(i, j) = qrs(i, j) + hrs
                        if (abs(hs(i, j) - (hr(i, j) - depth(i, j))) .lt. 0.00001d0) exit
                        hs_top = hs(i, j)
                        hr_top = hr(i, j) - depth(i, j)
                    end do
                    hr(i, j) = hs(i, j) + depth(i, j) ! 最終手段
                end if

                ! (Case b) : height > 0 and hr_top <= height and hr_top >= 0
                ! -> No exchange

           elseif (height_top .gt. 0.d0 .and. hs_top .le. height_top .and. hr_top .le. height_top .and. hr_top .ge. 0.d0) then
                qrs(i, j) = 0.d0
                continue

                ! (Case c) : hs <= hrt & hrt >= height
                ! -> From river to slope : overtopping (hrs : negative)
                ! (incl. hs = 0 and hrt > 0)

            elseif (hs_top .le. hr_top .and. hr_top .ge. height_top) then

                h1 = hr_top - height_top
                h2 = hs_top - height_top
                if (h1 .le. 0.d0) then
                    hrs = 0.d0
                elseif (h2/h1 .le. 2.d0/3.d0) then
                    !hrs = - mu2 * h1 * sqrt( 2.d0 * 9.81d0 * h1 ) * dt / length
                    hrs = -mu2*h1*sqrt(2.d0*9.81d0*h1)*rs_dt*len/area ! modified v1.4
                else
                    !hrs = - mu3 * h2 * sqrt( 2.d0 * 9.81d0 * (h1 - h2) ) * dt / length
                    hrs = -mu3*h2*sqrt(2.d0*9.81d0*(h1 - h2))*rs_dt*len/area ! modified v1.4
                end if

                call sec_h2b(hr(i, j), k, b)
                ar = len*b/area
                if (abs(hrs/ar) .gt. (hr_top - height_top)) hrs = -(hr_top - height_top)*ar
                !if( abs(hrs / area_ratio(i, j)) .gt. (hr_top-height(i,j)) ) hrs = - (hr_top-height(i,j)) * area_ratio(i, j)
                qrs(i, j) = hrs

                call distribute_riv_overtop_neighbors(i, j, hr_top, hrs, hs, hrs_self)
                if (hrs_self .ne. 0.d0) hs(i, j) = hs(i, j) - hrs_self
                !hr(i,j) = hr(i,j) + hrs / area_ratio(i,j)
                call hr_update(hr(i, j), hrs*area, k, hr_new) ! add v1.4
                hr(i, j) = hr_new ! add v1.4

                ! avoid the situation of hs_top > hr_top
                hs_top = hs(i, j)
                hr_top = hr(i, j) - depth(i, j)
                !if( hr_top .gt. -0.00001d0 .and. hs_top .gt. hr_top ) then ! modified from v1.4 (ここはどっちか？)
                if (hs_top .gt. hr_top) then
                    do count = 1, 10
                        call sec_h2b(hr(i, j), k, b)
                        ar = len*b/area
                        hrs = (hs_top - hr_top)/(1.d0 + 1.d0/ar)
                        hs(i, j) = hs(i, j) - hrs
                        !hr(i,j) = hr(i,j) + hrs / area_ratio(i,j)
                        call hr_update(hr(i, j), hrs*area, k, hr_new) ! add v1.4
                        hr(i, j) = hr_new ! add v1.4
                        qrs(i, j) = qrs(i, j) + hrs
                        if (abs(hs(i, j) - (hr(i, j) - depth(i, j))) .lt. 0.00001d0) exit
                        hs_top = hs(i, j)
                        hr_top = hr(i, j) - depth(i, j)
                    end do
                    hr(i, j) = hs(i, j) + depth(i, j) ! 最終手段
                end if

                ! (Case d) : hs > hrt & hs >= height
                ! -> From slope to river : overtopping (hrs : positive)
                ! (incl. hs = 0 and hrt > 0)

            elseif (hs_top .ge. hr_top .and. hs_top .ge. height_top) then

                h1 = hs_top - height_top
                h2 = hr_top - height_top
                if (h1 .le. 0.d0) then
                    hrs = 0.d0
                elseif (h2/h1 .le. 2.d0/3.d0) then
                    !hrs = mu2 * h1 * sqrt( 2.d0 * 9.81d0 * h1 ) * dt / length
                    hrs = mu2*h1*sqrt(2.d0*9.81d0*h1)*rs_dt*len/area ! v1.4
                else
                    !hrs = mu3 * h2 * sqrt( 2.d0 * 9.81d0 * (h1 - h2) ) * dt / length
                    hrs = mu3*h2*sqrt(2.d0*9.81d0*(h1 - h2))*rs_dt*len/area ! v1.4
                end if

                if (hrs .gt. (hs_top - height_top)) hrs = hs(i, j) - height_top
                qrs(i, j) = hrs

                hs(i, j) = hs(i, j) - hrs
                !hr(i,j) = hr(i,j) + hrs / area_ratio(i,j)
                call hr_update(hr(i, j), hrs*area, k, hr_new) ! add v1.4
                hr(i, j) = hr_new ! add v1.4

                ! avoid the situation of hr_top > hs_top
                hs_top = hs(i, j)
                hr_top = hr(i, j) - depth(i, j)
                if (hr_top .ge. -0.00001d0 .and. hr_top .gt. hs_top) then
                    do count = 1, 10
                        call sec_h2b(hr(i, j), k, b)
                        ar = len*b/area
                        hrs = (hs_top - hr_top)/(1.d0 + 1.d0/ar)
                        hs(i, j) = hs(i, j) - hrs
                        !hr(i,j) = hr(i,j) + hrs / area_ratio(i,j)
                        call hr_update(hr(i, j), hrs*area, k, hr_new) ! add v1.4
                        hr(i, j) = hr_new ! add v1.4
                        qrs(i, j) = qrs(i, j) + hrs
                        if (abs(hs(i, j) - (hr(i, j) - depth(i, j))) .lt. 0.00001d0) exit
                        hs_top = hs(i, j)
                        hr_top = hr(i, j) - depth(i, j)
                    end do
                    hr(i, j) = hs(i, j) + depth(i, j) ! 最終手段
                end if
            elseif(isnan(hr(i,j))) then!---added by Qin
               ! write(*,'(a,3i)') 'k/i/j=', k, i, j
                write(*,'(a,3i)') 'k/i/j=', k, j, ny+1-i ! convert to i,j in IRIC GUI 20250505
                write(*,'(a,4f12.5)') 'hs/hr/depth/height=', hs(i,j), hr(i,j),depth(i,j),height(i,j)
                stop "Error : RivSlo_hr=Nan"
            else

                ! Condition not considered above
                stop "Error : RivSlo"

            end if

            qrs(i, j) = qrs(i, j)/rs_dt ! [m/s]

        end do
    end do

end subroutine funcrs_dt

subroutine distribute_riv_overtop_neighbors(i, j, hr_top, hrs, hs, hrs_self)
    use globals
    use sediment_mod
    implicit none

    integer i, j
    real(8) hr_top, hrs, hs(ny, nx), hrs_self

    integer di(8), dj(8)
    integer ii, jj, n, count, src_sk, tgt_sk
    real(8) overflow_depth, source_level, target_lev, target_level
    real(8) weights(8), sum_weight, dist_fac
    integer cand_i(8), cand_j(8)

    hrs_self = hrs
    if (river_overtop_neighbor_switch .le. 0) return
    if (hrs .ge. 0.d0) return

    src_sk = slo_ij2idx(i, j)
    if (src_sk .le. 0) return

    di = (/ -1, -1, -1, 0, 0, 1, 1, 1 /)
    dj = (/ -1, 0, 1, -1, 1, -1, 0, 1 /)

    overflow_depth = -hrs
    source_level = zb_slo_idx(src_sk) + hr_top
    count = 0
    sum_weight = 0.d0

    do n = 1, 8
        ii = i + di(n)
        jj = j + dj(n)

        if (ii .lt. 1 .or. ii .gt. ny) cycle
        if (jj .lt. 1 .or. jj .gt. nx) cycle
        if (domain(ii, jj) .eq. 0) cycle
        if (riv(ii, jj) .ne. 0) cycle

        tgt_sk = slo_ij2idx(ii, jj)
        if (tgt_sk .le. 0) cycle

        call h2lev(hs(ii, jj), tgt_sk, target_lev)
        target_level = zb_slo_idx(tgt_sk) + target_lev
        dist_fac = max(source_level - target_level, 0.d0)
        if (dist_fac .le. 0.d0) cycle

        if (abs(di(n)) + abs(dj(n)) .eq. 2) dist_fac = dist_fac/sqrt(2.d0)

        count = count + 1
        cand_i(count) = ii
        cand_j(count) = jj
        weights(count) = dist_fac
        sum_weight = sum_weight + dist_fac
    end do

    if (count .eq. 0 .or. sum_weight .le. 0.d0) return

    do n = 1, count
        hs(cand_i(n), cand_j(n)) = hs(cand_i(n), cand_j(n)) + overflow_depth*weights(n)/sum_weight
    end do

    hrs_self = 0.d0

end subroutine distribute_riv_overtop_neighbors

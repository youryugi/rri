module RRI_iric

    use iric

    implicit none
    character(len=64):: cgns_name
    integer:: cgns_f, icount

    public iric_cgns_open, iric_cgns_close
    public iric_read_input_condition
    public iric_cgns_output_result
    public iric_write_cell_real, iric_write_cell_integer
    public iric_read_cell_attr_int, iric_read_cell_attr_real

contains
    subroutine iric_cgns_open()
        implicit none

        integer:: ierr

        !引数取得
        !icount = nargs()
        icount = iargc()
        !if (icount ==  2) then
        if (icount == 1) then
            !call getarg(1, cgns_name, ierr)
            call getarg(1, cgns_name)
        else
            write (*, "(a)") "You should specify an argument."
            stop
        end if

        call cg_iric_open(cgns_name, IRIC_MODE_MODIFY, cgns_f, ierr)
        if (ierr /= 0) stop "cg_iric_open failed"

        ! guiにcgnsファイルを読込みであることを知らせるファイルを生成
        call iric_initoption(IRIC_OPTION_CANCEL, ierr)

    end subroutine

    subroutine iric_cgns_close()
        implicit none
        integer:: ierr

        call cg_iric_close(cgns_f, ierr)
    end subroutine

    subroutine iric_read_input_condition()
        call iric_read_grid()
    end subroutine

    subroutine iric_read_cell_attr_int(name, v)
        use globals
        implicit none

        character(len=*), intent(in):: name
        integer, intent(out):: v(ny, nx)

        integer, dimension(:, :), allocatable:: tmpv
        integer:: i, j, jj, ierr

        allocate (tmpv(nx, ny))
        call cg_iric_read_grid_integer_cell(cgns_f, name, tmpv, ierr)

        do i = 1, nx
            do j = 1, ny
                jj = ny - j + 1
                v(j, i) = tmpv(i, jj)
            end do
        end do
    end subroutine

    subroutine iric_read_cell_attr_real(name, v)
        use globals
        implicit none

        character(len=*), intent(in):: name
        double precision, intent(out):: v(ny, nx)

        double precision, dimension(:, :), allocatable:: tmpv
        integer:: i, j, jj, ierr

        allocate (tmpv(nx, ny))
        call cg_iric_read_grid_real_cell(cgns_f, name, tmpv, ierr)

        do i = 1, nx
            do j = 1, ny
                jj = ny - j + 1
                v(j, i) = tmpv(i, jj)
            end do
        end do
    end subroutine

    subroutine iric_write_result_real(name, v)
        use globals
        implicit none

        character(len=*), intent(in):: name
        double precision, dimension(:, :), allocatable, intent(in):: v

        double precision, dimension(:, :), allocatable:: tmpv
        integer:: i, j, ierr

        allocate (tmpv(nx, ny))

        do i = 1, nx
            do j = 1, ny
                tmpv(i, j) = v(ny - j + 1, i)
                !tmpv(i, j) = v(j, i)
                !tmpv(i, j) = v(i, j)
            end do
        end do

        call cg_iric_write_sol_cell_real(cgns_f, name, tmpv, ierr)
    end subroutine

    subroutine iric_write_result_integer(name, v)
        use globals
        implicit none

        character(len=*), intent(in):: name
        integer, dimension(:, :), allocatable, intent(in):: v

        integer, dimension(:, :), allocatable:: tmpv
        integer:: i, j, ierr

        allocate (tmpv(nx, ny))

        do i = 1, nx
            do j = 1, ny
                tmpv(i, j) = v(ny - j + 1, i)
                !tmpv(i, j) = v(j, i)
                !tmpv(i, j) = v(i, j)
            end do
        end do

        call cg_iric_write_sol_cell_integer(cgns_f, name, tmpv, ierr)
    end subroutine

    subroutine iric_write_cell_real(name, v)
        use globals
        implicit none

        character(len=*), intent(in):: name
        double precision, dimension(:, :), allocatable, intent(in):: v

        double precision, dimension(:, :), allocatable:: tmpv
        integer:: i, j, ierr

        allocate (tmpv(nx, ny))

        do i = 1, nx
            do j = 1, ny
                tmpv(i, j) = v(ny - j + 1, i)
                !tmpv(i, j) = v(j, i)
                !tmpv(i, j) = v(i, j)
            end do
        end do

        call cg_iric_write_grid_real_cell(cgns_f, name, tmpv, ierr)
    end subroutine

    subroutine iric_write_cell_integer(name, v)
        use globals
        implicit none

        character(len=*), intent(in):: name
        integer, dimension(:, :), allocatable, intent(in):: v

        integer, dimension(:, :), allocatable:: tmpv
        integer:: i, j, ierr

        allocate (tmpv(nx, ny))

        do i = 1, nx
            do j = 1, ny
                tmpv(i, j) = v(ny - j + 1, i)
                !tmpv(i, j) = v(j, i)
                !tmpv(i, j) = v(i, j)
            end do
        end do

        call cg_iric_write_grid_integer_cell(cgns_f, name, tmpv, ierr)
    end subroutine

    subroutine iric_read_grid()
        use globals

        integer:: isize, jsize, ierr
        double precision, dimension(:, :), allocatable:: grid_x, grid_y

        call cg_iRIC_Read_Grid2d_Str_Size(cgns_f, isize, jsize, ierr)
        print *, "isize, jsize", isize, jsize
        if (ierr /= 0) stop "CGNS grid read error"
        allocate (grid_x(isize, jsize), grid_y(isize, jsize))
        print *, "grid_x, grid_y allocated"
        call cg_iRIC_Read_Grid2d_Coords(cgns_f, grid_x, grid_y, ierr)
        print *, "cg_iRIC_Read_Grid2d_Coords called"

        ! iRIC から読み込んだ格子はメートル単位の座標系
        ! utm = 1
        xllcorner = grid_x(1, 1)
        yllcorner = grid_y(1, 1)
        cellsize = grid_x(2, 1) - grid_x(1, 1)
        nx = isize - 1
        ny = jsize - 1

        allocate (zs(ny, nx), zb(ny, nx), zb_riv(ny, nx), domain(ny, nx))
        allocate (riv(ny, nx), acc(ny, nx))
        allocate (dir(ny, nx))
        allocate (land(ny, nx))

        print *, "iric_read_cell_attr_real Elevation"
        call iric_read_cell_attr_real("Elevation", zs)
        print *, "iric_read_cell_attr_int Acc"
        call iric_read_cell_attr_int("Acc", acc)
        print *, "iric_read_cell_attr_int Dir"
        call iric_read_cell_attr_int("Dir", dir)
        land = 1
        if (land_switch .eq. 1) then
            print *, "iric_read_cell_attr_int Land"
            call iric_read_cell_attr_int("Land", land)
        end if
        print *, "iric_read_cell_attr_real all ok"

        ! ONLY FOR DEBUGGING
        print *, "ISIZE, JSIZE = ", isize, jsize
        print *, "XLLCORNER, YLLCORNER = ", xllcorner, yllcorner
        print *, "CELLSIZE = ", cellsize

        deallocate (grid_x, grid_y)
    end subroutine

    subroutine iric_write_qu(qs_ave)
        use globals
        implicit none

        double precision, dimension(:, :, :), allocatable, intent(in):: qs_ave
        double precision, dimension(:, :), allocatable:: v
        integer:: i, j

        allocate (v(ny, nx))
        do i = 1, ny
            do j = 1, nx
                v(i, j) = ((qs_ave(1, i, j) + (qs_ave(3, i, j) - qs_ave(4, i, j))/2.d0)*area)
            end do
        end do
        call iric_write_result_real('qu', v)
        deallocate (v)
    end subroutine

    subroutine iric_write_qv(qs_ave)
        use globals
        implicit none

        double precision, dimension(:, :, :), allocatable, intent(in):: qs_ave
        double precision, dimension(:, :), allocatable:: v
        integer:: i, j

        allocate (v(ny, nx))
        do i = 1, ny
            do j = 1, nx
                v(i, j) = ((qs_ave(2, i, j) + (qs_ave(3, i, j) + qs_ave(4, i, j))/2.d0)*area)
            end do
        end do
        call iric_write_result_real('qv', v)
        deallocate (v)
    end subroutine

    subroutine iric_write_gu(qg_ave)
        use globals
        implicit none

        double precision, dimension(:, :, :), allocatable, intent(in):: qg_ave
        double precision, dimension(:, :), allocatable:: v
        integer:: i, j

        allocate (v(ny, nx))
        do i = 1, ny
            do j = 1, nx
                v(i, j) = ((qg_ave(1, i, j) + (qg_ave(3, i, j) - qg_ave(4, i, j))/2.d0)*area)
            end do
        end do
        call iric_write_result_real('gu', v)
        deallocate (v)
    end subroutine

    subroutine iric_write_gv(qg_ave)
        use globals
        implicit none

        double precision, dimension(:, :, :), allocatable, intent(in):: qg_ave
        double precision, dimension(:, :), allocatable:: v
        integer:: i, j

        allocate (v(ny, nx))
        do i = 1, ny
            do j = 1, nx
                v(i, j) = ((qg_ave(2, i, j) + (qg_ave(3, i, j) + qg_ave(4, i, j))/2.d0)*area)
            end do
        end do
        call iric_write_result_real('gv', v)
        deallocate (v)
    end subroutine

    subroutine iric_river_cell_center(k, x, y)
        use globals
        implicit none

        integer, intent(in):: k
        double precision, intent(out):: x, y
        integer:: ii, jj

        ii = riv_idx2i(k)
        jj = riv_idx2j(k)
        x = xllcorner + (dble(jj) - 0.5d0)*cellsize
        y = yllcorner + (dble(ny - ii) + 0.5d0)*cellsize
    end subroutine

    subroutine iric_write_river_polydata_values(id, k, qr_ave, hr, qsb, qss, sumdzb)
        use globals
        implicit none

        integer, intent(in):: id, k
        double precision, dimension(:, :), allocatable, intent(in):: qr_ave, hr
        double precision, dimension(:, :), allocatable, intent(in):: qsb, qss, sumdzb

        double precision:: v_qsb, v_qss, v_sumdzb
        integer:: ierr, ii, jj

        ii = riv_idx2i(k)
        jj = riv_idx2j(k)
        v_qsb = 0.d0
        v_qss = 0.d0
        v_sumdzb = 0.d0
        if (allocated(qsb)) v_qsb = qsb(ii, jj)
        if (allocated(qss)) v_qss = qss(ii, jj)
        if (allocated(sumdzb)) v_sumdzb = sumdzb(ii, jj)
        call cg_iric_write_sol_polydata_real(cgns_f, '01 River water level line[m]', hr(ii, jj), ierr)
        call cg_iric_write_sol_polydata_real(cgns_f, '02 River discharge line[m3_s]', qr_ave(ii, jj), ierr)
        call cg_iric_write_sol_polydata_real(cgns_f, '03 River bed change line[m]', v_sumdzb, ierr)
        call cg_iric_write_sol_polydata_real(cgns_f, '04 Bedload transport line[m3_s]', v_qsb, ierr)
        call cg_iric_write_sol_polydata_real(cgns_f, '05 Suspended sediment line[m3_s]', v_qss, ierr)
    end subroutine


    subroutine iric_simplify_polyline(x, y, npts)
        implicit none

        double precision, dimension(:), intent(inout):: x, y
        integer, intent(inout):: npts

        double precision:: dx1, dy1, dx2, dy2
        integer:: p, nout

        if (npts .le. 2) return

        nout = 1
        do p = 2, npts - 1
            dx1 = x(p) - x(p - 1)
            dy1 = y(p) - y(p - 1)
            dx2 = x(p + 1) - x(p)
            dy2 = y(p + 1) - y(p)
            if (abs(dx1 - dx2) .gt. 1.d-12 .or. abs(dy1 - dy2) .gt. 1.d-12) then
                nout = nout + 1
                x(nout) = x(p)
                y(nout) = y(p)
            end if
        end do

        nout = nout + 1
        x(nout) = x(npts)
        y(nout) = y(npts)
        npts = nout
    end subroutine
    subroutine iric_write_river_polydata(qr_ave, hr, qsb, qss, sumdzb)
        use globals
        use sediment_mod
        implicit none

        double precision, dimension(:, :), allocatable, intent(in):: qr_ave, hr
        double precision, dimension(:, :), allocatable, intent(in):: qsb, qss, sumdzb

        double precision, dimension(:), allocatable:: x, y
        integer:: ierr, l, k, kk, npts, nmax
        logical:: has_link

        if (riv_count .le. 0) return
        has_link = allocated(link_idx_k) .and. allocated(link_ups_k) .and. allocated(link_to_riv)

        call cg_iric_write_sol_polydata_groupbegin(cgns_f, 'River channel lines', ierr)

        if (has_link .and. link_count .gt. 0) then
            do l = 1, link_count
                if (link_idx_k(l) .le. 0) cycle
                nmax = 2
                if (allocated(link_cell_num)) nmax = max(2, link_cell_num(l) + 1)
                allocate (x(nmax), y(nmax))

                npts = 0
                k = link_ups_k(l)
                if (k .le. 0) k = link_idx_k(l)
                do
                    if (k .le. 0) exit
                    if (npts .ge. nmax) exit
                    npts = npts + 1
                    call iric_river_cell_center(k, x(npts), y(npts))
                    if (k .eq. link_idx_k(l)) exit
                    kk = down_riv_idx(k)
                    if (kk .le. 0) exit
                    if (link_to_riv(kk) .ne. l) exit
                    k = kk
                end do

                if (npts .eq. 1) then
                    kk = down_riv_idx(link_idx_k(l))
                    if (kk .gt. 0) then
                        npts = 2
                        call iric_river_cell_center(kk, x(npts), y(npts))
                    else
                        npts = 2
                        x(npts) = x(1) + 0.25d0*abs(cellsize)
                        y(npts) = y(1)
                    end if
                end if

                if (npts .ge. 2) then
                    call iric_simplify_polyline(x, y, npts)
                    call cg_iric_write_sol_polydata_polyline(cgns_f, npts, x, y, ierr)
                    call iric_write_river_polydata_values(l, link_idx_k(l), qr_ave, hr, qsb, qss, sumdzb)
                end if

                deallocate (x, y)
            end do
        else
            allocate (x(2), y(2))
            do k = 1, riv_count
                call iric_river_cell_center(k, x(1), y(1))
                kk = down_riv_idx(k)
                if (kk .gt. 0) then
                    call iric_river_cell_center(kk, x(2), y(2))
                else
                    x(2) = x(1) + 0.25d0*abs(cellsize)
                    y(2) = y(1)
                end if
                call cg_iric_write_sol_polydata_polyline(cgns_f, 2, x, y, ierr)
                call iric_write_river_polydata_values(k, k, qr_ave, hr, qsb, qss, sumdzb)
            end do
            deallocate (x, y)
        end if

        call cg_iric_write_sol_polydata_groupend(cgns_f, ierr)
    end subroutine
    subroutine iric_cgns_output_result( &
        sum_qp_t, qp_t, hs, hr, hg, qr_ave, qs_ave, qg_ave, &
        qsb, qss, sumdzb, sumqsb, sumqss) !added unitchannel display 20250328
        use globals
        use sediment_mod

        double precision, dimension(:, :), allocatable, intent(in):: &
            sum_qp_t, qp_t, hs, hr, hg, qr_ave
        double precision, dimension(:, :, :), allocatable, intent(in):: &
            qs_ave, qg_ave
        double precision, dimension(:, :), allocatable, intent(in):: &      !added for RSR 20240724
            qsb, qss, sumdzb, sumqsb, sumqss

        double precision, dimension(:, :), allocatable :: rain_rate, rain_vol
        double precision, dimension(:, :), allocatable :: tmpv
        integer :: i, j
        integer:: ierr
    !added unitchannel display 20250328    
        integer, dimension(:, :), allocatable :: tmpv_i
    !!!!!!!!!!!!!!!!!!    
        call cg_iric_write_sol_start(cgns_f, ierr)
        call cg_iric_write_sol_time(cgns_f, time, ierr)

        !node valueが出力されていないとグラフが表示できないための仮出力
        allocate (tmpv(1:ny + 1, 1:nx + 1))
        tmpv = 0.0d0
        call cg_iRIC_Write_Sol_Node_Real(cgns_f, "dummy", tmpv, ierr)
        deallocate (tmpv)
 !added unitchannel display 20250328 
        if(sed_switch==2)then
            allocate (tmpv_i(1:ny + 1, 1:nx + 1))
            tmpv_i = 0  
            call cg_iric_write_sol_node_integer (cgns_f,"dummy2", tmpv_i,ierr)
         deallocate (tmpv_i)
        endif 
 !!!!!!!!!!!!!        
        allocate (rain_rate(1:ny, 1:nx), rain_vol(1:ny, 1:nx))
        do i = 1, ny
            do j = 1, nx
                rain_rate(i, j) = qp_t(i, j)*3600.d0*1000.d0
                rain_vol(i, j) = sum_qp_t(i, j)*1000.d0
            end do
        end do

        call iric_write_result_real('total_qp_t[mm]', rain_vol)
        call iric_write_result_real('qp_t[mm_h]', rain_rate)
        call iric_write_result_real('hs[m]', hs)
        call iric_write_result_real('Surface depth[m]', h_surf)
        call iric_write_result_real('hr[m]', hr)
        call iric_write_result_real('qr[m3_s]', qr_ave)
        call iric_write_qu(qs_ave)
        call iric_write_qv(qs_ave)
        call iric_write_result_real('hg[m]', hg)
        call iric_write_gu(qg_ave)
        call iric_write_gv(qg_ave)
        call iric_write_result_real('gampt_ff', gampt_ff)
        call iric_write_river_polydata(qr_ave, hr, qsb, qss, sumdzb)

!-----For RSR model 20240724

        if(sed_switch .ne. 0)then
            call iric_write_result_real('Bedload transport rate[m3_s]', qsb)
            call iric_write_result_real('Suspended sediment transport rate[m3_s]', qss)
            call iric_write_result_real('Elevation change[m]', sumdzb)
            call iric_write_result_real('Total bedload transport[m3]', sumqsb)
            call iric_write_result_real('Total S.S. transport[m3]', sumqss)
           ! if(sed_switch==2) call iric_write_result_real('Suspended sediment concentration[m3/m3]', ssc_ij)
            call iric_write_result_real('Mean diameter [mm]', dmean_out)
            !for slope erosion
            if(slo_sedi_cal_switch>0)  call iric_write_result_real('Elevation change by slope erosion[m]', dzslo)
            if(slo_sedi_cal_switch>0)  call iric_write_result_real('Slope erosion volume[m3]', eroslovol)
            ! if(slo_sedi_cal_switch>0) call iric_write_result_real('Sediment concentration of slope cell[m3/m3]',ss_slope_ij)
            !  if(slo_sedi_cal_switch>0) call iric_write_result_real ('Sediment inflow discharge from slope to unit channel [m3/s]',inflow_sedi_ij)
             ! if(slo_sedi_cal_switch>0) call iric_write_result_real ('Sediment overflow discharge from unit channel [m3/s]', overflow_sedi_ij)            
            if(slo_sedi_cal_switch>0) call iric_write_result_real ('Total sediment supply volume from slope to unit channel[m3]', slo_supply)
            if(slo_sedi_cal_switch>0) call iric_write_result_real ('Total sediment overflow volume from unit channel[m3]', overflow_sed_sum) 
            if(debris_switch==1) call iric_write_result_real('Land slide occurence', LS)
            if(debris_switch==1) call iric_write_result_real('Debris flow path', hki_g_2d)  ! hki_g_2d must be real type
            if(debris_switch==1) call iric_write_result_real('Elevation change (debris flow) [m]', dzslo_mspnt)   !abs remove this later 20260322
            if(debris_switch==1) call iric_write_result_real('Cumulative elevation change (debris flow) [m]', dzslo_mspnt_cum)
            if(debris_switch==1) call iric_write_result_real('Total sediment supply from debris flow [m3]', debri_sup_sum_ij)            
            if(j_drf==1) call iric_write_result_real('Wood_deposition [m3_m2]', vw2d)
            if(j_drf==1) call iric_write_result_real('Wood_total', qwsum_2d)
            if(j_drf==1) call iric_write_result_real('Wood_concentration', cw2d)
            !if(sed_switch==2)
           ! call iric_write_result_real('Suspended sediment concentration[m3/m3]', ssc_ij)
            if(sed_switch==2) call iric_write_result_integer('Unit channel ID', link_ij) !added 20250328
        end if

!-----RSR until here
        deallocate (rain_rate)
        deallocate (rain_vol)

        !if (outswitch_hg /= 0) then
        !  call iric_write_result_real('hg', hg)
        !end if
        !if (outswitch_gu /= 0) then
        !  call iric_write_gu(qg_ave)
        !end if
        !if (outswitch_gv /= 0) then
        !  call iric_write_gv(qg_ave)
        !end if
        !if (outswitch_gampt_ff /= 0) then
        !  call iric_write_result_real('gampt_ff', gampt_ff)
        !end if

        call cg_iric_write_sol_end(cgns_f, ierr)

    end subroutine

end module

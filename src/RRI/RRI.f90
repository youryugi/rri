! RRI.f90
!
! coded by T.Sayama
!
! ver 1.4.2.3
!
program RRI
    use globals
    use runge_mod
    use dam_mod         !this line is modefied for RSR model 20240724
    use sediment_mod    !added for RSR model 20240724
    use tecout_mod
    use RRI_iric
    use iric
    use omp_lib
    implicit none

! variable definition

! topographic variable
    real(8) nodata
    real(8) x1, x2, y1, y2, d1, d2, d3, d4

! rainfall variable
    integer, allocatable :: rain_i(:), rain_j(:)
    integer tt_max_rain
    integer, allocatable :: t_rain(:)
    integer nx_rain, ny_rain
    real(8), allocatable :: qp(:, :, :), qp_t(:, :), qp_t_idx(:)
    real(8), allocatable :: sum_qp_t(:, :)

! calculation variable
    real(8) rho, total_area
    real(8) vr_out

    real(8), allocatable :: hs(:, :), hr(:, :), hg(:, :), inith(:, :)
    real(8), allocatable :: qs_ave(:, :, :), qg_ave(:, :, :), qr_ave(:, :)

    real(8), allocatable :: fs(:), hs_idx(:), fr(:), hr_idx(:), fg(:), hg_idx(:)
    real(8), allocatable :: qr_idx(:), qr_ave_idx(:), qr_ave_temp_idx(:)
    real(8), allocatable :: vr_idx(:)
    real(8), allocatable :: qs_idx(:, :), qs_ave_idx(:, :), qs_ave_temp_idx(:, :)
    real(8), allocatable :: qg_idx(:, :), qg_ave_idx(:, :), qg_ave_temp_idx(:, :)
    real(8), allocatable :: gampt_ff_idx(:), gampt_f_idx(:)
    real(8), allocatable :: min_hs_idx(:) !added 20250312
 
!real(8), allocatable :: rdummy_dim(:)

!-------added for RSR model 20240724
    real(8), allocatable :: qsur_ave(:,:,:) !added 20231226
	real(8), allocatable :: hr_idx2(:)
	real(8), allocatable :: hr_idxa(:)
    real(8), allocatable :: qsur_ave_temp_idx(:,:)
    real(8), allocatable :: ust(:,:), qsb(:,:), qss(:,:), qsw(:,:),ssc_ij(:,:) !added 20250405
    real(8), allocatable :: ust_idx(:), qsb_idx(:), qss_idx(:), qsw_idx(:),sumqsb_idx(:),sumqss_idx(:),ssc_idx(:)!---modified by Qin, added ssc_ij 20250405
    real(8), allocatable :: area_idx(:), water_v_idx(:)!----modified by Qin 2021/6/11
    real(8), allocatable :: dzb_temp(:), sumdzb_temp(:),sumdzb_idx(:)!---modified by Qin
    real(8), allocatable :: sumqsb(:,:), sumqss(:,:),sumdzb(:,:) !----modified by Qin 2021/7/28
		 
	type(sed_struct) :: sed_temp                          !***** added by Chilli	2015/07/01
	type(sed_struct), allocatable, save :: sed(:,:)
    type(sed_struct), allocatable, save :: sed_idx(:)
    type(sed_struct) :: sed_temp2                          !***** added by harada
    type(sed_struct), allocatable, save :: sed_lin(:)
    real(8), allocatable :: ust_lin(:), qsb_lin(:), qss_lin(:), qsw_lin(:), qss_b(:) !---modified by Qin 2021/5/27
    real(8), allocatable ::  water_v_lin(:)
    real(8), allocatable :: dzb_temp_lin(:), sumdzb_lin(:)
    integer maxloc_hs(2) !20230430
    real(8), allocatable :: hr_lin(:)

	real(8) slope11, leng_sum
    real(8) dhs,slo_Depo_sum,sum_fmslo,sum_fmgully ,d80 !20230430
    real(8) A1, AA !20220805
    real(8) lev_p,lev_n,I_sur,I_sur_x,I_sur_y !modified for slope erosion

	integer mmm, iexel
	integer m,mm,nn,kk, kk1, l, ll,f,n
    integer downstream_k
		 
	character*256 fi_temp
    character*4 t_char2
    character*256 hydro_hr_file
    character*256, allocatable :: location_name(:)
    parameter(hydro_hr_file = "hydro_hr.txt")

    character*256 qsb_file, qss_file
	parameter(qsb_file = "bedload.txt")
	parameter(qss_file = "susload.txt")
    integer kkk
    integer k1
    integer ku1, ku2
    integer kcc, ku, kd
    real(8) Fmall, dzb_cap
    real(8) slo_sedi_cal_duration !added 20231204
!------------------------------------------------------------

! other variable
    integer ni, nj
    integer i, j, t, k, ios, itemp, jtemp, tt, ii, jj,riv_k
    integer out_next
    real(8) out_dt
    real(8) rtemp
    real(8) ss, sr, si, sg, sinit, sout
    integer idummy
    real(8) rdummy
    real(8) rain_sum
    real(8) distance
    real(8) ddt_chk_riv, ddt_chk_slo
    character*256 ctemp
    character*6 t_char
    integer div_org_i, div_org_j, div_dest_i, div_dest_j
    integer :: tmp_idx(2), tmp_ii, tmp_jj, size, ier

! calcHydro
    character*256 hydro_file
    parameter(hydro_file="hydro.txt")
    integer, allocatable :: hydro_i(:), hydro_j(:)
    integer maxhydro

    integer :: ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! STEP 0: FILE NAME AND PARAMETER SETTING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call RRI_Read
    if (omp_num_threads .gt. 0) then
        call omp_set_num_threads(omp_num_threads)
        write (*, '("OpenMP threads : ", i5)') omp_num_threads
    else
        write (*, '("OpenMP threads : auto (max ", i5, ")")') omp_get_max_threads()
    end if
    if (run_type == 0) then
        call iric_cgns_close()
        write (*, "(a)") "***** check grid attributes *****"
        stop
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! STEP 1: FILE READING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! max timestep
    maxt = lasth*3600/dt

! load data from CGNS file
    call iric_cgns_open

!read grid
    call cg_iRIC_Read_Grid2d_Str_Size(cgns_f, ni, nj, ierr)
    allocate (gxx(ni, nj), gyy(ni, nj))
    call cg_iRIC_Read_Grid2d_Coords(cgns_f, gxx, gyy, ierr)
    xllcorner = gxx(1, 1); yllcorner = gyy(1, 1)
    cellsize = dabs(gxx(2, 1) - gxx(1, 1))
    nodata = -9999

!change nx, ny as cell numbers
    nx = ni - 1; ny = nj - 1

    allocate (zs(ny, nx), zb(ny, nx), zb_riv(ny, nx), domain(ny, nx))
    call iric_read_cell_attr_real('elevation_c', zs)

    allocate (riv(ny, nx), acc(ny, nx))
    call iric_read_cell_attr_int('acc_c', acc)

    allocate (dir(ny, nx))
    call iric_read_cell_attr_int('dir_c', dir)

    allocate (land(ny, nx))
    call iric_read_cell_attr_int('landuse_c', land)
    write (*, *) "num_of_landuse : ", num_of_landuse
    where (land .le. 0 .or. land .gt. num_of_landuse) land = num_of_landuse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! STEP 2: CALC PREPARATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! dx, dy calc
! d1: south side length
    x1 = xllcorner
    y1 = yllcorner
    x2 = xllcorner + nx*cellsize
    y2 = yllcorner
    if (utm .eq. 0) call hubeny_sub(x1, y1, x2, y2, d1)

! d2: north side length
    x1 = xllcorner
    y1 = yllcorner + ny*cellsize
    x2 = xllcorner + nx*cellsize
    y2 = yllcorner + ny*cellsize
    if (utm .eq. 0) call hubeny_sub(x1, y1, x2, y2, d2)

! d3: west side length
    x1 = xllcorner
    y1 = yllcorner
    x2 = xllcorner
    y2 = yllcorner + ny*cellsize
    if (utm .eq. 0) call hubeny_sub(x1, y1, x2, y2, d3)

! d4: east side length
    x1 = xllcorner + nx*cellsize
    y1 = yllcorner
    x2 = xllcorner + nx*cellsize
    y2 = yllcorner + ny*cellsize
    if (utm .eq. 0) call hubeny_sub(x1, y1, x2, y2, d4)

    if (utm .eq. 1) then
        dx = cellsize
        dy = cellsize
    else
        dx = (d1 + d2)/2.d0/real(nx)
        dy = (d3 + d4)/2.d0/real(ny)
    end if
    write (*, *) "dx [m] : ", dx, "dy [m] : ", dy

! length and area of each cell
    length = sqrt(dx*dy)
    area = dx*dy

! river widhth, depth, leavy height, river length, river area ratio
    allocate (width(ny, nx), depth(ny, nx), height(ny, nx), len_riv(ny, nx), area_ratio(ny, nx))

!-------added for RSR model 20240724
    allocate (riv_0th(ny,nx)) ! 0th order channel added by Qin 2021/6/23
    allocate (zb_riv0(ny, nx),dzb_riv(ny, nx),zb_riv_slope(ny, nx),zs_slope(ny, nx),riv_ij2idx(ny, nx))
    riv_0th = 0 ! not 0th order channel added by Qin 2021/6/23
    where(acc.gt.riv_thresh.and.acc.lt.max_acc_0th_riv.and.riv.eq.1) riv_0th =  1 ! 0th order channel
!-------RSR until here   

    width = 0.d0
    depth = 0.d0
    height = 0.d0
    len_riv = 0.d0

    area_ratio = 0.d0

    call iric_read_cell_attr_real('width_c', width) ! set the width value -9999 to eliminate the incorrect river cell generated by regime law
    call iric_read_cell_attr_real('depth_c', depth)
    call iric_read_cell_attr_real('height_c', height)

    riv = 0 ! slope cell
    if (riv_thresh .gt. 0) then
        where (acc .gt. riv_thresh .and. width .gt. 0) riv = 1 ! river cell
    end if

! river data is replaced by the information in files
    if (rivfile_switch .ge. 1) then
        riv = 0
        riv_thresh = 1
        where (width .gt. 0) riv = 1 ! river cells (if width >= 0.)
        where (height(:, :) .lt. 0.d0) height(:, :) = 0.d0
    end if
    where (riv .eq. 1) len_riv = length

! river cross section is set by section files
    allocate (sec_map(ny, nx))
    sec_map = 0
    if (sec_switch .eq. 1) then
        call read_gis_int(sec_map_file, sec_map)
        sec_id_max = maxval(sec_map(:, :))
        call set_section
    end if

! river length is set by input file
    allocate (sec_length(ny, nx))
    sec_length = 0
    if (sec_length_switch .eq. 1) then
        call read_gis_real(sec_length_file, sec_length)
        where (sec_length .gt. 0.d0) len_riv = sec_length
    end if

    if (rivfile_switch .eq. 2) then
        ! levee on both river and slope grid cells : zs is increased with height
        where (height .gt. 0.d0) zs = zs + height
    else
        ! levee on only slope grid cells : zs is increased with height
        where (height .gt. 0.d0 .and. riv .eq. 0) zs = zs + height
    end if

    where (riv .eq. 1) area_ratio = width/len_riv

    zb_riv = zs
    zb_riv0=zs   !this line is added for RSR 20240724
    do i = 1, ny
        do j = 1, nx
            zb(i, j) = zs(i, j) - soildepth(land(i, j))
            if (riv(i, j) .eq. 1) zb_riv(i, j) = zs(i, j) - depth(i, j)
        end do
    end do

! domain setting

! domain = 0 : outside the domain
! domain = 1 : inside the domain
! domain = 2 : outlet point (where dir(i,j) = 0 or dir(i,j) = -1),
!              and cells located at edges
    domain = 0
    num_of_cell = 0
    do i = 1, ny
        do j = 1, nx
            if (zs(i, j) .gt. -100.d0) then
                domain(i, j) = 1
                if (dir(i, j) .eq. 0) domain(i, j) = 2
                if (dir(i, j) .eq. -1) domain(i, j) = 2
                num_of_cell = num_of_cell + 1
            end if
        end do
    end do
    write (*, *) "num_of_cell : ", num_of_cell
    write (*, *) "total area [km2] : ", num_of_cell*area/(10.d0**6.0d0)

!-------added for RSR model 20240724
 ! reading dam file modified by Qin 2021/4/28
    if (dam_switch ==1) then
          call dam_read
          write(*,*) "for check01"
          do i  = 1, dam_num
          area_ratio(dam_iy(i), dam_ix(i)) = dam_reserv_area(i) / area
          if(reserv_elev(i)> 0.) zs(dam_iy(i), dam_ix(i)) = reserv_elev(i)
          if(dam_discharge_switch(i)==1) call dam_discharge_operation(i)
           write(*,'(a,i)') "for check02, i =", i
          enddo
          write(*,'(a,i)') "End time of dam outflow discharge(hour) : ",max_dam_discharge_t
    endif 


    do i = 1, ny
        do j = 1, nx
            if (riv(i, j) .eq. 1) zb_riv0(i,j) = zb_riv(i,j)
        end do
    end do
    kill_n = 0

do
!-------RSR model until here    

! river index setting (modified for RSR)
    call riv_idx_setting

!moved to here activate the diversion in RSR     
 ! div file
!div_id_max = 0
    if (div_switch .eq. 1) then
        allocate (div_org_idx(div_id_max), div_dest_idx(div_id_max), div_rate(div_id_max))
        !rewind(20)

        do k = 1, div_id_max
            call cg_iric_read_bc_indicessize(cgns_f, "div", k, size, ier)
            if (size /= 1) then
                write (*, *) "Error:Number of indicessize setting div must be one."
                stop
            end if

            call cg_iric_read_bc_indices(cgns_f, "div", k, tmp_idx, ier)
            div_org_i = ny - tmp_idx(2) + 1
            div_org_j = tmp_idx(1)
            call cg_iric_read_bc_integer(cgns_f, "div", k, "div_dest_i", tmp_ii, ier)
            call cg_iric_read_bc_integer(cgns_f, "div", k, "div_dest_j", tmp_jj, ier)
            call cg_iric_read_bc_real(cgns_f, "div", k, "div_rate", div_rate(k), ier)

            div_org_idx(k) = riv_ij2idx(div_org_i, div_org_j)
            !modified 20250313
            if(tmp_jj>0 .and. tmp_ii>0 ) then
            div_dest_i = ny - tmp_jj + 1
            div_dest_j = tmp_ii
            div_dest_idx(k) = riv_ij2idx(div_dest_i, div_dest_j)
            if(div_dest_idx(k)<1) div_dest_idx(k) = 0
            else
            div_dest_idx(k) = 0
            div_switch=2
            endif
            write (*, *) "done: reading div dest idx", div_dest_idx(k)
        end do
        write (*, *) "done: reading div file"
        !close(20)
    end if

!-------added for RSR model 20240724
    if(sed_switch.ge.1)call riv_set4sedi
    if(kill_n .le. 0) exit
end do

    if(sed_switch.ge.1)then
        open( 22222, file = 'kij_file.txt', status = 'unknown' )
        write(22222,'(a)') '    k rivi rivj'
        do k = 1, riv_count
            write(22222,'(3i5)') k, riv_idx2i(k), riv_idx2j(k)
        enddo
        close(22222)
    end if
!-------RSR model until here     
! slope index setting
    call slo_idx_setting

! added for RSR model 20240724  connect slope to river
    if (sed_switch.eq.2) call slo_to_riv

! reading dam file
!    if (sed_switch.eq.0) call dam_read  !this line is modified for RSR model 20240724 !not work skip later

! initial condition
    allocate (hs(ny, nx), hr(ny, nx), hg(ny, nx), gampt_ff(ny, nx))
    allocate (gampt_f(ny, nx), qrs(ny, nx))
    allocate (min_hs_idx(slo_count)) ! added 20250312

    hr = -0.1d0
    hs = -0.1d0
    hg = -0.1d0
    gampt_ff = 0.d0
    gampt_f = 0.d0
    qrs = 0.d0

!-------added for RSR model 20240724 modified 20250312
        where(riv.eq.1) hr = hr0
    do k= 1, slo_count
           min_hs_idx(k)= da_idx(k)*wc0
            i = slo_idx2i(k)
            j = slo_idx2j(k)
            hs(i,j)= min_hs_idx(k)
    enddo        
!--------RSR model until here 
! for slope cells
! if init_slo_switch = 1 => read from file

    if (init_slo_switch .eq. 1) then

        allocate (inith(ny, nx))
        inith = 0.d0
        open (13, file=initfile_slo, status="old")
        do i = 1, ny
            read (13, *) (inith(i, j), j=1, nx)
        end do
        where (inith .le. 0.d0) inith = 0.d0
        where (domain .eq. 1 .and. inith .ge. 0.d0) hs = inith
        deallocate (inith)
        close (13)

    end if

! for river cells
! if init_riv_switch = 1 => read from file

    if (init_riv_switch .eq. 1) then

        allocate (inith(ny, nx))
        inith = 0.d0
        open (13, file=initfile_riv, status="old")
        do i = 1, ny
            read (13, *) (inith(i, j), j=1, nx)
        end do
        where (inith .le. 0.d0) inith = 0.d0
        where (riv .eq. 1 .and. inith .ge. 0.d0) hr = inith
        deallocate (inith)
        close (13)

    end if

! for slope cells
! if init_gw_switch = 1 => read from file

    if (init_gw_switch .eq. 1) then

        allocate (inith(ny, nx))
        inith = 0.d0
        open (13, file=initfile_gw, status="old")
        do i = 1, ny
            read (13, *) (inith(i, j), j=1, nx)
        end do
        where (inith .le. 0.d0) inith = 0.d0
        where (domain .eq. 1 .and. inith .ge. 0.d0) hg = inith
        deallocate (inith)
        close (13)

    end if

! if init_gampt_ff_switch = 1 => read from file

    if (init_gampt_ff_switch .eq. 1) then

        allocate (inith(ny, nx))
        inith = 0.d0
        open (13, file=initfile_gampt_ff, status="old")
        do i = 1, ny
            read (13, *) (inith(i, j), j=1, nx)
        end do
        where (inith .le. 0.d0) inith = 0.d0
        where (domain .eq. 1) gampt_ff = inith
        deallocate (inith)
        close (13)

    end if

! boundary conditions
    call read_bound



! emb file
!if( emb_switch.eq.1 ) then
! allocate (emb_r(ny, nx), emb_b(ny, nx))
! call read_gis_real(embrfile, emb_r)
! call read_gis_real(embbfile, emb_b)

! call sub_slo_ij2idx(emb_r, emb_r_idx)
! call sub_slo_ij2idx(emb_b, emb_b_idx)
!endif

!-------added for RSR model 20240724    !sediment condition setting
    if(sed_switch.ge.1)then
        allocate (sed(ny,nx))
        call sed_distr2 (sed)
        call cal_zs_slope    !moved to RRI_sed2.f90
        if( hydro_switch .eq. 1 )open(1013, file = hydro_hr_file)
    end if

!--------RSR model until here 

! hydro file
    if (hydro_switch .eq. 1) then
        open (5, file=location_file, status="old")
        open (1012, file=hydro_file)
        i = 1
        do
            read (5, *, iostat=ios) ctemp, itemp, itemp
            if (ios .ne. 0) exit
            i = i + 1
        end do
        maxhydro = i - 1
        rewind (5)
        allocate (hydro_i(maxhydro), hydro_j(maxhydro))
        do i = 1, maxhydro
            read (5, *) ctemp, hydro_i(i), hydro_j(i)
        end do
        close (5)
    end if

! dynamic allocation
    allocate (qs_ave(i4, ny, nx), qr_ave(ny, nx), qg_ave(i4, ny, nx))

    allocate (qr_idx(riv_count), qr_ave_idx(riv_count), qr_ave_temp_idx(riv_count), hr_idx(riv_count))
    allocate (fr(riv_count), vr_temp(riv_count), hr_err(riv_count), vr_err(riv_count))
    allocate (vr_idx(riv_count))
    allocate (kr2(riv_count), kr3(riv_count), kr4(riv_count), kr5(riv_count), kr6(riv_count))

    allocate (qs_idx(i4, slo_count), qs_ave_idx(i4, slo_count), qs_ave_temp_idx(i4, slo_count), hs_idx(slo_count))
    allocate (qp_t_idx(slo_count))
    allocate (fs(slo_count), hs_temp(slo_count), hs_err(slo_count))
    allocate (ks2(slo_count), ks3(slo_count), ks4(slo_count), ks5(slo_count), ks6(slo_count))

    allocate (qg_idx(i4, slo_count), qg_ave_idx(i4, slo_count), qg_ave_temp_idx(i4, slo_count), hg_idx(slo_count))
    allocate (fg(slo_count), hg_temp(slo_count), hg_err(slo_count))
    allocate (kg2(slo_count), kg3(slo_count), kg4(slo_count), kg5(slo_count), kg6(slo_count))
    allocate (gampt_ff_idx(slo_count), gampt_f_idx(slo_count))

    allocate (rain_i(ny), rain_j(nx))
    allocate (qe_t_idx(slo_count))
    allocate (evp_i(ny), evp_j(nx))
    allocate (aevp(ny, nx), aevp_tsas(slo_count), exfilt_hs_tsas(slo_count), rech_hs_tsas(slo_count))

!----dynamic allocation for RSR model !moved 20240724
        allocate (hr_idx2(riv_count), hr_idxa(riv_count))   !moved for RSR 20240724
        allocate ( qsur_ave_temp_idx(i4,slo_count)) 
         allocate (ust(ny, nx), qsb(ny, nx), qss(ny, nx), qsw(ny,nx),ssc_ij(ny, nx)) !added 20250405
         allocate ( sumqsb(ny, nx), sumqss(ny, nx),sumdzb(ny,nx)) !----added by Qin 2021/07/28
         allocate (Emb(ny,nx), Et(ny,nx), Nb(ny,nx))
         allocate (Emb_idx(riv_count), Et_idx(riv_count), zb_roc_idx(riv_count))!, zb_air_idx(riv_count))
         allocate (ust_idx(riv_count), qsb_idx(riv_count), qss_idx(riv_count), qsw_idx(riv_count),ssc_idx(riv_count)) ! added ssc_idx 20250405
         allocate (zb_riv_slope_idx(riv_count))
	     allocate (main_riv(riv_count))
	     allocate (Nb_idx(riv_count))
	     allocate (zb_riv0_idx(riv_count))
         allocate (sed_idx(riv_count))
         allocate (dzb_temp(riv_count), sumdzb_temp(riv_count))
         allocate (sumqsb_idx(riv_count), sumqss_idx(riv_count),sumdzb_idx(riv_count)) !----added by Qin check 2021/5/20
         allocate ( area_idx(riv_count), water_v_idx(riv_count))!----modified by Qin 2021/6/11
         allocate (slo_s_dsum(riv_count,Np),slo_s_sum(riv_count))
! dynamic allocation for link model (Added by harada 2020_11_19)
         allocate (Emb_lin(link_count), Et_lin(link_count), zb_roc_lin(link_count))
	     allocate (Nb_lin(link_count))
         allocate (ust_lin(link_count), qsb_lin(link_count), qss_lin(link_count), qsw_lin(link_count),qss_b(link_count))  !---added by Qin 2021/5/27
         allocate (zb_riv_slope_lin(link_count)) !moved zb_riv_slope0_lin to subroutine riv_set4sedi;20240601
	     allocate (zb_riv0_lin(link_count))
         allocate (sed_lin(link_count))
         allocate (dzb_temp_lin(link_count), sumdzb_lin(link_count)) 
         allocate (ss_lin(link_count))
         allocate (hr_lin(link_count))
         allocate ( water_v_lin(link_count))
         allocate (dam_sedi_total(dam_num),dam_sedi_total_b(dam_num),dam_sedi_total_s(dam_num),dam_sedi_total_w(dam_num),dam_sedi_totalV(dam_num),dam_sedi_qsi(dam_num,Npmax))  !added by Qin 
         allocate(dam_outflow_qb(dam_num),dam_outflow_qs(dam_num),dam_outflow_cs(dam_num))
!pause'before dynamic allocation2'         
! dynamic allocation for slope erosion (Added 2021/11/1)
         allocate (slo_sur_zb(slo_count),slo_sur_zb_before(slo_count)) 
         allocate (fmslo(slo_count,Np),fgully(slo_count,Np)) !modified 20231226
         allocate (Ero_mslo(slo_count,Np,maxt+1),Ero_kslo(slo_count,maxt+1),dzslo_idx(slo_count), dzslo(ny, nx), Ero_slo_vol(slo_count), eroslovol(ny, nx),dzslo_fp_idx(slo_count,Np),inum_sed_dm(slo_count))
         allocate (slo_to_lin_sed(link_count,Np),slo_to_lin_sed_sum(link_count),slo_vol_remain(link_count,Np),slo_to_lin_sum_di(link_count,Np),slo_to_lin_sed_sum_idx(riv_count),slo_supply(ny,nx),slo_to_lin_idx_di(riv_count,Np),lin_to_slo_idx_di(riv_count,Np))!added 20250503
         allocate (lin_to_slo_sed_sum(link_count),lin_to_slo_sum_di(link_count,Np),lin_to_slo_sed_sum_idx(riv_count),overflow_sed_sum(ny,nx))
         allocate (B_chan(slo_count),D_chan(slo_count), h_chan(slo_count),q_chan(slo_count), I_chan(slo_count),area_chan(slo_count),gully_sedi_dep(slo_count, Np)) ! 2022/07/21 20231226
         allocate (infil_w_depth(slo_count)) !added 20250618
         allocate (ss_slope(slo_count),qss_slope(slo_count))
        allocate (slo_qsisum(slo_count,Np),slo_qsi(slo_count,Np),slo_ssi(slo_count,Np),slo_Dsi(slo_count,Np),slo_Esi(slo_count,Np))
        allocate (inflow_sedi(riv_count), overflow_sedi(riv_count),overflow_sedi_di(link_count,Np),overdepo_sedi_di(link_count,Np)) !modified 20230924
        allocate (h_surf(ny,nx),ss_slope_ij(ny,nx) ) !added 20231101
        allocate (inflow_sedi_ij(ny,nx),overflow_sedi_ij(ny,nx)) !for slope erosion
! allocation for debris (2021/11/17)
         allocate (c_dash(slo_count), hsc(slo_count), pw(slo_count), sf(slo_count))
         allocate (vol(slo_count), vcc(slo_count), vcf(slo_count),LS_idx(slo_count),LS(ny, nx), Qg(slo_count),hki_g(slo_count), hki_g_2d(ny,nx), water_v_cell(slo_count))
         allocate (dzslo_mspnt_idx(slo_count), dzslo_mspnt(ny, nx), soildepth_idx_deb(slo_count))
         allocate (vo_total(riv_count), vo_total_river(link_count), hki_area(riv_count), vo_total_l(link_count))
         allocate (debri_sup_sum(link_count), debri_sup_sum_di(link_count, Np),debri_sup_sum_ij(ny,nx)) !added 20240424
         allocate (cw(link_count), qw(link_count), vw(link_count), qwsum(link_count), qwsum_total(link_count))
         allocate (vw_idx(riv_count), qwsum_total_idx(riv_count), vw2d(ny,nx), dmean_out(ny,nx),cw_idx(riv_count), cw2d(ny,nx), qwsum_2d(ny,nx))  !20240315  !20240724 added !20241125 cw_idx,cw2d added !20251002 qwsum_2d added
         allocate (n_link_depth(link_count))
         allocate (depth_idx_ini(riv_count),width_idx_ini(riv_count))
! allocation for sedput
         allocate (fm_sedput(link_count,Np),put_depth(link_count),put_depth_remain(link_count),l_put(link_count))
!added for river width adjustment; 20240419
!----RSR until here

! array initialization
    qr_ave(:, :) = 0.d0
    qr_idx(:) = 0.d0
    qr_ave_idx(:) = 0.d0
    qr_ave_temp_idx(:) = 0.d0

    hr_idx(:) = 0.d0
    vr_idx(:) = 0.d0
    fr(:) = 0.d0
    hr_err(:) = 0.d0
    vr_temp(:) = 0.d0
    vr_err(:) = 0.d0
    kr2(:) = 0.d0
    kr3(:) = 0.d0
    kr4(:) = 0.d0
    kr5(:) = 0.d0
    kr6(:) = 0.d0

    qs_ave(:, :, :) = 0.d0
    qs_idx(:, :) = 0.d0
    qs_ave_idx(:, :) = 0.d0
    qs_ave_temp_idx(:, :) = 0.d0
    hs_idx(:) = 0.d0
    qp_t_idx(:) = 0.d0
    fs(:) = 0.d0
    hs_temp(:) = 0.d0
    hs_err(:) = 0.d0
    ks2(:) = 0.d0
    ks3(:) = 0.d0
    ks4(:) = 0.d0
    ks5(:) = 0.d0
    ks6(:) = 0.d0

    qg_ave(:, :, :) = 0.d0
    qg_idx(:, :) = 0.d0
    qg_ave_idx(:, :) = 0.d0
    qg_ave_temp_idx(:, :) = 0.d0
    hg_idx(:) = 0.d0
    fg(:) = 0.d0
    hg_temp(:) = 0.d0
    hg_err(:) = 0.d0
    gampt_ff_idx(:) = 0.d0
    gampt_f_idx(:) = 0.d0

    rain_i(:) = 0
    rain_j(:) = 0
    ksv(:) = 0.d0
    faif(:) = 0.d0
    aevp(:, :) = 0.d0
    hg_idx(:) = 0.d0
    evp_i(:) = 0
    evp_j(:) = 0
    aevp_tsas(:) = 0.d0
    exfilt_hs_tsas(:) = 0.d0
    rech_hs_tsas(:) = 0.d0

!-------added for RSR model 20240724
         hr_idx2(:) = 0.d0
         hr_idxa(:) = 0.d0
         !added 20231226
         qsur_ave(:,:,:) =0.d0
         !qsur_idx(:,:) =0.d0
         !qsur_ave_idx(:,:) = 0.d0
         qsur_ave_temp_idx(:,:) =0.d0
         ! Sediment variables initialization (Added by Yorozuya 2015_06_19)--
	     ust(:,:) = 0.d0
         qsb(:,:) = 0.d0
         qss(:,:) = 0.d0
         qsw(:,:) = 0.d0
	     Emb(:,:) = 0.d0
	     Et(:,:) = 0.d0
	     Nb(:,:) = 0
         ssc_ij(:,:) = 0.d0 !added 20250405
	
	     Emb_idx(:) = 0.d0
	     Et_idx(:) = 0.d0
	     zb_roc_idx(:) = 0.d0
	     !zb_air_idx(:) = 0.d0
	     Nb_idx(:) = 0
         zb_riv0_idx(:) = 0.d0
!---added by Qin 2021/6/12
         water_v_idx(:) = 0.d0
         area_idx(:) = 0.d0
         !ssc_idx(:) =0.d0
!---added by Qin 2021/7/28
         sumqsb(:,:) = 0.d0
         sumqss(:,:) = 0.d0         

! Sediment variables initialization (Added by harada 2020_11_19)--
	     Emb_lin(:) = 0.d0
	     Et_lin(:) = 0.d0
	     Nb_lin(:) = 0
	     zb_roc_lin(:) = 0.d0
         zb_riv0_lin(:) = 0.d0

         hr_lin(:) = 0.d0
         water_v_lin(:) =0.d0 
         
         ust_lin(:) =0.d0 
         qsb_lin(:) =0.d0 
         qss_lin(:) =0.d0 
         qsw_lin(:) =0.d0 
         qss_b(:) = 0.d0    !added by Qin 2021/5/27
         ss_lin(:) = 0.d0
!---added for slope erosion (2021/11/1)
         slo_sur_zb(:) =0.d0
         slo_sur_zb_before(:) = 0.d0
         fmslo(:,:) = 0.d0
         fgully(:,:) = 0.d0 !20231226
         Ero_mslo(:,:,:) = 0.d0
         Ero_kslo(:,:) = 0.d0 !added by Qin 2021/11/7
         Dzslo_idx(:) = 0.d0
         dzslo(:,:) = 0.d0
         Ero_slo_vol(:) = 0.d0
         eroslovol(:,:) = 0.d0
         slo_to_lin_sed(:,:) = 0.d0
         slo_to_lin_sed_sum(:)=0.d0   
         slo_vol_remain(:,:) = 0.d0  
         slo_to_lin_sum_di(:,:) = 0.d0  
         slo_s_dsum(:,:) =0.d0 !added by Qin 2021/11/8
         slo_s_sum(:) = 0.d0
         slo_to_lin_sed_sum_idx(:) = 0.d0
         slo_supply(:,:) = 0.d0
         dzslo_fp_idx(:,:) = 0.d0
         fm_sedput(:,:) =0.d0
         slope_erosion_total = 0.d0
         lin_to_slo_sed_sum(:)=0.d0
         lin_to_slo_sum_di(:,:)= 0.d0
         lin_to_slo_sed_sum_idx(:) = 0.d0
         overflow_sed_sum(:,:) = 0.d0
         overflow_sedi_di(:,:) = 0.d0
         overdepo_sedi_di(:,:) = 0.d0
         slo_to_lin_idx_di(:,:)=0.d0!added 20250503
         lin_to_slo_idx_di(:,:)=0.d0
!2022/07/21
         B_chan(:) = 0.d0
         D_chan(:) = 0.d0
         h_chan(:) = 0.d0
        q_chan(:) = 0.d0
         I_chan(:) = 0.d0
         area_chan(:) =0.d0
         ss_slope(:) = 0.d0
         qss_slope(:) = 0.d0 
         slo_qsisum(:,:)=0.d0
         slo_qsi(:,:)=0.d0
         slo_ssi(:,:)=0.d0
         slo_Dsi(:,:)=0.d0
         slo_Esi(:,:)=0.d0
         h_surf(:,:) = 0.d0 !added 20231101
         ss_slope_ij(:,:) =0.d0 !added 20231101
         gully_sedi_dep(:,:) =0.d0
         inflow_sedi_ij(:,:) =0.d0 !for slope erosion
         overflow_sedi_ij(:,:) = 0.d0 !for slope erosion
         infil_w_depth(:)  =0.d0 !added 20250618
!---added for debris (2021/11/17)
         c_dash(:) = 0.d0
         hsc(:) = 0.d0
         pw(:) = 0.d0
         sf(:) = 0.d0
         vol(:) = 0.d0
         vcc(:) = 0.d0 
         vcf(:) = 0.d0
         dzslo_mspnt_idx(:) = 0.d0
         dzslo_mspnt(:,:) = 0.d0
         soildepth_idx_deb(:) = 0.d0
         vo_total(:) = 0.d0
         vo_total_river(:) = 0.d0
         debris_total = 0.d0
         wood_total = 0.d0
         kkk1 = 0
         hki_total = 0.d0
         vo_total_l(:) = 0.d0
         LS_idx(:) = 0.d0
         LS(:,:) = 0.d0
!---added 20240424         
         debri_sup_sum(:) = 0.d0
         debri_sup_sum_di(:,:) = 0.d0
         debri_sup_sum_ij(:,:) = 0.d0
         Qg(:) = 0.d0
        water_v_cell(:) = 0.d0
         hki_g(:) = 0.
         hki_g_2d(:,:) = 0.
         n_link_depth(:) = 0
!---added for driftwood (2022/3/1)
         cw(:) = 0.d0
         qw(:) = 0.d0
         vw(:) = 0.d0
         qwsum(:) = 0.d0
         hki_area(:) = 0.d0
         vw_idx(:) = 0.d0
         qwsum_total(:) = 0.d0
         qwsum_total_idx(:) = 0.d0
         qwsum_2d(:,:) = 0.d0
         vw2d(:,:) = 0.d0
         dmean_out(:,:) = 0.d0
         cw_idx(:) = 0.d0
         cw2d(:,:) = 0.d0
!---added for sedput (20251002)
         put_depth(:) = 0.d0
         put_depth_remain(:) = 0.d0
         l_put(:) = 0
         depth_idx_ini(:) = 0.d0

	do k = 1, riv_count
      zb_riv0_idx(k) = zb_riv_idx(k)
      zb_roc_idx(k) = zb_riv_idx(k) - perosion
      if(sed_switch==2)then!modified for river width adjustment; 20240422
       wl4wid_expan(k) = chan_capa_decre_ratio*(depth_idx(k)+height_idx(k))+zb_riv0_idx(k) 
      endif
    enddo
! sediment layer setting
    Ed = Em
    do i = 1, ny
	   do j = 1, nx
		Emb(i,j) = Emc
		Et(i,j) = Em
	   end do
	end do
    do k = 1, riv_count
	  Nb_idx(k) = Nl
      ! define downstream k
      if(domain_riv_idx(k) == 2 ) downstream_k = k 
	enddo
    do l = 1, link_count
      Nb_lin(l) = Nl
    enddo
!need to modify
if( dam_switch .eq. 1) then!for dam
    do i = 1, dam_num
    if(outswitch_dam(i)>0)then
        open(900+i, file = trim(outfile_dam(i)))
        if(sed_switch > 0)then
            write(900+i,'(a)')" Time   Water_and_sediment_volume_in_dam    Total_sedimentation_volume    Total_bedload_sedimentation_volume    Total_suspenedload_sedimentation_volume    Inflow_discharge    Outflow_discharge     bedload_outflow    suspended_outflow   suspended_concentration      Deposited_volume_of_each_size_class "
        else
            write(900+i,'(a)')" Time   Water_volume_in_dam   Inflow_discharge    Outflow_discharge" 
        endif 
    endif    
    enddo      
endif 
!time output !---for yazagyodam
      if(hydro_switch .eq. 1 .and. sed_switch.ge.1)then
        do i = 1, maxhydro
            open(1100+i, file = trim(time_output) // trim(location_name(i)) //".txt") 
            if(sed_switch==1)then
            write(1100+i, '(a)')' time     qr     qb     qs        zb          slope         Ust              Fr       dzb   sumdzb      sumqsb      sumqss       Esim                         Dsim                                          Csi       fm'!----for check 2021/5/27    
            else !modified 20240125
                if(j_drf==1)then
                write(1100+i, '(a)')' time     qr     qb     qs   slo_sedi_supply_volume_to_unit_channel  sedi_over_flow  Num_of_unstable_mesh   dmean    cw   zb          slope         Ust              Fr       dzb   sumdzb      sumqsb      sumqss                 Esim                         Dsim                                          Csi       fm               slope_erosion_supply_di'!----for check 2021/5/27
                else
                write(1100+i, '(a)')' time     qr     qb     qs   slo_sedi_supply_volume_to_unit_channel  sedi_over_flow    dmean    zb          slope         Ust              Fr       dzb   sumdzb      sumqsb      sumqss                 Esim                         Dsim                                          Csi       fm               slope_erosion_supply_di'!----for check 2021/5/27
                endif
            endif
        enddo
    endif  
!RSR until here ---------------------------------------------       

! gw initial setting
    if (init_gw_switch .ne. 1) then
        call hg_init(hg_idx)
        call sub_slo_idx2ij(hg_idx, hg)
    end if

! initial storage calculation
    rain_sum = 0.d0
    aevp_sum = 0.d0
    pevp_sum = 0.d0
    sout = 0.d0
    si = 0.d0
    sg = 0.d0
    if (outswitch_storage == 1) open (1000, file=outfile_storage)
    call storage_calc(hs, hr, hg, ss, sr, si, sg)
    sinit = ss + sr + si + sg
    if (outswitch_storage == 1) write (1000, '(1000e15.7)') rain_sum, pevp_sum, aevp_sum, sout, ss + sr + si + sg, &
        (rain_sum - aevp_sum - sout - (ss + sr + si + sg) + sinit), ss, sr, si, sg

! reading rainfall data
    open (11, file=rainfile, status='old')

    tt = 0
    do
        read (11, *, iostat=ios) t, nx_rain, ny_rain
        do i = 1, ny_rain
            read (11, *, iostat=ios) (rdummy, j=1, nx_rain)
        end do
        if (ios .lt. 0) exit
        if (ios /= 0) exit
        tt = tt + 1

    end do
    tt_max_rain = tt - 1

    allocate (t_rain(0:tt_max_rain), qp(0:tt_max_rain, ny_rain, nx_rain), qp_t(ny, nx))
    allocate (sum_qp_t(ny, nx))
    sum_qp_t = 0.0d0
    rewind (11)

    qp = 0.d0
    do tt = 0, tt_max_rain
        read (11, *) t_rain(tt), nx_rain, ny_rain
        do i = 1, ny_rain
            read (11, *) (qp(tt, i, j), j=1, nx_rain)
        end do
    end do
! unit convert from (mm/h) to (m/s)
    qp = qp/3600.d0/1000.d0

    do j = 1, nx
        rain_j(j) = int((xllcorner + (dble(j) - 0.5d0)*cellsize - xllcorner_rain)/cellsize_rain_x) + 1
    end do
    do i = 1, ny
        rain_i(i) = ny_rain - int((yllcorner + (dble(ny) - dble(i) + 0.5d0)*cellsize - yllcorner_rain)/cellsize_rain_y)
    end do
    close (11)

    write (*, *) "done: reading rain file"

! reading evp data
    if (evp_switch .ne. 0) then
        open (11, file=evpfile, status='old')

        tt = 0
        do
            read (11, *, iostat=ios) t, nx_evp, ny_evp
            do i = 1, ny_evp
                read (11, *, iostat=ios) (rdummy, j=1, nx_evp)
            end do
            if (ios .lt. 0) exit
            tt = tt + 1
        end do
        tt_max_evp = tt - 1

        allocate (t_evp(0:tt_max_evp), qe(0:tt_max_evp, ny_evp, nx_evp), qe_t(ny, nx))
        rewind (11)

        write (*, *) "done: reading evp file"

        qe = 0.d0
        do tt = 0, tt_max_evp
            read (11, *) t_evp(tt), nx_evp, ny_evp
            do i = 1, ny_evp
                read (11, *) (qe(tt, i, j), j=1, nx_evp)
            end do
        end do
        ! unit convert from (mm/h) to (m/s)
        qe = qe/3600.d0/1000.d0

        do j = 1, nx
            evp_j(j) = int((xllcorner + (dble(j) - 0.5d0)*cellsize - xllcorner_evp)/cellsize_evp_x) + 1
        end do
        do i = 1, ny
            evp_i(i) = ny_evp - int((yllcorner + (dble(ny) - dble(i) + 0.5d0)*cellsize - yllcorner_evp)/cellsize_evp_y)
        end do
        close (11)
    end if

! For TSAS Output (Initial Condition)
    call sub_slo_ij2idx(hs, hs_idx)
    call sub_riv_ij2idx(hr, hr_idx)
!call RRI_TSAS(0, hs_idx, hr_idx, hg_idx, qs_ave_idx, &
!              qr_ave_idx, qg_ave_idx, qp_t_idx)

!-------added for RSR model 20240724    !RSR initial setting
call sub_riv_ij2idx( depth,  depth_idx )
call sub_riv_ij2idx( depth,  depth_idx_ini )   !20260221
if (sed_switch.ne.0) then
    call sub_sed_ij2idx( sed, sed_idx )
    call sub_riv_ij2idx( Emb, Emb_idx )
    call sub_riv_ij2idx( Et, Et_idx )
    !---added by Qin
    sumdzb_idx(:)  = 0.d0
    sumqsb_idx(:)  = 0.d0
    sumqss_idx(:)  = 0.d0
    slo_sedi_cal_duration = 0.d0 !added 20231204
    if (sed_switch==2)then
        call sub_sed_idxlin( sed_idx, sed_lin )
        call sub_riv_idxlin( Emb_idx, Emb_lin )  !this subroutine was moved to RRI_sed2.f90 !20240724
        call sub_riv_idxlin( Et_idx, Et_lin )
        do l = 1, link_count
            Nb_lin(l) = Nb_idx(link_idx_k(l))
        end do
        sumdzb_lin(:) = 0.d0 !added by Qin
        qsb_total = 0.d0
        qss_total = 0.d0
        qsw_total = 0.d0
        if(j_drf==1) qwood_total = 0.d0
        if(dam_switch.eq.1)then !added by Qin
            dam_sedi_total(:) = 0.d0
            dam_sedi_total_b(:) = 0.d0
            dam_sedi_total_s(:) = 0.d0
            dam_sedi_total_w(:) = 0.d0
            dam_sedi_totalV(:) = 0.d0
            dam_sedi_qsi(:,:) = 0.d0

            dam_outflow_qb(:) =0.d0
            dam_outflow_qs(:) =0.d0
            dam_outflow_cs(:) =0.d0
        endif
            do k = 1, slo_count
                i = slo_idx2i(k)
                j = slo_idx2j(k)
        !----sediment size distribution allocation at slope        
                do m = 1, Np
                    fmslo(k,m) = sed(i,j)%fms(m)
                    	if(fmslo(k,m)>1.d0) then
	                    !write(*,*) 'initial fmslo >1 ',i,j, k, m, dsi(m), fmslo(k,m),sed(i,j)%fms(m),sed(i,j)%fms(1),sed(i,j)%fms(2),sed(i,j)%fms(10)!modified 20240509
		                write(*,*) 'initial fmslo >1 ',j,ny+1-i, k, m, dsi(m), fmslo(k,m),sed(i,j)%fms(m),sed(i,j)%fms(1),sed(i,j)%fms(2),sed(i,j)%fms(10)!convert i,j in IRIC GUI 20250505
                        stop
	                    endif
                end do
             !---modified for slope erosion
             if(slope_ero_switch>0)then
             B_chan(k) = len_slo_1d_idx(k)*B_gully_r(sdt_s(i,j))               
             D_chan(k) = D_gully(sdt_s(i,j))
             I_chan(k) = slo_grad(k)
             area_chan(k) = B_chan(k)*l_chan(k)
             slo_sur_zb(k) = zb_slo_idx(k)+soildepth_idx(k)
             ! added 20250618
             infil_w_depth(k) = infil_s_depth(sdt_s(i,j))*gammaa_idx(k)
             endif
            end do
            deallocate(sdt_s) 
        !----Debris flow prepareation
        if(debris_switch .eq.1 )then
            call debris_setting    !this subroutine was moved to RRI_sed2.f90 !20240724
            LS_num = 0
            call cal_Landslide ( hs_idx )  !In some locations, landslide occurs even in the initial timestep. This subroutine is to prevent that.
            LS_num = 0
        end if
        !----Sediment_put
        if(j_sedput .eq.1 )then
        nn = ubound(iput, 1)
        allocate(iput2(nn), jput2(nn))
        iput2(:) = 0
        jput2(:) = 0
            do ii = 1, nn   ! Caution i and j are opposite in case of iRIC interface
                if(iput(ii) == 0 .and. jput(ii) == 0) cycle
                if(iput(ii) < 1 .or. iput(ii) > nx .or. &
                   jput(ii) < 1 .or. jput(ii) > ny)then
                    write(*,*) 'Invalid sediment supply location:', ii, iput(ii), jput(ii)
                    stop 'Sediment supply location is outside the calculation grid'
                end if
                iput2(ii) = ny - jput(ii) + 1
                jput2(ii) = iput(ii)
                kk = riv_ij2idx(iput2(ii),jput2(ii))
                if(kk <= 0)then
                    write(*,*) 'Sediment supply location is not a river cell:', ii, iput(ii), jput(ii)
                    stop 'Sediment supply must be specified at a river cell'
                end if
                write(*,*)iput2(ii),jput2(ii),kk
            end do

            do k = 1, riv_count
                l = link_to_riv(k)
                i = riv_idx2i(k) 
                j = riv_idx2j(k)                
                do ii = 1, nn
                    if(iput2(ii) == 0 .or. jput2(ii) == 0) cycle
                    if(i == iput2(ii).and.j == jput2(ii))then
                        l_put(l) = 1
                        ! put_v is bulk sediment volume including voids [m3].
                        ! Store the remaining supply as an equivalent bed thickness [m].
                        put_depth_remain(l) = put_depth_remain(l) + put_v(ii)/area_lin(l)
                        do m = 1, Np
                            fm_sedput(l,m) = sed(i,j)%fms(m)
                        end do

                        if(j_woodput==1)then   !added 20251031
                            vw(l) = vw(l) + put_w(ii)/area_lin(l)                            
                        end if
        
                        write(*,*)'sedput, i=, j=, k=, l=, put_v, put_depth, fm_put(1), put_wood', i, j, k, l, put_v(ii), put_depth_remain(l),fm_sedput(l,1), put_w(ii)
                    end if
                end do
            end do
            !added 20260318
            do k = 1, riv_count
                l = link_to_riv(k)
                if(l_put(l)==1)then
                    zb_roc_idx(k) = zb_roc_idx(k) - put_depth_remain(l)
                end if
            end do
        end if
    endif
endif
!-------sedout setting-----
if (sd_out_switch == 1)then
	i = 1
	do 
		if(itar2(i)==0.and.jtar2(i)==0)then
			exit
		else
            itar(i) = ny - jtar2(i) + 1
            jtar(i) = itar2(i)
            k = riv_ij2idx(itar(i),jtar(i))
            l = link_to_riv(k)
            write(*,*) 'i_out_location, j_out_location, k, l', itar(i) , jtar(i), k, l
			i = i + 1
		end if
	end do
end if

!-------RSR until here

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! STEP 3: CALCULATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    rain_sum = 0.d0
    aevp_sum = 0.d0
    sout = 0.d0

! output timestep
    out_dt = dble(maxt)/dble(outnum)
    out_dt = max(1.d0, out_dt)
    out_next = nint(out_dt)
    tt = 0

! rainfall for initial value
    do i = 1, ny
        if (rain_i(i) .lt. 1 .or. rain_i(i) .gt. ny_rain) cycle
        do j = 1, nx
            if (rain_j(j) .lt. 1 .or. rain_j(j) .gt. nx_rain) cycle
            qp_t(i, j) = qp(0, rain_i(i), rain_j(j))
        end do
    end do

!初期値出劁E    call iric_cgns_output_result(sum_qp_t, qp_t, hs, hr, hg, qr_ave, qs_ave, qg_ave,qsb, qss, sumdzb, sumqsb, sumqss) !this line was modified for RSR 20240724 for slope erosion

    do t = 1, maxt

        if (mod(t, 1) .eq. 0) write (*, *) t, "/", maxt

        !******* Comunication with iRIC GUI  ******************************

        call iric_check_cancel(ierr)
        if (ierr == 1) then
            write (*, *) "Solver is stopped because the STOP button was clicked."
            call iric_cgns_close()
            stop
        end if

        !******* RIVER CALCULATION ******************************
        if (riv_thresh .lt. 0) go to 2

        ! from time = (t - 1) * dt to t * dt
        time = (t - 1)*dt ! (current time)
        ! time step is initially set to be "dt_riv"
        ddt = dt_riv
        ddt_chk_riv = dt_riv

        qr_ave = 0.d0
        qr_ave_idx = 0.d0
        if (dam_switch .eq. 1) dam_vol_temp(:) = 0.d0

        ! hr -> hr_idx
        ! Memo: riv_ij2idx must be here.
        ! hr_idx cannot be replaced within the following do loop.
        call sub_riv_ij2idx(hr, hr_idx)

!-------added for RSR model 20240724
        do k = 1, riv_count
        		if(isnan(hr_idx(k))) then !check 2021/6/3
			write(*,*) hr_idx(k)
            stop 'hr is nan check1'
            endif  
            hr_idx2(k) = hr_idx(k)
        enddo
!-------RSR until here

        ! from hr_idx -> vr_idx
        do k = 1, riv_count
            call hr2vr(hr_idx(k), k, vr_idx(k))
        end do

        do

            ! "time + ddt" should be less than "t * dt"
            if (time + ddt .gt. t*dt) ddt = t*dt - time

            ! boundary condition for river (water depth boundary)
            if (bound_riv_wlev_switch .ge. 1) then
                itemp = -1
                do jtemp = 1, tt_max_bound_riv_wlev
                    if (t_bound_riv_wlev(jtemp - 1) .lt. (time + ddt) .and. (time + ddt) .le. t_bound_riv_wlev(jtemp)) itemp = jtemp
                end do
                do k = 1, riv_count
                    if (bound_riv_wlev_idx(itemp, k) .le. -100.0) cycle ! not boundary
                    hr_idx(k) = bound_riv_wlev_idx(itemp, k)
                    call hr2vr(hr_idx(k), k, vr_idx(k))
                end do
            end if

1           continue
            qr_ave_temp_idx(:) = 0.d0

            ! Adaptive Runge-Kutta
            ! (1)
            call funcr(vr_idx, fr, qr_idx)
            vr_temp = vr_idx + b21*ddt*fr
            where (vr_temp .lt. 0) vr_temp = 0.d0
            qr_ave_temp_idx = qr_ave_temp_idx + qr_idx*ddt

            ! (2)
            call funcr(vr_temp, kr2, qr_idx)
            vr_temp = vr_idx + ddt*(b31*fr + b32*kr2)
            where (vr_temp .lt. 0) vr_temp = 0.d0
            qr_ave_temp_idx = qr_ave_temp_idx + qr_idx*ddt

            ! (3)
            call funcr(vr_temp, kr3, qr_idx)
            vr_temp = vr_idx + ddt*(b41*fr + b42*kr2 + b43*kr3)
            where (vr_temp .lt. 0) vr_temp = 0.d0
            qr_ave_temp_idx = qr_ave_temp_idx + qr_idx*ddt

            ! (4)
            call funcr(vr_temp, kr4, qr_idx)
            vr_temp = vr_idx + ddt*(b51*fr + b52*kr2 + b53*kr3 + b54*kr4)
            where (vr_temp .lt. 0) vr_temp = 0.d0
            qr_ave_temp_idx = qr_ave_temp_idx + qr_idx*ddt

            ! (5)
            call funcr(vr_temp, kr5, qr_idx)
            vr_temp = vr_idx + ddt*(b61*fr + b62*kr2 + b63*kr3 + b64*kr4 + b65*kr5)
            where (vr_temp .lt. 0) vr_temp = 0.d0
            qr_ave_temp_idx = qr_ave_temp_idx + qr_idx*ddt

            ! (6)
            call funcr(vr_temp, kr6, qr_idx)
            vr_temp = vr_idx + ddt*(c1*fr + c3*kr3 + c4*kr4 + c6*kr6)
            where (vr_temp .lt. 0) vr_temp = 0.d0
            qr_ave_temp_idx = qr_ave_temp_idx + qr_idx*ddt

            ! (e)
            vr_err = ddt*(dc1*fr + dc3*kr3 + dc4*kr4 + dc5*kr5 + dc6*kr6)

            hr_err(:) = vr_err(:)/(area*area_ratio_idx(:))

            ! error evaluation
            where (domain_riv_idx .eq. 0) hr_err = 0
            errmax = maxval(hr_err)/eps
            !if(ddt.le. 4.*ddt_min_riv .and. errmax*eps.le.0.02 ) errmax= 1.d0 !modified 20231208   !tentative disable RSR 20240724
            if (errmax .gt. 1.d0 .and. ddt .ge. ddt_min_riv) then
                ! try smaller ddt
                ddt = max(safety*ddt*(errmax**pshrnk), 0.5d0*ddt)
                ddt_chk_riv = ddt
                write (*, *) "shrink (riv): ", ddt, errmax, maxloc(vr_err)                
                if (ddt .eq. 0) stop 'stepsize underflow'
                if (dam_switch .eq. 1) dam_vol_temp(:) = 0.d0
                go to 1
            else
                ! "time + ddt" should be less than "t * dt"
                if (time + ddt .gt. t*dt) ddt = t*dt - time
                time = time + ddt
                vr_idx = vr_temp
                qr_ave_idx = qr_ave_idx + qr_ave_temp_idx
            end if
            if (time .ge. t*dt) exit ! finish for this timestep
        end do
        qr_ave_idx = qr_ave_idx/dble(dt)/6.d0

        do k = 1, riv_count
            call vr2hr(vr_idx(k), k, hr_idx(k))
            if(isnan(hr_idx(k))) then !check 2021/6/3
			write(*,*) hr_idx(k),vr_idx(k),qr_ave_idx(k)
            stop 'hr is nan check2'
            endif  
        end do

!-------added for RSR model 20240724
        do k = 1, riv_count
            if(damflg(k)>0 .and. hr_idx(k) .ge. 0.95*depth_idx(k)) hr_idx(k) = 0.95 * depth_idx(k)
        end do
        call sub_riv_ij2idx(zb_riv,zb_riv_idx)!return the bed elevation to the value of present time step; modified 20240419

        if(sed_switch.ne.0)then !note 20250312
           do k = 1, riv_count
			kk = down_riv_idx(k)
			kk1 = down_riv_idx(kk)

			if(domain_riv_idx(kk1).eq.2) then
			 hr_idx2(kk) = hr_idx2(k)
			 hr_idx2(kk1) = hr_idx2(k)

			 hr_idx(kk) = hr_idx(k)
 			 hr_idx(kk1) = hr_idx(k)

			 qr_ave_idx(kk) = qr_ave_idx(k)
			 qr_ave_idx(kk1) = qr_ave_idx(k)
			endif			

           enddo
        endif
!-------RSR until here

        ! hr_idx -> hr, qr_ave_idx -> qr_ave
        call sub_riv_idx2ij(hr_idx, hr)
        call sub_riv_idx2ij(qr_ave_idx, qr_ave)

        if (dam_switch .eq. 1) call dam_checkstate(qr_ave)

!******* SEDIMENT CALCULATION RSR model main program ****** 20240724 moved 
if(t*dt.ge.t_beddeform_start-0.00001.and.sed_switch .ne. 0)then
!*************************************************************************
 
if(sed_switch == 1) then         

	      dlambda = 1./(1.- lambda)

          time = (t - 1) * dt ! (current time)
          ddt = dt/real(iidt)

	call det_rivebedslope

do 
       do k = 1, riv_count
         dzb_temp(k) = 0.d0
        do m = 1, Np
         sed_idx(k)%dzbpr(m) = 0.d0
        end do
       end do

         ! "time + ddt" should be less than "t * dt"
         if(time + ddt .gt. t * dt ) ddt = t * dt - time

	! Determine change in bed elevation 'zb_riv_idx'
        call funcd( sed_idx, hr_idx, hr_idx2,hr_idxa,  qr_ave_idx, ust_idx, qsb_idx, qss_idx, qsw_idx, t,area_idx,water_v_idx )
	
  do k = 1, riv_count
    if(hr_idx(k).le.0.07) then
	 do m = 1, Np
       sed_idx(k)%ffd(m) = 0.d0
     enddo
    endif
	 
	 do m = 1, Np
	  sed_idx(k)%dzbpr(m) = - ddt * dlambda * sed_idx(k)%ffd(m)
	  dzb_temp(k) = dzb_temp(k) + sed_idx(k)%dzbpr(m)
	 enddo
!-----added by Qin
     sumqsb_idx(k) = sumqsb_idx(k) + qsb_idx(k) * ddt 
     sumqss_idx(k) = sumqss_idx(k) + qss_idx(k) * ddt 
  enddo		
		
	 call mean_diameter(dzb_temp, sed_idx, qsb_idx,qss_idx )
    !$omp parallel do 
      do k = 1, riv_count
       zb_riv_idx(k) = zb_riv_idx(k) + dzb_temp(k)
       sumdzb_idx(k) = sumdzb_idx(k) + dzb_temp(k)!----added by Qin
      !modified 20240419 for setting depth_idx2 as the channel depth of previous time step 
       depth_idx(k) = depth_idx(k) - dzb_temp(k) !---addedby Qin 20230924
      enddo
      
         ! "time + ddt" should be less than "t * dt"
            if(time + ddt .gt. t * dt ) ddt = t * dt - time
              time = time + ddt
            if(time.ge.t * dt) exit ! finish for this timestep

enddo 

do k = 1, riv_count
!---modified 2021/6/25
if(qsb_idx(k).le.0.d0) qsb_idx(k) = 0.d0
if(qss_idx(k).le.0.d0) qss_idx(k) = 0.d0
if(qsw_idx(k).le.0.d0) qsw_idx(k) = 0.d0
enddo

		 ! pause'sub_riv_idx2ij'
         ! zb_riv_idx -> zb_riv
         ! qsb_idx -> qsb, qss_idx -> qss, qsw_idx -> qsw
         ! ust_idx -> ust, sed_idx -> sed
             call sub_riv_idx2ij( depth_idx, depth )!added by Qin 20230924; modified 20240419
             call sub_riv_idx2ij( zb_riv_idx, zb_riv )
             call sub_riv_idx2ij( qsb_idx, qsb )
             call sub_riv_idx2ij( qss_idx, qss )
             call sub_riv_idx2ij( qsw_idx, qsw )
             call sub_riv_idx2ij( ust_idx, ust )
		     call sub_sed_idx2ij( sed_idx, sed )
		     call sub_riv_idx2ij( Emb_idx, Emb )
             call sub_riv_idx2ij(sumqsb_idx,sumqsb) !added by Qin 2021/07/28
             call sub_riv_idx2ij(sumqss_idx,sumqss) !added by Qin 2021/07/28
             call sub_riv_idx2ij(sumdzb_idx,sumdzb)
		     call sub_riv_idx2ij( Et_idx, Et ) !reviesed by Qin 20210924

!--------------------------
!   Added by Harada
!--------------------------
elseif(sed_switch == 2) then

    dlambda = 1./(1.- lambda)

    time = (t - 1) * dt ! (current time)
    ddt = dt/real(iidt)

!    write(*,*) 'sediment computation   t=', time

!added by Qin 2021/08/02
    do 
        do l = 1, link_count
            dzb_temp_lin(l) = 0.d0  
         do m = 1, Np
          sed_lin(l)%dzbpr(m) = 0.d0
         end do    
        end do

        call det_rivebedslope2(hr_idx, hr_lin,water_v_lin, hr_idxa, qr_ave_idx)

        if(time + ddt .gt. t * dt ) ddt = t * dt - time

        call funcd2( sed_lin, hr_idx, hr_idx2,hr_idxa, hr_lin, qr_ave_idx, ust_idx, ust_lin, qsb_lin, qss_lin, qsw_idx, qsw_lin, t,water_v_lin, dzb_temp_lin,qss_b,sumdzb_lin)

        call sub_riv_linidx(dzb_temp_lin, dzb_temp)
        call sub_riv_linidx(Emb_lin, Emb_idx)
        call sub_riv_linidx(qsb_lin, qsb_idx)
        call sub_riv_linidx(qss_lin, qss_idx)
        call sub_riv_linidx(qsw_lin, qsw_idx)
        call sub_riv_linidx(ust_lin, ust_idx)
        call sub_riv_linidx(ss_lin, ssc_idx) !added 20250405

        !---Mean diameter computation on link variables
        call mean_diameter_lin(dzb_temp_lin, sed_lin, qsb_lin, qss_lin)
!added 20240422
        call sub_riv_linidx(width_lin, width_idx)   

        !---change riverbed elevation for all river cells
!!$omp parallel do      
        do k = 1, riv_count
            if(damflg(k).gt.0) dzb_temp(k) = 0.d0 !---- There is no bed elevation change of dam cell; added by Qin         
            if(domain_riv_idx(k).eq.2 .and. dzb_temp(k).gt.0.d0) dzb_temp(k) = 0.d0
            zb_riv_idx(k) = zb_riv_idx(k) + dzb_temp(k)
            !modified 20240419 for setting depth_idx2 as the channel depth of previous time step 
            depth_idx(k) = depth_idx(k) - dzb_temp(k) !---addedby Qin 20230924        
            sumdzb_idx(k) = sumdzb_idx(k) + dzb_temp(k)
            sumqsb_idx(k) = sumqsb_idx(k) + qsb_idx(k)*ddt
            sumqss_idx(k) = sumqss_idx(k) + qss_idx(k)*ddt
            if(isnan(sumdzb_idx(k))) then 
			    write(*,*) k, sumdzb_idx(k), dzb_temp(k)
			stop "sumdzb_idx(k) is nan"
		endif
        enddo

!!$omp parallel do private(k,i,j)    
        do l = 1, link_count!---- There is no bed elevation change of dam cell; added by Qin
            k = link_idx_k(l)
            sumdzb_lin(l) = sumdzb_lin(l) + dzb_temp(k)
            !added 20240424
            if (debris_switch >0) then !modified 20250216
             i = riv_idx2i(k) 
             j = riv_idx2j(k)           
             debri_sup_sum_ij(i,j) = debri_sup_sum(l)
            endif
        end do   

        if(time + ddt .gt. t * dt ) ddt = t * dt - time
        time = time + ddt
        if(time.ge.t * dt) exit ! finish for this timestep		

    end do

    !----link -> 1D sync for output and 2D conversion
    call sub_sed_linidx(sed_lin, sed_idx)
    call sub_riv_linidx(Et_lin, Et_idx)
    do k = 1, riv_count
        Nb_idx(k) = Nb_lin(link_to_riv(k))
    end do

!$omp parallel do 
    do k = 1, riv_count       
        if(qsb_idx(k).le.0.d0) qsb_idx(k) = 0.d0
        if(qss_idx(k).le.0.d0) qss_idx(k) = 0.d0
        if(qsw_idx(k).le.0.d0) qsw_idx(k) = 0.d0      
        !added 20240422
        area_ratio_idx(k) = width_idx(k)/len_riv_idx(k)
    enddo

    !----1D -> 2D 
    call sub_riv_idx2ij( depth_idx, depth ) !added by Qin 20230924
    call sub_riv_idx2ij( zb_riv_idx, zb_riv )
    call sub_riv_idx2ij( qsb_idx, qsb )
    call sub_riv_idx2ij( qss_idx, qss )
    call sub_riv_idx2ij( qsw_idx, qsw )
    call sub_riv_idx2ij( ust_idx, ust )
    call sub_sed_idx2ij( sed_idx, sed )
    call sub_riv_idx2ij( Emb_idx, Emb )
    call sub_riv_idx2ij( Et_idx, Et )
    call sub_riv_idx2ij(ssc_idx, ssc_ij) !added 20250405    

    call sub_riv_idx2ij(sumqsb_idx,sumqsb) !added by Qin 2021/07/28
    call sub_riv_idx2ij(sumqss_idx,sumqss) !added by Qin 2021/07/28
    call sub_riv_idx2ij(sumdzb_idx,sumdzb)
 !added 20240422       
    call sub_riv_idx2ij(width_idx,width)
    call sub_riv_idx2ij(area_ratio_idx, area_ratio)
    call sub_riv_idx2ij(hr_idx,hr)    
endif

    if(sed_switch == 2) then
       if(detail_console==1)then
        write(*,'(a)') '    l       k   n_link_depth    put_depth_remain(l)   depth(k)  hr(k)   hr(l)  qsb   qss    Emb     sumqb   sumqs     Ust      zb     sumdzb  link_0th  slo(deg)  ini_slo  width   depth   link_len  Slo_sed_sup(m3)   Inun_sed(m3)  Deb_sup(m3) Deb_remai(m3)'
        do l = 1, link_count
            k = link_idx_k(l)
            write(*,'(i5, i9, i5, 2f10.4,5f7.4, 2e10.2, 3f8.3, i5, 2f8.3,f10.2,f10.3,f10.2,2f12.5, 2e15.2)') l, k,n_link_depth(l),put_depth_remain(l),depth_idx(k),hr_idx(k),hr_lin(l), qsb_lin(l), qss_lin(l), Emb_lin(l), sumqsb_idx(k), sumqss_idx(k) , ust_lin(l), zb_riv_idx(k), sumdzb_lin(l), link_0th_order(l),zb_riv_slope_lin(l), zb_riv_slope0_lin(l), width_lin(l),depth_idx(k)+height_idx(k),Link_len(l), slo_to_lin_sed_sum(l), lin_to_slo_sed_sum(l), debri_sup_sum(l), vo_total_l(l)
        end do
       end if
        write(*,'(a,f10.2)') 'Discharge_downstream= ', qr_ave_idx(downstream_k)
        write(*,'(a,f10.2)') 'qb_downstream= ', qsb_total
        write(*,'(a,f10.2)') 'qs_downstream= ', qss_total
        if(j_drf==1)write(*,'(a,f10.2)') 'qwood_downstream= ', qwood_total
        !convert to i,j in IRIC GUI 20250505
        tmp_idx=maxloc(sumqsb)
        write(*,'(a,f10.2,a,2i7)') 'max_sumqsb= ', maxval(sumqsb), " loc : ", tmp_idx(2) , ny+1-tmp_idx(1) 
        tmp_idx=maxloc(sumqss)
        write(*,'(a,f10.2,a,2i7)') 'max_sumqss= ', maxval(sumqss), " loc : ", tmp_idx(2) , ny+1-tmp_idx(1)  
        !write(*,'(a,f10.2,a,2i7)') 'max_sumqsb= ', maxval(sumqsb), " loc : ", maxloc(sumqsb)
        !write(*,'(a,f10.2,a,2i7)') 'max_sumqss= ', maxval(sumqss), " loc : ", maxloc(sumqss)
        tmp_idx(:)=0
        if (dam_switch .eq. 1)then
            do f = 1, dam_num
                k=dam_loc(f)
            write(*,'(a,2f10.2)') 'sediment volume in the dam =', dam_sedi_totalV(f)
            write(*,'(a,2f10.2)') 'sediment  in the dam =', dam_sedi_total(f)
            write(*,'(a,2f10.2)') 'bedload  in the dam =', dam_sedi_total_b(f)
            write(*,'(a,2f10.2)') 'suspended load  in the dam =', dam_sedi_total_s(f)
            enddo
        endif    
    else
        if (dam_switch .eq. 1)then
            do f = 1, dam_num
                k=dam_loc(f)
            write(*,'(a,2f10.2)') 'sediment volume in the dam =', dam_sedi_totalV(f)
            write(*,'(a,2f10.2)') 'sediment  in the dam =', dam_sedi_total(f)
            write(*,'(a,2f10.2)') 'bedload  in the dam =', dam_sedi_total_b(f)
            write(*,'(a,2f10.2)') 'suspended load  in the dam =', dam_sedi_total_s(f)

            enddo
        endif    
    endif  
!----------------------------------------------------------------------------
endif     !endif for (t*dt.ge.t_beddeform_start.and.sed_switch .ne. 0)

!******* Until here RSR model 20240724  ******

        !******* SLOPE CALCULATION ******************************
2       continue

        ! from time = (t - 1) * dt to t * dt
        time = (t - 1)*dt  ! (current time)
        ! time step is initially set to be "dt"
        ddt = dt
        ddt_chk_slo = dt

        qs_ave = 0.d0
        qs_ave_idx = 0.d0

        ! hs -> hs_idx
        ! Memo: slo_ij2idx must be here.
        ! hs_idx cannot be replaced within the following do loop.
        call sub_slo_ij2idx(hs, hs_idx)
        call sub_slo_ij2idx(gampt_ff, gampt_ff_idx) ! modified by T.Sayama on June 10, 2017

        do

            if (time + ddt .gt. t*dt) ddt = t*dt - time

            ! rainfall
            itemp = -1
            do jtemp = 1, tt_max_rain
                if (t_rain(jtemp - 1) .lt. (time + ddt) .and. (time + ddt) .le. t_rain(jtemp)) itemp = jtemp
            end do
            do i = 1, ny
                if (rain_i(i) .lt. 1 .or. rain_i(i) .gt. ny_rain) cycle
                do j = 1, nx
                    if (rain_j(j) .lt. 1 .or. rain_j(j) .gt. nx_rain) cycle
                    qp_t(i, j) = qp(itemp, rain_i(i), rain_j(j))
                end do
            end do
            call sub_slo_ij2idx(qp_t, qp_t_idx)

            ! boundary condition for slope (water depth boundary)
            if (bound_slo_wlev_switch .ge. 1) then
                itemp = -1
                do jtemp = 1, tt_max_bound_slo_wlev
                    if (t_bound_slo_wlev(jtemp - 1) .lt. (time + ddt) .and. (time + ddt) .le. t_bound_slo_wlev(jtemp)) itemp = jtemp
                end do
                do k = 1, slo_count
                    if (bound_slo_wlev_idx(itemp, k) .le. -100.0) cycle ! not boundary
                    hs_idx(k) = bound_slo_wlev_idx(itemp, k)
                end do
            end if


3           continue
            qs_ave_temp_idx(:, :) = 0.d0
            ! for RSR model slope sediment calculation, funcs is not modified yet 20240724 
            ! Adaptive Runge-Kutta
            ! (1)
            call funcs(hs_idx, qp_t_idx, fs, qs_idx)
            hs_temp = hs_idx + b21*ddt*fs
            where (hs_temp .lt. 0) hs_temp = 0.d0
            qs_ave_temp_idx = qs_ave_temp_idx + qs_idx*ddt

            ! (2)
            call funcs(hs_temp, qp_t_idx, ks2, qs_idx)
            hs_temp = hs_idx + ddt*(b31*fs + b32*ks2)
            where (hs_temp .lt. 0) hs_temp = 0.d0
            qs_ave_temp_idx = qs_ave_temp_idx + qs_idx*ddt

            ! (3)
            call funcs(hs_temp, qp_t_idx, ks3, qs_idx)
            hs_temp = hs_idx + ddt*(b41*fs + b42*ks2 + b43*ks3)
            where (hs_temp .lt. 0) hs_temp = 0.d0
            qs_ave_temp_idx = qs_ave_temp_idx + qs_idx*ddt

            ! (4)
            call funcs(hs_temp, qp_t_idx, ks4, qs_idx)
            hs_temp = hs_idx + ddt*(b51*fs + b52*ks2 + b53*ks3 + b54*ks4)
            where (hs_temp .lt. 0) hs_temp = 0.d0
            qs_ave_temp_idx = qs_ave_temp_idx + qs_idx*ddt

            ! (5)
            call funcs(hs_temp, qp_t_idx, ks5, qs_idx)
            hs_temp = hs_idx + ddt*(b61*fs + b62*ks2 + b63*ks3 + b64*ks4 + b65*ks5)
            where (hs_temp .lt. 0) hs_temp = 0.d0
            qs_ave_temp_idx = qs_ave_temp_idx + qs_idx*ddt

            ! (6)
            call funcs(hs_temp, qp_t_idx, ks6, qs_idx)
            hs_temp = hs_idx + ddt*(c1*fs + c3*ks3 + c4*ks4 + c6*ks6)
            where (hs_temp .lt. 0) hs_temp = 0.d0
            qs_ave_temp_idx = qs_ave_temp_idx + qs_idx*ddt

            ! (e)
            hs_err = ddt*(dc1*fs + dc3*ks3 + dc4*ks4 + dc5*ks5 + dc6*ks6)

            ! error evaluation
            where (domain_slo_idx .eq. 0) hs_err = 0.d0
            errmax = maxval(hs_err)/eps

            if (errmax .gt. 1.d0 .and. ddt .ge. ddt_min_slo) then
                ! try smaller ddt
                ddt = max(safety*ddt*(errmax**pshrnk), 0.5d0*ddt)
                ddt_chk_slo = ddt
                write (*, *) "shrink (slo): ", ddt, errmax, maxloc(hs_err)                 
                if (ddt .eq. 0) stop 'stepsize underflow'
                go to 3
            else
                ! "time + ddt" should be less than "t * dt"
                if (time + ddt .gt. t*dt) ddt = t*dt - time
                time = time + ddt
                hs_idx = hs_temp
                qs_ave_idx = qs_ave_idx + qs_ave_temp_idx
            end if  
            ! cumulative rainfall
            do i = 1, ny
                do j = 1, nx
                    if (domain(i, j) .ne. 0) rain_sum = rain_sum + dble(qp_t(i, j)*area*ddt)
                end do
            end do

            if (time .ge. t*dt) exit ! finish for this timestep
        end do
        qs_ave_idx = qs_ave_idx/dble(dt)/6.d0 ! modified on ver 1.4.1
!modified for slope erosion
                            if ((t-1)*dt.ge.t_beddeform_start-0.00001.and.sed_switch==2 .and.slo_sedi_cal_switch>0)then 
                                qsur_ave = 0.d0
                                qsur_ave_temp_idx = 0.d0
                               slo_sedi_cal_duration =dt
                                !$omp parallel do private(i,j,m)
                                    do k = 1, slo_count
                                        i = slo_idx2i(k)
                                        j = slo_idx2j(k)
                                        if( domain(i,j) == 0 ) cycle
                                        if(soildepth_idx(k).le.0.d0) soildepth_idx(k)=0.d0
                                        slo_sur_zb(k) = zb_slo_idx(k)+soildepth_idx(k)
                                        if(soildepth_idx(k).le.0.)then
                                            soildepth_idx(k)=0.
                                            da_idx(k)= 0.
                                            dm_idx(k) = 0.
                                            slo_sur_zb(k) = zb_slo_idx(k)
                                            do m = 1, Np
                                            fmslo(k,m)=0.d0
                                            enddo
                                        else
                                            da_idx(k)= (slo_sur_zb(k)-zb_slo_idx(k))*gammaa_idx(k)
                                            dm_idx(k) =  (slo_sur_zb(k)-zb_slo_idx(k)) * gammam_idx(k)
                                        endif
                                        !modified 20250618 this surface water depth is for hillslope erosion calculation; it will be recaluated as the surface water depth of RRI later
			                            if(infil_w_depth(k)>da_idx(k))then
                                            h_surf(i,j) = hs_idx(k)-da_idx(k)
                                        else
                                            h_surf(i,j) = hs_idx(k)-infil_w_depth(k) !added 20250618
                                        endif

                                        if(h_surf(i,j).le.0.d0) h_surf(i,j) = 0.d0
                                        if(isnan(h_surf(i,j)))then 
                                        write(*,*) k,h_surf(i,j),hs_idx(k),da_idx(k), dm_idx(k), zb_slo_idx(k), soildepth_idx(k), slo_sur_zb(k)
                                        stop 'surface water depth at slope is nan check3 '
                                        endif
                                    enddo    

                                        			!modified for slope erosion
                                    if(maxval(h_surf).ge.surflowdepth)then                
                                    !$omp parallel do private(i,j,kk,l,I_sur,I_sur_x,I_sur_y,lev_p,lev_n)                
		                            do k = 1, slo_count
                                        i = slo_idx2i(k)
                                        j = slo_idx2j(k)
                                     if(h_surf(i,j)>surflowdepth)then
                                        if(dif_slo_idx(k)==0)then
                                            kk = down_slo_1d_idx(k)
                                            if(kk == -1) cycle
                                ! 1-direction : kinematic wave
                                            I_sur = max( (zb_slo_idx(k) - zb_slo_idx(kk)) / dis_slo_1d_idx(k), 0.001 )
                                            qsur_ave_temp_idx(1,k) = 1./ns_slo_idx(k)*h_surf(i,j)** (5./3.)*sqrt(I_sur)*len_slo_1d_idx(k)/area
                                            qsur_ave_temp_idx(2,k) = 0.
                                        else
                                            do l = 1, lmax
                                            kk = down_slo_idx(l,k)
                                            if(kk == -1) cycle
                                            call  h2lev(hs_idx(k),k,lev_p)
                                            call  h2lev(hs_idx(kk),kk,lev_n)
                                ! diffusion wave
                                                I_sur = ((zb_slo_idx(k) + lev_p) - (zb_slo_idx(kk) + lev_n)) / dis_slo_idx(l,k)
                                !-----added to compute resultant
                                                !to avoid very small value
                                                if(I_sur >= 0.) I_sur = max(I_sur , 0.00001)
                                                if(I_sur < 0.) I_sur = min(I_sur , -0.00001)

                                                if( l.eq.1 ) then
                                                    I_sur_x = I_sur
                                                else     ! l=2
                                                    I_sur_y = I_sur
                                                end if
                                            enddo

                                            qsur_ave_temp_idx(1,k) = 1./ns_slo_idx(k)*h_surf(i,j)** (5./3.)* sqrt((I_sur_x *I_sur_x)/sqrt(I_sur_x**2.+I_sur_y**2.))*len_slo_idx(1,k)/area
                                            if(I_sur_x<0.)qsur_ave_temp_idx(1,k) = -qsur_ave_temp_idx(1,k)
                                            qsur_ave_temp_idx(2,k) = 1./ns_slo_idx(k)*h_surf(i,j)** (5./3.)* sqrt((I_sur_y *I_sur_y)/sqrt(I_sur_x**2.+I_sur_y**2.))*len_slo_idx(2,k)/area
                                            if(I_sur_y<0.)qsur_ave_temp_idx(2,k) = -qsur_ave_temp_idx(2,k)
                                        endif
                                       endif
                                    enddo   	
                                        call slo_sedi_cal ( sed_lin, hs_idx,qs_ave_idx,qsur_ave_temp_idx,slo_sedi_cal_duration) !modified 20231226  
                                    else
                !$omp parallel do private(i,j,riv_k,m,slo_Depo_sum,sum_fmslo,nn,d80)	 
                                        do k = 1, slo_count
                                            i = slo_idx2i(k)
                                            j = slo_idx2j(k)
                                            if( domain(i,j) == 0 ) cycle	
                                            riv_k=riv_ij2idx(i,j)
                            !modified 20240115				
                                        if(debris_switch.ne.0)then		
                                        !if(slo_grad(k).ge.0.25.and.zb(i,j)>100.d0 .and. riv(i,j).ne.1) cycle !modified for gofukuya river 20231101;20240121	
                                       ! if(slo_grad(k).ge.0.25) cycle !20240126
                                        !if(slo_grad(k).ge.0.36397)cycle
                                        if(hki_g(k)>0.) cycle	! no suspended sediment transportation taking place in debris flow grid cell; modified for gofukuya river 20240109	
                                        endif                    
                                         water_v_cell(k) = hs_idx(k)*area !20250505
                                        !Modified 20250503   ! turnded off sediment inundation 20250513  
                                        !if(riv(i,j)==1 .and. qrs(i,j).ge.0.d0)then
                                         !   do m =1,Np
                                          !   slo_to_lin_idx_di(riv_k,m)=slo_to_lin_idx_di(riv_k,m)+slo_ssi(k,m)*water_v_cell(k)
                                         !    slo_ssi(k,m) = 0.
                                         !    slo_Dsi(k,m)=0.
                                         !    slo_Esi(k,m) = 0.  
                                          !   enddo
                                            ! qss_slope(k)=0.
                                            ! ss_slope(k) = 0.
                                        !elseif (riv(i,j)==1 .and. qrs(i,j)<0.d0)then 
                                         !   do m =1, Np
                                         !   slo_qsi(k,m) = -1.d0*lin_to_slo_idx_di(riv_k,m)
                                          !  slo_ssi(k,m)= slo_qsi(k,m)/(sqrt(qs_ave_idx(1,k) ** 2.d0 + qs_ave_idx(2,k) ** 2.d0) * area)
                                         !   enddo     
                                        !else
                                            if(ss_slope(k)>0.d0)then
                                            slo_Depo_sum=0.d0
                                                do m = 1, Np
                                                slo_Dsi(k,m) = slo_ssi(k,m)*water_v_cell(k)/area/slo_sedi_cal_duration !revised 20231226
                                                if(slo_Dsi(k,m).le.0.d0) slo_Dsi(k,m) = 0.d0
                                                slo_qsi(k,m)  = 0.d0 
                                                slo_ssi(k,m) = 0.d0 
                                                slo_Depo_sum=slo_Depo_sum+slo_Dsi(k,m)
                                                dzslo_fp_idx(k,m) = slo_Dsi(k,m)*slo_sedi_cal_duration/(1.d0-gammaa_idx(k))+ dzslo_fp_idx(k,m) !revised 20231101
                                                Ero_slo_vol(k) =  Ero_slo_vol(k) -slo_Dsi(k,m)*slo_sedi_cal_duration/(1.d0-gammaa_idx(k))*area !negative value means deposition revised 20231101
                                                enddo                          
                                                sum_fmslo =0.d0
                                                if((slo_sur_zb(k)+slo_Depo_sum).gt.zb_slo_idx(k)) then	                        
                                                    do m = 1,Np		
                                                            !fmslo(k,m) = ((slo_sur_zb(k) -zb_slo_idx(k))*fmslo(k,m)+slo_Dsi(k,m) *dt_slo_sed/(1.d0-gammaa_idx(k)))/(slo_sur_zb(k)+slo_Dsi(k,m)-zb_slo_idx(k))
                                                        fmslo(k,m) = ((slo_sur_zb(k) -zb_slo_idx(k))*fmslo(k,m)+slo_Dsi(k,m) *slo_sedi_cal_duration/(1.d0-gammaa_idx(k)))/(slo_sur_zb(k)+slo_Dsi(k,m)-zb_slo_idx(k))   !revised 20231101
                                                            if(fmslo(k,m).lt.0.d0) fmslo(k,m) = 0.d0
                                                        sum_fmslo =sum_fmslo + fmslo(k,m)
                                                        if(isnan(fmslo(k,m)).or.fmslo(k,m)<0.d0)then
                                                        write(*,*) k,m, fmslo(k,m),dzslo_idx(k),dzslo_fp_idx(k,m),(slo_sur_zb(k)-zb_slo_idx(k))*fmslo(k,m),soildepth_idx(k)
                                                        stop "fm of slope is incorrect" 
                                                        endif
                                                    end do 	
                                                else
                                                    do m =1, Np
                                                        fmslo(k,m) = 0.d0
                                                    enddo
                                                        slo_sur_zb(k) = zb_slo_idx(k)
                                                        sum_fmslo=0.d0
                                                endif
                                                    if(slope_ero_switch==1)then	
                                                        d80=0.d0
                                                        nn = 1
                                                		if(sum_fmslo.gt.0.d0)then
                                                            do m = 1,np
                                                                d80 = d80 + fmslo(k,m)
                                                                if(d80<=0.8d0) then
                                                                nn=nn+1
                                                                if(nn>np) nn = np
                                                                endif
                                                            enddo
                                                            ns_MS(k) = (2.d0*dsi(nn))**(1./6.)/(7.66*grav**0.5)
                                                                if(ns_MS(k)<0.02)then 
                                                                    ns_MS(K)=0.02
                                                                elseif(ns_MS(k)>0.03d0)then
                                                                    ns_MS(k) = 0.03d0
                                                                endif
                                                        else    
                                                            ns_MS(K)=0.02	
                                                        endif
                                                    endif    
                                            qss_slope(k)=0.d0
                                            else !added 20240101
                                             qss_slope(k)=0.
                                             ss_slope(k) = 0.
                                             do m = 1, Np
                                             slo_ssi(k,m) = 0.
                                             slo_Dsi(k,m)=0.
                                             slo_Esi(k,m) = 0.    
                                             enddo            
                                            endif !endif for qss_slope>0.d0 
                                        !endif !20250503
                                        enddo !endo slope k  
                                    endif !end if dhs.ge.surflowdepth
                               ! 7 continue            
                                endif                                 
                            !endif                                                     

        !******* GW CALCULATION ******************************
        if (gw_switch .eq. 0) go to 6

        ! from time = (t - 1) * dt to t * dt
        time = (t - 1)*dt  ! (current time)
        ! time step is initially set to be "dt"
        ddt = dt
        ddt_chk_slo = dt

        qg_ave = 0.d0
        qg_ave_idx = 0.d0

        ! hg -> hg_idx
        ! Memo: slo_ij2idx must be here.
        ! hg_idx cannot be replaced within the following do loop.
        call sub_slo_ij2idx(hg, hg_idx)
        !call sub_slo_ij2idx( gampt_ff, gampt_ff_idx ) ! modified by T.Sayama on June 10, 2017

        ! GW Recharge
        call gw_recharge(hs_idx, gampt_ff_idx, hg_idx)

        ! GW Lose
        call gw_lose(hg_idx)

        do

            if (time + ddt .gt. t*dt) ddt = t*dt - time

5           continue
            qg_ave_temp_idx(:, :) = 0.d0

            ! Adaptive Runge-Kutta
            ! (1)
            call funcg(hg_idx, fg, qg_idx)
            hg_temp = hg_idx + b21*ddt*fg
            qg_ave_temp_idx = qg_ave_temp_idx + qg_idx*ddt

            ! (2)
            call funcg(hg_temp, kg2, qg_idx)
            hg_temp = hg_idx + ddt*(b31*fg + b32*kg2)
            qg_ave_temp_idx = qg_ave_temp_idx + qg_idx*ddt

            ! (3)
            call funcg(hg_temp, kg3, qg_idx)
            hg_temp = hg_idx + ddt*(b41*fg + b42*kg2 + b43*kg3)
            qg_ave_temp_idx = qg_ave_temp_idx + qg_idx*ddt

            ! (4)
            call funcg(hg_temp, kg4, qg_idx)
            hg_temp = hg_idx + ddt*(b51*fg + b52*kg2 + b53*kg3 + b54*kg4)
            qg_ave_temp_idx = qg_ave_temp_idx + qg_idx*ddt

            ! (5)
            call funcg(hg_temp, kg5, qg_idx)
            hg_temp = hg_idx + ddt*(b61*fg + b62*kg2 + b63*kg3 + b64*kg4 + b65*kg5)
            qg_ave_temp_idx = qg_ave_temp_idx + qg_idx*ddt

            ! (6)
            call funcg(hg_temp, kg6, qg_idx)
            hg_temp = hg_idx + ddt*(c1*fg + c3*kg3 + c4*kg4 + c6*kg6)
            qg_ave_temp_idx = qg_ave_temp_idx + qg_idx*ddt

            ! (e)
            hg_err = ddt*(dc1*fg + dc3*kg3 + dc4*kg4 + dc5*kg5 + dc6*kg6)

            ! error evaluation
            where (domain_slo_idx .eq. 0) hg_err = 0.d0
            errmax = maxval(hg_err)/eps

            if (errmax .gt. 1.d0 .and. ddt .ge. ddt_min_slo) then
                ! try smaller ddt
                ddt = max(safety*ddt*(errmax**pshrnk), 0.5d0*ddt)
                ddt_chk_slo = ddt
                write (*, *) "shrink (gw): ", ddt, errmax, maxloc(hg_err)
                if (ddt .eq. 0) stop 'stepsize underflow'
                go to 5
            else
                ! "time + ddt" should be less than "t * dt"
                if (time + ddt .gt. t*dt) ddt = t*dt - time
                time = time + ddt
                hg_idx = hg_temp
                qg_ave_idx = qg_ave_idx + qg_ave_temp_idx
            end if

            if (time .ge. t*dt) exit ! finish for this timestep
        end do
        qg_ave_idx = qg_ave_idx/dble(dt)/6.d0

        time = t*dt

        !******* GW Exfiltration ********************************
        call gw_exfilt(hs_idx, gampt_ff_idx, hg_idx)

6       continue

        !******* Evapotranspiration *****************************
        if (evp_switch .ne. 0) call evp(hs_idx, gampt_ff_idx)

        ! hs_idx -> hs
        call sub_slo_idx2ij(hs_idx, hs)
        call sub_slo_idx2ij4(qs_ave_idx, qs_ave)
        !call sub_slo_idx2ij4 (qsur_ave_temp_idx,qsur_ave)
        call sub_slo_idx2ij(hg_idx, hg)
        call sub_slo_idx2ij4(qg_ave_idx, qg_ave)
        call sub_slo_idx2ij(gampt_ff_idx, gampt_ff)

        !******* LEVEE BREAK ************************************
        !call levee_break(t, hr, hs, xllcorner, yllcorner, cellsize)

        !******* RIVER-SLOPE INTERACTIONS ***********************
        if (riv_thresh .ge. 0) call funcrs(hr, hs)
        call sub_riv_ij2idx(hr, hr_idx)
        call sub_slo_ij2idx(hs, hs_idx)

        if(t*dt.ge.t_beddeform_start-0.00001.and.sed_switch==2 ) then
            call sed_exchange(sed_lin, hs_idx) !modified for slope erosion
        end if

        !******* INFILTRATION (Green Ampt) **********************

        call infilt(hs_idx, gampt_ff_idx, gampt_f_idx)
        call sub_slo_idx2ij(hs_idx, hs)
        call sub_slo_idx2ij(gampt_ff_idx, gampt_ff)
        call sub_slo_idx2ij(gampt_f_idx, gampt_f)

        !******* SET WATER DEPTH 0 AT DOMAIN = 2 ****************
        do i = 1, ny
            do j = 1, nx
                if (domain(i, j) .eq. 2) then
                    sout = sout + hs(i, j)*area
                    hs(i, j) = 0.d0
                    if (riv(i, j) .eq. 1) then
                        call hr2vr(hr(i, j), riv_ij2idx(i, j), vr_out)
                        sout = sout + vr_out
                        hr(i, j) = 0.d0
                    end if
                end if
            end do
        end do

        ! hs -> hs_idx, hr -> hr_idx, hg -> hg_idx
        call sub_riv_ij2idx(hr, hr_idx)
        call sub_slo_ij2idx(hs, hs_idx)
        call sub_slo_ij2idx(hg, hg_idx)
        !convert to i,j in IRIC GUI 20250505
        tmp_idx=maxloc(hr)
        write (*, *) "max hr: ", maxval(hr), "loc : ", tmp_idx(2) , ny+1-tmp_idx(1) 
        tmp_idx=maxloc(hs)
        write (*, *) "max hs: ", maxval(hs), "loc : ", tmp_idx(2) , ny+1-tmp_idx(1) 
        if (gw_switch .eq. 1) then
        tmp_idx=maxloc(hg)
        write (*, *) "max hg: ", maxval(hg), "loc : ", tmp_idx(2) , ny+1-tmp_idx(1) 
        endif
        tmp_idx(:)=0
        !write (*, *) "max hr: ", maxval(hr), "loc : ", maxloc(hr)
        !write (*, *) "max hs: ", maxval(hs), "loc : ", maxloc(hs)
        !if (gw_switch .eq. 1) write (*, *) "max hg: ", maxval(hg), "loc : ", maxloc(hg)

        !***** count total rain volume *****

        sum_qp_t = sum_qp_t + qp_t*dt

!-----RSR model: Landslide, debris flow, Slope erosion, Driftwood------

        !***** Debris flow calculation (added RSR model) 20240724 *****
                     if (t.gt.T_bebris_off/dt .and. debris_end_switch==1 .and. debris_end_switch==1)then
                      debris_switch= -1 !skip debris flow computation after t>T_bebris_off
                      j_drf =0
                      write(*,*) "Landslide Debris flow and drift wood computations have been switch off at time =", T_bebris_off,"s", "t/maxT=", t,"/",maxt
                     endif
!added 20240115
                    
                      if(debris_switch==-1)then
                        write(*,*) 'Number_of_unstable_mesh=', LS_num
                        write(*,*) "debris_total = ", debris_total
                        write(*,*) "hki_total = ", hki_total
                      endif
                     if(debris_switch == 1 .and. t*dt.ge.t_beddeform_start-0.00001) then

                        !Initialize vo_total_l
		                do l = 1, link_count
			                vo_total_l(l) = 0.d0
		                end do
                        call cal_Landslide ( hs_idx )
                        n_LS = 0
!$omp parallel do private(k) reduction(+ : n_LS) 
                        do k = 1, slo_count
                            if(sf(k)>=1.0d0)then
                                 n_LS = n_LS + 1
                            end if
                        end do

                        if(tt.gt.1)then
                        !added 20240717
!                        !$omp parallel do
                            do k = 1, slo_count
                             LS(slo_idx2i(k), slo_idx2j(k))=0.
                            if(LS_idx(k)==1.) LS(slo_idx2i(k), slo_idx2j(k)) = LS_idx(k)
                            enddo
                        else
                            call sub_slo_idx2ij( LS_idx, LS )
                        endif
                        if(n_LS >= 1)then
                            call cal_mspnt ( hs_idx )
                            call sub_slo_idx2ij( dzslo_mspnt_idx, dzslo_mspnt )
                            do k = 1, riv_count
                                l = link_to_riv(k)
                                vo_total_l(l) = vo_total_l(l) + vo_total(k)  !consider remaining sediment in previous timestep
                                vo_total(k) = 0.d0
                            end do
                            if(j_drf == 1)then
                                do k = 1, riv_count
                                    l = link_to_riv(k)
                                    vw(l) = vw(l) + Wood_density*hki_area(k)/area_lin(l)
                                    wood_total = wood_total + Wood_density*hki_area(k)
                                    hki_area(k) = 0.d0
                                end do
                            end if
                        end if
                        write(*,*) "debris_total = ", debris_total
                        write(*,*) "hki_total = ", hki_total
                      if(detail_console==1)then
                        if(j_drf == 1)then
                            write(*,*) "wood_total = ", wood_total
                            write(*,'(a)') '    l    k    qr     hr      area        vo_total_l(l)         cw(l)       vw(l)      qw(l)'
                            do l = 1, link_count
                                k = link_idx_k(l)
                                write(*,'(2i5, 3f10.4, 5f14.5)') l, k, qr_ave_idx(k),hr_idxa(k), area_lin(l),vo_total_l(l),cw(l),vw(l),qw(l)
                            end do
                        end if
                      end if
                     end if
        !**********************************************

!******* renew the information of slope area & the sediment supply from slope to river channel during dt duration time 203231204
                    if(sed_switch == 2 .and. t*dt.ge.t_beddeform_start-0.00001) then
                    
                    !20240328 moved here: width modification to avoid channel closing
!                        do k = 1, riv_count
!                          l = link_to_riv(k)
!                          if(sumdzb_idx(k)>0.)then
!                            if(hr_idx(k)>0.)then
!                                if(depth_idx_ini(k)-sumdzb_idx(k)-hr_idx(k) < depth_idx_ini(k)*0.4)then   !check this 0.4 later
!                                    width_idx(k) = width_idx(k)*1.2    !check this ratio later
                                    !write(*,'(a, 3i5, 4f14.5)') "k, l, n_link_depth, modified width, hr, depth, depth_ini", k, l, n_link_depth(l), width_idx(k), hr_idx(k), depth_idx(k), depth_idx_ini(k)
!                                end if
!                            end if
!                          end if
!                          if(width_idx(k) > max_width*2. ) width_idx(k) = max_width*2.   !check this later
!                          if(hr_idx(k)<-20.)then
!                            write(*,'(a, 2i5, 10f14.5)')'l, k, width_idx(k) ,hr_idx(k),depth_idx(k),depth_idx_ini(k),ust_lin(l),slope,zb_riv_idx(k),dzb_temp(k),sumdzb_idx(k)',l, k, width_idx(k) ,hr_idx(k),depth_idx(k),depth_idx_ini(k),ust_lin(l),tan(zb_riv_slope_lin(l)*3.14159/180.),zb_riv_idx(k),dzb_temp(k),sumdzb_idx(k)
!                          end if
!                        end do

                        do l = 1, link_count    !calculate width_lin again 20240327
		                    width_lin(l) =0.
		                    k = link_idx_k(l)
                            link_cell_num(l) = 0
		                    do
			                    width_lin(l) = width_lin(l) + width_idx(k)
                                link_cell_num(l) = link_cell_num(l) + 1
			                    if(node_ups(k).ge.1.or.up_riv_idx(k,1).eq.0) then
				                    width_lin(l) = width_lin(l)/real(link_cell_num(l))
				                    area_lin(l) = width_lin(l)*link_len(l)
				                    exit
			                    end if
			                    k = up_riv_idx(k,1)
		                    end do
	                    end do

!$omp parallel do private(i,j,m)
                            do k = 1, slo_count
                                i = slo_idx2i(k)
                                j = slo_idx2j(k)
                                if(soildepth_idx(k).le.0.d0) soildepth_idx(k)=0.d0
                                slo_sur_zb(k) = zb_slo_idx(k)+soildepth_idx(k)
                                if(soildepth_idx(k).le.0.)then
                                    soildepth_idx(k)=0.
                                    da_idx(k)= 0.
                                    dm_idx(k) = 0.
                                    slo_sur_zb(k) = zb_slo_idx(k)
                                    do m =1, Np
                                        fmslo(k,m)=0.d0
                                    enddo 
                                else
!0801                                    da_idx(k)= (slo_sur_zb(k)-zb_slo_idx(k))*gammaa_idx(k)
!0801                                    dm_idx(k) =  (slo_sur_zb(k)-zb_slo_idx(k)) * gammam_idx(k)
                                endif
                                dzslo_idx(k) = -1.d0*Ero_slo_vol(k)/area !elevation change of due to sediment deposition and erosion on slope area; modified 20231105  
                                if(dzslo_idx(k)>0.d0)then
                                    sed(i,j)%dmean_slo=0.d0
                                    do m =1, Np
                                        if(dzslo_fp_idx(k,m).le.0.d0)cycle
                                        sed(i,j)%dmean_slo=sed(i,j)%dmean_slo+dzslo_fp_idx(k,m)/dzslo_idx(k)*dsi(m)
                                    enddo
                                endif
                            enddo  	   

                            slope_erosion_total= 0.d0	               

!$omp parallel do private(k) reduction(+ : slope_erosion_total)            
                        do k = 1,slo_count
                        slope_erosion_total = slope_erosion_total + Ero_slo_vol(k)
                        enddo	

!$omp parallel do private(k,m)	
                        do l = 1, link_count
                            k = link_idx_k(l)	            
                            if(water_v_lin(l).lt.0.01.or.qr_ave_idx(k).le.0.d0.or.hr_lin(l).le.0.d0) then
                                do m = 1,Np
                                    slo_vol_remain(l,m) = slo_vol_remain(l,m)+slo_to_lin_sed(l,m)
                                    slo_to_lin_sed(l,m)=0.
                                enddo
                            else 
                                do m = 1 , Np
                                    sed_lin(l)%qsi(m) = sed_lin(l)%qsi(m) + slo_vol_remain(l,m)/water_v_lin(l)*qr_ave_idx(k)+ slo_to_lin_sed(l,m)/dt !m3/s    !actual supply is for the next timestep modified 20230924
                            		if (isnan(sed_lin(l)%qsi(m)).or.sed_lin(l)%qsi(m) .lt.0.) then !check 2021/6//3
		                            	write(*,*) l, k,sed_lin(l)%qsi(m),qr_ave_idx(k),slo_to_lin_sed(l,m)/dt,slo_vol_remain(l,m)
			                            stop "qsi is not correct"
		                            endif	

                                    if(damflg(k).gt.0)then
                                        sed_lin(l)%qsi(m) = 0.d0 !modified 2021/11/29
                                        slo_s_dsum(k,m) = 0.d0
                                        slo_s_sum(k) = 0.d0
                                        sed_lin(l)%ssi(m)=0.d0
                                        slo_s_dsum(k,m) = slo_s_dsum(k,m) - slo_to_lin_sed(l,m)/area_lin(l)/dt ![m/s]all supply sediment deposit into the link at next timestep; note that the value is minus 
                                        slo_s_sum(k)=slo_s_sum(k)+slo_s_dsum(k,m)
                                    endif
                                    slo_to_lin_sed(l,m)=0.d0
                                    slo_vol_remain(l,m)=0.d0    
                                end do  
                            endif !endif for hr_lin(l)<=0                                    
                        end do                    
 !                       call sub_riv_idx2ij( dzslo_idx, dzslo )
                        call sub_slo_idx2ij( dzslo_idx, dzslo ) !modified by Qin; 2021/11/29
                        call sub_slo_idx2ij( Ero_slo_vol, eroslovol ) !modified by Qin;
                        call sub_slo_idx2ij(ss_slope, ss_slope_ij) !added 20231101
                        call sub_slo_idx2ij(hki_g, hki_g_2d) !added 20260321                        
                        !call sub_riv_linidx (slo_to_lin_sed_sum,slo_to_lin_sed_sum_idx)
                        call sub_riv_idx2ij(slo_to_lin_sed_sum_idx,slo_supply)
                        call sub_riv_idx2ij(lin_to_slo_sed_sum_idx,overflow_sed_sum)
                        !for slope erosion
                        call sub_riv_idx2ij(inflow_sedi, inflow_sedi_ij)
                        call sub_riv_idx2ij(overflow_sedi, overflow_sedi_ij)
                        if(slo_sedi_cal_switch>0)then
                         !convert to i,j in IRIC GUI 20250505
                         tmp_idx=maxloc(-1.d0*dzslo)
                         write(*,*) "max slope erosion depth: ", maxval(-1.d0*dzslo), " loc: ", tmp_idx(2) , ny+1-tmp_idx(1)   
                         tmp_idx=maxloc(eroslovol)
                         write(*,*) "max slope erosion volume: ", maxval(eroslovol), " loc : ", tmp_idx(2) , ny+1-tmp_idx(1) 
!                        write(*,*) "max slope erosion depth: ", maxval(-1.d0*dzslo), " loc(i,j) : ", maxloc(-1.d0*dzslo)                                                                      
!                        write(*,*) "max slope erosion volume: ", maxval(eroslovol), " loc(i,j) : ", maxloc(eroslovol)
                         write(*,*) "total slope erosion volume = ", slope_erosion_total
                         tmp_idx(:)=0
                        end if
                     end if

        !***** surface water level calculation (added RSR model) 20240724 *****
        do k = 1, slo_count
            i = slo_idx2i(k)
            j = slo_idx2j(k)
            !h_surf(i,j) = hs_idx(k) - da_idx(k) ! modified 20231204
            if(riv(i,j).eq.0 .and. domain(i,j).eq.1) h_surf(i,j) = hs_idx(k) - da_idx(k)
			if(h_surf(i,j).le.0.) h_surf(i,j) = 0.d0
        end do
!-------RSR until here

        !******* OUTPUT *****************************************
!$omp single
        ! For TSAS Output
        !call RRI_TSAS(t, hs_idx, hr_idx, hg_idx, qs_ave_idx, &
        !              qr_ave_idx, qg_ave_idx, qp_t_idx)

        if (hydro_switch .eq. 1 .and. mod(int(time), 3600) .eq. 0) write (1012, '(f12.2, 10000f14.5)') time, &
            (qr_ave(hydro_i(k), hydro_j(k)), k=1, maxhydro)

        ! open output files
        if (t .eq. out_next) then

            write (*, *) "OUTPUT :", t, time
            if(dam_switch == 1) call dam_write !this line was added for RSR 20240724 !need to modify

        !-------added for RSR model 20240724
        if(sed_switch.ne.0)then  !20240806 modified
            do k = 1, slo_count
                i = slo_idx2i(k)
                j = slo_idx2j(k)
                if(riv(i,j)==1)dmean_out(i,j) = sed(i,j)%dmean*1000.
            end do
            if(j_drf==1)then
                call sub_riv_linidx(vw, vw_idx)
                call sub_riv_idx2ij( vw_idx, vw2d )
                call sub_riv_linidx(cw, cw_idx)
                call sub_riv_idx2ij( cw_idx, cw2d )
                call sub_riv_linidx(qwsum_total, qwsum_total_idx)
                call sub_riv_idx2ij(qwsum_total_idx, qwsum_2d)
            end if
        end if

        !------RSR until here

            tt = tt + 1
            out_next = nint((tt + 1)*out_dt)
            call int2char(tt, t_char)

            where (domain .eq. 0) hs = -0.1d0
            if (riv_thresh .ge. 0) where (domain .eq. 0) hr = -0.1d0
            if (riv_thresh .ge. 0) where (domain .eq. 0) qr_ave = -0.1d0
            where (domain .eq. 0) gampt_ff = -0.1d0
            where (domain .eq. 0) aevp = -0.1d0
            if (evp_switch .ne. 0) where (domain .eq. 0) qe_t = -0.1d0
            where (domain .eq. 0) hg = -0.1d0

            if (outswitch_hs .eq. 1) ofile_hs = trim(outfile_hs)//trim(t_char)//".out"
            if (outswitch_hs .eq. 2) ofile_hs = trim(outfile_hs)//trim(t_char)//".bin"
            if (outswitch_hr .eq. 1) ofile_hr = trim(outfile_hr)//trim(t_char)//".out"
            if (outswitch_hr .eq. 2) ofile_hr = trim(outfile_hr)//trim(t_char)//".bin"
            if (outswitch_hg .eq. 1) ofile_hg = trim(outfile_hg)//trim(t_char)//".out"
            if (outswitch_hg .eq. 2) ofile_hg = trim(outfile_hg)//trim(t_char)//".bin"
            if (outswitch_qr .eq. 1) ofile_qr = trim(outfile_qr)//trim(t_char)//".out"
            if (outswitch_qr .eq. 2) ofile_qr = trim(outfile_qr)//trim(t_char)//".bin"
            if (outswitch_qu .eq. 1) ofile_qu = trim(outfile_qu)//trim(t_char)//".out"
            if (outswitch_qu .eq. 2) ofile_qu = trim(outfile_qu)//trim(t_char)//".bin"
            if (outswitch_qv .eq. 1) ofile_qv = trim(outfile_qv)//trim(t_char)//".out"
            if (outswitch_qv .eq. 2) ofile_qv = trim(outfile_qv)//trim(t_char)//".bin"
            if (outswitch_gu .eq. 1) ofile_gu = trim(outfile_gu)//trim(t_char)//".out"
            if (outswitch_gu .eq. 2) ofile_gu = trim(outfile_gu)//trim(t_char)//".bin"
            if (outswitch_gv .eq. 1) ofile_gv = trim(outfile_gv)//trim(t_char)//".out"
            if (outswitch_gv .eq. 2) ofile_gv = trim(outfile_gv)//trim(t_char)//".bin"
            if (outswitch_gampt_ff .eq. 1) ofile_gampt_ff = trim(outfile_gampt_ff)//trim(t_char)//".out"
            if (outswitch_gampt_ff .eq. 2) ofile_gampt_ff = trim(outfile_gampt_ff)//trim(t_char)//".bin"

            if (outswitch_hs .eq. 1) open (100, file=ofile_hs)
            if (outswitch_hr .eq. 1) open (101, file=ofile_hr)
            if (outswitch_hg .eq. 1) open (102, file=ofile_hg)
            if (outswitch_qr .eq. 1) open (103, file=ofile_qr)
            if (outswitch_qu .eq. 1) open (104, file=ofile_qu)
            if (outswitch_qv .eq. 1) open (105, file=ofile_qv)
            if (outswitch_gu .eq. 1) open (106, file=ofile_gu)
            if (outswitch_gv .eq. 1) open (107, file=ofile_gv)
            if (outswitch_gampt_ff .eq. 1) open (108, file=ofile_gampt_ff)

            if (outswitch_hs .eq. 2) open (100, file=ofile_hs, form='unformatted', access='direct', recl=nx*ny*4)
            if (outswitch_hr .eq. 2) open (101, file=ofile_hr, form='unformatted', access='direct', recl=nx*ny*4)
            if (outswitch_hg .eq. 2) open (102, file=ofile_hr, form='unformatted', access='direct', recl=nx*ny*4)
            if (outswitch_qr .eq. 2) open (103, file=ofile_qr, form='unformatted', access='direct', recl=nx*ny*4)
            if (outswitch_qu .eq. 2) open (104, file=ofile_qu, form='unformatted', access='direct', recl=nx*ny*4)
            if (outswitch_qv .eq. 2) open (105, file=ofile_qv, form='unformatted', access='direct', recl=nx*ny*4)
            if (outswitch_gu .eq. 2) open (106, file=ofile_gu, form='unformatted', access='direct', recl=nx*ny*4)
            if (outswitch_gv .eq. 2) open (107, file=ofile_gv, form='unformatted', access='direct', recl=nx*ny*4)
            if (outswitch_gampt_ff .eq. 2) open (108, file=ofile_gampt_ff, form='unformatted', access='direct', recl=nx*ny*4)

!added for sediment output 20250405
            if(outswitch_dzb ==1) then
                ofile_dzb = trim(outfile_dzb)//trim(t_char)//".out"
                open(109, file = ofile_dzb)            
            endif    
            if(outswitch_mspnt ==1) then
                ofile_mspnt = trim(outfile_mspnt)//trim(t_char)//".out"
                open(110, file = ofile_mspnt)            
            endif   
            if(outswitch_LS ==1) then
                ofile_LS = trim(outfile_LS)//trim(t_char)//".out"
                open(111, file = ofile_LS)            
            endif                   
            if(outswitch_dzslo ==1) then
                ofile_dzslo = trim(outfile_dzslo)//trim(t_char)//".out"
                open(112, file = ofile_dzslo)            
            endif                   
            if(outswitch_eroslovol ==1) then
                ofile_eroslovol = trim(outfile_eroslovol)//trim(t_char)//".out"
                open(113, file = ofile_eroslovol)         
            endif               
            if(outswitch_erosupply ==1) then
                ofile_erosupply = trim(outfile_erosupply)//trim(t_char)//".out"
                open(114, file = ofile_erosupply)            
            endif   
            
            if(outswitch_h_surf ==1) then
                ofile_h_surf = trim(outfile_h_surf)//trim(t_char)//".out"
                open(116, file = ofile_h_surf)            
            endif

            if (sd_out_switch==1) then
            !call sed_output (sed,sed_idx ,t_char)!--check 5/20
!               --> RRI_sed2.f90 line 745
            !if(sed_switch==1)then!--check 5/20
                !call sed_output2(sed,t_char, dzb_temp,sumdzb_idx,hr_idxa,qr_ave_idx,ust_idx,qsb_idx,qss_idx,sumqsb_idx,sumqss_idx)
            !elseif(sed_switch == 2)then
                call sed_output3(sed,t_char, dzb_temp,sumdzb_lin,hr_idxa,qr_ave_idx,ust_idx,qsb_idx,qss_idx,sumqsb_idx,sumqss_idx,area_lin)
            !endif 
	        endif   

            ! output (ascii)
            do i = 1, ny
                if (outswitch_hs .eq. 1) write (100, '(10000f14.5)') (hs(i, j), j=1, nx)
                if (outswitch_hr .eq. 1) write (101, '(10000f14.5)') (hr(i, j), j=1, nx)
                !if(outswitch_hr .eq. 1) write(101,'(10000f14.5)') ((hr(i, j) + zb_riv(i, j)), j = 1, nx)
                if (outswitch_hg .eq. 1) write (102, '(10000f14.5)') (hg(i, j), j=1, nx)
                if (outswitch_qr .eq. 1) write (103, '(10000f14.5)') ((qr_ave(i, j)), j=1, nx) ! [m3/s]
                if (outswitch_qu .eq. 1) write (104, '(10000f14.5)') &
                    (((qs_ave(1, i, j) + (qs_ave(3, i, j) - qs_ave(4, i, j))/2.d0)*area), j=1, nx)
                !if(outswitch_qv .eq. 1) write(105,'(10000f14.5)') &
                ! (((qs_ave(2, i, j) + (qs_ave(3, i, j) + qs_ave(4, i, j)) / 2.d0) * area), j = 1, nx)
                if (outswitch_qv .eq. 1) write (105, '(10000e14.5)') &
                    (((qs_ave(2, i, j) + (qs_ave(3, i, j) + qs_ave(4, i, j))/2.d0)*area), j=1, nx)
                if (outswitch_gu .eq. 1) write (106, '(10000f14.8)') &
                    (((qg_ave(1, i, j) + (qg_ave(3, i, j) - qg_ave(4, i, j))/2.d0)*area), j=1, nx)
                if (outswitch_gv .eq. 1) write (107, '(10000f14.5)') &
                    (((qg_ave(2, i, j) + (qg_ave(3, i, j) + qg_ave(4, i, j))/2.d0)*area), j=1, nx)
                if (outswitch_gampt_ff .eq. 1) write (108, '(10000f14.5)') (gampt_ff(i, j), j=1, nx)
 !added sediment ascii 20250405
                if (outswitch_dzb .eq. 1) write (109, '(10000f14.5)') (sumdzb(i, j), j=1, nx)  
                if (outswitch_mspnt .eq. 1) write (110, '(10000f14.5)') (dzslo_mspnt(i, j), j=1, nx) 
                if (outswitch_LS .eq. 1) write (111, '(10000f14.5)') (LS(i, j), j=1, nx)  
                if (outswitch_dzslo .eq. 1) write (112, '(10000f14.5)') (dzslo(i, j), j=1, nx) 
                if(outswitch_eroslovol .eq. 1) write (113, '(10000f14.5)') (eroslovol(i, j), j=1, nx) 
                if(outswitch_erosupply .eq. 1) write (114, '(10000f14.5)') (slo_supply(i, j), j=1, nx)                
                if(outswitch_h_surf .eq. 1) write(116,'(10000e14.5)') (h_surf(i,j), j = 1, nx)           
                 !added sediment ascii 20250405                                       
            end do
                if (outswitch_dzb > 0) close(109)
                if (outswitch_mspnt > 0) close(110)
                if (outswitch_LS> 0) close(111) 
                if (outswitch_dzslo > 0) close(112)
                if(outswitch_eroslovol > 0)close(113)
                if(outswitch_erosupply > 0) close(114)                
                if(outswitch_h_surf > 0) close(116)            

        !********TEST FILE FOR SEDIMENT CALCULATION output (txt) <-- test (./out/test/test_)
    if(outswitch_test .eq. 1) then

	 ofile_test = trim(outfile_test) // trim(t_char) // ".txt"
	  open( 115, file = ofile_test )

      if(sed_switch.eq.1)then
	    write(115,'(a)') '    k   kk       distance     zb_riv_idx          sumdz      width_idx        slope_idx        hr_idxa         ust_idx        qr_idx     qsb_idx     qss_idx      sedmean   Nb      dzb           Frn      Sumqsb_idx   Sumqss_idx      Exchange_layer_depth             transition_layer_depth      area         Csi               f_qbi                f_qsi               ft                                  fsur'
	     leng_sum = 0.0
            do k = 1, mriv_count
             kk = main_riv(k)
	         leng_sum = leng_sum + dis_riv_idx(kk)
             write(115,'(2i5,8f15.5,2es12.3,f12.5,i5,e15.5,f12.5,3es12.3,3f20.3,10000f7.4)') k, kk, leng_sum, zb_riv_idx(kk), zb_riv_idx(kk)-zb_riv0_idx(kk), width_idx(kk), tan(zb_riv_slope_idx(kk)*3.14159/180.),hr_idxa(kk), ust_idx(kk), qr_ave_idx(kk), qsb_idx(kk), qss_idx(kk),sed_idx(kk)%dmean, Nb_idx(kk), dzb_temp(kk),qr_ave_idx(kk)/width_idx(kk)/hr_idxa(kk)/sqrt(9.81*hr_idxa(kk)), sumqsb_idx(kk), sumqss_idx(kk), Emb_idx(kk),Et_idx(kk),width_idx(kk)*dis_riv_idx(kk), (sed_idx(kk)%ssi(m), m = 1, Np), (sed_idx(kk)%fqbi(m), m = 1, Np),(sed_idx(kk)%fqsi(m), m = 1, Np), (sed_idx(kk)%ft(m), m =1, Np),(sed_idx(k)%fsur(m), m =1, Np)
            enddo
      elseif(sed_switch.eq.2)then     
        !modified 20250405
        write(115,'(a)') '    l   qr     qsb      qss       hr        ust        slope        zb      dzb       sumdzb       sumqsb                 sumqss     Exchange_layer_depth             transition_layer_depth       area        length              dmean     fr            fm               Csi               f_qbi                f_qsi                  ft                            fsur                        Esi                   Dsi                 Esisum                Dsisum                       slope_erosion_supply_sum              slope_erosion_supply_di      inun_sedi_sum             inun_sedi_sum_di           total_storaged_sed_volum        total_storaged_sed_volum_di           width_expansion_rate            Total_sedi_sup_from_debri_vol             Total_sedi_sup_from_debri_vol_di'
            do l = 1, link_count
                k = link_idx_k(l)
                write(115,'(i9,6f10.5,f15.4,f15.8,f15.4,6f20.3,2f10.5,100000f20.8)') l, qr_ave_idx(k), qsb_lin(l), qss_lin(l),hr_idx(k), ust_lin(l),tan(zb_riv_slope_lin(l)*3.14159/180.),zb_riv_idx(k),dzb_temp_lin(l),sumdzb_lin(l), sumqsb_idx(k),sumqss_idx(k),Emb_idx(k),Et_idx(k),area_lin(l), Link_len(l), sed_idx(k)%dmean,qr_ave_idx(k)/width_idx(k)/hr_lin(l)/sqrt(9.81*hr_lin(l)),(sed_idx(k)%fm(m), m = 1, Np ), (sed_idx(k)%ssi(m), m = 1, Np), (sed_idx(k)%fqbi(m), m = 1, Np), (sed_idx(k)%fqsi(m),m = 1, Np),(sed_idx(k)%ft(m), m =1, Np),(sed_idx(k)%fsur(m), m =1, Np),(sed_lin(l)%Esi(m), m = 1, Np ), (sed_lin(l)%Dsi(m), m = 1, Np ),(sed_lin(l)%Esisum(m), m = 1, Np ), & 
                 (sed_lin(l)%Dsisum(m), m = 1, Np ), slo_to_lin_sed_sum(l),(slo_to_lin_sum_di(l,m), m = 1,Np),lin_to_slo_sed_sum(l),(lin_to_slo_sum_di(l,m), m = 1,Np),sed_lin(l)%stor_sed_sum, (sed_lin(l)%stor_sed_sum_i(m), m = 1, Np), width_lin(l)/ width_lin_0(l), debri_sup_sum(l), (debri_sup_sum_di(l,m), m = 1,Np)
            end do    
        write(115,'(a,f10.2)') 'Discharge_downstream= ', qr_ave_idx(downstream_k)
        write(115,'(a,f10.2)') 'qb_downstream= ', qsb_total
        write(115,'(a,f10.2)') 'qs_downstream= ', qss_total
        if(j_drf==1)write(115,'(a,f10.2)') 'qwood_downstream= ', qwood_total 
!modified 20240115    
        if(debris_switch.ne.0)write(115,'(a,i10)') 'Number_of_unstable_mesh=', LS_num    
        if(debris_switch.ne.0)write(115,'(a,f10.2)') 'debris_total= ', debris_total
        if(debris_switch.ne.0)write(115,'(a,f10.2)') 'hki_total= ', hki_total 
        if(j_drf==1)write(115,'(a,f10.2)') 'wood total= ', wood_total
        if(slope_ero_switch.gt.0 .and. outswitch_slope>0)then      
	     ofile_slope = trim(outfile_slope) // trim(t_char) // ".txt"        
         open(1151, file = ofile_slope )
         write (1151, '(a)') ' riv_k     i   j   slo_k     l            exchange_flow_from_slo_to_riv(m3/s)       inflow_sedi(m3/s)    overflow_sedi(m3/s)   Water_lev_slo(m)       Water_lev_riv(m)       SSC_of_slope(m3/m3)           SSC_of_river(m3/m3)     Total_Ele_change_slo(m)           Total_Ele_change_riv(m)'
         do  k = 1, riv_count
         !modified for test 20240311        
           l = link_to_riv(k)  
           i = riv_idx2i(k)
           j = riv_idx2j(k)
           kk = slo_ij2idx(i,j)
           write(1151, '(5i9, f10.3, 2f10.5, 2f8.2,4f15.8)') k,j, ny+1-i, kk, l, qrs(i,j)*area, inflow_sedi(k), overflow_sedi(k), hs(i,j)+zb(i,j), hr(i,j)+zb_riv(i,j), ss_slope(kk), ss_lin(l), dzslo_idx(kk), sumdzb_lin(l)
         enddo
        
            write(1151,'(a)') 'Info_for_gully'
            write(1151,'(a)') '   k  i  j  Total_ele_change(m)   Erosion_volume(m3)  Slope_Elevation(m)   Surface_flow_Dep(m)      qss_slope          ss_slope     Grid_discharge(m3/s)  Grid_cell_sur_water_vol width_(regime)(m)    Grid(=Gully)_slope    total_channel_width(m)  channel_length(m)      channel_depth(m)   flow_depth(m)     discharge_in_chanel(m3/s)    slope_of_channel       sin_alpha       cos_alpha         Roughness_channel    Fr     Esi    Dsi '
            do k = 1, slo_count
                i = slo_idx2i(k)
                j = slo_idx2j(k)
                write(1151,'(3i7,200f20.8)') k, j, ny+1-i,  dzslo_idx(k), Ero_slo_vol(k), slo_sur_zb(k), h_surf(i,j), qss_slope(k),ss_slope(k), Qg(k), water_v_cell(k), width(i,j), slo_grad(k), B_chan(k),l_chan(k),D_chan(k), h_chan(k),q_chan(k), I_chan(k),slo_ang_sin(k),slo_ang_cos(k),ns_MS(k), 1./ns_MS(k)*I_chan(k)**0.5*h_chan(k)**(2./3.)/(grav*h_chan(k))**0.5, (slo_Esi(k,m), m = 1, Np), (slo_Dsi(k,m), m = 1, Np)
            end do
        end if 
        close(1151)                  
      end if
	  close(115)
    endif
!-----RSR until here

            ! output (binary)
            if (outswitch_hs .eq. 2) write (100, rec=1) ((hs(i, j), j=1, nx), i=ny, 1, -1)
            if (outswitch_hr .eq. 2) write (101, rec=1) ((hr(i, j), j=1, nx), i=ny, 1, -1)
            if (outswitch_hg .eq. 2) write (102, rec=1) ((hg(i, j), j=1, nx), i=ny, 1, -1)
            if (outswitch_qr .eq. 2) write (103, rec=1) ((qr_ave(i, j), j=1, nx), i=ny, 1, -1) ! [m3/s]
            if (outswitch_qu .eq. 2) write (104, rec=1) (((qs_ave(1, i, j) + (qs_ave(3, i, j) - qs_ave(4, i, j))/2.d0)*area), &
                                                         i=ny, 1, -1)
            if (outswitch_qv .eq. 2) write (105, rec=1) (((qs_ave(2, i, j) + (qs_ave(3, i, j) + qs_ave(4, i, j))/2.d0)*area), &
                                                         i=ny, 1, -1)
            if (outswitch_gu .eq. 2) write (106, rec=1) (((qg_ave(1, i, j) + (qg_ave(3, i, j) - qg_ave(4, i, j))/2.d0)*area), &
                                                         i=ny, 1, -1)
            if (outswitch_gv .eq. 2) write (107, rec=1) (((qg_ave(2, i, j) + (qg_ave(3, i, j) + qg_ave(4, i, j))/2.d0)*area), &
                                                         i=ny, 1, -1)
            if (outswitch_gampt_ff .eq. 2) write (108, rec=1) ((gampt_ff(i, j), j=1, nx), i=ny, 1, -1)

            if (outswitch_hs .ne. 0) close (100)
            if (outswitch_hr .ne. 0) close (101)
            if (outswitch_hg .ne. 0) close (102)
            if (outswitch_qr .ne. 0) close (103)
            if (outswitch_qu .ne. 0) close (104)
            if (outswitch_qv .ne. 0) close (105)
            if (outswitch_gu .ne. 0) close (106)
            if (outswitch_gv .ne. 0) close (107)
            if (outswitch_gampt_ff .ne. 0) close (108)

            !iRIC Output  !modified for RSR 20240724
            call iric_cgns_output_result(sum_qp_t, qp_t, hs, hr, hg, qr_ave, qs_ave, qg_ave, &
                                       qsb, qss, sumdzb, sumqsb, sumqss)
            if (tec_switch .eq. 1) then
                if (tt .eq. 1) then
                    call Tecout_alloc(nx, ny, 4)
                    call Tecout_mkGrid(dx, dy, zs)
                    call Tecout_write_initialize(tt, width, depth, height, area_ratio)
                end if
                call Tecout_write(tt, qp_t, hr, qr_ave, hs, area)
            end if

            ! For dt_check
            !call dt_check_riv(hr_idx, tt, ddt_chk_riv)
            !call dt_check_slo(hs_idx, tt, ddt_chk_slo)

        end if

        ! check water balance
        if (mod(t, 1) .eq. 0) then
            call storage_calc(hs, hr, hg, ss, sr, si, sg)
            write (*, '(6e15.3)') rain_sum, pevp_sum, aevp_sum, sout, ss + sr + si + sg, &
                (rain_sum - aevp_sum - sout - (ss + sr + si + sg) + sinit)
            if (outswitch_storage == 1) write (1000, '(1000e15.7)') rain_sum, pevp_sum, aevp_sum, sout, ss + sr + si + sg, &
                (rain_sum - aevp_sum - sout - (ss + sr + si + sg) + sinit), ss, sr, si, sg
        end if

        !iRIC Cancel check and Flush
        call iric_check_cancel(ierr)
        if (ierr == 1) then
            write (*, *) "Solver is stopped because the STOP button was clicked."
            call iric_cgns_close()
            stop
        end if
!$omp end single
    end do

    if (sed_switch .ne. 0) then
    end if

    call iric_cgns_close

!pause

end program RRI

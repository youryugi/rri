subroutine RRI_Read

    use globals
    use dam_mod
    use tecout_mod
    use sediment_mod    !added for RSR model 20240724
    use runge_mod       !added for RSR model 20240724
    use RRI_iric
    use iric

    implicit none

    integer i, j,num_of_bound_point
    character*256 format_version

    integer :: ier, tmp,mm,k,n
!--------------------------------------------------
!CGNSファイルを開く
!--------------------------------------------------
    icount = iargc()
    if (icount == 1) then
        call getarg(1, cgns_name)
    else
        write (*, "(a)") "You should specify an argument."
        stop
    end if
    call cg_iric_open(cgns_name, IRIC_MODE_MODIFY, cgns_f, ier)
    if (ier /= 0) stop "cg_iric_open failed"

!--------------------------------------------------
!RRI バージョン情報
!--------------------------------------------------
!read(1,'(a)') format_version
!call cg_iric_read_string(cgns_f, "rri_ver", format_version, ier)
!
!write(*,'("format_version : ", a)') trim(adjustl(format_version))
!if( format_version .ne. "Ver1_4_2 for iRIC" ) stop "This RRI model requires RRI_Input_Format_Ver1_4_2"
!write(*,*)

!Run Type
    call cg_iric_read_integer(cgns_f, "run_type", run_type, ier)

!--------------------------------------------------
!外部ファイル　→　格子属性として設定
!--------------------------------------------------
!read(1,'(a)') rainfile
    call cg_iric_read_string(cgns_f, "rainfile", rainfile, ier)

!read(1,'(a)') demfile
    call cg_iric_read_string(cgns_f, "demfile", demfile, ier)

!read(1,'(a)') accfile
    call cg_iric_read_string(cgns_f, "accfile", accfile, ier)

!read(1,'(a)') dirfile
    call cg_iric_read_string(cgns_f, "dirfile", dirfile, ier)

    write (*, '("rainfile : ", a)') trim(adjustl(rainfile))
    write (*, '("demfile : ", a)') trim(adjustl(demfile))
    write (*, '("accfile : ", a)') trim(adjustl(accfile))
    write (*, '("dirfile : ", a)') trim(adjustl(dirfile))

!
!read(1,*)
    write (*, *)

!--------------------------------------------------
!基本オプション
!--------------------------------------------------
!read(1,*) utm
call cg_iric_read_integer(cgns_f, "utm", utm, ier)
!    utm = 0

!read(1,*) eight_dir
call cg_iric_read_integer(cgns_f, "eight_dir", eight_dir, ier)
   ! eight_dir = 1
!    eight_dir = 0 !using 4 direction 20240808
    omp_num_threads = 0
    call cg_iric_read_integer(cgns_f, "omp_num_threads", omp_num_threads, ier)
    if (ier /= 0) omp_num_threads = 0
    write (*, '("utm : ", i5)') utm
    write (*, '("eight_dir : ", i5)') eight_dir
    write (*, '("omp_num_threads : ", i5)') omp_num_threads
    write (*, *)

!--------------------------------------------------
!時間条件
!--------------------------------------------------
!read(1,*) lasth
    call cg_iric_read_integer(cgns_f, "lasth", lasth, ier)

!read(1,*) dt
    call cg_iric_read_integer(cgns_f, "dt", dt, ier)

!read(1,*) dt_riv
    call cg_iric_read_integer(cgns_f, "dt_riv", dt_riv, ier)

!read(1,*) outnum
    call cg_iric_read_integer(cgns_f, "outnum", outnum, ier)

    write (*, '("lasth : ", i8)') lasth
    write (*, '("dt : ", i12)') dt
    write (*, '("dt_riv : ", i8)') dt_riv
    write (*, '("outnum : ", i8)') outnum
    write (*, *)

!--------------------------------------------------
!降雨条件　格子属性として設定
!--------------------------------------------------
!read(1,*) xllcorner_rain
!read(1,*) yllcorner_rain
!read(1,*) cellsize_rain_x, cellsize_rain_y

    call cg_iric_read_real(cgns_f, "xllcorner_rain", xllcorner_rain, ier)
    call cg_iric_read_real(cgns_f, "yllcorner_rain", yllcorner_rain, ier)
    call cg_iric_read_real(cgns_f, "cellsize_rain_x", cellsize_rain_x, ier)
    call cg_iric_read_real(cgns_f, "cellsize_rain_y", cellsize_rain_y, ier)

    write (*, '("xllcorner_rain : ", f15.5)') xllcorner_rain
    write (*, '("yllcorner_rain : ", f15.5)') yllcorner_rain
    write (*, '("cellsize_rain_x : ", f15.5, "  cellsize_rain_y : ", f15.5)') cellsize_rain_x, cellsize_rain_y

    write (*, *)

!--------------------------------------------------
!slopeパラメータ　→　格子セル属性　複合条件として設定
!--------------------------------------------------
!read(1,*) num_of_landuse
!call cg_iric_read_complex_count(cgns_f, "landuse_c", num_of_landuse, ier)
!RRI for iRIC version
!Number of Land use type is 3 as fixed.
    num_of_landuse = 5
!
    allocate (dif(num_of_landuse))
    allocate (ns_slope(num_of_landuse), soildepth(num_of_landuse))
    allocate (gammaa(num_of_landuse))
!
!dif
    call cg_iric_read_integer(cgns_f, "dif_1", dif(1), ier)
    call cg_iric_read_integer(cgns_f, "dif_2", dif(2), ier)
    call cg_iric_read_integer(cgns_f, "dif_3", dif(3), ier)
    call cg_iric_read_integer(cgns_f, "dif_4", dif(4), ier)
    call cg_iric_read_integer(cgns_f, "dif_5", dif(5), ier)

!ns_slope
    call cg_iric_read_real(cgns_f, "ns_slope_1", ns_slope(1), ier)
    call cg_iric_read_real(cgns_f, "ns_slope_2", ns_slope(2), ier)
    call cg_iric_read_real(cgns_f, "ns_slope_3", ns_slope(3), ier)
    call cg_iric_read_real(cgns_f, "ns_slope_4", ns_slope(4), ier)
    call cg_iric_read_real(cgns_f, "ns_slope_5", ns_slope(5), ier)

!solid depth
    call cg_iric_read_real(cgns_f, "soildepth_1", soildepth(1), ier)
    call cg_iric_read_real(cgns_f, "soildepth_2", soildepth(2), ier)
    call cg_iric_read_real(cgns_f, "soildepth_3", soildepth(3), ier)
    call cg_iric_read_real(cgns_f, "soildepth_4", soildepth(4), ier)
    call cg_iric_read_real(cgns_f, "soildepth_5", soildepth(5), ier)

!Porosity
    call cg_iric_read_real(cgns_f, "gammaa_1", gammaa(1), ier)
    call cg_iric_read_real(cgns_f, "gammaa_2", gammaa(2), ier)
    call cg_iric_read_real(cgns_f, "gammaa_3", gammaa(3), ier)
    call cg_iric_read_real(cgns_f, "gammaa_4", gammaa(4), ier)
    call cg_iric_read_real(cgns_f, "gammaa_5", gammaa(5), ier)

!

    write (*, '("num_of_landuse : ", i5)') num_of_landuse
    write (*, '("dif : ", 100i5)') (dif(i), i=1, num_of_landuse)
    write (*, '("ns_slope : ", 100f12.3)') (ns_slope(i), i=1, num_of_landuse)
    write (*, '("soildepth : ", 100f12.3)') (soildepth(i), i=1, num_of_landuse)
    write (*, '("gammaa : ", 100f12.3)') (gammaa(i), i=1, num_of_landuse)
!
!
!read(1,*)
    write (*, *)
!
    allocate (ksv(num_of_landuse), faif(num_of_landuse))
!
!ksv
    call cg_iric_read_real(cgns_f, "ksv_1", ksv(1), ier)
    call cg_iric_read_real(cgns_f, "ksv_2", ksv(2), ier)
    call cg_iric_read_real(cgns_f, "ksv_3", ksv(3), ier)
    call cg_iric_read_real(cgns_f, "ksv_4", ksv(4), ier)
    call cg_iric_read_real(cgns_f, "ksv_5", ksv(5), ier)

!Sf
    call cg_iric_read_real(cgns_f, "faif_1", faif(1), ier)
    call cg_iric_read_real(cgns_f, "faif_2", faif(2), ier)
    call cg_iric_read_real(cgns_f, "faif_3", faif(3), ier)
    call cg_iric_read_real(cgns_f, "faif_4", faif(4), ier)
    call cg_iric_read_real(cgns_f, "faif_5", faif(5), ier)

!
    write (*, '("ksv : ", 100e12.3)') (ksv(i), i=1, num_of_landuse)
    write (*, '("faif : ", 100f12.3)') (faif(i), i=1, num_of_landuse)
!
!read(1,*)
    write (*, *)
!
    allocate (ka(num_of_landuse), gammam(num_of_landuse), beta(num_of_landuse),min_wc4latflow(num_of_landuse))
    min_wc4latflow(:)=0.d0
!
!ka
    call cg_iric_read_real(cgns_f, "ka_1", ka(1), ier)
    call cg_iric_read_real(cgns_f, "ka_2", ka(2), ier)
    call cg_iric_read_real(cgns_f, "ka_3", ka(3), ier)
    call cg_iric_read_real(cgns_f, "ka_4", ka(4), ier)
    call cg_iric_read_real(cgns_f, "ka_5", ka(5), ier)

!Unsat. porosity
    call cg_iric_read_real(cgns_f, "gammam_1", gammam(1), ier)
    call cg_iric_read_real(cgns_f, "gammam_2", gammam(2), ier)
    call cg_iric_read_real(cgns_f, "gammam_3", gammam(3), ier)
    call cg_iric_read_real(cgns_f, "gammam_4", gammam(4), ier)
    call cg_iric_read_real(cgns_f, "gammam_5", gammam(5), ier)

!beta
    call cg_iric_read_real(cgns_f, "beta_1", beta(1), ier)
    call cg_iric_read_real(cgns_f, "beta_2", beta(2), ier)
    call cg_iric_read_real(cgns_f, "beta_3", beta(3), ier)
    call cg_iric_read_real(cgns_f, "beta_4", beta(4), ier)
    call cg_iric_read_real(cgns_f, "beta_5", beta(5), ier)
!min_wc
    !call cg_iric_read_real(cgns_f, "min_wc_1", min_wc4latflow(1), ier)
    !call cg_iric_read_real(cgns_f, "min_wc_2", min_wc4latflow(2), ier)
    !call cg_iric_read_real(cgns_f, "min_wc_3", min_wc4latflow(3), ier)
    !call cg_iric_read_real(cgns_f, "min_wc_4", min_wc4latflow(4), ier)
    !call cg_iric_read_real(cgns_f, "min_wc_5", min_wc4latflow(5), ier)    

    write (*, '("ka : ", 100e12.3)') (ka(i), i=1, num_of_landuse)
    write (*, '("gammam : ", 100f12.3)') (gammam(i), i=1, num_of_landuse)
    write (*, '("beta : ", 100f12.3)') (beta(i), i=1, num_of_landuse)
   ! write (*, '("Minimum water content for forming lateral flow : ", 100f12.3)') (min_wc4latflow(i), i=1, num_of_landuse)

    do i = 1, num_of_landuse
        if (gammam(i) .gt. gammaa(i)) stop "gammag must be smaller than gammaa"
    end do
!
!read(1,*)
    write (*, *)
!
    allocate (ksg(num_of_landuse), gammag(num_of_landuse), kg0(num_of_landuse), fpg(num_of_landuse), rgl(num_of_landuse))
!
! ksg = 0.0 means that this model has not been implemented into iRIC version. ****
!
    ksg = 0.0; gammag = 0.0; kg0 = 0.0; fpg = 0.0; rgl = 0.0

!
    write (*, '("ksg : ", 100e12.3)') (ksg(i), i=1, num_of_landuse)
    write (*, '("gammag : ", 100f12.3)') (gammag(i), i=1, num_of_landuse)
    write (*, '("kg0 : ", 100e12.3)') (kg0(i), i=1, num_of_landuse)
    write (*, '("fpg : ", 100f12.3)') (fpg(i), i=1, num_of_landuse)
    write (*, '("rgl : ", 100e12.3)') (rgl(i), i=1, num_of_landuse)
!read(1,*)
    write (*, *)

!--------------------------------------------------
!河道パラメータ　→　wc,ws,dc,dsで指定する
!--------------------------------------------------
!read(1,*) ns_river
    call cg_iric_read_real(cgns_f, "ns_river", ns_river, ier)
    write (*, '("ns_river : ", f12.3)') ns_river
    write (*, *)

!read(1,*) riv_thresh
    call cg_iric_read_integer(cgns_f, "riv_thresh", riv_thresh, ier)
    write (*, '("riv_thresh : ", i7)') riv_thresh
    write (*, *)

!read(1,*) width_param_c
!read(1,*) width_param_s
    call cg_iric_read_real(cgns_f, "width_param_c", width_param_c, ier)
    call cg_iric_read_real(cgns_f, "width_param_s", width_param_s, ier)

!read(1,*) depth_param_c
!read(1,*) depth_param_s
    call cg_iric_read_real(cgns_f, "depth_param_c", depth_param_c, ier)
    call cg_iric_read_real(cgns_f, "depth_param_s", depth_param_s, ier)

!read(1,*) height_param
!read(1,*) height_limit_param
    call cg_iric_read_real(cgns_f, "height_param", height_param, ier)
    call cg_iric_read_integer(cgns_f, "height_limit_param", height_limit_param, ier)

!read(1,*) rivfile_switch
!call cg_iric_read_integer(cgns_f, "rivfile_switch", rivfile_switch, ier)
!rivfile_switchは利用しない
    rivfile_switch = 0

!read(1,'(a)') widthfile
!read(1,'(a)') depthfile
!read(1,'(a)') heightfile

    if (rivfile_switch .eq. 0) then
        write (*, '("width_param_c : ", f12.2)') width_param_c
        write (*, '("width_param_s : ", f12.2)') width_param_s
        write (*, '("depth_param_c : ", f12.2)') depth_param_c
        write (*, '("depth_param_s : ", f12.2)') depth_param_s
        write (*, '("height_param : ", f12.2)') height_param
        write (*, '("height_limit_param : ", i10)') height_limit_param
    else
        write (*, '("widthfile : ", a)') trim(adjustl(widthfile))
        write (*, '("depthfile : ", a)') trim(adjustl(depthfile))
        write (*, '("heightfile : ", a)') trim(adjustl(heightfile))
    end if

!read(1,*)
    write (*, *)

!--------------------------------------------------
!-------------------
!hotstart用　初期条件
!--------------------------------------------------
!read(1,*) init_slo_switch, init_riv_switch, init_gw_switch, init_gampt_ff_switch
    call cg_iric_read_integer(cgns_f, "init_slo_switch", init_slo_switch, ier)
    call cg_iric_read_integer(cgns_f, "init_riv_switch", init_riv_switch, ier)
    call cg_iric_read_integer(cgns_f, "init_gw_switch", init_gw_switch, ier)
    call cg_iric_read_integer(cgns_f, "init_gampt_ff_switch", init_gampt_ff_switch, ier)

!read(1,"(a)") initfile_slo
!read(1,'(a)') initfile_riv
!read(1,'(a)') initfile_gw
!read(1,'(a)') initfile_gampt_ff

    initfile_slo = ''; initfile_riv = ''; initfile_gw = ''; initfile_gampt_ff = ''
    if (init_slo_switch == 1) call cg_iric_read_string(cgns_f, "initfile_slo", initfile_slo, ier)
    if (init_riv_switch == 1) call cg_iric_read_string(cgns_f, "initfile_riv", initfile_riv, ier)
    if (init_gw_switch == 1) call cg_iric_read_string(cgns_f, "initfile_gw", initfile_gw, ier)
    if (init_gampt_ff_switch == 1) call cg_iric_read_string(cgns_f, "initfile_gampt_ff", initfile_gampt_ff, ier)

    if (init_slo_switch .ne. 0) write (*, '("initfile_slo : ", a)') trim(adjustl(initfile_slo))
    if (init_riv_switch .ne. 0) write (*, '("initfile_riv : ", a)') trim(adjustl(initfile_riv))
    if (init_gw_switch .ne. 0) write (*, '("initfile_gw : ", a)') trim(adjustl(initfile_gw))
    if (init_gampt_ff_switch .ne. 0) write (*, '("initfile_gampt_ff : ", a)') trim(adjustl(initfile_gampt_ff))

!read(1,*)
    write (*, *)

!--------------------------------------------------
!hs, hr境界条件　→　境界条件設定に実装
!--------------------------------------------------
    bound_slo_wlev_switch = 0; bound_riv_wlev_switch = 0
    call cg_iric_read_bc_count(cgns_f, "bound_hs", num_of_bound_point)
    if (num_of_bound_point > 0) bound_slo_wlev_switch = 1

    call cg_iric_read_bc_count(cgns_f, "bound_hr", num_of_bound_point)
    if (num_of_bound_point > 0) bound_riv_wlev_switch = 1

!read(1,*) bound_slo_wlev_switch, bound_riv_wlev_switch
!read(1,'(a)') boundfile_slo_wlev
!read(1,'(a)') boundfile_riv_wlev
!if(bound_slo_wlev_switch.ne.0) write(*,'("boundfile_slo_wlev : ", a)') trim(adjustl(boundfile_slo_wlev))
!if(bound_riv_wlev_switch.ne.0) write(*,'("boundfile_riv_wlev : ", a)') trim(adjustl(boundfile_riv_wlev))
!read(1,*)
!write(*,*)

!--------------------------------------------------
!qs, qr境界条件　→　境界条件設定に実装
!--------------------------------------------------
    bound_slo_disc_switch = 0; bound_riv_disc_switch = 0
    call cg_iric_read_bc_count(cgns_f, "bound_qs", num_of_bound_point)
    if (num_of_bound_point > 0) bound_slo_disc_switch = 1

    call cg_iric_read_bc_count(cgns_f, "bound_qr", num_of_bound_point)
    if (num_of_bound_point > 0) bound_riv_disc_switch = 1

!read(1,*) bound_slo_disc_switch, bound_riv_disc_switch
!read(1,'(a)') boundfile_slo_disc
!read(1,'(a)') boundfile_riv_disc
!if(bound_slo_disc_switch.ne.0) write(*,'("boundfile_slo_disc : ", a)') trim(adjustl(boundfile_slo_disc))
!if(bound_riv_disc_switch.ne.0) write(*,'("boundfile_riv_disc : ", a)') trim(adjustl(boundfile_riv_disc))
!read(1,*)
!write(*,*)

!--------------------------------------------------
!土地利用条件　→　格子属性として与える
!--------------------------------------------------
!read(1,*) land_switch
!read(1,'(a)') landfile
!if(land_switch.eq.1) write(*,'("landfile : ", a)') trim(adjustl(landfile))
!
!read(1,*)
!write(*,*)

!--------------------------------------------------
!Dam条件　→　境界条件設定に実装
!--------------------------------------------------
    dam_switch = 0
    call cg_iric_read_bc_count(cgns_f, "dam", dam_num)
    if (dam_num > 0) dam_switch = 1
    !for dam
    if(dam_switch>0)then
    write(*,'("Number of dam : ", i7)') dam_num
    endif
!read(1,*) dam_switch
!read(1,'(a)') damfile
!if(dam_switch.eq.1) write(*,'("damfile : ", a)') trim(adjustl(damfile))
!read(1,*)
!write(*,*)

!--------------------------------------------------
!div条件　→　境界条件設定に実装
!--------------------------------------------------
    div_id_max = 0
    call cg_iric_read_bc_count(cgns_f, "div", div_id_max)
    if (div_id_max > 0) div_switch = 1
!read(1,*) div_switch
!call cg_iric_read_integer(cgns_f, "div_switch", div_switch, ier)
!read(1,'(a)') divfile
!if(div_switch.eq.1) write(*,'("divfile : ", a)') trim(adjustl(divfile))
!read(1,*)
!write(*,"(a)") "div condition has not implemented yet."
!write(*,*)

!--------------------------------------------------
!蒸発条件　→
!--------------------------------------------------
!read(1,*) evp_switch
!read(1,'(a)') evpfile
!read(1,*) xllcorner_evp
!read(1,*) yllcorner_evp
!read(1,*) cellsize_evp_x, cellsize_evp_y
!if( evp_switch .ne. 0 ) then
! write(*,'("evpfile : ", a)') trim(adjustl(evpfile))
! write(*,'("xllcorner_evp : ", f15.5)') xllcorner_evp
! write(*,'("yllcorner_evp : ", f15.5)') yllcorner_evp
! write(*,'("cellsize_evp_x : ", f15.5, " cellsize_evp_y : ", f15.5)') cellsize_evp_x, cellsize_evp_y
!endif
!read(1,*)
!write(*,"(a)") "evp condition can be set as the grid attribute."
!write(*,*)

!--------------------------------------------------
!河道断面データ　→　実装しない
!--------------------------------------------------
!read(1,*) sec_length_switch
!read(1,'(a)') sec_length_file        !これはセルごとに河道長を設定するファイル

!if(sec_length_switch.eq.1) write(*,'("sec_length : ", a)') trim(adjustl(sec_length_file))
!read(1,*)
!write(*,*)
!read(1,*) sec_switch
!read(1,'(a)') sec_map_file        !これはセルごとにsecfileのidを指定するファイル
!read(1,'(a)') sec_file                !これはsection形状のファイル　なければwc,ws,dc,dsで指定された値となる
!if(sec_switch.eq.1) write(*,'("sec_map_file : ", a)') trim(adjustl(sec_map_file))
!if(sec_switch.eq.1) write(*,'("sec_file : ", a)') trim(adjustl(sec_file))
!read(1,*)
!write(*,"(a)") "Cross section condition has not implemented yet."
!write(*,*)

!--------------------------------------------------
!??? 　→　要確認
!--------------------------------------------------
!read(1,*) emb_switch
!read(1,'(a)') embrfile
!read(1,'(a)') embbfile
!if(emb_switch.eq.1) write(*,'("embrfile : ", a)') trim(adjustl(embrfile))
!if(emb_switch.eq.1) write(*,'("embbfile : ", a)') trim(adjustl(embbfile))
!write(*,*)

!--------------------------------------------------
!計算結果出力・ファイル
!--------------------------------------------------
!read(1,*) outswitch_hs, outswitch_hr, outswitch_hg, outswitch_qr, outswitch_qu, outswitch_qv, &
!          outswitch_gu, outswitch_gv, outswitch_gampt_ff, outswitch_storage

!read(1,'(a)') outfile_hs
    outfile_hs = ''
    call cg_iric_read_integer(cgns_f, "outswitch_hs", outswitch_hs, ier)
    if (outswitch_hs == 1) call cg_iric_read_string(cgns_f, "outfile_hs", outfile_hs, ier)

!read(1,'(a)') outfile_hr
    outfile_hr = ''
    call cg_iric_read_integer(cgns_f, "outswitch_hr", outswitch_hr, ier)
    if (outswitch_hr == 1) call cg_iric_read_string(cgns_f, "outfile_hr", outfile_hr, ier)

!read(1,'(a)') outfile_hg
    outfile_hg = ''
    call cg_iric_read_integer(cgns_f, "outswitch_hg", outswitch_hg, ier)
    if (outswitch_hg == 1) call cg_iric_read_string(cgns_f, "outfile_hg", outfile_hg, ier)

!read(1,'(a)') outfile_qr
    outfile_qr = ''
    call cg_iric_read_integer(cgns_f, "outswitch_qr", outswitch_qr, ier)
    if (outswitch_qr == 1) call cg_iric_read_string(cgns_f, "outfile_qr", outfile_qr, ier)

!read(1,'(a)') outfile_qu
    outfile_qu = ''
    call cg_iric_read_integer(cgns_f, "outswitch_qu", outswitch_qu, ier)
    if (outswitch_qu == 1) call cg_iric_read_string(cgns_f, "outfile_qu", outfile_qu, ier)

!read(1,'(a)') outfile_qv
    outfile_qv = ''
    call cg_iric_read_integer(cgns_f, "outswitch_qv", outswitch_qv, ier)
    if (outswitch_qv == 1) call cg_iric_read_string(cgns_f, "outfile_qv", outfile_qv, ier)

!read(1,'(a)') outfile_gu
    outfile_gu = ''
    call cg_iric_read_integer(cgns_f, "outswitch_gu", outswitch_gu, ier)
    if (outswitch_gu == 1) call cg_iric_read_string(cgns_f, "outfile_gu", outfile_gu, ier)

!read(1,'(a)') outfile_gv
    outfile_gv = ''
    call cg_iric_read_integer(cgns_f, "outswitch_gv", outswitch_gv, ier)
    if (outswitch_gv == 1) call cg_iric_read_string(cgns_f, "outfile_gv", outfile_gv, ier)

!read(1,'(a)') outfile_gampt_ff
    outfile_gampt_ff = ''
    call cg_iric_read_integer(cgns_f, "outswitch_gampt_ff", outswitch_gampt_ff, ier)
    if (outswitch_gampt_ff == 1) call cg_iric_read_string(cgns_f, "outfile_gampt_ff", outfile_gampt_ff, ier)

!read(1,'(a)') outfile_storage
    outfile_storage = ''
    call cg_iric_read_integer(cgns_f, "outswitch_storage", outswitch_storage, ier)
    if (outswitch_storage == 1) call cg_iric_read_string(cgns_f, "outfile_storage", outfile_storage, ier)

    if (outswitch_hs .ne. 0) write (*, '("outfile_hs : ", a)') trim(adjustl(outfile_hs))
    if (outswitch_hr .ne. 0) write (*, '("outfile_hr : ", a)') trim(adjustl(outfile_hr))
    if (outswitch_hg .ne. 0) write (*, '("outfile_hg : ", a)') trim(adjustl(outfile_hg))
    if (outswitch_qr .ne. 0) write (*, '("outfile_qr : ", a)') trim(adjustl(outfile_qr))
    if (outswitch_qu .ne. 0) write (*, '("outfile_qu : ", a)') trim(adjustl(outfile_qu))
    if (outswitch_qv .ne. 0) write (*, '("outfile_qv : ", a)') trim(adjustl(outfile_qv))
    if (outswitch_gu .ne. 0) write (*, '("outfile_gu : ", a)') trim(adjustl(outfile_gu))
    if (outswitch_gv .ne. 0) write (*, '("outfile_gv : ", a)') trim(adjustl(outfile_gv))
    if (outswitch_gampt_ff .ne. 0) write (*, '("outfile_gampt_ff : ", a)') trim(adjustl(outfile_gampt_ff))
    if (outswitch_storage .ne. 0) write (*, '("outfile_storage : ", a)') trim(adjustl(outfile_storage))

!read(1,*)
    write (*, *)

!----read advanced settings--modified 20250317
 !------Advanced setting eps 20240724
        call cg_iric_read_real(cgns_f, "eps", eps, ier)
        write (*, '("eps: ", f12.5)') eps
        call cg_iric_read_real(cgns_f, "ddt_min_riv", ddt_min_riv, ier)
        write (*, '("ddt_min_riv: ", f12.5)') ddt_min_riv
        call cg_iric_read_real(cgns_f, "ddt_min_slo", ddt_min_slo, ier)
        write (*, '("ddt_min_slo: ", f12.5)') ddt_min_slo     

!------------Initial conditions of water contents of surface soil and river flow depth
!moved to RRI setting 20250312
    call cg_iric_read_real(cgns_f, "hr0", hr0, ier)
    call cg_iric_read_real(cgns_f, "wc0", wc0, ier)        
        write (*, '("Initial water content of surface soil: ", f12.5)') wc0
        write (*, '("Initial flow depth : ", f12.5)') hr0    

!----setting of downstream boundary condition    
    call cg_iric_read_integer(cgns_f, "DBC_switch", DBC_switch, ier)
    if (DBC_switch==0)then 
    write(*,*) "Downstream boundary condition is same with the original RRI"
    elseif(DBC_switch==1)then
        write(*,*) "Downstream boundary condition is setted as free flow"
    endif    
!------------------------------------

 !------added for RSR model 20240724
    allocate( kgv(num_of_landuse), tg(num_of_landuse))

    call cg_iric_read_integer(cgns_f, "sed_switch", sed_switch, ier)
    write (*, '("sed_switch : ", i7)') sed_switch
    write (*, *)
    channel_blockage_switch = 0
    channel_blockage_leakage = 0.d0
    channel_blockage_exponent = 5.d0/3.d0
    channel_capacity_qr_switch = 0
    river_overtop_neighbor_switch = 0
    river_slope_preexchange_switch = 0
    !In case sediment computation, read the following
    if (sed_switch >= 1)then
        call cg_iric_read_real(cgns_f, "t_beddeform_start", t_beddeform_start, ier)
        write (*, '("t_beddeform_start(hour): ", f12.1)') t_beddeform_start
        t_beddeform_start = t_beddeform_start*3600.
        t_bed_freeze = 0.d0
        call cg_iric_read_real(cgns_f, "t_bed_freeze", t_bed_freeze, ier)
        if(ier/=0) t_bed_freeze = 0.d0
        write (*, '("t_bed_freeze(hour): ", f12.1)') t_bed_freeze
        t_bed_freeze = t_bed_freeze*3600.   ! hours -> seconds; during time<t_bed_freeze: transport & fm update, but bed elevation frozen
        call cg_iric_read_real(cgns_f, "s", s, ier)
        write (*, '("s : ", f12.2)') s
        grav = 9.81d0
        call cg_iric_read_real(cgns_f, "grav", grav, ier)
        if(ier/=0) grav = 9.81d0
        write (*, '("grav : ", f12.2)') grav
        call cg_iric_read_real(cgns_f, "t_Crit", t_Crit, ier)
        write (*, '("t_Crit : ", f12.3)') t_Crit
        kin_visc = 1.004d-6
        call cg_iric_read_real(cgns_f, "kin_visc", kin_visc, ier)
        if(ier/=0) kin_visc = 1.004d-6
        write (*, '("kin_visc : ", f12.3)') kin_visc
        karmans = 0.4d0
        call cg_iric_read_real(cgns_f, "karmans", karmans, ier)
        if(ier/=0) karmans = 0.4d0
        write (*, '("karmans : ", f12.1)') karmans
        call cg_iric_read_real(cgns_f, "lambda", lambda, ier)
        write (*, '("lambda : ", f12.1)') lambda
        write (*, *)
        call cg_iric_read_integer(cgns_f, "sed_type_switch", sed_type_switch, ier)
        write (*, '("sed_type_switch : ", i7)') sed_type_switch
        call cg_iric_read_real(cgns_f, "ds_river", ds_river, ier)
        write (*, '("ds_river : ", f12.5)') ds_river
        call cg_iric_read_integer(cgns_f, "iidt", iidt, ier)
        write (*, '("iidt : ", i7)') iidt
        call cg_iric_read_integer(cgns_f, "isedeq", isedeq, ier)
        write (*, '("isedeq : ", i7)') isedeq
        call cg_iric_read_integer(cgns_f, "isuseq", isuseq, ier)
        write (*, '("isuseq : ", i7)') isuseq
        write (*, *)
        !Regulations
!        call cg_iric_read_integer(cgns_f, "max_acc_0th_riv", max_acc_0th_riv, ier)
!        write (*, '("max_acc_0th_riv : ", i7)') max_acc_0th_riv
        call cg_iric_read_integer(cgns_f, "min_num_cell_link", min_num_cell_link, ier)
        write (*, '("min_num_cell_link : ", i7)') min_num_cell_link
        call cg_iric_read_real(cgns_f, "perosion", perosion, ier)
        write (*, '("perosion : ", f12.5)') perosion
        put_rate_fac = 1.d0
        call cg_iric_read_real(cgns_f, "put_rate_fac", put_rate_fac, ier)
        if(ier/=0 .or. put_rate_fac<=0.d0) put_rate_fac = 1.d0
        write (*, '("put_rate_fac : ", f12.5)') put_rate_fac
        d_wash = 1.d-4
        call cg_iric_read_real(cgns_f, "d_wash", d_wash, ier)
        if(ier/=0 .or. d_wash<=0.d0) d_wash = 1.d-4
        write (*, '("d_wash : ", es12.4)') d_wash
        call cg_iric_read_real(cgns_f, "min_slope", min_slope, ier)
        write (*, '("min_slope : ", f12.5)') min_slope
        call cg_iric_read_real(cgns_f, "max_slope", max_slope, ier)
        write (*, '("max_slope : ", f12.2)') max_slope
        call cg_iric_read_real(cgns_f, "min_hr", min_hr, ier)
        write (*, '("min_hr : ", f12.5)') min_hr
        alpha_ss1 = 0.5d0
        call cg_iric_read_real(cgns_f, "alpha_ss1", alpha_ss1, ier)
        if(ier/=0) alpha_ss1 = 0.5d0
        write (*, '("alpha_ss1 : ", f12.5)') alpha_ss1
        alpha_ss2 = 0.1d0
        call cg_iric_read_real(cgns_f, "alpha_ss2", alpha_ss2, ier)
        if(ier/=0) alpha_ss2 = 0.1d0
        write (*, '("alpha_ss2 : ", f12.5)') alpha_ss2
        thresh_ss = 0.3d0
        call cg_iric_read_real(cgns_f, "thresh_ss", thresh_ss, ier)
        if(ier/=0) thresh_ss = 0.3d0
        write (*, '("thresh_ss : ", f12.5)') thresh_ss
       !moved to rri setting 20250312 
      !  call cg_iric_read_real(cgns_f, "hr0", hr0, ier)
      !  write (*, '("Initial flow depth : ", f12.5)') hr0
       ! call cg_iric_read_real(cgns_f, "wc0", wc0, ier)
      !  write (*, '("Initial water depth of slope cells: ", f12.5)') wc0
        call cg_iric_read_integer(cgns_f, "cut_overdepo_switch", cut_overdepo_switch, ier) 
         write (*, '("Enforcing sediment overflow: ", i7)') cut_overdepo_switch
        pass_bedload_switch = 0
        call cg_iric_read_integer(cgns_f, "pass_bedload_switch", pass_bedload_switch, ier)
        if(ier/=0) pass_bedload_switch = 0
         write (*, '("Pass bedload downstream when channel capacity is reduced: ", i7)') pass_bedload_switch
        pass_bedload_depth_thresh = 0.05d0
        call cg_iric_read_real(cgns_f, "pass_bedload_depth_thresh", pass_bedload_depth_thresh, ier)
        if(ier/=0 .or. pass_bedload_depth_thresh<0.d0) pass_bedload_depth_thresh = 0.05d0
         write (*, '("Remaining channel depth threshold for passing bedload downstream(m): ", f12.5)') pass_bedload_depth_thresh
        call cg_iric_read_integer(cgns_f, "channel_blockage_switch", channel_blockage_switch, ier)
        if(ier/=0) channel_blockage_switch = 0
         write (*, '("Reduce river conveyance according to bed aggradation: ", i7)') channel_blockage_switch
        call cg_iric_read_real(cgns_f, "channel_blockage_leakage", channel_blockage_leakage, ier)
        if(ier/=0) channel_blockage_leakage = 0.d0
        if(channel_blockage_leakage<0.d0) channel_blockage_leakage = 0.d0
        if(channel_blockage_leakage>1.d0) channel_blockage_leakage = 1.d0
         write (*, '("Residual conveyance ratio at complete blockage: ", f12.5)') channel_blockage_leakage
        call cg_iric_read_real(cgns_f, "channel_blockage_exponent", channel_blockage_exponent, ier)
        if(ier/=0 .or. channel_blockage_exponent<=0.d0) channel_blockage_exponent = 5.d0/3.d0
         write (*, '("Exponent for conveyance reduction by remaining channel depth ratio: ", f12.5)') channel_blockage_exponent
        channel_capacity_qr_switch = 0
        call cg_iric_read_integer(cgns_f, "channel_capacity_qr_switch", channel_capacity_qr_switch, ier)
        if(ier/=0) channel_capacity_qr_switch = 0
         write (*, '("Limit river discharge by remaining channel capacity: ", i7)') channel_capacity_qr_switch
        river_overtop_neighbor_switch = 0
        call cg_iric_read_integer(cgns_f, "river_overtop_neighbor_switch", river_overtop_neighbor_switch, ier)
        if(ier/=0) river_overtop_neighbor_switch = 0
         write (*, '("Distribute river overtopping water to adjacent non-river cells: ", i7)') river_overtop_neighbor_switch
        river_slope_preexchange_switch = 0
        call cg_iric_read_integer(cgns_f, "river_slope_preexchange_switch", river_slope_preexchange_switch, ier)
        if(ier/=0) river_slope_preexchange_switch = 0
         write (*, '("Perform river-slope exchange before slope routing: ", i7)') river_slope_preexchange_switch
        !Non-uniform
        if(sed_type_switch ==2)then
            call cg_iric_read_real(cgns_f, "Em", Em, ier)
            write (*, '("Em : ", f12.5)') Em
            call cg_iric_read_integer(cgns_f, "no_of_layers", no_of_layers, ier)
            write (*, '("no_of_layers : ", i7)') no_of_layers
            call cg_iric_read_integer(cgns_f, "Nl", Nl, ier)
            write (*, '("Nl : ", i7)') Nl
            call cg_iric_read_real(cgns_f, "Emc", Emc, ier)
            write (*, '("Emc : ", f12.5)') Emc
        end if

        call cg_iric_read_integer(cgns_f, "detail_console", detail_console, ier)

    end if

 !------Landslide and debris flow 20240724
    call cg_iric_read_integer(cgns_f, "debris_switch", debris_switch, ier)
    write (*, '("debris_switch : ", i7)') debris_switch
    if (debris_switch >= 1)then
        debris_rate_fac = 1.d0
        call cg_iric_read_real(cgns_f, "debris_rate_fac", debris_rate_fac, ier)
        if(ier/=0 .or. debris_rate_fac<=0.d0) debris_rate_fac = 1.d0
        write (*, '("debris_rate_fac : ", f12.5)') debris_rate_fac
        call cg_iric_read_real(cgns_f, "cohe", cohe, ier)
        write (*, '("cohesion: ", f12.2)') cohe
        call cg_iric_read_real(cgns_f, "phi", phi, ier)
        write (*, '("internal friction angle: ", f12.2)') phi
        if (phi <= 0.d0 .or. phi >= 90.d0) then
            write (*, *) 'Invalid landslide parameters.'
            write (*, *) 'phi must satisfy 0 < phi < 90 degrees.'
            write (*, *) 'phi =', phi
            stop 'Invalid phi'
        endif
        call cg_iric_read_real(cgns_f, "pwc", pwc, ier)
        write (*, '("pwc: ", f12.2)') pwc
        if (pwc < 0.d0 .or. pwc >= lambda) then
            write (*, *) 'Invalid landslide parameters.'
            write (*, *) 'pwc must satisfy 0 <= pwc < lambda.'
            write (*, *) 'pwc =', pwc, 'lambda =', lambda
            stop 'Invalid pwc and lambda'
        end if
        call cg_iric_read_real(cgns_f, "pf", pf, ier)
        write (*, '("pf: ", f12.2)') pf
        call cg_iric_read_real(cgns_f, "b_mp", b_mp, ier)
        write (*, '("b_mp: ", f12.2)') b_mp
        call cg_iric_read_real(cgns_f, "d_mp_ini", d_mp_ini, ier)
        write (*, '("d_mp_ini: ", f12.2)') d_mp_ini
        call cg_iric_read_real(cgns_f, "L_rain_ini", L_rain_ini, ier)
        write (*, '("L_rain_ini: ", f12.2)') L_rain_ini
        call cg_iric_read_integer(cgns_f, "debris_end_switch", debris_end_switch, ier)
        write (*, '("debris_end_switch: ", i7)') debris_end_switch
        if(debris_end_switch==1)then
            call cg_iric_read_real(cgns_f, "T_bebris_off", T_bebris_off, ier)
            write (*, '("T_bebris_off(hour): ", f12.2)') T_bebris_off
            T_bebris_off = T_bebris_off*3600.
        else
            T_bebris_off = 1000000000.
        end if            
    end if

 !------Slope erosion 20240724
 !---modified for slope erosion
    call cg_iric_read_integer(cgns_f, "slo_sedi_cal_switch", slo_sedi_cal_switch, ier)
    write (*, '("slo_sedi_cal_switch : ", i7)') slo_sedi_cal_switch
    if(slo_sedi_cal_switch>0) then
        call cg_iric_read_integer(cgns_f, "slope_ero_switch", slope_ero_switch, ier)
        write (*, '("slope_ero_switch : ", i7)') slope_ero_switch
        if(slope_ero_switch>0)then
        nm_cell =9 !tentative
        allocate (B_gully_r(nm_cell+1), D_gully(nm_cell+1),infil_s_depth(nm_cell+1))
        do i = 0, nm_cell
        write(cm,'(i1)') i
        gullyB_label = 'B_gully_r_'//trim(cm)  
        gullyD_label = 'D_gully_'//trim(cm)  
        infil_s_d_label='infil_s_d_'//trim(cm)
        call cg_iric_read_real(cgns_f, gullyB_label, B_gully_r(i+1), ier)      
        call cg_iric_read_real(cgns_f, gullyD_label, D_gully(i+1), ier)     
        call cg_iric_read_real(cgns_f, infil_s_d_label, infil_s_depth(i+1), ier) 
        enddo
    endif

        call cg_iric_read_real(cgns_f, "modirate_to_slo_dt", modirate_to_slo_dt, ier)
        write (*, '("modirate_to_slo_dt: ", f12.2)') modirate_to_slo_dt
        dt_slo_sed = dt/modirate_to_slo_dt 
        call cg_iric_read_real(cgns_f, "surflowdepth", surflowdepth, ier)
        write (*, '("surflowdepth: ", f12.2)') surflowdepth
    else
    slope_ero_switch= 0    
    end if

 !------Driftwood 20240724
    call cg_iric_read_integer(cgns_f, "j_drf", j_drf, ier)
    write (*, '("j_drf : ", i7)') j_drf
    call cg_iric_read_real(cgns_f, "Wood_density", Wood_density, ier)
    write (*, '("Wood_density : ", f12.5)') Wood_density
    call cg_iric_read_real(cgns_f, "C_wood_kd", C_wood_kd, ier)
    write (*, '("Wood_deposition_rate : ", f12.5)') C_wood_kd

 !------Advanced setting eps 20240724
        !call cg_iric_read_real(cgns_f, "eps", eps, ier)
       ! write (*, '("eps: ", f12.5)') eps
       ! call cg_iric_read_real(cgns_f, "ddt_min_riv", ddt_min_riv, ier)
       ! write (*, '("ddt_min_riv: ", f12.5)') ddt_min_riv
       ! call cg_iric_read_real(cgns_f, "ddt_min_slo", ddt_min_slo, ier)
       ! write (*, '("ddt_min_slo: ", f12.5)') ddt_min_slo
!--Settings for division of unit channels;added 20250304
        call cg_iric_read_integer(cgns_f, "link_divi_switch", link_divi_switch, ier)
        if(link_divi_switch==1) then
        write (*, *) 'Enabled the division of unit channels' 
        call cg_iric_read_integer(cgns_f, "merg_cell_num", merg_cell_num, ier)
        write (*, *) 'The number of intervel cells for division is', merg_cell_num 
        call CG_IRIC_READ_FUNCTIONALSIZE(cgns_f,"division_sec", tmp, ier)
            if(ier==0)then
                sele_l_num = tmp
                mm = sele_l_num
            endif
            if(sele_l_num.le.0) then
                write(*,*) 'No division section has been specified '
                stop
            endif    
            write(*,*) 'Number of division sections:', sele_l_num
             allocate(xtmp(tmp),ytmp(tmp))
            allocate(sele_loc(4, sele_l_num),divi_sec_No(sele_l_num)) 
            call cg_iric_read_functionalwithname(cgns_f, 'division_sec', 'sec_ij', xtmp, ier)
            do k=1, sele_l_num
               divi_sec_No(k) = xtmp(k)
            enddo

            do n=1,4
             write(cm,'(i1)') n
             cordi_label = 'Cordi_'//trim(cm)          
             call cg_iric_read_functionalwithname(cgns_f, 'division_sec', cordi_label, ytmp, ier)          
             do k=1,sele_l_num
                sele_loc(n,k) = ytmp(k)  
             end do                    
            end do
            do k=1,sele_l_num
             write(*,*)"ij coordinate of the upstream cell of No.",divi_sec_No(k) ,"is", sele_loc(1,k), sele_loc(2,k)
             write(*,*)"ij coordinate of the downstream cell of No.",divi_sec_No(k) ,"is", sele_loc(3,k), sele_loc(4,k)
            enddo
            DEALLOCATE(xtmp, STAT = ier)
             DEALLOCATE(ytmp, STAT = ier)

        endif

 !---modified for sedput 20251002  !caution i adn j are opposite in case of iRIC interface
    ! put_v is the bulk sediment volume including voids [m3].
    call cg_iric_read_integer(cgns_f, "j_sedput", j_sedput, ier)
    call cg_iric_read_integer(cgns_f, "j_woodput", j_woodput, ier)
    write (*, '("j_sedput : ", i7)') j_sedput
    if(j_sedput>0) then
        nm_cell =9 !tentative
        allocate (iput(nm_cell+1), jput(nm_cell+1),put_v(nm_cell+1),put_w(nm_cell+1))
        do i = 0, nm_cell
            write(cm,'(i1)') i
            iput_label = 'put_i_'//trim(cm)  
            jput_label = 'put_j_'//trim(cm)  
            put_v_label ='put_v_'//trim(cm)
            put_w_label ='put_w_'//trim(cm)
            call cg_iric_read_integer(cgns_f, iput_label, iput(i+1), ier)      
            call cg_iric_read_integer(cgns_f, jput_label, jput(i+1), ier)     
            call cg_iric_read_real(cgns_f, put_v_label, put_v(i+1), ier)
            call cg_iric_read_real(cgns_f, put_w_label, put_w(i+1), ier)  
        enddo
    endif


 !------Advanced output setting 20240724
    outfile_test = ''
    call cg_iric_read_integer(cgns_f, "outswitch_test", outswitch_test, ier)
    if (outswitch_test == 1)then
        call cg_iric_read_string(cgns_f, "outfile_test", outfile_test, ier)
        write(*,'("outfile_test : ", a)') trim(adjustl(outfile_test))
    end if
!added 20250405    
    outfile_dzb = ''
    call cg_iric_read_integer(cgns_f, "outswitch_dzb", outswitch_dzb, ier)
    if (outswitch_dzb == 1)then
        call cg_iric_read_string(cgns_f, "outfile_dzb", outfile_dzb, ier)
        write(*,'("outfile_dzb : ", a)') trim(adjustl(outfile_dzb))
    end if    

!for slope erosion
    outfile_slope = ''
    call cg_iric_read_integer(cgns_f, "outswitch_slope", outswitch_slope, ier)
    if (outswitch_slope == 1)then
        call cg_iric_read_string(cgns_f, "outfile_slope", outfile_slope, ier)
        write(*,'("outfile_slope : ", a)') trim(adjustl(outfile_slope))
    end if
!added 20250405    
    outfile_dzslo = ''
    call cg_iric_read_integer(cgns_f, "outswitch_dzslo", outswitch_dzslo, ier)
    if (outswitch_dzslo == 1)then
        call cg_iric_read_string(cgns_f, "outfile_dzslo", outfile_dzslo, ier)
        write(*,'("outfile_dzslo : ", a)') trim(adjustl(outfile_dzslo))
    end if  

    outfile_eroslovol = ''
    call cg_iric_read_integer(cgns_f, "outswitch_eroslovol", outswitch_eroslovol, ier)
    if (outswitch_eroslovol == 1)then
        call cg_iric_read_string(cgns_f, "outfile_eroslovol", outfile_eroslovol, ier)
        write(*,'("outfile_eroslovol: ", a)') trim(adjustl(outfile_eroslovol))
    end if    

    outfile_erosupply = ''
    call cg_iric_read_integer(cgns_f, "outswitch_erosupply", outswitch_erosupply, ier)
    if (outswitch_erosupply == 1)then
        call cg_iric_read_string(cgns_f, "outfile_erosupply", outfile_erosupply, ier)
        write(*,'("outfile_erosupply: ", a)') trim(adjustl(outfile_erosupply))
    end if
! for landslide and debris flows
    outfile_mspnt = ''
    call cg_iric_read_integer(cgns_f, "outswitch_mspnt", outswitch_mspnt, ier)
    if (outswitch_mspnt == 1)then
        call cg_iric_read_string(cgns_f, "outfile_mspnt", outfile_mspnt, ier)
        write(*,'("outfile_mspnt: ", a)') trim(adjustl(outfile_mspnt))
    end if

    outfile_LS = ''
    call cg_iric_read_integer(cgns_f, "outswitch_LS", outswitch_LS, ier)
    if (outswitch_LS == 1)then
        call cg_iric_read_string(cgns_f, "outfile_LS", outfile_LS, ier)
        write(*,'("outfile_LS: ", a)') trim(adjustl(outfile_LS))
    end if
       
    outfile_h_surf = ''
    call cg_iric_read_integer(cgns_f, "outswitch_h_surf", outswitch_h_surf, ier)
    if (outswitch_h_surf == 1)then
        call cg_iric_read_string(cgns_f, "outfile_h_surf", outfile_h_surf, ier)
        write(*,'("outfile_h_surf: ", a)') trim(adjustl(outfile_h_surf))
    end if

    debug_inundation_switch = 0
    outfile_debug_inundation = 'debug_inundation.csv'
    outfile_debug_sed_budget = 'debug_sed_budget.csv'
    outfile_debug_spread = 'debug_inundation_spread.csv'
    debug_inundation_i = 114
    debug_inundation_j = 92
    debug_inundation_radius = 2
    debug_inundation_header_written = 0
    debug_sed_budget_header_written = 0
    debug_spread_header_written = 0
    ! Distribution build: keep the point-debug machinery disabled so ordinary
    ! runs do not create debug_inundation/progress/heartbeat CSV files.
    debug_inundation_switch = 0

!added 20260208
    outfile_sdout = ''
    call cg_iric_read_integer(cgns_f, "sd_out_switch", sd_out_switch, ier)
    if (sd_out_switch == 1)then
        call cg_iric_read_string(cgns_f, "outfile_sdout", outfile_sdout, ier)
        write(*,'("outfile_sdout : ", a)') trim(adjustl(outfile_sdout))
        nm_cell =9 !tentative
        allocate (itar(nm_cell+1), jtar(nm_cell+1),itar2(nm_cell+1), jtar2(nm_cell+1))
        do i = 0, nm_cell
            write(cm,'(i1)') i
            i_tar_label = 'tar_i_'//trim(cm)  
            j_tar_label = 'tar_j_'//trim(cm)  
            call cg_iric_read_integer(cgns_f, i_tar_label, itar2(i+1), ier)      
            call cg_iric_read_integer(cgns_f, j_tar_label, jtar2(i+1), ier) 
            ! Caution i and j are opposite in case of iRIC interface
        enddo
    endif

!---tentatively turn off the channel width expansion and over deposition cutting 20250122
    riv_wid_expan_switch = 0 
    !cut_overdepo_switch = 0
!--------------------------------------------------
!iRICの基本機能で対応
!--------------------------------------------------
!read(1,*) hydro_switch
!read(1,'(a)') location_file
!!location_file = trim(rri_dir)//location_file(3:len(location_file))
!
!if(hydro_switch .eq. 1) write(*,'("location_file : ", a)') trim(adjustl(location_file))
!
    write (*, *)

!close(1)
    call cg_iric_close(cgns_f, ier)

! Parameter Check
    do i = 1, num_of_landuse
        if (ksv(i) .gt. 0.d0 .and. ka(i) .gt. 0.d0) &
            stop "Error: both ksv and ka are non-zero."
        if (gammam(i) .gt. gammaa(i)) &
            stop "Error: gammam must be smaller than gammaa."
    end do

! Set da, dm and infilt_limit
    allocate (da(num_of_landuse), dm(num_of_landuse), infilt_limit(num_of_landuse))
    da(:) = 0.d0
    dm(:) = 0.d0
    infilt_limit(:) = 0.d0
    do i = 1, num_of_landuse
        if (soildepth(i) .gt. 0.d0 .and. ksv(i) .gt. 0.d0) infilt_limit(i) = soildepth(i)*gammaa(i)
        if (soildepth(i) .gt. 0.d0 .and. ka(i) .gt. 0.d0) da(i) = soildepth(i)*gammaa(i)
        if (soildepth(i) .gt. 0.d0 .and. ka(i) .gt. 0.d0 .and. gammam(i) .gt. 0.d0) &
            dm(i) = soildepth(i)*gammam(i)
    end do

! if ksg(i) = 0.d0 -> no gw calculation
    gw_switch = 0
    do i = 1, num_of_landuse
        if (ksg(i) .gt. 0.d0) then
            gw_switch = 1
        else
            gammag(i) = 0.d0
            kg0(i) = 0.d0
            fpg(i) = 0.d0
            rgl(i) = 0.d0
        end if
    end do

end subroutine RRI_Read

! Variables used for "sediment parts"

      module sediment_mod

         implicit none
        ! move from RRI_mod.f90 20240724
        real(8), allocatable, save :: kgv(:), tg(:), init_cond_gw(:)
        real(8), allocatable, save :: kgv_idx(:), tg_idx(:), init_cond_gw_idx(:)

!move3 20240422
! added by Qin 2021/6/9
         real(8), allocatable, save :: zb_riv0(:,:), dzb_riv(:,:), zb_riv_slope(:,:), zs_slope(:,:)
         real(8), allocatable, save :: zb_riv0_idx(:), zb_riv_slope_idx(:)
         integer,save:: ave_slope_switch 
         integer,save:: ave_cell_num

         !------parameters added by harada
         integer, save :: link_count,kill_n!20240601
         integer, allocatable, save :: up_riv_idx(:,:)
         integer, allocatable, save :: link_idx(:,:), node_ups(:), link_idx_k(:), link_ups_k(:), link_0th_order(:), link_to_riv(:)
         integer, allocatable, save :: up_riv_lin(:,:),link_ij(:,:)
         real(8), allocatable, save :: width_lin(:), depth_lin(:), height_lin(:), area_ratio_lin(:),width_lin_0(:)! added intial width of nunit channel 20240422
         real(8), allocatable, save :: Link_len(:),area_lin(:)
         real(8), allocatable, save :: zb_riv0_lin(:), zb_riv_slope_lin(:),zb_riv_slope0_lin(:)
         !added by Qin 2021/6/23  
         integer, allocatable, save :: riv_0th(:,:)   !20240724 integer
         integer, allocatable, save :: riv_0th_idx(:) !20240724 integer
         integer, save ::  max_acc_0th_riv  
         !added by Qin 2021/8/18
         integer, save:: min_num_cell_link
         integer, allocatable,save:: link_cell_num(:)  
         !added 2021/11/5
         integer, allocatable, save :: slo_riv_idx(:)
         integer, allocatable, save :: up_slo_idx(:,:) !20220805
!monve3 20240422         
         integer, allocatable, save :: main_riv(:)
         integer, save :: mriv_count
!	 integer, mmm
	     real(8) t_Crit, s, grav, kin_visc, karmans
         real(8) a_const, Ar_const
         real(8) lambda, dlambda
         real(8) t_beddeform_start   !added 20240724
		 
	     integer no_of_sedfiles, no_of_layers , Np, Nl
	     integer, parameter :: Npmax = 20
	     integer, parameter :: Nbmax = 40
	     character*256 infile_sed
	       
         integer outswitch_qsb
         integer outswitch_qss
         integer outswitch_qsw
	     integer outswitch_qdd               !***** added by Chilli	2015/06/18
	     integer outswitch_dzb                !***** added by Chilli	2015/06/18
	     integer outswitch_test              !***** added by Chilli	2015/06/18
	     integer outswitch_sdm
	     integer outswitch_nb
       integer outswitch_dzslo
       integer outswitch_mspnt
       integer outswitch_LS
        integer outswitch_eroslovol
        integer outswitch_erosupply
        integer outswitch_qrs, outswitch_overflowsed !added 202309024
        integer outswitch_slope !for slope erosion
        integer outswitch_h_surf
      

         character*256 outfile_qsb
         character*256 outfile_qss
         character*256 outfile_qsw
         character*256 outfile_qdd           !***** added by Chilli	2015/06/18
         character*256 outfile_dzb            !***** added by Chilli	2015/06/18
         character*256 outfile_test          !***** added by Chilli	2015/06/18
         character*256 outfile_sdm
         character*256 outfile_nb
         character*256 outfile_sumqsb         !***** added by Qin	2021/07/28
         character*256 outfile_sumqss         !***** added by Qin	2021/07/28
         character*256 outfile_dzslo          !added 2021/11/5
         character*256 outfile_mspnt          !added 2021/11/5
         character*256 outfile_mspnt2
         character*256 outfile_LS
         character*256 outfile_eroslovol 
         character*256 outfile_erosupply
         character*256 outfile_qrs
         character*256 outfile_slope
         character*256 outfile_overflow_sed_vol
         character*256 outfile_h_surf !added 20250405
         character*256 outfile_sdout  !added 20260208
      

         character*256 outfile_qsb1
         character*256 outfile_qsb2
         character*256 outfile_qsb3
         character*256 outfile_qsb4
         character*256 outfile_qsb5
         character*256 outfile_qsb6
         character*256 outfile_qsb7
         character*256 outfile_qsb8
         character*256 outfile_qsb9
         character*256 outfile_qsb10
         character*256 outfile_qsb11
         character*256 outfile_qsb12
         character*256 outfile_qsb13
         character*256 outfile_qsb14
         character*256 outfile_qsb15
         character*256 outfile_qsb16
         character*256 outfile_qsb17
         character*256 outfile_qsb18

         character*256 outfile_qss1
         character*256 outfile_qss2
         character*256 outfile_qss3
         character*256 outfile_qss4
         character*256 outfile_qss5
         character*256 outfile_qss6
         character*256 outfile_qss7
         character*256 outfile_qss8
         character*256 outfile_qss9
         character*256 outfile_qss10
         character*256 outfile_qss11
         character*256 outfile_qss12
         character*256 outfile_qss13
         character*256 outfile_qss14
         character*256 outfile_qss15
         character*256 outfile_qss16
         character*256 outfile_qss17
         character*256 outfile_qss18

		 
         integer sed_switch                  !***** added by Chilli 2015/06/18
         character*256 sediment_file	     !***** added by Chilli	2015/06/18

	     integer sed_type_switch
	     real(8) ds_river
	     real(8) Em, Emc, Ed
       real(8) Dmax_relea_sedi
		 
         character*256 ofile_qsb, ofile_qss, ofile_qsw, ofile_qdd, ofile_dzb, ofile_sdm, ofile_nb, ofile_dzslo, ofile_mspnt,ofile_mspnt2, ofile_LS, ofile_eroslovol,ofile_erosupply,ofile_qrs,ofile_overflowsedVol !modified 20240424
         character*256 ofile_qsb1, ofile_qsb2, ofile_qsb3, ofile_qsb4, ofile_qsb5, ofile_qsb6, ofile_qsb7, ofile_qsb8, ofile_qsb9, ofile_qsb10
         character*256 ofile_qsb11, ofile_qsb12, ofile_qsb13, ofile_qsb14, ofile_qsb15, ofile_qsb16, ofile_qsb17, ofile_qsb18
         character*256 ofile_qss1, ofile_qss2, ofile_qss3, ofile_qss4, ofile_qss5, ofile_qss6, ofile_qss7, ofile_qss8, ofile_qss9, ofile_qss10
         character*256 ofile_qss11, ofile_qss12, ofile_qss13, ofile_qss14, ofile_qss15, ofile_qss16, ofile_qss17, ofile_qss18
         character*256 ofile_test
         character*256 ofile_sumqsb, ofile_sumqss, ofile_slope, ofile_h_surf!--added by Qin 2021/07/28
		 
       integer sdt_switch
       integer isedeq, isuseq                     !***** added by Harada	2021/04/20
	     character*256 sdtfile
	     integer, allocatable, save :: sdt_c(:,:),sdt_s(:,:) ! modified for slope erosion

       !----added to read GSD 20240724
       real(8), allocatable, save :: ddist_mm(:)
       REAL(8), DIMENSION(:), ALLOCATABLE :: xtmp, ytmp
       real(8),dimension(:,:),allocatable :: pdist_m_100, pmk0
       real(8),dimension(:,:),allocatable :: pdist_d_100, pdk0
       real(8),dimension(:,:),allocatable :: pdist_slosoil_100
       		 
	     integer sd_out_switch
	     character*256 sdt_location, sdout_file
       character*256 time_output!for yazagyodam
       character(50) :: i_tar_label, j_tar_label !added 20260208
       integer, allocatable,save :: itar(:), jtar(:), itar2(:), jtar2(:)   !added 20260208

	     real(8), allocatable, save :: zb_rock(:,:), Emb(:,:), Et(:,:)
	     real(8), allocatable, save :: zb_roc_idx(:), Emb_idx(:), Et_idx(:)!, zb_air_idx(:)
	     integer, allocatable, save :: Nb_idx(:), Nb(:,:)
       real(8), allocatable, save :: dmean_out(:,:)  ! added 20240724
		 
         type sed_struct                                       !***** added by Chilli	2015/07/01
		 real(8), dimension(Npmax) :: dsed, fm, ft, dzbpr, dzbtem, ffd, fm2, fm1, fms
		 real(8), dimension(Npmax) :: qbi, qsi, qswi, qdi, qdsum, ssi, qsisum, swi, qswisum
		 real(8), dimension(Npmax) :: Esi, Dsi, Es, Ds
		 real(8), dimension(Npmax) :: Ewi, Dwi, Ew, Dw
		 real(8) :: dmean, dmean_slo
		 real(8), dimension(Npmax, Nbmax) :: fd
     real(8), dimension(Npmax) :: fqbi, fqsi,fsur, Esisum, Dsisum ! bedload, suspended load GSD; added by Qin 2021/6/25
     real(8), dimension(Npmax) :: stor_sed_sum_i ! added 20240424
     real(8):: Bed_D90, Bed_D60, samp_dep, stor_sed_sum !added bed material sampling depth 20240419; 20240422, 20240424
	     end type
		 
!--------------------local information added by yorozuya 2016/07/14
	     real(8) perosion
	     real(8) th_em, zm_re, zm_ss
	     integer(8) ibedpro1,ibedpro2,ibedpro3,ibedpro4,ibedpro5,ibedpro6,ibedpro7,ibedpro8
	     real(8) raint, mzbt
	     integer(8) ismooth
	     integer iidt
	     integer(8) isus
           real(8) qtsum
           
!-------------------added by harada for link model
       real(8) qsb_total, qss_total, qsw_total
       real(8),allocatable,save:: dam_sedi_total(:),dam_sedi_total_b(:),dam_sedi_total_s(:),dam_sedi_total_w(:),dam_sedi_totalV(:)!added by Qin
       real(8),allocatable, save:: dam_outflow_qb(:),dam_outflow_qs(:),dam_outflow_cs(:)
       real(8),allocatable,save:: dam_sedi_qsi(:,:)
       !real(8),allocatable,save:: dam_sedi_total2(:),dam_sedi_total_s2(:),dam_sedi_total_w2(:),dam_sedi_totalV2(:)  !---check       
       real(8), allocatable, save :: zb_roc_lin(:), Emb_lin(:), Et_lin(:)!, zb_air_lin(:)
	     integer, allocatable, save :: Nb_lin(:)
       real(8), allocatable, save:: ss_lin(:)
       integer cut_overdepo_switch !added for preventing the over-deposition in channels 20240304

         type sed_struct2                                     
           real(8), dimension(Npmax) :: dsed, fm, ft, dzbpr, dzbtem, ffd, fm2, fm1, fms
           real(8), dimension(Npmax) :: qbi, qsi, qswi, qdi, qdsum, ssi, qsisum, swi, qswisum
           real(8), dimension(Npmax) :: Esi, Dsi, Es, Ds
           real(8), dimension(Npmax) :: Ewi, Dwi, Ew, Dw
           real(8) :: dmean, dmean_slo
           real(8), dimension(Npmax, Nbmax) :: fd
           real(8), dimension(Npmax) :: fqbi, fqsi, fsur, Esisum, Dsisum ! bedload, suspended load GSD; added by Qin 2021/6/25
           real(8), dimension(Npmax) :: stor_sed_sum_i ! added 20240424
           real(8):: Bed_D90, Bed_D60, samp_dep,stor_sed_sum  !added bed material sampling depth 20240419; 20240422, 20240424
         end type
           real(8) samp_dep_min!added minimum sampling depth 20240419 
           
!-------------------added by harada for sediment from slope
         real(8),allocatable,save:: slo_grad(:),slograd(:,:) ,slo_ang_sin(:),slo_ang_cos(:)    !slopeの最急勾酁E         integer,save:: slo_inflow_num     
!--------------------added for slope erosion
         real(8),allocatable,save:: slo_sur_zb(:),slo_sur_zb_before(:)
         real(8),allocatable,save:: fmslo(:,:),fgully(:,:)!, dsi_slo(:,:),Dmslo(:) 20231226
         real(8),allocatable,save:: Ero_mslo(:,:,:), Ero_kslo(:,:),dzslo_idx(:), dzslo(:,:),Ero_slo_vol(:), eroslovol(:,:),dzslo_fp_idx(:,:),inum_sed_dm(:)
         real(8),allocatable,save:: slo_to_lin_sed(:,:),slo_s_dsum(:,:),slo_s_sum(:),slo_to_lin_sed_sum(:),slo_vol_remain(:,:),slo_to_lin_sum_di(:,:),slo_to_lin_sed_sum_idx(:), slo_supply(:,:),lin_to_slo_sed_sum(:),lin_to_slo_sum_di(:,:),lin_to_slo_sed_sum_idx(:),overflow_sed_sum(:,:),slo_to_lin_idx_di(:,:),lin_to_slo_idx_di(:,:)!added 20250503
         real(8),allocatable,save:: B_chan(:),D_chan(:),h_chan(:),  q_chan(:),I_chan(:),l_chan(:),area_chan(:), gully_sedi_dep(:,:) !modified for rill and gully erosion on slope 2022/07/21 20231226
         real(8), allocatable, save:: infil_s_depth(:), infil_w_depth(:) !maximum the infitration depth of the rainfall; addedd 2025/06/18
         real(8),allocatable,save:: qss_slope(:),ss_slope(:)
         real(8),allocatable,save:: dsi(:),w0(:),ns_MS(:)
         real(8),allocatable,save:: slo_qsisum(:,:),slo_qsi(:,:),slo_ssi(:,:),slo_Dsi(:,:),slo_Esi(:,:)
         real(8) slope_erosion_total
         real(8),allocatable,save:: inflow_sedi(:), overflow_sedi(:),overflow_sedi_di(:,:), overdepo_sedi_di(:,:) !modified 20230924
         integer slo_sedi_cal_switch,slope_ero_switch, Ouput_TSAS_switch! modified 20230924
         real(8) dt_slo_sed,surflowdepth,modirate_to_slo_dt  !added 20230427/20230430
 !        real(8) infil_s_depth, infil_w_depth !maximum the infitration depth of the rainfall; addedd 2025/06/18
         !integer,allocatable,save:: landslide(:,:)
          real(8),allocatable,save:: h_surf(:,:),ss_slope_ij(:,:)    !surface flow depth;  added 20231101
          real(8),allocatable,save:: inflow_sedi_ij(:,:), overflow_sedi_ij(:,:)!for slope erosion
!--------------------added for debris
         integer debris_switch, n_LS, LS_num, detail_console
         real(8),allocatable,save:: c_dash(:), hsc(:), pw(:), sf(:)
         real(8),allocatable,save:: vol(:), vcc(:), vcf(:)
         real(8),allocatable,save:: dzslo_mspnt_idx(:), dzslo_mspnt(:,:), soildepth_idx_deb(:)  !added 20240219
         real(8),allocatable,save:: vo_total(:), vo_total_river(:), vo_total_l(:)
         real(8),allocatable,save:: debri_sup_sum(:), debri_sup_sum_di(:,:), debri_sup_sum_ij(:,:) !added 20240424
         real(8) cohe, phi, pwc, pc, pf, b_mp, d_mp_ini
         integer, allocatable, save :: k_LS(:)
         real(8),allocatable,save:: LS_idx(:), LS(:,:), Qg(:), water_v_cell(:), hki_g(:), hki_g_2d(:,:)
         real(8) debris_total, L_rain_ini, wood_total, hki_total
         real(8) T_bebris_off !added for setting the time to switch off debris flow and lanslide computation
         integer debris_end_switch   !20240724
         integer, allocatable, save :: n_link_depth(:)    !0 or 1   added 20240229
         real(8) max_width    !added 20240229
         real(8),allocatable,save:: depth_idx_ini(:),width_idx_ini(:) !added 20240229
!---added for past Landslide event 20240717
         integer past_LS_switch
         character*256 past_ls_file, past_debris_file
         real(8), allocatable,save:: past_ls(:,:), past_Debri_dzb(:,:)
!--------------------added for driftwood
         real(8) Wood_density, qwood_total, C_wood_kd
         integer j_drf
         real(8),allocatable,save:: cw(:), qw(:), vw(:), qwsum(:), hki_area(:), qwsum_total(:)
         real(8),allocatable,save:: vw_idx(:), qwsum_total_idx(:), vw2d(:,:),cw_idx(:), cw2d(:,:), qwsum_2d(:,:)
         integer kkk1
         !--------------------added for sediment put
         integer j_sedput, l_sedput, j_woodput
         real(8) sedput_depth
         real(8),allocatable,save:: fm_sedput(:,:)
         integer, allocatable,save :: iput(:), jput(:), l_put(:), iput2(:), jput2(:)   !------20251002
         ! put_v is bulk sediment volume including voids [m3]; put_w is wood volume [m3].
         ! put_depth and put_depth_remain are equivalent bed thicknesses [m].
         real(8),allocatable,save:: put_v(:), put_depth(:), put_depth_remain(:), put_w(:)  !------20260220
         character(50) :: iput_label,jput_label, put_v_label, put_w_label

!---------------added for division of unit channel
        character*256 infile_divi
        integer link_divi_switch, sele_l_num, merg_cell_num
        !integer, allocatable,save :: sele_up_loci(:),sele_up_locj(:), sele_down_loci(:),sele_down_locj(:), divi_cell(:)
        integer, allocatable,save :: sele_loc(:,:),divi_sec_No(:), divi_cell(:) !modified for IRIC GUI
        character(50):: cordi_label
!-------added for setting some threshhold value for sediment transportation computation
        real(8)  min_slope, max_slope, min_hr,alpha_ss1, alpha_ss2, thresh_ss    
!-------added for river width adjustment; 20240419
         integer riv_wid_expan_switch
         real(8) chan_capa_decre_ratio, max_expan_rate
         real(8), allocatable, save:: wl4wid_expan(:),wid_expa_rate(:)
         !real(8), allocatable, save:: chan_capa4wid_expan(:), depth_uni_chan(:)
!-----added for initial hs hr setting; 20240521
         !real(8) hr0, wc0 !moved 20250314
!-----modified for slope erosion
        integer nm_cell
      	character(50) :: mix_label, gullyB_label,gullyD_label,cm,infil_s_d_label     
        real(8),allocatable,save:: B_gully_r(:), D_gully(:) 
      end module sediment_mod

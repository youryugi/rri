! RRI_Sediment
! calculation of sediment
      
!      subroutine funcd( sed_idx, hr_idxa,  hr_idxa2, qr_ave_idx, ust_idx, qsb_idx, qss_idx, qsw_idx, width_idx, len_riv_idx, zb_riv_idx )
      subroutine funcd( sed_idx, hr_idx, hr_idx2,hr_idxa, qr_ave_idx, ust_idx, qsb_idx, qss_idx, qsw_idx, t,area_idx,water_v_idx )
         use globals
         use sediment_mod
		 use dam_mod!, only:dam_switch,dam_num, dam_loc,damflg ! added by Qin 
         implicit none

         type(sed_struct) sed_idx(riv_count)

	     real(8) hr_idx(riv_count), hr_idx2(riv_count), qr_ave_idx(riv_count), ust_idx(riv_count)
         real(8) hr_idxa(riv_count)
	     real(8) qsb_idx(riv_count), qss_idx(riv_count), qsw_idx(riv_count)
         real(8) qd_sum_idx(riv_count)
	     real(8) dzb_temp(riv_count)
         real(8) area_idx(riv_count),water_v_idx(riv_count) !added by Qin 2021/6/11
         real dzb_cap, ffd_total
         integer t

! 2015/12/29 added by yorozuya
!         real(8) Es_idx(riv_count), Ds_idx(riv_count), Ew_idx(riv_count), Dw_idx(riv_count)
         
         integer k, kk, kkk, m, n, f,nn
         integer kk_div
         real(8) :: u_flax, Frn
         real(8) :: Emb_min

         qd_sum_idx(:) = 0.d0
         Emb_min = Emc

!write(*,*) 't=', t
!pause'in funcd'

!do k = 1, riv_count
!write(*,*) 'k hr_idx=', k, hr_idx(k)
!enddo
!pause

        do k = 1, riv_count
           dzb_temp(k) = 0.d0
	     do m = 1, Np
	       sed_idx(k)%ffd(m) = 0.d0
	       sed_idx(k)%qdi(m) = 0.d0
           sed_idx(k)%qdsum(m) = 0.d0
	     end do
		enddo

!pause'before call bedload'
!         call bedload(sed_idx, hr_idxa, qr_ave_idx , ust_idx, qsb_idx, width_idx, len_riv_idx, zb_riv_idx )
         call bedload(sed_idx, hr_idxa, qr_ave_idx , ust_idx, qsb_idx )

!---calculation the water volume and area of each cell  added by Qin 2021/6/11         
         do k = 1, riv_count
            if (damflg(k).gt.0)then 
                water_v_idx(k) = dam_w_vol(damflg(k))
                area_idx(k) = dam_reserv_area(damflg(k))
            else
                area_idx(k) = dis_riv_idx(k) * width_idx(k)
                water_v_idx(k) = dis_riv_idx(k) * width_idx(k) * hr_idxa(k)
            endif    
            if(isnan(water_v_idx(k))) then !check 2021/6/3
                write(*,*) water_v_idx(k), hr_idxa(k),k
                stop
            endif   
    
        end do	


!pause'before call susload'
!         call susload(sed_idx, hr_idxa, hr_idxa2, ust_idx, qss_idx, qr_ave_idx, width_idx, len_riv_idx)
        call susload(sed_idx, hr_idxa,  ust_idx, qss_idx, qr_ave_idx, t,area_idx,water_v_idx )

    do k = 1, riv_count
	  if(hr_idx(k).eq.0.0) then
       u_flax = 0.0
	   Frn = 0.0
	  else
	   u_flax = qr_ave_idx(k)/width_idx(k)/hr_idxa(k)
	   Frn = u_flax/sqrt(9.81*hr_idxa(k))
	  endif

	if(Frn.ge.100.0) then
!	 if(Frn.ge.1.0) then !Fr?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ»èÔøΩÃèÍçá

!	 pause'Frn>1.0'

	  do m = 1, Np
            sed_idx(k)%qdsum(m) = sed_idx(k)%qdsum(m) - sed_idx(k)%qbi(m) !m3/s
            kk = down_riv_idx(k)
            ! qd_sum minus (flowing into) discharge at the destination cell
            sed_idx(k)%qdsum(m) = sed_idx(k)%qdsum(m) + sed_idx(kk)%qbi(m) !m3/s
	  enddo

	 else !Frn?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ»âÔøΩ?ÅEΩÅEΩÃèÍçá

!	 pause'Frn<1.0'
 	 if (riv_0th_idx(k)== 1)then ! set the sediment supply condition from the 0th order channel;modified by Qin 2021/6/23		
	  do m = 1, Np
!	    do n = 1, 8
 !	     kkk = up_riv_idx(k,n)
!	     if(domain_riv_idx(kkk).eq.0) then
        sed_idx(k)%qdsum(m) = sed_idx(k)%qdsum(m) !unlimit supply (no beddeform at upstream end)
 !       sed_idx(k)%qdsum(m) = sed_idx(k)%qdsum(m)+sed_idx(k)%qbi(m)
		!enddo	
		enddo			  
      else
		do m = 1, Np
			do n = 1,8
				kkk = up_riv_idx(k,n)	
			if(domain_riv_idx(kkk).ne.0) sed_idx(k)%qdsum(m) = sed_idx(k)%qdsum(m) - sed_idx(kkk)%qbi(m) !m3/s
!	     endif
            enddo
             sed_idx(k)%qdsum(m) = sed_idx(k)%qdsum(m) + sed_idx(k)%qbi(m) !m3/s
	  enddo

	  endif
	endif

    enddo !do k = 1, riv_count

!pause'after sediment equation'
!----?ÅEΩÅEΩ|?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩÃåv?ÅEΩÅEΩZ
    do k = 1, riv_count
	  do m = 1, Np
!	   sed_idx(k)%ffd(m) = sed_idx(k)%qdsum(m)/width_idx(k)/dis_riv_idx(k) !?ÅEΩÅEΩP?ÅEΩÅEΩ ÇÔøΩ[m/s]
	   sed_idx(k)%ffd(m) = sed_idx(k)%qdsum(m)/area_idx(k) !?ÅEΩÅEΩP?ÅEΩÅEΩ ÇÔøΩ[m/s]
	  enddo
	enddo

    do k = 1, riv_count
!----?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩV?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩÃëÔøΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩZ
		kk = down_riv_idx(k)
!	if(domain_riv_idx(k)==2.and.domain_riv_idx(kk).ne.0)then ! modified by Qin 2021/6/18
	if (riv_0th_idx(k)==1)then ! no bed deformation at 0th order channel +modified by Qin 2021/6/23	
		do m = 1, Np
		sed_idx(k)%ffd(m) = 0.d0!unlimit supply (no beddeform at upstream end) 2021/6/14
!		sed_idx(k)%ffd(m) = sed_idx(k)%ffd(m)	!unlimted supply of suspended sediment; 2021/7/16
		enddo
	else
		do m = 1, Np
			sed_idx(k)%ffd(m) = sed_idx(k)%ffd(m) + (sed_idx(k)%Esi(m) - sed_idx(k)%Dsi(m)) !?ÅEΩÅEΩP?ÅEΩÅEΩ ÇÔøΩ[m/s]
		enddo	
	endif	

	  do m = 1, Np
	   sed_idx(k)%dzbpr(m) = - ddt * sed_idx(k)%ffd(m) * dlambda !?ÅEΩÅEΩP?ÅEΩÅEΩ ÇÔøΩ[m]
	   dzb_temp(k) = dzb_temp(k) + sed_idx(k)%dzbpr(m) ![m]
	  enddo

!dzb_cap = 0.000002
!if(dzb_temp(k).gt.dzb_cap) then
!do m = 1, Np
!sed_idx(k)%dzbpr(m) = sed_idx(k)%dzbpr(m) * dzb_cap / dzb_temp(k) * 0.01
!sed_idx(k)%ffd(m) = sed_idx(k)%ffd(m) * dzb_cap / dzb_temp(k) * 0.01
!enddo
!dzb_temp(k) = dzb_cap * 0.01
!endif

!?ÅEΩÅEΩÕèÔøΩ?ÅEΩÅEΩœìÔøΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩÂÇ´?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩÅEΩ??ÅEΩÅEΩC?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ}?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩD(?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩÍÇº?ÅEΩÅEΩ?ÅEΩÅEΩÃóÔøΩ?ÅEΩÅEΩa?ÅEΩÅEΩ…ÇÔøΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩÕèÔøΩ?ÅEΩÅEΩœìÔøΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩl?ÅEΩÅEΩ…âÔøΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩK?ÅEΩÅEΩv?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩD)

!	if(dzb_temp(k).gt.mzbt) then
!	  dzb_temp(k) = mzbt
!	endif
!	if(dzb_temp(k).lt.-1.* mzbt) then
!	  dzb_temp(k) = - mzbt
!	endif

    enddo

!pause'before call washload'

 !         call washload(sed_idx, hr_idxa,  ust_idx, qsw_idx, qr_ave_idx, dzb_temp )

  !  do k = 1, riv_count
!----washload??????Z

!	if(zb_riv_idx(k).gt.zm_ss) then
!	  do m = 1, Np
!	   sed_idx(k)%ffd(m) = sed_idx(k)%ffd(m)
!	  enddo
!	else
!	  do m = 1, Np
!	   sed_idx(k)%ffd(m) = sed_idx(k)%ffd(m) + (sed_idx(k)%Ewi(m) - sed_idx(k)%Dwi(m)) !?P???[m/s]
!	  enddo
!	endif

!	 if(zb_riv_idx(k).gt.zm_ss) then !?|???????????l??????B
!	  do m = 1, Np
!	   sed_idx(k)%ffd(m) = sed_idx(k)%qdsum(m)
!	  enddo
!	 else !?|?????{???V???{washload???l??????B
!	  do m = 1, Np
!	   sed_idx(k)%ffd(m) = sed_idx(k)%ffd(m)
!	  enddo
!	 endif

!	enddo

!pause'after call washload'
!pause'endof funcd'
! it was made sure that 'width_idx(k), len_riv_idx(k), dlambda' have good number

!----sum of sediment volume in the dam; Added by Qin 2021/6/12
	if(dam_switch == 1) then
		do f = 1, dam_num
			k = dam_loc(f)
			do n = 1, 8
                kkk = up_riv_idx(k,n)
                if(domain_riv_idx(kkk).eq.0) cycle
				dam_sedi_total_b(f)= dam_sedi_total_b (f) + (qsb_idx(kkk))*ddt 
			!	dam_sedi_total_s2(f) = dam_sedi_total_s2(f) + qss_b(lll)*ddt !---check
				do m = 1, Np ! Get the GSD of the deposited sediment in the dam reservior  modified by Qin 2021/6/22
					dam_sedi_qsi(f,m) = dam_sedi_qsi(f,m) + sed_idx(kkk)%qbi(m)*ddt 
				enddo
			end do	
				
			dam_sedi_total_b(f)= dam_sedi_total_b (f) - qsb_idx(k)*ddt
!			dam_sedi_total_s2(f) = dam_sedi_total_s2(f) - qss_lin(l)*ddt!---check
			do m = 1, Np ! Get the GSD of the deposited sediment in the dam reservior  modified by Qin 2021/6/22
				dam_sedi_qsi(f,m) = dam_sedi_qsi(f,m) - sed_idx(k)%qbi(m)*ddt 
			enddo
			do m = 1, Np
			dam_sedi_total_s(f)= dam_sedi_total_s(f) + (sed_idx(k)%Dsi(m)-sed_idx(k)%Esi(m))*area_idx(k)*ddt 
!			dam_sedi_total_w (f)= dam_sedi_total_w (f) +  (sed_lin(l)%Dwi(m)-sed_lin(l)%Ewi(m))*area_lin(l)*ddt
			dam_sedi_qsi(f,m) = dam_sedi_qsi(f,m) + (sed_idx(k)%Dsi(m)-sed_idx(k)%Esi(m))*area_idx(k)*ddt 
			if(isnan(dam_sedi_total_s(f)))then
				write(*,'(a,f)') "dam_sedi_total_s=", dam_sedi_total_s (f)
				write(*,'(a,f)') "Dsi=", sed_idx(k)%Dsi(m)
				write(*,'(a,f)') "Esi=", sed_idx(k)%Esi(m)
				stop
			endif
			enddo
			dam_sedi_total(f) = dam_sedi_total_b(f) + dam_sedi_total_s(f)
			dam_sedi_totalV(f) = dam_sedi_total(f)* dlambda
!			dam_sedi_total2(f) = dam_sedi_total_b(f) + dam_sedi_total_s2(f) !---check
	!		dam_sedi_total2(f) = dam_sedi_total_b(f) + dam_sedi_total_s2(k) !---check2021/5/28
	!		dam_sedi_totalV2(f) = dam_sedi_total2(f)*dlambda!---check
		enddo	

	endif		

      end subroutine funcd

! >>>>>>>>>>>>>>>>>>>>>>> CALCULATING wash laod
!      subroutine washload(sed_idx, hr_idxa, hr_idxa2, ust_idx, qsw_idx, qr_ave_idx, dzb_temp, width_idx, len_riv_idx)
      subroutine washload(sed_idx, hr_idxa, ust_idx, qsw_idx, qr_ave_idx, dzb_temp )
         use globals
         use sediment_mod
         implicit none
		 
         type(sed_struct) sed_idx(riv_count)
         real(8) hr_idx(riv_count), hr_idx2(riv_count), ust_idx(riv_count)
         real(8) hr_idxa(riv_count)
         real(8) qsw_idx(riv_count), qr_ave_idx(riv_count)
	     real(8) ss_idx(riv_count)
         real(8) qsw_sum_idx(riv_count)
         real(8) dzb_temp(riv_count)
         integer k, kk, kkk, n, mm, mmax, m
	     real(8) :: ss11, ss1, ss2, ss3, ss22
	     real(8) :: ss4, ss5, DDD
	     real(8) :: u_flax, Frn
	     integer mmm

		DDD = 0.01
!pause'in washload'

	    do k = 1, riv_count
!		 E_idx(k) = 0.0
!		 D_idx(k) = 0.0
		 qsw_idx(k) = 0.0
		 ss_idx(k) = 0.0
	 	enddo

    do k = 1, riv_count
		do m = 1, Np
		 if(qr_ave_idx(k).le.0.0) then
		  sed_idx(k)%swi(m) = 0.0
		 else
		  sed_idx(k)%swi(m) = sed_idx(k)%qswi(m) / qr_ave_idx(k)
		 endif
		enddo

		do m = 1, Np
		  ss_idx(k) = ss_idx(k) + sed_idx(k)%swi(m)
		enddo

		 if(ss_idx(k).gt.0.1) then
		  do m = 1, Np
		    sed_idx(k)%swi(m) = sed_idx(k)%swi(m) / ss_idx(k)* 0.1
		  enddo
		 endif
    enddo

    do k = 1, riv_count
         do m = 1, Np
		   sed_idx(k)%qswisum(m) = 0.0
		 enddo
    enddo

!pause'before call washsedi'
	call washsedicom(sed_idx, dzb_temp, ust_idx)
!pause'after call washsedi'

    do k = 1, riv_count

		   u_flax = qr_ave_idx(k)/width_idx(k)/hr_idxa(k)
		   Frn = u_flax/sqrt(9.81*hr_idxa(k))

		 if(Frn.ge.100.0) then 
!		 if(Frn.lt.1.0) then !!Fr?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ»èÔøΩÃèÍçá

		 do m = 1, Np
		   sed_idx(k)%qswisum(m) = sed_idx(k)%qswisum(m) - sed_idx(k)%qswi(m)
            kk = down_riv_idx(k)
           	! qd_sum minus (flowing into) discharge at the destination cell
		   sed_idx(k)%qswisum(m) = sed_idx(k)%qswisum(m) + sed_idx(kk)%qswi(m)
		 enddo

		 else !Fr?????????

		 do m = 1, Np
	      do n = 1, 8
 	        kkk = up_riv_idx(k,n)
	        if(domain_riv_idx(kkk).eq.0) then
		     sed_idx(k)%qswisum(m) = sed_idx(k)%qswisum(m)
            else
		     sed_idx(k)%qswisum(m) = sed_idx(k)%qswisum(m) - sed_idx(kkk)%qswi(m)
		    endif
          enddo
		   sed_idx(k)%qswisum(m) = sed_idx(k)%qswisum(m) + sed_idx(k)%qswi(m)
		 enddo

		endif
    enddo !end of do k = 1, riv_count

    do k = 1, riv_count
	
		 do m = 1, Np
		  ss11 = sed_idx(k)%swi(m)
!		  ss1 = ss11 * hr_idxa2(k)
!		  ss2 = -1.* sed_idx(k)%qswisum(m)
      	 if(domain_riv_idx(k).eq.2) sed_idx(k)%qswisum(m) = 0.0
		  ss22 = ss2/width_idx(k)/dis_riv_idx(k)
		  ss3 = sed_idx(k)%Ewi(m) - sed_idx(k)%Dwi(m)
         if(hr_idxa(k).le.min_hr) then
		  sed_idx(k)%swi(m) = 0.0
		 else
		  sed_idx(k)%swi(m) = sed_idx(k)%swi(m) + (-1.* sed_idx(k)%qswisum(m)/width_idx(k)/dis_riv_idx(k) + sed_idx(k)%Ewi(m) - sed_idx(k)%Dwi(m)) * ddt/hr_idxa(k)
		  if(sed_idx(k)%swi(m).le.0.0) then
		   sed_idx(k)%Dwi(m) = -1.* sed_idx(k)%qswisum(m)/width_idx(k)/dis_riv_idx(k) + sed_idx(k)%Ewi(m)
           sed_idx(k)%swi(m) = 0.0
		  endif
		 endif

!         if(sed_idx(k)%swi(m).lt.0.0) then
!		  sed_idx(k)%swi(m) = 0.0
!		 endif

if(abs(ust_idx(k)).gt.1000.001.and.m.eq.1) then
write(*,'(a)') '    k    m        ss11        ss22 sed_idx%Ewi sed_idx%Dwi         ss3 sed_idx%swi'
write(*,'(2i5,6e12.3)') k, m, ss11, ss22, sed_idx(k)%Ewi(m), sed_idx(k)%Dwi(m), ss3, sed_idx(k)%swi(m)
pause
endif
		 enddo
!write(*,'(a)') 'k, a2(k), hr_idxa(k)'
!pause
		 
!		if(ust_idx(k).ge.0.1) then
		  do m = 1, Np
		   sed_idx(k)%qswi(m) = sed_idx(k)%swi(m) * qr_ave_idx(k) !m3/s
		   qsw_idx(k) = qsw_idx(k) + sed_idx(k)%qswi(m) !m3/s
		  enddo
!		else
!		  do m = 1, Np
!		   sed_idx(k)%qswi(m) = 0.0
!		   sed_idx(k)%swi(m) = 0.0
!		  enddo
!		endif		  
		  
		  
		  
		  
if(abs(ust_idx(k)).gt.10000.0) then
write(*,'(a)') '    m         ss1        ss22         ss3 sed_idx%swi'
write(*,'(a,i)') 'k=           ', k
write(*,'(a,f)') 'ddt=         ', ddt
write(*,'(a,f)') 'dzb_temp(k)= ', dzb_temp(k)
write(*,'(a,f)') 'hr_idxa(k)=  ', hr_idxa(k)
write(*,'(a,f)') 'ust_idx(k)=  ', ust_idx(k)
write(*,'(a,f)') 'qr_ave_idx(k)=  ', qr_ave_idx(k)
write(*,'(a,f)') 'qsw_idx(k)=  ', qsw_idx(k)
write(*,'(a)') '  m      dsi(m)          fp      swi(m)      Dwi(m)      Ewi(m)   sedqsw(m)'
do m = 1, Np
write(*,'(i3,6f12.5)') m, sed_idx(k)%dsed(m), sed_idx(k)%fm(m), sed_idx(k)%swi(m), sed_idx(k)%Dwi(m), sed_idx(k)%Ewi(m), sed_idx(k)%qswi(m)
enddo
pause'in washsedicom'
endif

    enddo

!pause'end of washload'

end subroutine washload

! >>>>>>>>>>>>>>>>>>>>>>> CALCULATING washload sediment
      subroutine washsedicom(sed_idx, dzb_temp, ust_idx)
     use globals
     use sediment_mod
     implicit none
		 
	 type(sed_struct) sed_idx(riv_count)
	 real(8) dzb_temp(riv_count), ust_idx(riv_count)
     real(8) ss
     real(8) dsm, fp, ramdaw
	 integer k, m, iwcheck
     real(8) Caew, dzb


	 ramdaw = 0.4

!pause'in washsedi'

         do k = 1, riv_count
! For non-uniform sized sediment
!         	if(ust_idx(k).gt.0.0) then
		     do m = 1, Np
                dzb = sed_idx(k)%dzbpr(m)
		      if (dsi(m).le.0.0001) then ! di=< 0.1mm??wash load?????????

			    fp = sed_idx(k)%fm(m)
			    ss = sed_idx(k)%swi(m)

		       if(dzb_temp(k).le.0.d0) then
!   		       if(dzb_temp(k).le.0.0.and.ust_idx(k).ge.0.12) then
		        Caew = -(1.-ramdaw)*fp*dzb_temp(k)/ddt
		        sed_idx(k)%Ewi(m) = Caew
		       else
		        sed_idx(k)%Ewi(m) = 0.d0
               endif
			  sed_idx(k)%Dwi(m) = ss * w0(m)

		      else !middle of if (dsi.le.0.0001) then
			    sed_idx(k)%Ewi(m) = 0.d0
			    sed_idx(k)%Dwi(m) = 0.d0
		      end if !end of if (dsi.le.0.0001) then

		     end do

if(abs(dzb_temp(k)).gt.100000000000.0) then
write(*,'(a,i)') 'k=             ', k
write(*,'(a,f)') 'ddt=           ', ddt
write(*,'(a,f)') 'dzb_temp(k)=  ', dzb_temp(k)
write(*,'(a,f)') 'Emb_idx(k)=   ', Emb_idx(k)
write(*,'(a)') '  m      dsi(m)          fp      swi(m)      Dwi(m)      Ewi(m)'
do m = 1, Np
write(*,'(i3,5f12.5)') m, sed_idx(k)%dsed(m), sed_idx(k)%fm(m), sed_idx(k)%swi(m), sed_idx(k)%Dwi(m), sed_idx(k)%Ewi(m)
enddo
pause'in washsedicom'
endif

 
        end do ! end of do k = 1, riv_count

!pause'after washsedi'

      end subroutine washsedicom
	  
! >>>>>>>>>>>>>>>>>>>>>>> CALCULATING suspended laod
!      subroutine susload(sed_idx, hr_idxa, hr_idxa2, ust_idx, qss_idx, qr_ave_idx, width_idx, len_riv_idx)
     subroutine susload(sed_idx, hr_idxa,ust_idx, qss_idx, qr_ave_idx, t,area_idx,water_v_idx)
     use globals
     use sediment_mod
     use dam_mod
     implicit none

	 type(sed_struct) sed_idx(riv_count)	      
     real(8) hr_idx(riv_count), hr_idx2(riv_count), ust_idx(riv_count)
     real(8) hr_idxa(riv_count)
     real(8) qss_idx(riv_count), qr_ave_idx(riv_count)
	 real(8) ss_idx(riv_count)
!	 real(8) width_idx(riv_count), len_riv_idx(riv_count)
     real(8) area_idx(riv_count), water_v_idx(riv_count) !Added by Qin 2021/6/12
     integer k, kk, kkk, n, m
	 real(8) :: ss1, ss11, ss2, ss3, ss22
!    real(8) :: ss4, ss5, DDD
!	 real(8) :: qsss1
	 integer ktest
	 real(8) :: u_flax, Frn

     real(8) :: bbeta, ss
     real(8) :: Cae, Cae1, Csb,Ce
     real(8) :: fp
	 real(8) :: z_star
     integer kkmax, i, j, t
     real(8) :: alfa2  !for bedlock
     real(8) :: ristar, Emb_min
	 real(8) :: thet !2021/7/18
 
     
     real(8), parameter :: pai = 3.141592
     real(8), parameter :: csbar = 0.2
     real(8), parameter :: kappa = 0.0015
     real(8), parameter :: rw = 1.0
     real(8), parameter :: rs = 2.65
 
 
 !	 DDD = 0.01
 Emb_min = Emc  
       do k = 1, riv_count
		 qss_idx(k) = 0.0
		 ss_idx(k) = 0.0
        do m = 1, Np
		 sed_idx(k)%qsisum(m) = 0.0
		 sed_idx(k)%ssi(m) = 0.0
		enddo
	   enddo
!pause'1'
    do k = 1, riv_count
		do m = 1, Np
!		 if(qr_ave_idx(k).le.0.0) then
         if(qr_ave_idx(k).le.1e-3) then !modified by Qin 2021/6/12          
		  sed_idx(k)%ssi(m) = 0.0
		 else
		  sed_idx(k)%ssi(m) = sed_idx(k)%qsi(m) / qr_ave_idx(k)
		 endif
		enddo
!pause'2'
		do m = 1, Np
		  ss_idx(k) = ss_idx(k) + sed_idx(k)%ssi(m)
		enddo
!pause'3'
	    if(ss_idx(k).gt.alpha_ss2) then
		  do m = 1, Np
		   sed_idx(k)%ssi(m) = sed_idx(k)%ssi(m) / ss_idx(k) * alpha_ss2
		  enddo
		endif

	enddo ! end of do k = 1, riv_count

!pause'before deposition yoro-work'

!------------------calculation of Erosion and Deposition of suspended sediment
!-----------------Calculation of Deposition and Erosion
    do k = 1, riv_count
		thet = zb_riv_slope_idx(k) !modified by Qin 2021/7/18
     do m = 1, Np


!write(*,'(a,i,f)') 'm/dsi', m, dsi

		     fp = sed_idx(k)%fm(m)
		     ss = sed_idx(k)%ssi(m)

!	    if (dsi.gt.0.0001) then
!		di=> 0.1mm??sus load?????????

		   if(ust_idx(k).eq.0.0) then
!		   if(hr_idxa(k).le.0.001) then
		     sed_idx(k)%Dsi(m) = ss * w0(m)
		   else
  		     bbeta = 6.0*w0(m)/karmans/ust_idx(k)
!		     bbeta = w0 * hr_idxa(k) * 1.0E-4
		     Csb = ss*bbeta/(1.-exp(-1. * bbeta))
!			if(Csb.gt.3.0*ss) Csb = ss * 3.0 !??[???????????z??????????????}????B?????????R?????
!			if(Csb.lt.0.0) Csb = 0.0
		     sed_idx(k)%Dsi(m) = Csb * w0(m)
		   endif

         if(isuseq==1) then
            !Lane-Kalinske formula  
           if (w0(m) .lt. ust_idx(k)) then
!           if (w0 .lt. ust_idx(k).and.ust_idx(k).ge.0.12) then
		      z_star = ust_idx(k)/w0(m)
!		      Cae1 = 5.55*(0.5*z_star*exp(-(1./z_star)))**1.61 * 1.0E-4
		      Cae1 = 5.55*(0.5*z_star*exp(-(1./z_star)**2.))**1.61 * 1.0E-4 

				if(Cae1.gt.0.1d0) Cae1 = 0.1d0
		      Cae = Cae1 *fp
           else
	          Cae = 0.d0
           endif
		     sed_idx(k)%Esi(m) = Cae * w0(m)
        
        elseif(isuseq == 2) then 
            ! Density straitified flow
            u_flax = qr_ave_idx(k)/width_idx(k)/hr_idxa(k)			
			if(u_flax > 10.) u_flax = 10.   !check this reguration	
			if (w0(m) .ge. ust_idx(k)) then
				sed_idx(k)%Esi(m) = 0.d0
			elseif(ss >= 0.5*csbar) then
				sed_idx(k)%Esi(m) = 0.d0
			elseif(hr_idxa(k) < min_hr .and. (u_flax**3./hr_idxa(k))>1.)then
				sed_idx(k)%Esi(m) = fp*kappa/((rs/rw-1.)*grav)
			else
!				ristar = (rs/rw-1.)*(csbar-ss)*grav*hr_idxa(k)/(u_flax**2.)
!				if( ristar<0.01 ) ristar = 0.01
!				sed_idx(k)%Esi(m) = fp*kappa/ristar*u_flax*(csbar-ss)
!Richardson number is calculated under the assumption as Csbar>>ss; modified by Qin 2021/7/9 
				ristar = (rs/rw-1.)*csbar*grav*hr_idxa(k)/(u_flax**2.)
				Ce = kappa/ristar*u_flax*csbar/w0(m)
					if(Ce.gt.0.1d0) Ce = 0.1d0

				sed_idx(k)%Esi(m) = fp*Ce*w0(m)
				!deposition modified 
			    sed_idx(k)%Dsi(m) = ss * w0(m)  
			endif

        endif ! end for for formula selection  
        !bedlock part
		  if(zb_riv_idx(k) - zb_roc_idx(k) .lt. Emb_idx(k))then 
			alfa2 = (zb_riv_idx(k) - zb_roc_idx(k))/Emb_idx(k)
			if(alfa2.ge.1.)alfa2=1.
			if(alfa2 <= Emb_min)alfa2 = 0.
			sed_idx(k)%Esi(m) = sed_idx(k)%Esi(m) * alfa2
!			if(zb_riv_idx(k) - zb_roc_idx(k) < 0.)then
!				write(*,'(2i5,3f17.10)')l, m, zb_riv_idx(k) - zb_roc_idx(k), alfa2, sed_idx(k)%Esi(m)
!			endif 
	   	  end if

!		else
!		  sed_idx(k)%Dsi(m) = 0.d0
!		  sed_idx(k)%Esi(m) = 0.d0
!		endif

if(sed_idx(k)%Esi(2).gt.1.0E+2) then
!if(ust_idx(k).gt.0.01) then
write(*,'(i3,4f9.5,4e15.5)') m, dsi(m), fp, w0(m), sed_idx(k)%ssi(m), Csb, sed_idx(k)%Dsi(m), Cae, sed_idx(k)%Esi(m)
endif
 
 enddo !end of do m = 1, Np

if(sed_idx(k)%Esi(2).gt.1.0E+2) then
!if(ust_idx(k).gt.0.01) then
write(*,'(a)') '  m      dsi       fp       w0      ssi            Csb           Dsi            Cae            Esi'
write(*,'(a,f)') 'ust=', ust_idx(k)
write(*,'(a,i)') 't=', t
!kkmax = int(real(k)/real(idn))
!i = riv_idx2i(kkmax)
!j = riv_idx2j(kkmax)
write(*,'(a,i7)') 'k=', k
pause'check deposition'
endif

enddo !end of do k = 1, riv_count

!pause'after deposition'

!--------------------trasportation of suspended sediment
   do k = 1, riv_count
!pause'in do roope'        	
		 do m = 1, Np
		  sed_idx(k)%qsisum(m) = 0.0
		 enddo
		  ss_idx(k) = 0.0
!pause'before write'
!write(*,*) k, qr_ave_idx(k), width_idx(k), hr_idxa(k)

          u_flax = qr_ave_idx(k)/width_idx(k)/hr_idxa(k)
	      Frn = u_flax/sqrt(9.81*hr_idxa(k))
!pause'after Frn'
	    if(Frn.ge.100.0) then !Frn????????
!	    if(Frn.lt.1.0) then !Frn????????

		 do m = 1, Np
		   sed_idx(k)%qsisum(m) = sed_idx(k)%qsisum(m) - sed_idx(k)%qsi(m)
   	       if(domain_riv_idx(k).eq.2) cycle
       		kk = down_riv_idx(k)
           	! qd_sum minus (flowing into) discharge at the destination cell
		   sed_idx(k)%qsisum(m) = sed_idx(k)%qsisum(m) + sed_idx(kk)%qsi(m)
		 enddo

	    else !Frn?????????
		if(riv_0th_idx(k)==1)then !modified by Qin 2021/7/18
			do m =1, Np
				sed_idx(k)%qsisum(m) = sed_idx(k)%qsisum(m)
			enddo
		else			

		 do m = 1, Np
	        do n = 1, 8
 	         kkk = up_riv_idx(k,n)
	         if(domain_riv_idx(kkk).eq.0) then
		       sed_idx(k)%qsisum(m) = sed_idx(k)%qsisum(m)
             else
		       sed_idx(k)%qsisum(m) = sed_idx(k)%qsisum(m) - sed_idx(kkk)%qsi(m)
             endif
		    enddo
		     sed_idx(k)%qsisum(m) = sed_idx(k)%qsisum(m) + sed_idx(k)%qsi(m)
		 enddo
		endif

		endif !Frn??????if??

   enddo ! end of do k = 1, riv_count

!pause'sum of ssload'

   do k = 1, riv_count
	kk = down_riv_idx(k)
		 do m = 1, Np
!	     if(domain_riv_idx(k).eq.2) cycle
		ss11 = sed_idx(k)%ssi(m)
!		if(domain_riv_idx(k).eq.2)  sed_idx(k)%qsisum(m) = 0.d0
		if(domain_riv_idx(k).eq.2.and.domain_riv_idx(kk)==0)  sed_idx(k)%qsisum(m) = 0.d0 !no suspended load at downstream end check 2021/6/18
		if(water_v_idx(k).lt.0.01.or.qr_ave_idx(k).le.1e-3) then    !check this regulation
!		if(qr_ave_idx(k).le.1e-3) then !modified by Qin 2021/6/3
			sed_idx(k)%ssi(m) = 0.d0
			sed_idx(k)%Esi(m) = 0.d0 !modified by Qin 2021/6/3
			sed_idx(k)%Dsi(m) =-1.*sed_idx(k)%qsisum(m)/area_idx(k)
			if(sed_idx(k)%Dsi(m).lt.0.d0) then
				sed_idx(k)%Dsi(m) = 0.d0
				sed_idx(k)%qsisum(m) = 0.d0
			endif
			if (isnan(sed_idx(k)%Dsi(m))) then
				write(*,*) sed_idx(k)%qsisum(m), area_idx(k)
				stop
			endif
		else !-----revised by Harada 2021/5/21	
			sed_idx(k)%ssi(m) = sed_idx(k)%ssi(m) + (-1.* sed_idx(k)%qsisum(m) + sed_idx(k)%Esi(m)*area_idx(k) - sed_idx(k)%Dsi(m)*area_idx(k)) * ddt/water_v_idx(k)
			if(sed_idx(k)%ssi(m).le.0.0) then
			sed_idx(k)%ssi(m) = sed_idx(k)%ssi(m) - (-1.* sed_idx(k)%qsisum(m) + sed_idx(k)%Esi(m)*area_idx(k) - sed_idx(k)%Dsi(m)*area_idx(k)) * ddt/water_v_idx(k) !?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ…ñﬂÇÔøΩ?ÅEΩÅEΩ?ÅEΩÅEΩ
			sed_idx(k)%Dsi(m) = -1.* sed_idx(k)%qsisum(m)/area_idx(k) + sed_idx(k)%Esi(m)      !D?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩC?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ
			if(sed_idx(k)%Dsi(m) < 0.) sed_idx(k)%Dsi(m) = 0.d0 !?ÅEΩÅEΩ‹ÇÔøΩD?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩÃéÔøΩ?ÅEΩÅEΩÕè„â∫?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩÃÉA?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩo?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩX?ÅEΩÅEΩ…ÇÔøΩ?ÅEΩÅEΩÃÇÔøΩD?ÅEΩÅEΩ?ÅEΩÅEΩ0?ÅEΩÅEΩ∆Ç»ÇÔøΩ
			sed_idx(k)%ssi(m) = sed_idx(k)%ssi(m) + (-1.* sed_idx(k)%qsisum(m) + sed_idx(k)%Esi(m)*area_idx(k) - sed_idx(k)%Dsi(m)*area_idx(k)) * ddt/water_v_idx(k) !?ÅEΩÅEΩV?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩD?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩg?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩƒÇÔøΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩxssi?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩv?ÅEΩÅEΩZ
		 !----modified by Qin 2021/5/30
			 if(sed_idx(k)%ssi(m).lt.0.d0) then
				sed_idx(k)%ssi(m) = sed_idx(k)%ssi(m) - (-1.* sed_idx(k)%qsisum(m) + sed_idx(k)%Esi(m)*area_idx(k) - sed_idx(k)%Dsi(m)*area_idx(k)) * ddt/water_v_idx(k) !?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ…ñﬂÇÔøΩ?ÅEΩÅEΩ?ÅEΩÅEΩ
				sed_idx(k)%Esi(m) = sed_idx(k)%qsisum(m)/area_idx(k)!E?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩC?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ
!				sed_idx(k)%ssi(m) = sed_idx(k)%ssi(m) + (-1.* sed_idx(k)%qsisum(m) + sed_idx(k)%Esi(m)*area_idx(k) - sed_idx(k)%Dsi(m)*area_idx(k)) * ddt/water_v_idx(k) !?ÅEΩÅEΩV?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩE?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩg?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩƒÇÔøΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩxssi?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩv?ÅEΩÅEΩZ	 
				sed_idx(k)%ssi(m) = 0.d0
			endif	
			  if (sed_idx(k)%ssi(m)<0.d0) then
				write(*,'(a,f)') "k=", k
				write(*,'(a,f)') "m=",m
				write(*,'(a,f)') "ssi=", sed_idx(k)%ssi(m)
				stop"suspended load concentration < 0"	
			  endif
		   endif	  
		if (isnan(sed_idx(k)%ssi(m))) then
			write(*,*) k, m, sed_idx(k)%ssi(m),water_v_idx(k)
			stop "ssi is Nan"
		endif
		endif
		 
if(sed_idx(k)%Esi(m).gt.1.0E+8.and.k.eq.20) then
!if(t.ge.83) then
write(*,'(i5,7e12.3)') m, ss11, -1.* sed_idx(k)%qsisum(m), sed_idx(k)%Esi(m), sed_idx(k)%Dsi(m), sed_idx(k)%Esi(m) - sed_idx(k)%Dsi(m), sed_idx(k)%ssi(m), sed_idx(k)%ssi(m)- ss11
endif

		enddo !end of do m = 1, Np
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(sed_idx(k)%Esi(2).gt.1.0E+8.and.k.eq.20) then
!if(t.ge.83) then
write(*, '(a)') '    m        ss11        ss22         Esi         Dsi         ss3 sed_%ssi(m)        dssi'
write(*,'(a,2f)') 'hr_idxa(k)=',  hr_idxa(k)
write(*,'(a,i7)') 'k=', k
pause
endif		  

!		if(ust_idx(k).ge.0.1) then
		  do m = 1, Np
		   sed_idx(k)%qsi(m) = sed_idx(k)%ssi(m) * qr_ave_idx(k) !m3/s
!			if(thet>12.) sed_idx(k)%qsi(m) = 0.d0 !modified by Qin 2021/7/18
		   qss_idx(k) = qss_idx(k) + sed_idx(k)%qsi(m) !m3/s
		  enddo
!		else
!		  do m = 1, Np
!		   sed_idx(k)%qsi(m) = 0.d0
!		   sed_idx(k)%ssi(m) = 0.d0
!		  enddo
!		endif
		  
if(sed_idx(k)%Esi(2).gt.1.0E+8) then
write(*,'(a)') '    ust_idx    hr_idxa      qraveidx     qss_idx      ss_idx'
write(*,'(6e12.5)') ust_idx(k), hr_idxa(k),  qr_ave_idx(k), qss_idx(k), ss_idx(k)
write(*,'(a,i)') 't=', t
write(*,'(a,i7)') 'k=', k
write(*,*) 'ddt =', ddt
write(*,*)
write(*,'(a)') '    m    sed_idx%qsm    sed_idx%Esi    sed_idx%Dsi    sed_idx%ssi    sed_idx%qsi     sed_idx%fm'
do m = 1, Np
write(*,'(i5,6e15.5)') m, sed_idx(k)%qsisum(m), sed_idx(k)%Esi(m), sed_idx(k)%Dsi(m), sed_idx(k)%ssi(m), sed_idx(k)%qsi(m), sed_idx(k)%fm(m)
enddo
pause
endif

    enddo ! end of do k = 1, riv_count

!pause'end of suspended sediment'

end subroutine susload

! >>>>>>>>>>>>>>>>>>>>>>> CALCULATING BED LOAD
!      subroutine bedload( sed_idx, hr_idxa, qr_ave_idx, ust_idx, qsb_idx, width_idx, len_riv_idx, zb_riv_idx )
      subroutine bedload( sed_idx, hr_idxa, qr_ave_idx, ust_idx, qsb_idx )
	     use globals
	     use sediment_mod
         use dam_mod
	     implicit none

	     type(sed_struct) sed_idx(riv_count)
	     real(8) ust_idx(riv_count)		 
	     real(8) hr_idx(riv_count), qr_ave_idx(riv_count)
	     real(8) hr_idxa(riv_count)
	     real(8) qsb_idx(riv_count)
!	     real(8) width_idx(riv_count), len_riv_idx(riv_count)
!	     real(8) zb_riv_idx(riv_count)
	     real(8) :: fp, u_flax, ks, lg
	     integer :: k, m, kk

	     real(8) :: uscm, usci, u_se, qbi
	     real(8) :: tscm,  dsm
	     real(8) :: t_star, t_starE, qbe, den
	     real(8) :: R_tcts
	     real(8) :: zb55, Frn

	     real(8) :: ust_idx1, ust_idx22, ust_idx33, slope, slope11, h_dam
	     integer kkk, kkkk, kzb55, n
	     real(8) :: slope22, slope33
!----Egashira Formula-----added by Qin 2021/6/11
         real(8) :: S_angle, C_K1, C_K2, C_fd, C_ff 
         real(8) :: Dcr, alfa, hs_l
         real(8) :: thet, thetr, Emb_min

         real(8), parameter :: C_kf		= 0.16
         real(8), parameter :: C_kd		= 0.0828
         real(8), parameter :: C_e		= 0.85
         real(8), parameter :: Csta		= 0.6
 
         Emb_min = Emc		
         S_angle = 34./180.d0*3.14159d0
 
		 
!do k = 10, riv_count1 + 9
!write(*,*) 'k /sed_idxdmean=', k, (k), sed_idx(k)%dmean
!write(*,*) 'k Emb_idx(k)=', k, Emb_idx(k)
!enddo
!pause'in bedload'

!-----------?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ ÇÃçƒåv?ÅEΩÅEΩZ?ÅEΩÅEΩn?ÅEΩÅEΩ?ÅEΩÅEΩ----------------

!write(*,*) 'tan(45.)=', tan(45.*3.14159/180.)

do k = 1, riv_count
 slope11 = tan(zb_riv_slope_idx(k)*3.14159/180.)
 if(slope11.le.min_slope) slope11 = min_slope
 !if(qr_ave_idx(k) .le. 1e-3) qr_ave_idx(k) = 1e-3
 if(qr_ave_idx(k) .le. 0.d0) qr_ave_idx(k) = 0.d0	!modified by Qin 2021/6/3
 hr_idxa(k) = ( ns_river * qr_ave_idx(k) / width_idx(k) / sqrt(slope11) )**(3./5.)
 if(hr_idxa(k).le.0.d0) hr_idxa(k) = 0.d0 !modified by Qin 2021/6/3
 if(isnan(hr_idxa(k)))then
	write(*,*) k, slope11, zb_riv_slope_idx(k),width_idx(k), qr_ave_idx(k)
	stop "Depth is NaN"
 endif	
enddo
!write(*,'(a)') '??????v?Z'
!write(*,'(a)') 'k, qr_ave_idx(kk), width_idx(kk), hr_idx(kk), (kk), zb_riv_slope_idx(kk)'
!do k = 1, 90
!kk = main_riv(k)
!write(*,'(i5,5f12.3)') k, qr_ave_idx(kk), width_idx(kk), hr_idx(kk), (kk), zb_riv_slope_idx(kk)
!enddo

!----set the water level of dam reservoir same as the average level of its neighbouring upstream cells; added by Qin 2021/6/11

do k = 1, riv_count
    if(damflg(k).gt.0)then 
		h_dam = 0.d0
		m = 0
		do n = 1,8
			kkk= up_riv_idx(k,n)
			if (domain_riv_idx(kkk).eq.0) cycle
			m = m+1
			h_dam = h_dam+ hr_idxa(kkk)+zb_riv_idx(kkk)
		enddo
		hr_idxa(k) = h_dam/real(m)-zb_riv_idx(k)
		if(hr_idxa(k).le.0.d0) hr_idxa(k) = 0.d0 
	endif		

enddo   


!pause'in bedload'
!-----------?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ ÇÃçƒåv?ÅEΩÅEΩZ?ÅEΩÅEΩI?ÅEΩÅEΩ?ÅEΩÅEΩ----------------
!pause'1'

do k = 1, riv_count
    
	! Determine Mean shear stress (Iwagaki Eqn)
				 dsm = sed_idx(k)%dmean*100.
				 if (dsm .ge. 0.303 ) then
					 uscm = sqrt(80.9*dsm)
				 elseif (dsm .ge. 0.118) then
					 uscm = sqrt(134.6*dsm**(31./22.))
				 elseif (dsm .ge. 0.0565) then
					 uscm = sqrt(55.0*dsm)
				 elseif (dsm .ge. 0.0065) then
					 uscm = sqrt(8.41*dsm**(11./32.))
				 elseif (dsm .lt. 0.0065) then
					 uscm = sqrt(226.0*dsm)
				 endif
				 dsm = dsm/100.
				 uscm = uscm/100.
				 tscm = uscm**2./((s-1.)*grav*dsm)
	
	!pause'2'
	
	!	   if(Emb_idx(k).ne.0.0) then
     !            if(damflg(k).gt.0)then !added by Qin
     !            u_flax = 0.d0
      !           else
				 u_flax = qr_ave_idx(k)/width_idx(k)/hr_idxa(k)
      !           endif    
	!?ÅEΩÅEΩ}?ÅEΩÅEΩj?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩO?ÅEΩÅEΩÃëe?ÅEΩÅEΩx?ÅEΩÅEΩW?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩp?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩC?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩx?ÅEΩÅEΩÃéZ?ÅEΩÅEΩ?ÅEΩÅEΩ
	!			if(hr_idxa(k).lt.0.001) then
	!		         ust_idx1 = 0.0
	!			else
				! if(damflg(k).gt.0.or. hr_idxa(k).le.0.01d0)then !added by Qin
				if( hr_idxa(k).le.0.01d0)then !added by Qin	
                    ust_idx1 = 0.d0
                    else
                    ust_idx1 = (ns_river*sqrt(grav)* u_flax)/(hr_idxa(k)**(1./6.))
                    endif
                    if(hr_idxa(k).le.0.01d0) then !modified by Qin 2021/6/3
                    Frn = 0.d0
                    else
                    Frn = u_flax/sqrt(9.81*hr_idxa(k))
                    endif
	
	!?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ åÔøΩ?ÅEΩÅEΩz?ÅEΩÅEΩ?ÅEΩÅEΩp?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩC?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩx?ÅEΩÅEΩÃéZ?ÅEΩÅEΩ?ÅEΩÅEΩ
	
	!if(Frn.ge.100.0) then !Fr >= 1.0
	!if(Frn.ge.1.0) then !Fr < 1.0
	
	!pause'3'
	
	if(ave_slope_switch.ne.1)then			 
	!		     if(domain_riv_idx(k).eq.2) cycle
				  kk = down_riv_idx(k)
	!			  slope = (zb_riv_idx(k)-zb_riv_idx(kk))/len_riv_idx(k)
				  slope22 = (zb_riv_idx(k)-zb_riv_idx(kk))/dis_riv_idx(k)
	 !			 if(slope.le.0.0) slope = 0.02
				  
	!else !Fr<1.0
	
	!pause'4'
	
					zb55 = 0.0
					kzb55 = 0
				   do n = 1, 8
					  kkk = up_riv_idx(k,n)
					if(domain_riv_idx(kkk).eq.0.and.kkk==0) cycle  !----modified by Qin
					  zb55 = zb55 + zb_riv_idx(kkk)
					  kzb55 = kzb55 + 1
					enddo
	!			      zb55 = zb55/real(kzb55)
	
	!			    if(domain_riv_idx(kkk).eq.0) then
					if(kzb55.lt.1.or.domain_riv_idx(k).eq.0) then	 !----modified by Qin			  
					  slope33 = tan(zb_riv_slope_idx(k)*3.14159/180.)
	!			      slope = tan(zb_riv_slope_idx(k)*3.14159/180.)
					else
	!			      slope = (zb55 - zb_riv_idx(k))/len_riv_idx(k)
						zb55 = zb55/real(kzb55)	
					   slope33 = (zb55 - zb_riv_idx(k))/dis_riv_idx(k)
					endif
					
	!			    if(slope.le.0.0) slope = 0.002
					
	
	!endif
	!if(Frn.ge.1.0) then
	!slope = slope22
	!else
	!slope = slope33
	!endif
	if(Frn.ge.1.1) then
	slope = slope22
	elseif(Frn.lt.1.1.and.Frn.ge.0.9) then
	slope = (slope22 + slope33)/2.0
	else
	slope = slope33
	endif
!----modified by Qin 2021/6/9
else
slope = tan(zb_riv_slope_idx(k)*3.14159/180.)
			   
endif

!pause'5'
				
			   if(slope.le.0.0) slope = 0.001 
 !              if(damflg(k).gt.0.or. hr_idxa(k).le.0.01d0)then !added by Qin
			   if(hr_idxa(k).le.0.01d0)then !added by Qin
                ust_idx22 = 0.d0
               else     
			   ust_idx22 = sqrt(grav*hr_idxa(k)*slope)
               endif

!			   else
!			     ust_idx22 = 0.0
!			   endif

!                if (hr_idxa(k).lt.1.0) then
!                if (hr_idxa(k).lt.100.0) then
!                  ust_idx33 = ust_idx22 !?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩl?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩÃóp?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩÅEΩ?
!			     else
!                 ust_idx33 = ust_idx1 !?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩl?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩÃóp?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩÅEΩ?
!			     endif

!                if (hr_idxa(k).gt.0.1) then
			   if(dam_switch==1)then !modifed by Qin 2021/6/14
				ust_idx(k) = ust_idx22
				if(damflg(k).gt.0) ust_idx(k) = ust_idx1
			   else
	              ust_idx(k) = ust_idx22
			   endif
!                  if(domain_riv_idx(k).eq.2) ust_idx(k) = ust_idx(k-1)			  
!	              ust_idx(k) = ust_idx33
!	              ust_idx(k) = ust_idx1
!		         else
!	              ust_idx(k) = 0.0
!	             end if
    enddo !end of do k=1, riv_count

    do k = 1, riv_count
!	   if(Emb_idx(k).ne.0.0) then


			  qsb_idx(k) = 0.d0

!------Exchange layer thickness added by Qin 2021/6/11
		   thet = zb_riv_slope_idx(k)
!		   if (thet > 4.) thet = 4.
		   thetr = thet/180.*3.14159	  
		  ! if (thetr.le. 0.) thetr = 0.001 
		    if (thetr.le.  min_slope) thetr =atan(min_slope) 
            Emb_idx(k) = ust_idx(k)**2./((s-1.)*csta/2.*grav*cos(thetr)	 &      !Exchange layer thickness?@e_mov
		   *(tan(S_angle)-tan(thetr)))
!		    if(Emb_idx(k) >= 0.5*hr_idxa(k) ) Emb_idx(k)  = 0.5*hr_idxa(k)
		   if(Emb_idx(k) >= hr_idxa(k)) Emb_idx(k) = hr_idxa(k) !modified by Qin;2021/7/16
!		   if(Emb_idx(k) >= 0.2d0* hr_idxa(k)) Emb_idx(k) = 0.2d0*hr_idxa(k) !modified by Qin;2021/7/16
			if(Emb_idx(k)  <= Emb_min) Emb_idx(k)  = Emb_min
			!bedrock_alfa
			alfa = (zb_riv_idx(k) - zb_roc_idx(k))/Emb_idx(k)
			if(alfa.ge.1.)alfa=1.
			if(alfa <= Emb_min) alfa= 0.d0 


			do m = 1, Np
!	if(dsi.ge.0.0001) then ! di=> 0.1mm?ÅEΩÅEΩ?ÅEΩÅEΩbedload ?ÅEΩÅEΩ?ÅEΩÅEΩ?ÅEΩÅEΩv?ÅEΩÅEΩZ
			  fp = sed_idx(k)%fm(m)
			  if (dsi(m)/dsm .le. 0.4) then
!			     R_tcts = 0.85
			     R_tcts = 0.85 * dsm / dsi(m)
			  else
!			     R_tcts = dsi/dsm * (1.28/(1.28+log10(dsi/dsm)))**2.
			     R_tcts = (1.28/(1.28+log10(dsi(m)/dsm)))**2.
			  endif

			     t_Crit = R_tcts * tscm
			     den = (s-1.)*grav*dsi(m)

			     usci = sqrt(t_Crit * den)

                 !if (hr_idx(k).gt.0.1) then
				 if (hr_idxa(k).gt.0.1) then !modified 2021/6/10
	               ks = 1.+2.*ust_idx(k)**2./(s-1)/grav/dsi(m)
	               lg = 6.0 + 2.5 * log(hr_idxa(k)/dsi(m)/ks)
	               u_se = u_flax/lg
		         else
	               u_se = 0.0 
	             end if

		         t_star = ust_idx(k)**2./den
     	         t_starE = u_se**2./den
		         !isedeq = 2 !1: Ashida-michiue, 2: PMP, 3: Egashira formula
		if(isedeq.eq.1) then
		!ashida-michiue equation
		           if (ust_idx(k).gt.usci) then
		             qbe=(17.*t_starE**1.5)*(1-t_Crit/t_star)*(1-usci/ust_idx(k))
		             sed_idx(k)%qbi(m) = qbe*sqrt((s-1)*grav*dsi(m)**3.)*width_idx(k)*fp*alfa
                   else
		             sed_idx(k)%qbi(m) = 0.0
	               end if  					 
		elseif(isedeq.eq.2)then
		!MPM equation
		           if (ust_idx(k).gt.usci) then
		             qbe = 8.0 * (t_star - t_Crit)**1.5
		             sed_idx(k)%qbi(m) = qbe*sqrt((s-1)*grav*dsi(m)**3.)*width_idx(k)*fp*alfa !?ÅEΩÅEΩP?ÅEΩÅEΩ ÇÔøΩm3/s
                   else
		             sed_idx(k)%qbi(m) = 0.0
	               end if 
                   
         elseif(isedeq.eq.3) then
	   !Egashira formula--added by Qin 2021/6/11
				Dcr = ust_idx(k)**2./(0.05*(s-1.)*grav)
				hs_l = Emb_idx(k)
				!if(hs_l > 0.5*hr_idxa(k)) hs_l = 0.5*hr_idxa(k)
				if(hs_l > hr_idxa(k)) then !modified by Qin 2021/7/16
					hs_l = hr_idxa(k) 
					Emb_idx(k) = hs_l
				endif	
				if (dsi(m) > Dcr) then
					sed_idx(k)%qbi(m) = 0.d0
				!elseif(hs_l.le.Emb_min)then
				!	sed_idx(k)%qbi(m) = 0.0
!				elseif(hs_l > hr_idxa(k)*0.2d0)then  !check this criteria
!					sed_idx(k)%qbi(m) = 0.d0
				elseif(hr_idxa(k) < min_hr)then       !check this criteria
					sed_idx(k)%qbi(m) = 0.d0	
				else
					C_K1 = 1./cos(thetr)/(tan(S_angle)-tan(thetr))
					C_K2 = 1./Csta*2.*(1.-hs_l/hr_idxa(k))**0.5
					!C_K2 = 1./Csta*2.
					C_fd = C_kd*(1.-C_e**2.)*(1.+(s-1.))*(Csta/2.)**(1./3.)
					C_ff = C_kf*(1.-Csta/2.)**(5./3.)*(Csta/2.)**(-2./3.)

					sed_idx(k)%qbi(m) = 4./15.*C_K1**2.*C_K2/sqrt(C_fd+C_ff)*t_star**2.5       &
										*sqrt((s-1.)*grav*dsi(m)**3.)*width_idx(k)*fp*alfa										

					if(sed_idx(k)%qbi(m) <= 0.) sed_idx(k)%qbi(m) = 0.d0
					if(sed_idx(k)%qbi(m).ge.0.01)then
						!write(*,'(2i4,5f10.4)')l,m,C_K1,C_K2,hs_l,hr_idxa(k),thetr
					endif
				endif
	   endif                 

!	else !if(dsi.ge.0.0001) then
!      sed_idx(k)%qbi(m) = 0.d0
!	endif !if(dsi.ge.0.0001) then
       if(damflg(k).gt.0)	sed_idx(k)%qbi(m) = 0.d0  
	  qsb_idx(k) = qsb_idx(k) + sed_idx(k)%qbi(m)
	! return the exchange layer thickness to the initial setted value; Qin 2021/6/18	
!	  Emb_idx(k) = Em	

!pause'6'

if(ust_idx(k).ge.1000000) then
write(*,'(i5,3f12.5,2f15.5)') m, dsi(m), fp, usci, t_Crit/tscm, sed_idx(k)%qbi(m)
endif
	end do ! end of do m = 1, Np

if(ust_idx(k).ge.1000000..or.isnan(ust_idx(k))) then
!write(*,'(a)') '    m         dsi          fp        usci         sed_idxqbi'
write(*,'(a,f)') 'zb_riv_idx(k)=', zb_riv_idx(k)
write(*,'(a,f)') 'length=', len_riv_idx(k)
write(*,'(a,f)') 'distance=', dis_riv_idx(k)
write(*,'(a,f)') 'zb_slope=', zb_riv_slope_idx(k)
write(*,'(a,f)') 'tan(zb_slope)=', tan(zb_riv_slope_idx(k)*3.14/180.)
write(*,'(a,f)') 'slope=', slope
write(*,*) 'hr_idx(k)=', hr_idx(k)
write(*,*) 'hr_idxa(k)=', hr_idxa(k)
write(*,'(a,f)') 'ust=', sqrt(9.81*hr_idxa(k)*slope)
write(*,'(a,i,f)') 'k, ust_idx(k)', k, ust_idx(k)
write(*,'(a,f)') 'ust_idx22=', ust_idx22
write(*,'(a,f)') 'dsm=      ', dsm

write(*,'(a)') '  m           dsi            fp         qbi(m)'
do m = 1, Np
write(*,'(i3,2f14.5,f15.5)') m, sed_idx(k)%dsed(m), sed_idx(k)%fm(m), sed_idx(k)%qbi(m)
enddo
pause
endif

    enddo !end of do k = 1, riv_count
!pause'last of bedload'


      end subroutine bedload

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! added by Harada
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

       subroutine funcd2( sed_lin, hr_idx, hr_idx2,hr_idxa, hr_lin, qr_ave_idx, ust_idx, ust_lin,  qsb_lin, qss_lin, qsw_idx, qsw_lin, t,water_v_lin, dzb_temp_lin,qss_b,sumdzb_lin)
		use globals
		use sediment_mod
		use dam_mod!, only:dam_switch,dam_num,dam_loc, damflg, dam_w_vol,dam_reserv_area ! added by Qin 
		use omp_lib
		implicit none

		type(sed_struct) sed_lin(link_count)

		real(8) hr_idx(riv_count), hr_idx2(riv_count), qr_ave_idx(riv_count), ust_idx(riv_count)
		real(8) hr_idxa(riv_count) 
		real(8) hr_lin(link_count)
		real(8) ust_lin(link_count)
		real(8) qsw_idx(riv_count)
		real(8) qsb_lin(link_count), qss_lin(link_count), qsw_lin(link_count)
		real(8) qss_b(link_count) !added by Qin 2021/5/27
		real(8) qd_sum_idx(riv_count)
		real(8) dzb_temp_lin(link_count),sumdzb_lin(link_count)
		real(8) water_v_lin(link_count)
		real dzb_cap, ffd_total, slope11, h, thetr
		integer t
		
		integer k, kk, kkk, m, n, l, ll,lll,f, h_dam, k_put
		integer kk_div
		real(8) :: u_flax, Frn,thet
		real(8) ::  Emb_min 
		real(8), parameter :: Csta = 0.6
		!Driftwood
		real(8), parameter :: C_Dw = 0.4d0
		real(8), parameter :: C_Dr = 0.8d0
		real(8) :: C_rt, D_depo, D_ero,chan_capa ,wid, expa_rate, divloss, S_D
		
		qd_sum_idx(:) = 0.d0
		Emb_min = Emc
		chan_capa=1.d0
		divloss=1.d0
!write(*,*) 'sediment computation   t=', t
!$omp parallel do private(m)
	do l = 1, link_count
		  dzb_temp_lin(l) = 0.d0
		  sed_lin(l)%stor_sed_sum = 0.d0 !20240424
		do m = 1, Np
		  sed_lin(l)%ffd(m) = 0.d0
		  sed_lin(l)%qdi(m) = 0.d0
		  sed_lin(l)%qdsum(m) = 0.d0
		end do
	end do
			
		call bedload2(sed_lin, hr_idxa, hr_lin, qr_ave_idx, ust_idx, ust_lin, qsb_lin) !--modified by Qin 21/5/20

        call susload2(sed_lin, hr_idxa, hr_lin, ust_idx, ust_lin, qss_lin, qr_ave_idx, t, water_v_lin, qss_b)

!----sum of sediment volume in the dam; Added by Qin
	if(dam_switch == 1) then
		do f = 1, dam_num
			k = dam_loc(f)
			l = link_to_riv(k)
			qss_lin(l) = 0.d0
			do n = 1, 8
				if(up_riv_lin(l,n) == 0) exit 
				lll = up_riv_lin(l, n)
				do m = 1, Np ! Get the GSD of the deposited in the dam reservior  modified by Qin 2021/6/22
!----suspended sediment size <= Dmax_relea_sedi will be discharged directly to the downstream of the dam without depostion in dam reservior  ;modified by Qin 20230420
			   	if (dsi(m).le.Dmax_relea_sedi.and.qr_ave_idx(k).gt.0.d0) then
				sed_lin(l)%qsi(m) = sed_lin(l)%qsi(m)+sed_lin(lll)%qsi(m)! modified 20230603
				dam_sedi_qsi(f,m) = dam_sedi_qsi(f,m) + sed_lin(lll)%qbi(m)*ddt! modified 20230603
				!if(dam_outflow_cs(f).ge.thresh_ss)then
				!dam_sedi_qsi(f,m) = dam_sedi_qsi(f,m) + (dam_outflow_cs(f)-thresh_ss)/dam_outflow_cs(f)*sed_lin(l)%qsi(m)*ddt !add the sediment as deposition sediment that cutted off by thresh concentration
				!sed_lin(l)%qsi(m) = thresh_ss/dam_outflow_cs(f)*sed_lin(l)%qsi(m)
				!endif
				else
!				sed_lin(l)%qsi(m) = 0.d0
				dam_sedi_qsi(f,m) = dam_sedi_qsi(f,m) + sed_lin(lll)%qbi(m)*ddt +sed_lin(lll)%qsi(m)*ddt!modified by Qin 2021/10/19
!				dam_sedi_total_b(f)= dam_sedi_total_b(f) + sed_lin(lll)%qbi(m)*ddt ! modfied 20230603
				dam_sedi_total_s(f)= dam_sedi_total_s(f)+sed_lin(lll)%qsi(m)*ddt
				endif
				dam_sedi_total_b(f)= dam_sedi_total_b(f) + sed_lin(lll)%qbi(m)*ddt ! modfied 20230603				
				enddo
			end do	!enddo n=1,8

			do m=1,np
			if(dsi(m).le.Dmax_relea_sedi.and. qr_ave_idx(k).gt.0.d0)then
			!if(slope_ero_switch.gt.0) sed_lin(l)%qsi(m) =	sed_lin(l)%qsi(m) -slo_s_dsum(k,m)*area_lin(l)
			sed_lin(l)%qsi(m) =	sed_lin(l)%qsi(m) -slo_s_dsum(k,m)*area_lin(l)
			qss_lin(l)=qss_lin(l) + sed_lin(l)%qsi(m) !modified 20230603
			else
		    !if(slope_ero_switch.gt.0) then
			dam_sedi_qsi(f,m)= dam_sedi_qsi(f,m) - slo_s_dsum(k,m)*area_lin(l)*ddt
			dam_sedi_total_s(f)= dam_sedi_total_s(f) -slo_s_dsum(k,m)*area_lin(l)*ddt
			!endif
			endif
			enddo
				
			dam_sedi_total(f) = dam_sedi_total_b(f) + dam_sedi_total_s(f)
			dam_sedi_totalV(f) = dam_sedi_total(f)* dlambda
			dam_outflow_qb(f) = qsb_lin(l)
			dam_outflow_qs(f) = qss_lin(l)
			if(qr_ave_idx(k).gt.0.d0) then !modified 20230603
			dam_outflow_cs(f) = qss_lin(l)/qr_ave_idx(k)
			else
			dam_outflow_cs(f) = 0.d0
			endif
		enddo	

	endif		

!  In case of channel capasity loss due to deposition, pass the bedload to the downstream link  !Machino 20240119
!if(debris_switch == 1)then !20241123
	do l = 1, link_count    !check remaining channel capacity for all river cells k
		n_link_depth(l) = 0
		k = link_idx_k(l)
		do
			if(node_ups(k).ge.1.or.up_riv_idx(k,1).eq.0) then
				if(hr_idx(k)> depth_idx_ini(k)*0.9.and.depth_idx(k)<depth_idx_ini(k))then
					n_link_depth(l) = 1
					!write(*,*) k, l, hr_idx(k), depth_idx_ini(k), depth_idx(k)
					!stop
				end if
				exit
			end if
			if(hr_idx(k)> depth_idx_ini(k)*0.9.and.depth_idx(k)<depth_idx_ini(k))then
				n_link_depth(l) = 1
				!write(*,*) k, l, hr_idx(k), depth_idx_ini(k), depth_idx(k)
				!stop
			end if
			k = up_riv_idx(k,1)
		end do
!		if(depth_idx(k) < min_hr*2. )then
!			n_link_depth(l) = 1
!			write(*,*) k, l, depth_idx(k),min_hr
!		end if
	end do

	do l = 1,link_count
		k = link_idx_k(l)
		if(n_link_depth(l)==1) then   !20240229
		!if(zb_riv_idx(k)-zb_riv0_idx(k) + depth_idx(k)*0.1 >  depth_idx(k) ) then !check the condition 'depth_idx(k)*0.2'
		!if(depth_idx(k) < min_hr*2. )then !check this condition
				do m = 1, Np
					do n = 1, 8
						if(up_riv_lin(l,n) == 0) exit 	
						lll = up_riv_lin(l, n)			
						sed_lin(l)%qdsum(m) = sed_lin(l)%qdsum(m) - sed_lin(lll)%qbi(m)    !qdsum and ffd; erosion positive
					end do    									!calculate input sediment from two unit channels at first 
					sed_lin(l)%qbi(m) = -sed_lin(l)%qdsum(m)   !modify the output sediment so that in = out
					sed_lin(l)%qdsum(m) = 0.                  !initialization for the next loop
				end do

			if(depth_idx(k) < 0.)then
				!write(*,*) 'depth(k)<0,  k ,l, up_iln ', k, l, up_riv_lin(l, 1), up_riv_lin(l, 2)
			end if
		endif	
	end do
!end if
!!$omp parallel do private(m,n,lll,thet,k)
	do l =1, link_count
		thet = zb_riv_slope0_lin(l)
	!added for diversion 20250329	
		divloss =1.
		if(div_switch==2)then
			k = link_idx_k(l)
			do f = 1, div_id_max
			if(div_org_idx(f)==k) divloss=1.-div_rate(f)
			enddo
		endif
			if(link_0th_order(l)== 0 .and. thet.gt.max_slope)then !modified by Qin 2022/12/06	
				do m = 1, Np
					sed_lin(l)%qdsum(m) = 0.d0
				enddo
			else
			if(link_0th_order(l) == 0)then
				do m = 1, Np
					sed_lin(l)%qdsum(m) = sed_lin(l)%qdsum(m) + sed_lin(l)%qbi(m) !m3/s	 !check this condition  !?????????!  
		!			sed_lin(l)%qdsum(m) = sed_lin(l)%qdsum(m)  !m3/s                         !0th order channel: unlimit supply (no beddeform at upstream end)
				end do
			else
				do m = 1, Np
					do n = 1, 8
		!				if(up_riv_lin(l,n)==0) cycle !revised by Qin 2021/5/27
						if(up_riv_lin(l,n) == 0) exit 	
						lll = up_riv_lin(l, n)
						sed_lin(l)%qdsum(m) = sed_lin(l)%qdsum(m) - sed_lin(lll)%qbi(m)*divloss !m3/s   !?????????!!!  ffd ??-???????ÅEΩÅEΩ?ÅEΩÅEΩA+?????N?H modified for diversion 20250329
					end do
					sed_lin(l)%qdsum(m) = sed_lin(l)%qdsum(m) + sed_lin(l)%qbi(m)
				end do
			end if
			endif 
	
	end do
!!$omp end parallel do

!pause'after sediment equation'
!----?|??????????v?Z only for bedload
!$omp parallel do private(m)
   do l = 1, link_count
	 do m = 1, Np	
	  sed_lin(l)%ffd(m) = sed_lin(l)%qdsum(m)/area_lin(l) !?P???[m/s]
	 enddo
   enddo

!$omp single
if(debris_switch == 1)then !20230924
!!$omp parallel do private(k,m)
  do l = 1, link_count
   vo_total_river(l) = 0.d0
   if(n_link_depth(l)==0)then		
  	 k = link_idx_k(l)
	!add sediment from debris flow
	if(area_lin(l) > 100. .and. Emb_lin(l)>0.01) then
	 if( vo_total_l(l)/area_lin(l)*dlambda > Emb_lin(l)*0.01 ) then  !avoid huge deposition in small link
		vo_total_l(l) = vo_total_l(l) - Emb_lin(l)*0.01*area_lin(l)/dlambda
		vo_total_river(l) = Emb_lin(l)*0.01*area_lin(l)/dlambda
		if ( vo_total_l(l) < 0.d0 ) then
			vo_total_l(l) = vo_total_l(l) + Emb_lin(l)*0.01*area_lin(l)/dlambda
			vo_total_river(l) = vo_total_l(l)
			vo_total_l(l) = 0.d0
		end if
	 else
		vo_total_river(l) = vo_total_l(l)
		vo_total_l(l) = 0.d0
	 end if
	else
	 if( vo_total_l(l)/area_lin(l)*dlambda > Emb_lin(l)*0.0005 ) then  !avoid huge deposition in small link
		vo_total_l(l) = vo_total_l(l) - Emb_lin(l)*0.0005*area_lin(l)/dlambda
		vo_total_river(l) = Emb_lin(l)*0.0005*area_lin(l)/dlambda
		if ( vo_total_l(l) < 0.d0 ) then
			vo_total_l(l) = vo_total_l(l) + Emb_lin(l)*0.0005*area_lin(l)/dlambda
			vo_total_river(l) = vo_total_l(l)
			vo_total_l(l) = 0.d0
		end if
	 else
		vo_total_river(l) = vo_total_l(l)
		vo_total_l(l) = 0.d0
	 end if
	endif 
   end if   !enf if for (n_link_depth(l)==0)
	if(vo_total_l(l)<0.)then
	  write(*,*) 'vo_total_l < 0', l, vo_total_l(l)
	  pause
	end if	
!-----
!-----added 20240424
	debri_sup_sum(l) = debri_sup_sum(l) + vo_total_river(l)
	do m = 1, NP
	debri_sup_sum_di(l,m) =	debri_sup_sum_di(l,m) + vo_total_river(l)*fmslo(k,m)
	if(fmslo(k,m)>1.d0) then
	    write(*,*) 'fmslo >1 ', k, m, dsi(m), fmslo(k,m) !modified 20240509
		stop
	endif
	enddo
  end do  
endif
!$omp end single 

!-----added 20260220 sedput gradually
!$omp single
if(j_sedput == 1)then
  do l = 1, link_count
   k = link_idx_k(l)
   put_depth(l) = 0.d0
	if(area_lin(l) > 100. .and. Emb_lin(l)>0.01) then
	 if( put_depth_remain(l) > Emb_lin(l)*0.01 ) then  !avoid huge deposition in small link
	    if(sumdzb_lin(l) < 40.)then !check this condition
	 		put_depth(l) = Emb_lin(l)*0.01
			put_depth_remain(l) = put_depth_remain(l) - put_depth(l)
		end if
	 else
		put_depth(l) = put_depth_remain(l)
		put_depth_remain(l) = 0.d0
	end if
	else   !for very small link
	 if( put_depth_remain(l) > Emb_lin(l)*0.0005 ) then
	 	put_depth(l) = Emb_lin(l)*0.0005
		put_depth_remain(l) = put_depth_remain(l) - put_depth(l)
	 else
		put_depth(l) = put_depth_remain(l)
		put_depth_remain(l) = 0.d0
	 end if
	end if

	if(put_depth_remain(l)<0.) put_depth_remain(l) =0.
  end do
  
  k_put = 0
  do l = 1, link_count
	if(put_depth(l) > 0. .or. put_depth_remain(l) > 0.)then
		k_put = k_put + 1
	end if
  end do
  if(k_put == 0) j_sedput = 0

end if
!$omp end single

!	if(slo_sedi_cal_switch==0)then    ! added 20240806
!	    do k = 1, slo_count
!			do m = 1, Np
!				fmslo(k,m) = 0.d0
!			end do
!		end do
!	end if

!!$omp parallel do private(k,m)
!turnoff the parallel computation 20240424
!$omp single
  do l = 1, link_count
	k = link_idx_k(l)
	do m = 1, Np
!		if(link_0th_order(l) == 0) then
!		sed_lin(l)%ffd(m) = sed_lin(l)%ffd(m) !!0th order channel: unlimit supply of suspended sediment;2021/7/16
!		else
		sed_lin(l)%ffd(m) = sed_lin(l)%ffd(m) + (sed_lin(l)%Esi(m) - sed_lin(l)%Dsi(m)) !?P???[m/s]
!		endif	
!			if(link_0th_order(l) == 0) sed_lin(l)%ffd(m) = 0.d0 !!0th order channel: unlimit supply (no beddeform at upstream end);2021/6/14
		!sed_lin(l)%dzbpr(m) = - ddt * sed_lin(l)%ffd(m) * dlambda !?P???[m]	
		sed_lin(l)%dzbpr(m) = - ddt * sed_lin(l)%ffd(m) * dlambda + fmslo(k,m)*vo_total_river(l)/area_lin(l)*dlambda   !debris supply was added
!-------Sediment put: put_depth is bulk bed thickness including voids [m].
		if(j_sedput == 1 .and. l_put(l) == 1 )then
			sed_lin(l)%dzbpr(m) = sed_lin(l)%dzbpr(m) + fm_sedput(l,m)*put_depth(l)
		end if
!---added 20211126		
		sed_lin(l)%Esisum(m) = sed_lin(l)%Esisum(m) + sed_lin(l)%Esi(m)*ddt
		sed_lin(l)%Dsisum(m) = sed_lin(l)%Dsisum(m) + sed_lin(l)%Dsi(m)*ddt			
!		if(vo_total_l(l)>0.d0)then
!			write(*,'(a,2i4,3f15.9)')'k,m,dzbpr,fmslo,vo_total_l=',k,m,sed_lin(l)%dzbpr(m), fmslo(k,m), vo_total_l(l)
!		end if	
		dzb_temp_lin(l) = dzb_temp_lin(l) + sed_lin(l)%dzbpr(m) ![m]
!added 20240424		
		sed_lin(l)%stor_sed_sum_i(m) = sed_lin(l)%stor_sed_sum_i(m) + sed_lin(l)%dzbpr(m) *area_lin(l)/dlambda
		sed_lin(l)%stor_sed_sum = sed_lin(l)%stor_sed_sum + sed_lin(l)%stor_sed_sum_i(m)	
	 enddo
 enddo

!$omp end single

   ! to prevent over depostion in unit channel 20231103; added for river width adjustment 20240419	 20240422
   ! adjust the overflow state of each river cell to make the interaction states of flow between slope and river are same with their unit channel ;modified 20240424
   if (cut_overdepo_switch == 0 .and. riv_wid_expan_switch > 0) then 
!!$omp parallel do private(k,l,kk, h, wid, expa_rate,thetr) 
	 do k = 1, riv_count
	 	l = link_to_riv(k)
		!if(l>40) cycle ! for gofukuya 20240523
		thetr = zb_riv_slope_lin(l) /180.*3.14159 !20240508
		!if(dzb_temp_lin(l).ge.0.d0)then
		if(width_lin(l)/width_lin_0(l) .ge. max_expan_rate) cycle		
			!k = link_idx_k(l)
			kk = down_riv_idx(k)
		expa_rate = 1.	
		if((hr_idx(k)+zb_riv_idx(k) + dzb_temp_lin(l)).le.wl4wid_expan(k) .and. (depth_idx(k)+height_idx(k) - dzb_temp_lin(l)).gt.0.d0) cycle
			!h = 0.05 * (s-1.d0) * sed_lin(l)%Bed_D60 / ws_slope_lin(l)
			!wid = ns_river * qr_ave_idx(k) / h**(5./3.) / sqrt(ws_slope_lin(l)) 
			!20240508
			h = 0.05 * (s-1.d0) * sed_lin(l)%Bed_D60 / tan(thetr)
			wid = ns_river * qr_ave_idx(k) / h**(5./3.) / sqrt(tan(thetr)) 
			if(wid .ge.width_idx(k))then
			! expa_rate =max( wid/ width_idx(k), 0.1)
			 expa_rate = wid / width_idx(k)
			 wid_expa_rate(l) = max(width_lin(l)*expa_rate/ width_lin_0(l), wid_expa_rate(l))
			 if (wid_expa_rate(l) .ge.max_expan_rate) wid_expa_rate(l)  =  max_expan_rate
			endif					
	 enddo
!!$omp parallel do private(m)
	 do l = 1, link_count
		dzb_temp_lin(l)=0.d0 
		width_lin(l) =  wid_expa_rate(l) *width_lin_0(l)
		do m = 1,Np 
		sed_lin(l)%dzbpr(m)=sed_lin(l)%dzbpr(m)/wid_expa_rate(l)
		dzb_temp_lin(l) = dzb_temp_lin(l) + sed_lin(l)%dzbpr(m) ![m]
		enddo  
		area_lin(l) = width_lin(l) * link_len(l)
	 enddo
   elseif(cut_overdepo_switch > 0 .and. riv_wid_expan_switch> 0)then
   ! width expansion
	!!$omp parallel do private(k,l,kk, h, wid, expa_rate,thetr) 
	 do k = 1, riv_count
	 	l = link_to_riv(k)
		!if(dzb_temp_lin(l).ge.0.d0)then
		!if(l>40) cycle ! for gofukuya 20240523
		if(width_lin(l)/width_lin_0(l) .ge. max_expan_rate) cycle		
			!k = link_idx_k(l)
			kk = down_riv_idx(k)
		expa_rate = 1.	
		if((hr_idx(k)+zb_riv_idx(k) + dzb_temp_lin(l)).le.wl4wid_expan(k) .and. (depth_idx(k)+height_idx(k) - dzb_temp_lin(l)).gt.0.d0) cycle
			!h = 0.05 * (s-1.d0) * sed_lin(l)%Bed_D60 / ws_slope_lin(l)
			!wid = ns_river * qr_ave_idx(k) / h**(5./3.) / sqrt(ws_slope_lin(l)) 
			!20240508
			h = 0.05 * (s-1.d0) * sed_lin(l)%Bed_D60 / tan(thetr)
			wid = ns_river * qr_ave_idx(k) / h**(5./3.) / sqrt(tan(thetr)) 
			if(wid .ge.width_idx(k))then
			 !expa_rate =max( wid/ width_idx(k), 0.1)
			 expa_rate = wid / width_idx(k)
			 wid_expa_rate(l) = max(width_lin(l)*expa_rate/ width_lin_0(l), wid_expa_rate(l))
			 if (wid_expa_rate(l) .ge.max_expan_rate) wid_expa_rate(l)  =  max_expan_rate
			endif					
	 enddo
!!$omp parallel do private(m)
	 do l = 1, link_count
		dzb_temp_lin(l)=0.d0 
		width_lin(l) =  wid_expa_rate(l) *width_lin_0(l)
		do m = 1,Np 
		sed_lin(l)%dzbpr(m)=sed_lin(l)%dzbpr(m)/wid_expa_rate(l)
		dzb_temp_lin(l) = dzb_temp_lin(l) + sed_lin(l)%dzbpr(m) ![m]
		enddo  
		area_lin(l) = width_lin(l) * link_len(l)
	 enddo
	! cutting over-deposition	
	!!$omp parallel do private(l,k,m,chan_capa)
		do l = 1, link_count
				k = link_idx_k(l)
		!do k = 1, riv_count
	 		!l = link_to_riv(k)				
			if((depth_idx(k)+height_idx(k)).le. 0.)then
				dzb_temp_lin(l)=0.d0 
				do m = 1,Np 
				if(sed_lin(l)%dzbpr(m).ge.0.d0)then
					overdepo_sedi_di(l,m) = overdepo_sedi_di(l,m) + sed_lin(l)%dzbpr(m)/dlambda
					sed_lin(l)%stor_sed_sum_i(m) = sed_lin(l)%stor_sed_sum_i(m) -sed_lin(l)%dzbpr(m)/dlambda*area_lin(l)
					sed_lin(l)%stor_sed_sum = sed_lin(l)%stor_sed_sum - sed_lin(l)%dzbpr(m)/dlambda*area_lin(l)
					sed_lin(l)%dzbpr(m)= 0.d0
					sed_lin(l)%Dsisum(m) = sed_lin(l)%Dsisum(m) - sed_lin(l)%Dsi(m)*ddt	
				endif
					dzb_temp_lin(l) = dzb_temp_lin(l) + sed_lin(l)%dzbpr(m) ![m]
				enddo	
			else
				chan_capa=depth_idx(k)+height_idx(k)-dzb_temp_lin(l)
				if(chan_capa .lt. 0.d0)then		
					chan_capa = chan_capa/dzb_temp_lin(l)
					if(dzb_temp_lin(l)==0.d0) chan_capa=0.d0
					dzb_temp_lin(l)=0.d0 
					do m = 1,Np 
						if(sed_lin(l)%dzbpr(m).ge.0.d0)then	
							overdepo_sedi_di(l,m) = overdepo_sedi_di(l,m) - (sed_lin(l)%dzbpr(m)*chan_capa)/dlambda
							sed_lin(l)%stor_sed_sum_i(m) = sed_lin(l)%stor_sed_sum_i(m) + (sed_lin(l)%dzbpr(m)*chan_capa)/dlambda * area_lin(l)
							sed_lin(l)%stor_sed_sum = sed_lin(l)%stor_sed_sum  + (sed_lin(l)%dzbpr(m)*chan_capa)/dlambda * area_lin(l)							
							sed_lin(l)%dzbpr(m)=sed_lin(l)%dzbpr(m)+sed_lin(l)%dzbpr(m)*chan_capa
							sed_lin(l)%Dsisum(m) = sed_lin(l)%Dsisum(m) + chan_capa*sed_lin(l)%Dsi(m)*ddt	
						endif	
							dzb_temp_lin(l) = dzb_temp_lin(l) + sed_lin(l)%dzbpr(m) ![m]		
					enddo  
				endif
			endif	
			!endif
		enddo
   elseif(cut_overdepo_switch > 0 .and. riv_wid_expan_switch== 0)then
!!$omp parallel do private(k,m,chan_capa)  
		do l = 1, link_count
				k = link_idx_k(l)
		!do k = 1, riv_count
	 		!l = link_to_riv(k)				
			if((depth_idx(k)+height_idx(k)) .le. 0.)then
				dzb_temp_lin(l)=0.d0 
				do m = 1,Np 
				if(sed_lin(l)%dzbpr(m).ge.0.d0)then
					overdepo_sedi_di(l,m) = overdepo_sedi_di(l,m) + sed_lin(l)%dzbpr(m)/dlambda
					sed_lin(l)%stor_sed_sum_i(m) = sed_lin(l)%stor_sed_sum_i(m) -sed_lin(l)%dzbpr(m)/dlambda*area_lin(l)
					sed_lin(l)%stor_sed_sum = sed_lin(l)%stor_sed_sum - sed_lin(l)%dzbpr(m)/dlambda*area_lin(l)
					sed_lin(l)%dzbpr(m)= 0.d0
					sed_lin(l)%Dsisum(m) = sed_lin(l)%Dsisum(m) - sed_lin(l)%Dsi(m)*ddt	
				endif
					dzb_temp_lin(l) = dzb_temp_lin(l) + sed_lin(l)%dzbpr(m) ![m]
				enddo	
			else
				chan_capa=depth_idx(k)+height_idx(k)-dzb_temp_lin(l)
				if(chan_capa .le. 0.)then		
					chan_capa = chan_capa/dzb_temp_lin(l)
					if(dzb_temp_lin(l)==0.d0) chan_capa=0.d0
					dzb_temp_lin(l)=0.d0 
					do m = 1,Np 
						if(sed_lin(l)%dzbpr(m).ge.0.d0)then	
							overdepo_sedi_di(l,m) = overdepo_sedi_di(l,m) - (sed_lin(l)%dzbpr(m)*chan_capa)/dlambda
							sed_lin(l)%stor_sed_sum_i(m) = sed_lin(l)%stor_sed_sum_i(m) + (sed_lin(l)%dzbpr(m)*chan_capa)/dlambda * area_lin(l)
							sed_lin(l)%stor_sed_sum = sed_lin(l)%stor_sed_sum  + (sed_lin(l)%dzbpr(m)*chan_capa)/dlambda * area_lin(l)							
							sed_lin(l)%dzbpr(m)=sed_lin(l)%dzbpr(m)+sed_lin(l)%dzbpr(m)*chan_capa
							sed_lin(l)%Dsisum(m) = sed_lin(l)%Dsisum(m) + chan_capa*sed_lin(l)%Dsi(m)*ddt	
						endif	
							dzb_temp_lin(l) = dzb_temp_lin(l) + sed_lin(l)%dzbpr(m) ![m]		
					enddo  
				endif
			endif	
			!endif
		enddo
   endif

!  Downstream end unit channel is the outlet boundary.
!  Positive bed change there is exported out of the basin instead of raising the outlet bed.
  do l = 1, link_count
    k = link_idx_k(l)
    if(domain_riv_idx(k).eq.2)then
        dzb_temp_lin(l)=0.d0
        do m = 1, Np
            if(overdepo_sedi_di(l,m).gt.0.d0)then
                qsb_total = qsb_total + overdepo_sedi_di(l,m)*area_lin(l)
                overdepo_sedi_di(l,m)=0.d0
            endif
            if(sed_lin(l)%dzbpr(m).gt.0.d0)then
                qsb_total = qsb_total + sed_lin(l)%dzbpr(m)*area_lin(l)/dlambda
                sed_lin(l)%stor_sed_sum_i(m) = sed_lin(l)%stor_sed_sum_i(m) - sed_lin(l)%dzbpr(m)*area_lin(l)/dlambda
                sed_lin(l)%stor_sed_sum = sed_lin(l)%stor_sed_sum - sed_lin(l)%dzbpr(m)*area_lin(l)/dlambda
                sed_lin(l)%dzbpr(m)=0.d0
            endif
            dzb_temp_lin(l)=dzb_temp_lin(l)+sed_lin(l)%dzbpr(m)
        enddo
    endif
  enddo
  !j_sedput = 0 

if(j_drf == 1)then
!!$omp parallel do private(k,C_rt,D_depo,D_ero,S_D, kkk1)  
	do l = 1, link_count
		k = link_idx_k(l)
		if(hr_idxa(k) < min_hr)then
			cw(l) = 0.d0
		elseif(domain_riv_idx(k).eq.2) then
			!cw(l) = 0.d0 !no DW in downstream end
		else
		
		!Driftwood erosion and deposition
			if(dzb_temp_lin(l) > 0.d0) then   !Sediment Deposition
				if(hr_idxa(k)> 2.d0*C_Dw)then
					C_rt=0.
                else
                	C_rt=-0.5d0*hr_idxa(k)/C_Dw+1.d0
				end if
				cw(l) = cw(l) + (-1.*qwsum(l) - C_wood_kd*C_rt*cw(l)*(1.-Csta)*dzb_temp_lin(l)*area_lin(l)/ddt) * ddt/water_v_lin(l)   !C_wood_kd added 20251023
				vw(l) = vw(l) + ddt*C_wood_kd*C_rt*cw(l)*(1.-Csta)*dzb_temp_lin(l)/ddt                                                 !C_wood_kd added 20251023
				if(cw(l).lt.0.0) then 
					cw(l) = cw(l) - (-1.*qwsum(l) - C_wood_kd*C_rt*cw(l)*(1.-Csta)*dzb_temp_lin(l)*area_lin(l)/ddt) * ddt/water_v_lin(l)  !‰∏ÄÂõûÊàª„Åó„Å¶  !C_wood_kd added 20251023
					D_depo = cw(l)*hr_idxa(k) + qwsum(l)*ddt/area_lin(l)
					if(D_depo < 0.d0) D_depo = 0.d0
					cw(l) = 0.d0
					vw(l) = vw(l) + D_depo
				end if
			else 							  !Sediment Erosion
				if(hr_idxa(k)>= 2.*C_Dw)then
					C_rt=1.d0
				else
					C_rt=0.5d0*hr_idxa(k)/C_Dw
				endif
				S_D = vw(l)/C_Dr       !20260223 added S_D (=S/D) to avoid huge erosion 
				if(S_D > 1.) S_D = 1.
				if(cw(l) < alpha_ss1 )then    !20260223 added this if statement
					cw(l) = cw(l) + (-1.*qwsum(l) - C_rt*S_D*(1.-Csta)*dzb_temp_lin(l)*area_lin(l)/ddt) * ddt/water_v_lin(l)
					vw(l) = vw(l) + ddt*C_rt*S_D*(1.-Csta)*dzb_temp_lin(l)/ddt
				else      ! cw(l) > alpha_ss1     no wood entrainment, transport only
					cw(l) = cw(l) + (-1.*qwsum(l)) * ddt/water_v_lin(l)
				end if
				if(vw(l).lt.0.0) then
					vw(l) = vw(l) - ddt*C_rt*S_D*(1.-Csta)*dzb_temp_lin(l)/ddt
					D_ero = -vw(l)
					cw(l) = cw(l) + (-1.*qwsum(l) - D_ero*area_lin(l)/ddt) * ddt/water_v_lin(l)
					vw(l) = 0.d0
					if(cw(l).lt.0.d0) cw(l) = 0.d0
				end if
			end if
		end if
		qw(l) = cw(l) * qr_ave_idx(k)
		qwsum_total(l) = qwsum_total(l) + qw(l)*ddt
	enddo
end if

!pause'before call washload'
!          call washload(sed_idx, hr_idxa,  ust_idx, qsw_idx, qr_ave_idx, dzb_temp, width_idx, len_riv_idx)
		   !call washload2(sed_lin, hr_idxa,  ust_idx, ust_lin, qsw_idx, qsw_lin, qr_ave_idx, dzb_temp_lin,water_v_lin )

 !  do l = 1, link_count
	!----washload??????Z
	!	dzb_temp_lin(l) = 0.d0 
!		if (dam_switch==1)then !added by Qin
!			k = link_idx_k(l)
!			if(damflg(k).gt.0)then
!			dzb_temp_lin(l)= 0.d0
!			else	
!			do m = 1, Np
!			sed_lin(l)%ffd(m) = sed_lin(l)%ffd(m) + (sed_lin(l)%Ewi(m) - sed_lin(l)%Dwi(m)) !?P???[m/s]
!			if(link_0th_order(l) == 0) sed_lin(l)%ffd(m) = 0.         !0?????????ÅEΩÅEΩp??????????B???V???Awashload????????s????????????????u??ffd??0??
!			sed_lin(l)%dzbpr(m) = - ddt * sed_lin(l)%ffd(m) * dlambda !?P???[m]
!			dzb_temp_lin(l) = dzb_temp_lin(l) + sed_lin(l)%dzbpr(m) ![m]
!			enddo	
!			endif	
!		else	
!		do m = 1, Np
!			sed_lin(l)%ffd(m) = sed_lin(l)%ffd(m) + (sed_lin(l)%Ewi(m) - sed_lin(l)%Dwi(m)) !?P???[m/s]
!			if(link_0th_order(l) == 0) sed_lin(l)%ffd(m) = 0.         !0?????????ÅEΩÅEΩp??????????B???V???Awashload????????s????????????????u??ffd??0??
!			sed_lin(l)%dzbpr(m) = - ddt * sed_lin(l)%ffd(m) * dlambda !?P???[m]
!			dzb_temp_lin(l) = dzb_temp_lin(l) + sed_lin(l)%dzbpr(m) ![m]
!		enddo
!		endif
	
!	   enddo

!----sum of sediment transport at downstream boundary
    do l = 1, link_count
		k = link_idx_k(l)
!-----for check the qss budget 2021/5/28
		!do n = 1, 8
!			if(up_riv_lin(l,n) == 0) cycle !revised by Qin 2021/5/27
		!	if(up_riv_lin(l,n) == 0) exit 
		!	lll = up_riv_lin(l, n)
		!	dam_sedi_total_s2(k) = dam_sedi_total_s2(k) + qss_b(lll)*ddt !---check
		!enddo
		!dam_sedi_total_s2(k) = dam_sedi_total_s2(k) - qss_lin(l)*ddt!---check
		if( domain_riv_idx(k) == 2 )then
		!if(k==181)then   !chcek
			qsb_total =  qsb_total+ qsb_lin(l)*ddt
			qss_total = qss_total + qss_lin(l)*ddt
			qsw_total = qsw_total + qsw_lin(l)*ddt
			if(j_drf==1)qwood_total = qwood_total + qw(l)*ddt
		end if
    end do

!pause'after call washload'
!pause'endof funcd'

	 end subroutine funcd2

! >>>>>>>>>>>>>>>>>>>>>>> CALCULATING wash laod
!      subroutine washload(sed_idx, hr_idxa, ust_idx, qsw_idx, qr_ave_idx, dzb_temp, width_idx, len_riv_idx)
	 subroutine washload2(sed_lin, hr_idxa,  ust_idx, ust_lin, qsw_idx, qsw_lin, qr_ave_idx, dzb_temp_lin, water_v_lin )
		use globals
		use sediment_mod
		use dam_mod!, only: damflg !added by Qin
		implicit none
		
		type(sed_struct) sed_lin(link_count)
		real(8) hr_idx(riv_count), hr_idx2(riv_count), ust_idx(riv_count)
		real(8) hr_idxa(riv_count)
		real(8) qsw_idx(riv_count), qr_ave_idx(riv_count)
		real(8) ss_idx(riv_count)
	!	real(8) ss_lin(link_count)
		real(8) qsw_sum_idx(riv_count)
		real(8) dzb_temp(riv_count)
		real(8) ust_lin(link_count), qsw_lin(link_count)
		real(8) dzb_temp_lin(link_count), water_v_lin(link_count)
		integer k, kk, kkk, n, mm, mmax, m, l, lll
		real(8) :: ss11, ss1, ss2, ss3, ss22
		real(8) :: ss4, ss5, DDD
		real(8) :: u_flax, Frn
		integer mmm

	   DDD = 0.01
!pause'in washload'

	   do l = 1, link_count
		qsw_lin(l) = 0.d0
		ss_lin(l) = 0.d0
	   enddo

   do l = 1, link_count
	k = link_idx_k(l)
	   do m = 1, Np
		if(qr_ave_idx(k).le.0.0) then
		 sed_lin(l)%swi(m) = 0.0
		else
		 sed_lin(l)%swi(m) = sed_lin(l)%qswi(m) / qr_ave_idx(k)
		endif
	   enddo

	   do m = 1, Np
		 ss_lin(l) = ss_lin(l) + sed_lin(l)%swi(m)
	   enddo

		if(ss_lin(l).gt.0.1) then
		 do m = 1, Np
		   sed_lin(l)%swi(m) = sed_lin(l)%swi(m) / ss_lin(l)* 0.1
		 enddo
		endif
   enddo

   do l = 1, link_count
		do m = 1, Np
		  sed_lin(l)%qswisum(m) = 0.d0
		enddo
   enddo

!pause'before call washsedi'
   call washsedicom2(sed_lin, dzb_temp_lin, ust_lin)
!pause'after call washsedi'
    !?????S????????????
   do l = 1, link_count
	if(link_0th_order(l) == 0)then
		do m = 1, Np
			!sed_lin(l)%qswisum(m) = sed_lin(l)%qswisum(m) + sed_lin(l)%qsi(m)    !?????????
			sed_lin(l)%qswisum(m) = sed_lin(l)%qswisum(m)    !????????ÅEΩÅEΩ??
		end do
	else
			
		do m = 1, Np
			do n = 1, 8
!				if(up_riv_lin(l,n) == 0) cycle !revised by Qin 2021/5/27
			if(up_riv_lin(l,n) == 0) exit 
				lll = up_riv_lin(l, n)
				sed_lin(l)%qswisum(m) = sed_lin(l)%qswisum(m) - sed_lin(lll)%qswi(m) !m3/s   !?????????!!!  ffd ??-???????ÅEΩÅEΩ?ÅEΩÅEΩA+?????N?H
			end do
			sed_lin(l)%qswisum(m) = sed_lin(l)%qswisum(m) + sed_lin(l)%qswi(m)
		end do
	end if
   end do

   do l = 1, link_count
	k = link_idx_k(l)  
		do m = 1, Np
		 ss11 = sed_lin(l)%swi(m)
!		  ss1 = ss11 * hr_idxa2(k)
!		  ss2 = -1.* sed_idx(k)%qswisum(m)
		  if(domain_riv_idx(k).eq.2) sed_lin(l)%qswisum(m) = 0.0
		 ss22 = ss2/width_idx(k)/dis_riv_idx(k)
		 ss3 = sed_lin(l)%Ewi(m) - sed_lin(l)%Dwi(m)
		
		if(water_v_lin(l).le.0.1) then
			sed_lin(l)%swi(m) = 0.0
		else	
			sed_lin(l)%swi(m) = sed_lin(l)%swi(m) + (-1.* sed_lin(l)%qswisum(m) + sed_lin(l)%Ewi(m)*area_lin(l) - sed_lin(l)%Dwi(m)*area_lin(l)) * ddt/water_v_lin(l)
		 if(sed_lin(l)%swi(m).le.0.0) then
		  !if (damflg(k).gt.0)then !added by Qin	
			!sed_lin(l)%Dwi(m) = 0.0
		  !else	
		  sed_lin(l)%Dwi(m) = -1.* sed_lin(l)%qswisum(m)/area_lin(l) + sed_lin(l)%Ewi(m)
		  !endif
		  sed_lin(l)%swi(m) = 0.0
		 endif
		 !endif
		endif

		enddo
!write(*,'(a)') 'k, a2(k), hr_idxa(k)'
!pause
		
!		if(ust_idx(k).ge.0.1) then
		 do m = 1, Np
		  sed_lin(l)%qswi(m) = sed_lin(l)%swi(m) * qr_ave_idx(k) !m3/s
		  qsw_lin(l) = qsw_lin(l) + sed_lin(l)%qswi(m) !m3/s
		 enddo			 

   enddo

!pause'end of washload'

end subroutine washload2

! >>>>>>>>>>>>>>>>>>>>>>> CALCULATING washload sediment
	subroutine washsedicom2(sed_lin, dzb_temp_lin, ust_lin)
	use globals
	use sediment_mod
	use dam_mod!, only: damflg !added by Qin
	implicit none
		
	type(sed_struct) sed_lin(link_count)
	real(8) dzb_temp_lin(link_count), ust_lin(link_count)
	real(8) ss
	real(8) dsm, fp, ramdaw
	integer k, m, iwcheck, l
	real(8) Caew, dzb


	ramdaw = 0.4

!pause'in washsedi'

		do l = 1, link_count
! For non-uniform sized sediment
!         	if(ust_idx(k).gt.0.0) then
			k = link_idx_k(l) 
			do m = 1, Np
			   dzb = sed_lin(l)%dzbpr(m)
			 if (dsi(m).le.0.0001) then ! di=< 0.1mm??wash load?????????

			   fp = sed_lin(l)%fm(m)
			   ss = sed_lin(l)%swi(m)
			 if(dzb_temp_lin(l).le.0.d0) then	
			 ! if(dzb_temp_lin(l).le.0.0.and.damflg(k).eq.0) then!modifid by Qin
!   		       if(dzb_temp(k).le.0.0.and.ust_idx(k).ge.0.12) then
			   Caew = -(1.-ramdaw)*fp*dzb_temp_lin(l)/ddt
			   sed_lin(l)%Ewi(m) = Caew
			  else
			   sed_lin(l)%Ewi(m) = 0.d0
			  endif
			 sed_lin(l)%Dwi(m) = ss * w0(m)

			 else !middle of if (dsi.le.0.0001) then
			   sed_lin(l)%Ewi(m) = 0.d0
			   sed_lin(l)%Dwi(m) = 0.d0
			 end if !end of if (dsi.le.0.0001) then

			end do

	   end do ! end of do l = 1, link_count

!pause'after washsedi'

	 end subroutine washsedicom2
	 
! >>>>>>>>>>>>>>>>>>>>>>> CALCULATING suspended laod
!      subroutine susload(sed_idx, hr_idxa, ust_idx, qss_idx, qr_ave_idx, width_idx, len_riv_idx)
	subroutine susload2(sed_lin, hr_idxa, hr_lin, ust_idx, ust_lin, qss_lin, qr_ave_idx, t, water_v_lin,qss_b)
	use globals
	use sediment_mod
	use dam_mod!, only:damflg
	use omp_lib
	implicit none

	type(sed_struct) sed_lin(link_count)	      
	real(8) ust_idx(riv_count)
	real(8) hr_idxa(riv_count)
	real(8) qr_ave_idx(riv_count)
	real(8) qss_lin(link_count)
	real(8) ss_idx(riv_count)
!	real(8) ss_lin(link_count)   
	real(8) ust_lin(link_count)
	real(8) water_v_lin (link_count), hr_lin(link_count)
	real(8) qss_b(link_count) !added by Qin 2021/5/27
!	 real(8) width_idx(riv_count), len_riv_idx(riv_count)
	integer k, kk, kkk, n, m, l, lll,f
	real(8) :: ss1, ss11, ss2, ss3, ss22
!    real(8) :: ss4, ss5, DDD
!	 real(8) :: qsss1
	integer ktest
	real(8) :: u_flax, Frn

	real(8) :: bbeta, ss
	real(8) :: Cae, Cae1, Csb,Ce
	real(8) :: fp
	real(8) :: z_star
	integer kkmax, i, j, t
	real(8) :: alfa2 , sum_sus_dzb !for bedlock
	real(8) :: ristar, Emb_min
	real(8) :: thet !2021/7/18
	real(8) :: ss_c, ss_q, divloss

	
	real(8), parameter :: pai = 3.141592
	real(8), parameter :: csbar = 0.2
	real(8), parameter :: kappa = 0.0015
	real(8), parameter :: rw = 1.0
	real(8), parameter :: rs = 2.65


!	 DDD = 0.01
Emb_min = Emc
divloss =1.d0


!$omp parallel do private(k,m)
	  do l = 1, link_count
		qss_lin(l) = 0.d0
		ss_lin(l) = 0.d0  !????
		qss_b(l) = 0.d0 !adde by Qin 2021/5/27
	   	do m = 1, Np
			sed_lin(l)%qsisum(m) = 0.d0
!			sed_lin(l)%ssi(m) = 0.d0 
			sed_lin(l)%Dsi(m) = 0.d0 !20220511
			sed_lin(l)%Esi(m) = 0.d0 !20220924
	    enddo
		k = link_idx_k(l)
		if(damflg(k)>0) cycle
!modified 20230924	
		!if (hr_idxa(k).le.0.d0 .or. water_v_lin(l).lt.0.01 .or. qr_ave_idx(k).le.0.) then !modified 20230924
		if (hr_idxa(k).le.min_hr .or. water_v_lin(l).lt.0.01 .or. qr_ave_idx(k).le.0.) then !modified 202406160924
	   		do m = 1, Np
				if (sed_lin(l)%qsi(m).gt.0.d0) sed_lin(l)%Dsi(m) = sed_lin(l)%qsi(m)/area_lin(l)
				sed_lin(l)%ssi(m) = 0.d0			
				sed_lin(l)%qsi(m) =0.d0 
				ss_lin(l) = ss_lin(l) + sed_lin(l)%ssi(m)
		    	qss_b(l) = qss_b(l)+sed_lin(l)%qsi(m) !added by Qin 2021/5/30
			enddo
		else
	 		do m = 1, Np
				sed_lin(l)%ssi(m) = sed_lin(l)%qsi(m) / qr_ave_idx(k) + slo_vol_remain(l,m)/water_v_lin(l)! add the remained slope erosion sediment to suspended sediment 2022/3/30
				if (isnan(sed_lin(l)%ssi(m))) then !check 2022/4/7
					write(*,*) l, k,m,sed_lin(l)%qsi(m),qr_ave_idx(k),water_v_lin(l),slo_vol_remain(l,m)
					stop "ssi is Nan"
				endif	
				sed_lin(l)%qsi(m) = sed_lin(l)%ssi(m)*qr_ave_idx(k)	
				slo_vol_remain(l,m) = 0.d0
				ss_lin(l) = ss_lin(l) + sed_lin(l)%ssi(m)
				qss_b(l) = qss_b(l)+sed_lin(l)%qsi(m) !added by Qin 2021/5/30
	   		enddo
		end if
	end do  !20260223
!pause'2'
!$omp single
    do l = 1, link_count    !20260223 added this loop to evaluate ss_lin(l) properly
		k = link_idx_k(l)
		if(damflg(k)>0) cycle
	    if(ss_lin(l) > alpha_ss2) then
	    	ss_c = 0.d0
	    	ss_q = 0.d0
			do m = 1, Np
				!write(*,*)'l, m, ss_lin(l),slo_vol_remain(l,m) ',l, m, ss_lin(l),slo_vol_remain(l,m)
				slo_vol_remain(l,m) = sed_lin(l)%ssi(m) / ss_lin(l) * (ss_lin(l)-alpha_ss2)*water_v_lin(l) ! To avoid the high concentration suspended sediment flow due to the large sediment supply from slope erosion 2022/3/30
	!			slo_vol_remain(l,m) = sed_lin(l)%ssi(m) / ss_lin(l) * (ss_lin(l)-alpha_ss2)*qr_ave_idx(k)*ddt! To avoid the high concentration suspended sediment flow due to the large sediment supply from slope erosion 20230924
				if (isnan(slo_vol_remain(l,m))) then !check 2021/6//3
					write(*,*) l, k,m,sed_lin(l)%qsi(m),sed_lin(l)%ssi(m),ss_lin(l),alpha_ss2,qr_ave_idx(k),water_v_lin(l), hr_lin(l)
					stop "slo_vol_remain is not correct check2"   
				endif
				if (slo_vol_remain(l,m).lt.0.0001d0) slo_vol_remain(l,m) = 0.d0 !check this regulation 2022/3/30
				sed_lin(l)%ssi(m) = sed_lin(l)%ssi(m) / ss_lin(l) * alpha_ss2
				sed_lin(l)%qsi(m) = sed_lin(l)%ssi(m)*qr_ave_idx(k) ! To avoid the high concentration suspended sediment flow due to the large sediment supply from slope erosion 2022/3/30
				ss_c = ss_c + sed_lin(l)%ssi(m)
				ss_q = ss_q+sed_lin(l)%qsi(m)
				!write(*,*)'l, m, slo_vol_remain(l,m) ',l, m, slo_vol_remain(l,m)
		 	enddo
		 	ss_lin(l) =  ss_c
		 	qss_b(l) = ss_q
			!write(*,*)'l, ss_lin(l) ',l, ss_lin(l) 
			!stop
	    endif
	
		!----Driftwood setting  
		if(j_drf ==1)then
			if(qr_ave_idx(k).le.1e-2) cw(l) = 0.d0
			cw(l) = qw(l)/qr_ave_idx(k)
!			if(cw(l).gt.alpha_ss1) cw(l) = alpha_ss1   !removed 20250225
		end if
    enddo
!$omp end single
!pause'before deposition yoro-work'

!------------------calculation of Erosion and Deposition of suspended sediment
!$omp parallel do private(k,m,fp,ss,z_star,Cae1,Cae,u_flax,ristar,Ce,sum_sus_dzb,alfa2)
    do l = 1, link_count
		k = link_idx_k(l)
		if(damflg(k)>0) cycle !There is no bed erosion in dam reservior; added by Qin 
		sum_sus_dzb = 0.d0	
		!if (hr_idxa(k).le.0.d0 .or. water_v_lin(l).lt.0.01 .or. qr_ave_idx(k).le.0.) then  !modified 20231025
		if (hr_idxa(k).le.min_hr .or. water_v_lin(l).lt.0.01 .or. qr_ave_idx(k).le.0.) then !modified 20240616
			do m =1 , Np		
					sed_lin(l)%Esi(m) = 0.d0
			enddo		
	 	else	
				do m = 1, Np
					fp = sed_lin(l)%fm(m)
					ss = sed_lin(l)%ssi(m)
					sed_lin(l)%Dsi(m) = ss * w0(m)+ sed_lin(l)%Dsi(m)
					if (isnan(sed_lin(l)%Dsi(m))) then !check 2021/6//3
						!write(*,*) l, k, Csb, ss, w0(m),bbeta
						write(*,*) l, k, ss, w0(m)
						stop "Dsi is Nan"
					endif	
					if (sed_lin(l)%Dsi(m) < 0.d0) then !check 20240429
						!write(*,*) l, k, Csb, ss, w0(m),bbeta
						write(*,*) l, k,sed_lin(l)%Dsi(m), w0(m), qr_ave_idx(k), Emb_lin(l), hr_idxa(k), zb_riv_slope_lin(l)
						stop "Dsi < 0"
					endif	
			
					if(isuseq == 1)then		  
						!Lane-Kalinske formula
						if (w0(m) .lt. ust_lin(l)) then
			!           	if (w0 .lt. ust_idx(k).and.ust_idx(k).ge.0.12) then
							z_star = ust_lin(l)/w0(m)
							Cae1 = 5.55*(0.5*z_star*exp(-(1./z_star)**2.))**1.61 * 1.0E-4 
							if(Cae1.ge.0.1d0) Cae1 = 0.1d0
							Cae = Cae1 *fp
						else
							Cae = 0.d0
						endif
						if (isnan(Cae)) then
							write(*,*) k, m, w0(m), ust_lin(l), z_star,fp
							stop "Cae is Nan"
						endif
						sed_lin(l)%Esi(m) = Cae * w0(m)
						if(Emb_lin(l) .le. Emb_min) sed_lin(l)%Esi(m) = 0.d0 
						sed_lin(l)%Dsi(m) = Csb * w0(m) + sed_lin(l)%Dsi(m)
							
					elseif(isuseq == 2)then
						! density stratified flow
			!			u_flax = qr_ave_idx(k)/width_idx(k)/hr_idxa(k)
						u_flax = qr_ave_idx(k)/width_lin(l)/hr_idxa(k)	
						if(hr_idxa(k).le.min_hr) u_flax = 0.d0 !modified by Qin 2021/11/06		
			!			if(u_flax > 10.) u_flax = 10.   !check this reguration	
						if (w0(m) .ge. ust_lin(l).or.ust_lin(l).le.0.d0) then
							sed_lin(l)%Esi(m) = 0.d0
						elseif(ss >= 0.5*csbar.or.u_flax.le.0.d0) then
							sed_lin(l)%Esi(m) = 0.d0
						elseif(Emb_lin(l) .le. Emb_min .or. hr_idxa(k).le.min_hr) then !no erosion when bedload layer less than minmum bedload layer depth or water depth less than 0.1 m 20211006
							sed_lin(l)%Esi(m) = 0.d0					
						else
			!Richardson number 
							ristar = (rs/rw-1.)*csbar*grav*hr_lin(l)/(u_flax**2.) 
							Ce = kappa/ristar*u_flax*csbar/w0(m)
			!revised by Qin 20211111
							if(Ce.ge.0.1d0) then
								Ce = 0.1d0	
								ristar=kappa/Ce*u_flax*csbar/w0(m)
							endif
			!				sed_lin(l)%Esi(m) = fp*Ce*w0
							sed_lin(l)%Esi(m) = fp*kappa/ristar*u_flax*csbar
							if (isnan(sed_lin(l)%Esi(m))) then
								write(*,*) k, m, w0(m), ust_lin(l), u_flax, Ce, fp
								stop "Esi(m) is Nan"
							endif

						endif
					endif
				!bedlock part	  
					if(zb_riv_idx(k)- zb_roc_idx(k) .lt. Emb_lin(l))then 
						alfa2 = (zb_riv_idx(k) - zb_roc_idx(k))/Emb_lin(l)
						if(alfa2.ge.1.)alfa2=1.
						if(alfa2 <= 0.d0 .or.Emb_lin(l).le.Emb_min) alfa2 = 0. !modified 20221206
						sed_lin(l)%Esi(m) = sed_lin(l)%Esi(m) * alfa2
					end if
!added 20240429					
				sum_sus_dzb = sum_sus_dzb + (sed_lin(l)%Dsi(m) - sed_lin(l)%Esi(m) )*ddt*dlambda
				enddo !end of do m = 1, Np
				if((zb_riv_idx(k)- zb_roc_idx(k) + sum_sus_dzb) .lt. 0.d0 .and. sum_sus_dzb < 0.d0) then
				alfa2 = (zb_riv_idx(k)- zb_roc_idx(k) + sum_sus_dzb)/ sum_sus_dzb 
				do m = 1, Np
				 sed_lin(l)%Esi(m) = sed_lin(l)%Esi(m) * (1.d0-alfa2)
				enddo 
				endif
		endif !endif for hr<0 water_v_lin <0.01 qr<0
		if(sed_lin(l)%Esi(2).gt.1.0E+2) then
			!if(ust_idx(k).gt.0.01) then
			write(*,'(a)') '  m      dsi       fp       w0      ssi            Csb           Dsi            Cae            Esi'
			write(*,'(a,f)') 'ust=', ust_lin(l)
			write(*,'(a,i)') 't=', t
			write(*,'(a,i7)') 'l=', l
			pause'check deposition'
		endif	

    enddo !end of do l = 1, link_count
!pause'after deposition'

!  In case of channel capasity loss due to deposition, pass sus_load to the downstream link  !Machino 20240119
	do l = 1, link_count
		k = link_idx_k(l)
	 	thet = zb_riv_slope0_lin(l)
		!if(zb_riv_idx(k)-zb_riv0_idx(k) + depth_idx(k)*0.1 >  depth_idx(k) ) then !check the condition 'depth_idx(k)*0.1'
		!if(depth_idx(k) < min_hr*2. )then !check this condition
		if(n_link_depth(l)==1) then  !20240229
			if (thet.gt.max_slope) exit		 
			if (link_0th_order(l) == 0) exit			
			do m = 1, Np
				do n = 1, 8
					if(up_riv_lin(l,n) == 0) exit 	
					lll = up_riv_lin(l, n)	
					sed_lin(l)%qsisum(m) = sed_lin(l)%qsisum(m) - sed_lin(lll)%qsi(m) !m3/s	 erosion negative
				end do										!calculate input sediment from two unit channels at first 
				sed_lin(l)%qsi(m) = -sed_lin(l)%qsisum(m)    !modify the output sediment so that in = out
				sed_lin(l)%qsisum(m) = 0.					!initialization for the next loop
			end do
			if(j_drf == 1)then
				do n = 1, 8
					if(up_riv_lin(l,n) == 0) exit 	
					lll = up_riv_lin(l, n)	
!					qwsum(l) = qwsum(l) - qw(lll) !m3/s	 
				end do										!calculate input sediment from two unit channels at first 
!				qw(l) = -qwsum(l)    !modify the output sediment so that in = out
!				qwsum(l) = 0.					!initialization for the next loop		
			end if
		end if
	end do

!!$omp parallel do private(k,m,n,lll,divloss,f)
!--------------------trasportation of suspended sediment
  	do l = 1, link_count
		k = link_idx_k(l)
		divloss =1.
		if(damflg(k)>0)cycle	
		!added for diversion 20250329	
		if(div_switch==2)then
			do f = 1, div_id_max
			if(div_org_idx(f)==k) divloss=1.-div_rate(f)
			enddo
		endif
		do m = 1, Np
			do n = 1, 8
	!				if(up_riv_lin(l,n) == 0) cycle !revised by Qin 2021/5/27
					if(up_riv_lin(l,n) == 0) exit 
					lll = up_riv_lin(l, n)
					sed_lin(l)%qsisum(m) = sed_lin(l)%qsisum(m) - sed_lin(lll)%qsi(m)*divloss !m3/s   !Á¨¶Âè∑„Å´Ê≥®ÅEΩ?!!!  ffd ÅEΩ?-„ÅÆ„Å®„ÅçÔøΩ??Á©çÔøΩ?+„ÅÆ„Å®„Åç‰æµÅEΩ? modified for diversion 20250329
				end do
					sed_lin(l)%qsisum(m) = sed_lin(l)%qsisum(m) + sed_lin(l)%qsi(m)
		end do
		!-----Driftwood Transport
		if(j_drf == 1)then
			qwsum(l) = 0.d0
			if(link_0th_order(l) == 0)then
				qwsum(l) = qwsum(l) + qw(l)
			else
				do n = 1, 8
!					if(up_riv_lin(l,n) == 0) cycle !revised by Qin 2021/5/27
					if(up_riv_lin(l,n) == 0) exit 
					lll = up_riv_lin(l, n)
					qwsum(l) = qwsum(l) - qw(lll)*divloss !m3/s   !Á¨¶Âè∑„Å´Ê≥®ÅEΩ?!!!  ffd ÅEΩ?-„ÅÆ„Å®„ÅçÔøΩ??Á©çÔøΩ?+„ÅÆ„Å®„Åç‰æµÅEΩ? modified for diversion 20250329
				end do
				qwsum(l) = qwsum(l) + qw(l)
			end if
		end if
  	end do
!pause'sum of ssload'

!!$omp parallel do private(k,m,alfa2)	  
	do l = 1, link_count
		k = link_idx_k(l) 	
		ss_lin(l) = 0.d0 !20220512 
		qss_lin(l) = 0.d0
		if(damflg(k)>0) then !
			do m = 1, Np
				sed_lin(l)%qsi(m) = sed_lin(l)%ssi(m) * qr_ave_idx(k) !m3/s
				qss_lin(l) = qss_lin(l) + sed_lin(l)%qsi(m) !m3/s
			enddo
			cycle
		endif
		!if (hr_idxa(k).le.0.d0 .or. water_v_lin(l).lt.0.01 .or. qr_ave_idx(k).le.0.) then  !modified 20231025
		if (hr_idxa(k).le.min_hr .or. water_v_lin(l).lt.0.01 .or. qr_ave_idx(k).le.0.) then !modified 20240616
				do m = 1, Np
					sed_lin(l)%ssi(m) = 0.d0
					sed_lin(l)%Esi(m) = 0.d0 !modified by Qin 2021/6/3
					sed_lin(l)%Dsi(m) =sed_lin(l)%Dsi(m)-1.*sed_lin(l)%qsisum(m)/area_lin(l) !modified 20230924
					if(sed_lin(l)%Dsi(m).lt.0.d0) then
						sed_lin(l)%qsisum(m) = -1.*sed_lin(l)%Dsi(m)*area_lin(l) !modified 20230924
						sed_lin(l)%Dsi(m) = 0.d0
					endif
					if (isnan(sed_lin(l)%Dsi(m))) then
						write(*,*)l,m, sed_lin(l)%qsisum(m), area_lin(l), sed_lin(l)%qsi(m)
						stop "Check the suspended sediment budget"
					endif
				enddo
		elseif(n_link_depth(l)==1) then    !  In case of channel capasity loss due to deposition, pass sus_load to the downstream link  !20240229
				do m = 1, Np
					sed_lin(l)%Esi(m) = 0.d0 
					sed_lin(l)%Dsi(m) = 0.d0 
				enddo			
		else !-----revised by Harada 2021/5/21	
				do m =1,Np
					sed_lin(l)%ssi(m) = sed_lin(l)%ssi(m) + (-1.* sed_lin(l)%qsisum(m) + sed_lin(l)%Esi(m)*area_lin(l) - sed_lin(l)%Dsi(m)*area_lin(l)) * ddt/water_v_lin(l)
					if(sed_lin(l)%ssi(m).lt.0.d0) then
						sed_lin(l)%ssi(m) = sed_lin(l)%ssi(m) - (-1.* sed_lin(l)%qsisum(m) + sed_lin(l)%Esi(m)*area_lin(l) - sed_lin(l)%Dsi(m)*area_lin(l)) * ddt/water_v_lin(l) !‰∏ÄÂõûÔøΩ??„Å´Êàª„Åó„Å¶
						sed_lin(l)%Dsi(m) = -1.* sed_lin(l)%qsisum(m)/area_lin(l) + sed_lin(l)%Esi(m)      !D„Çí‰øÆÊ≠£„Åó„Å¶
						if(sed_lin(l)%Dsi(m) < 0.) sed_lin(l)%Dsi(m) = 0.d0 !„ÅæÅEΩ?D„ÅåÔøΩ?„ÅÆÊôÇÔøΩ?ÅEΩ‰∏ä‰∏ãÊµÅEøΩ?ÅEΩ„Ç¢„É≥„Éê„É©„É≥„Çπ„Å´„Çà„Çã„ÅÆ„ÅßD„ÅØ0„Å®„Å™ÅEΩ?
						sed_lin(l)%ssi(m) = sed_lin(l)%ssi(m) + (-1.* sed_lin(l)%qsisum(m) + sed_lin(l)%Esi(m)*area_lin(l) - sed_lin(l)%Dsi(m)*area_lin(l)) * ddt/water_v_lin(l) !Êñ∞„Åó„ÅÑD„Çí‰Ωø„Å£„Å¶„ÇÇ„ÅÜ‰∏ÄÂ∫¶ssi„ÇíË®àÔøΩ?
						if(sed_lin(l)%ssi(m).lt.0.d0) sed_lin(l)%ssi(m) = 0.d0 !modified 20240101	
						if (isnan(sed_lin(l)%ssi(m))) then
							write(*,*) k, m, sed_lin(l)%ssi(m),water_v_lin(l),sed_lin(l)%Dsi(m),sed_lin(l)%Esi(m),sed_lin(l)%qsisum(m)
							stop "suspended load concentration is Nan"
						endif
					endif
					ss_lin(l) = ss_lin(l) + sed_lin(l)%ssi(m)
				enddo
		endif !hr_idxa(k).le.0.d0.or. water_v_lin(l).lt.0.01 .or. qr_ave_idx(k).le.0.d0	 
			
		if (ss_lin(l).gt.alpha_ss2) then
			do m =1, Np
				slo_vol_remain(l,m) = slo_vol_remain(l,m)+ sed_lin(l)%ssi(m) / ss_lin(l) * (ss_lin(l)-alpha_ss2)*water_v_lin(l) ! no bed variation at upstream end unit channel (>8 degree) 20221206
				if (isnan(slo_vol_remain(l,m))) then !check 2021/6//3
					write(*,*) l, k,m,sed_lin(l)%qsi(m),sed_lin(l)%ssi(m),ss_lin(l),alpha_ss2,qr_ave_idx(k),water_v_lin(l), hr_lin(l)
					stop "slo_vol_remain is not correct check3"   
				endif
				if (slo_vol_remain(l,m).le.0.d0) slo_vol_remain(l,m) = 0.d0 
				sed_lin(l)%ssi(m) = sed_lin(l)%ssi(m)/ss_lin(l)*alpha_ss2
			enddo
			ss_lin(l) = alpha_ss2 !modified 20240315
		endif	

		do m = 1, Np
			sed_lin(l)%qsi(m) = sed_lin(l)%ssi(m) * qr_ave_idx(k) !m3/s
			qss_lin(l) = qss_lin(l) + sed_lin(l)%qsi(m) !m3/s
		enddo
	enddo
!pause'end of suspended sediment'

end subroutine susload2

! >>>>>>>>>>>>>>>>>>>>>>> CALCULATING BED LOAD
!      subroutine bedload( sed_idx, hr_idxa, qr_ave_idx, ust_idx, qsb_idx, width_idx, len_riv_idx, zb_riv_idx )
	 subroutine bedload2( sed_lin, hr_idxa, hr_lin, qr_ave_idx, ust_idx, ust_lin,qsb_lin)
		use globals
		use sediment_mod
		use dam_mod!, only:dam_switch, damflg
		implicit none

		type(sed_struct) sed_lin(link_count)
		real(8) ust_idx(riv_count)		 
		real(8) qr_ave_idx(riv_count)
		real(8) hr_lin(link_count), ust_lin(link_count), qsb_lin(link_count)
		real(8) hr_idxa(riv_count)
!	     real(8) width_idx(riv_count), len_riv_idx(riv_count)
!	     real(8) zb_riv_idx(riv_count)
		real(8) :: fp, u_flax, ks, lg
		integer :: k, m, l
		real(8) :: uscm, usci, u_se, qbi
		real(8) :: tscm, dsm
		real(8) :: t_star, t_starE, qbe, den
		real(8) :: R_tcts
		real(8) :: zb55, Frn

		real(8) :: ust_idx1, ust_idx22, ust_idx33, slope
		integer :: kkk, kkkk, kzb55, n
		real(8) :: S_angle, C_K1, C_K2, C_fd, C_ff 
		real(8) :: Dcr, alfa, hs_l
		real(8) :: thet, thetr, Emb_min,thet2
		real(8), parameter :: C_kf		= 0.16
		real(8), parameter :: C_kd		= 0.0828
		real(8), parameter :: C_e		= 0.85
		real(8), parameter :: Csta		= 0.6

		Emb_min = Emc		
		S_angle = 34./180.d0*3.14159d0



!pause'in bedload'
!-----------??????v?Z?I??----------------
!pause'1'

!$omp parallel do private(k, m,dsm,uscm,tscm,u_flax,ust_idx1,Frn,thet2,thetr,alfa, &
!$omp                     fp,R_tcts,t_Crit,den,usci,ks,lg,u_se,t_star,t_starE,qbe,Dcr,hs_l,C_K1,C_K2,C_fd,C_ff)
    do l = 1, link_count   
   		k = link_idx_k(l)
! Determine Mean shear stress (Iwagaki Eqn)
		dsm = sed_lin(l)%dmean*100.
		if (dsm .ge. 0.303 ) then
			uscm = sqrt(80.9*dsm)
		elseif (dsm .ge. 0.118) then
			uscm = sqrt(134.6*dsm**(31./22.))
		elseif (dsm .ge. 0.0565) then
			uscm = sqrt(55.0*dsm)
		elseif (dsm .ge. 0.0065) then
			uscm = sqrt(8.41*dsm**(11./32.))
		elseif (dsm .lt. 0.0065) then
			uscm = sqrt(226.0*dsm)
		endif
		dsm = dsm/100.
		uscm = uscm/100.
		tscm = uscm**2./((s-1.)*grav*dsm)

!pause'2'

!---use the width of link	modified by Qin	
		thet2 = zb_riv_slope_lin(l) !unit: degree; 20240508
		!thet2  = ws_slope_lin(l) !unit: gradient;use water surface slope ;modified 20240424 
		thetr = thet2/180.*3.14159 !20240508
		!thetr = atan(thet2) !20240424
		u_flax = qr_ave_idx(k)/width_lin(l)/hr_idxa(k)	
		if(hr_idxa(k).le.min_hr .or. u_flax .le.0.d0) u_flax = 0.d0 !2021/11/6 20240428
		if(hr_idxa(k).le.min_hr)then !added by Qin			
			ust_idx1 = 0.d0
		else
			ust_idx1 = (ns_river*sqrt(grav)* u_flax)/(hr_idxa(k)**(1./6.))
		endif
		if(hr_idxa(k).le.min_hr) then !modified by Qin 2021/6/3
			Frn = 0.d0
		else
			Frn = u_flax/sqrt(9.81*hr_idxa(k))
		endif
		ust_lin(l) = ust_idx1 
		qsb_lin(l) = 0.d0

			!if(thetr.le.min_slope) thetr = min_slope
			! when water surface slope <= min_slope, then no sediment transportation;modified 20240424
			if(tan(thetr).le.min_slope)then ! modified 20240427
			Emb_lin(l) = Emb_min
			else
		    Emb_lin(l) = ust_lin(l)**2./((s-1.)*csta/2.*grav*cos(thetr)	 &      !Exchange layer thickness?@e_mov
		    *(tan(S_angle)-tan(thetr)))
			endif
		    if(Emb_lin(l) >= hr_lin(l)) Emb_lin(l) = hr_lin(l) !modified by Qin;2021/7/16
			if(Emb_lin(l) <= Emb_min) 	Emb_lin(l) = Emb_min
			hs_l = Emb_lin(l)

			!bedrock_alfa
			alfa = (zb_riv_idx(k) - zb_roc_idx(k))/Emb_lin(l) !check the top of bedloadlayer has been defined as the theoritical riverbed surface  20231028
			if(alfa.ge.1.) alfa=1.
			if(Emb_lin(l)<= Emb_min.or.alfa<=0.d0) alfa= 0.d0 
		    do m = 1, Np
				fp = sed_lin(l)%fm(m)
				if (dsi(m)/dsm .le. 0.4) then
			    	R_tcts = 0.85 * dsm / dsi(m)
				else
			    	R_tcts = (1.28/(1.28+log10(dsi(m)/dsm)))**2.
				endif
				t_Crit = R_tcts * tscm
				den = (s-1.)*grav*dsi(m)
				usci = sqrt(t_Crit * den)

				if (hr_idxa(k).gt.min_hr.and.ust_lin(l)>0.d0) then	!--modified by Qin 2021/6/3
					ks = 1.+2.*ust_lin(l)**2./(s-1)/grav/dsi(m)
					lg = 6.0 + 2.5 * log(hr_idxa(k)/dsi(m)/ks)
					u_se = u_flax/lg
				else
					u_se = 0.d0 
				end if

				t_star = ust_lin(l)**2./den				 
				t_starE = u_se**2./den
				   
				!isedeq = 1 !1: Ashida-michiue, 2: MPM
				if(isedeq.eq.1) then
				!ashida-michiue equation
					if(hs_l.gt.0.5d0*hr_lin(l).and.hr_lin(l).gt.0.d0) hs_l= 0.5d0*hr_lin(l)
							if (ust_lin(l).gt.usci) then
								qbe=(17.*t_starE**1.5)*(1-t_Crit/t_star)*(1-usci/ust_lin(l))
								sed_lin(l)%qbi(m) = qbe*sqrt((s-1)*grav*dsi(m)**3.)*width_lin(l)*fp*alfa	
								if(hs_l/dsi(m).lt.0.18d0.or.hr_lin(l).le.min_hr)then !modified by Qin; 20210829
									sed_lin(l)%qbi(m) = 0.d0															
								elseif(Emb_lin(l).le.Emb_min)then
								sed_lin(l)%qbi(m) = 0.d0	
								endif		
							else
								sed_lin(l)%qbi(m) = 0.d0
							end if  					 
				elseif(isedeq.eq.2) then
				!MPM equation
					if(hs_l.gt.0.5d0*hr_lin(l).and.hr_lin(l).gt.0.d0) hs_l = 0.5d0*hr_lin(l)
						if (ust_lin(l).gt.usci) then
							qbe = 8.0 * (t_star - t_Crit)**1.5
							sed_lin(l)%qbi(m) = qbe*sqrt((s-1)*grav*dsi(m)**3.)*width_lin(l)*fp*alfa !?P???m3/s	
							if(hs_l/dsi(m).lt.0.18d0.or.hr_lin(l).le.min_hr)then !modified by Qin; 20210829
								sed_lin(l)%qbi(m) = 0.d0					
			!------------------modified by Qin 2021/8/3
							elseif(Emb_lin(l).le.Emb_min)then
								sed_lin(l)%qbi(m) = 0.d0
							endif							
						else
							sed_lin(l)%qbi(m) = 0.d0
						end if
				elseif(isedeq.eq.3) then
				!Egashira formula
					Dcr = ust_lin(l)**2./(tscm*(s-1.)*grav) !revised by Qin; 20210829
					if(hs_l > 0.5d0*hr_lin(l).and.hr_lin(l).gt.0.d0) hs_l = 0.5d0*hr_lin(l)	
					if (dsi(m) > Dcr) then
							sed_lin(l)%qbi(m) = 0.d0
					elseif(hs_l/dsi(m).lt.0.18d0.or.hr_lin(l).le.min_hr)then !modified by Qin; 20210829
							sed_lin(l)%qbi(m) = 0.d0	
					elseif(Emb_lin(l).le.Emb_min)then !no erosion when bedload layer less than minmum bedload layer depth 
							sed_lin(l)%qbi(m) = 0.d0
					else
						C_K1 = 1./cos(thetr)/(tan(S_angle)-tan(thetr))
						C_K2 = 1./Csta*2.*(1.-hs_l/hr_lin(l))**0.5
						C_fd = C_kd*(1.-C_e**2.)*(1.+(s-1.))*(Csta/2.)**(1./3.)
						C_ff = C_kf*(1.-Csta/2.)**(5./3.)*(Csta/2.)**(-2./3.)
						sed_lin(l)%qbi(m) = 4./15.*C_K1**2.*C_K2/sqrt(C_fd+C_ff)*t_star**2.5       &
													*sqrt((s-1.)*grav*dsi(m)**3.)*width_lin(l)*fp*alfa										

						if(sed_lin(l)%qbi(m) <= 0.) sed_lin(l)%qbi(m) = 0.d0
					endif
				endif

				if (dam_switch== 1)then !added by Qin 2021/4/22
					k = link_idx_k(l)
					if(damflg(k).gt.0)	sed_lin(l)%qbi(m) = 0.d0 
				endif   
				qsb_lin(l) = qsb_lin(l) + sed_lin(l)%qbi(m)

			!pause'6'
	    	end do ! end of do m = 1, Np

		if(ust_lin(l).ge.1000000.d0.or.isnan(ust_lin(l))) then	   
		write(*,'(a)') '    m         dsi          fp        usci    t_Crit/tscm     sed_idxqbi'
		write(*,'(a,f,a,f,a,i,a,i)') 'zb_riv_idx(k)=', zb_riv_idx(k), 'ust_lin(l)=',ust_lin(l), 'l=',l,'k=', k
		write(*,'(a,f)') 'length=', len_riv_idx(k)
		write(*,'(a,f)') 'distance=', dis_riv_idx(k)
		!write(*,'(a,f)') 'zb_slope=', zb_riv_slope_idx(k)
		!write(*,'(a,f)') 'tan(zb_slope)=', tan(zb_riv_slope_idx(k)*3.14/180.)
		write(*,'(a,f)') 'zb_slope_lin=',zb_riv_slope_lin(l)
		write(*,'(a,f)') 'Frn=', Frn
		write(*,'(a,f)') 'slope=', slope
		write(*,'(a,f)') 'zb55=', zb55
		write(*,'(a,f)') 'kzb55=', kzb55
		write(*,*) 'hr_idxa(k)=', hr_idxa(k)
		write(*,'(a,f)') 'ust=', sqrt(9.81*hr_idxa(k)*slope)
		write(*,*)
		write(*,'(a,f)') 'ust_idx22=', ust_idx22
		write(*,'(a,f)') 'dsm=      ', dsm

		write(*,'(a)') '  m           dsi            fp         qbi(m)'
		do m = 1, Np
		write(*,'(i3,2f14.5,f15.5)') m, sed_lin(l)%dsed(m), sed_lin(l)%fm(m), sed_lin(l)%qbi(m)
		enddo
		stop
		endif

	enddo !end of do l = 1, link_count
!!pause'last of bedload'

	 end subroutine bedload2


!-----Sediment transportation calculation in slope area
	subroutine slo_sedi_cal ( sed_lin, hs_idx,qs_ave_idx,qsur_ave_temp_idx,slo_sedi_cal_duration) 
		use globals
        use sediment_mod
        implicit none
        type(sed_struct) sed_lin(link_count)
		real(8) hs_idx(slo_count),qs_ave_idx(i4,slo_count),qsur_ave_temp_idx(i4,slo_count) !revised
		!real(8) qss_slope(slo_count),slo_Dsi(slo_count,Np),slo_qsi(slo_count,Np),slo_ssi(slo_count,Np),slo_sur_zb(slo_count),fmslo(slo_count,Np),ns_MS(slo_count)
		real(8) :: V, Frn, V_ave
		real(8) :: bbeta, ss,cc
		real(8) :: Cae, Cae1, Csb,Ce
		real(8) :: fp,sedi
		real(8) :: z_star, ristar
		real(8) :: I_slo, ns,hl,emb_gully,S_angle,thetr,h_dash,d80
		real(8) :: usta, E,dh,D,dzslo_before,dzslo_fp_before,sum_fmslo,soildepth_fp_idx,soildepth_fp_before,ql,total_w_v,x_up_inflow ,y_up_inflow, sum_fmgully,x_outflow,y_outflow !modified 20231206 20231226 20240101
		real(8) :: E_v,E_speci,h_p,E_min
		real(8):: wl !flood propagate velocity under kinematic wave approximate
		real(8) :: RC,DR,Ke_DT,Ke_LD,dsm, d50,TC, uscm,omega,omega_c2,param_c, param_eta,param_beta,param_b,gammaa_k,dgammaa
		real(8) :: AA,BB,h_c
		real(8) :: ddt_slo,	cal_time, slo_sedi_cal_duration,grad_H,over_ss,D_inf
		integer:: lt !delay time = t : time step
		integer :: depo_t ! depostion time
		real(8), parameter :: pai = 3.141592
		real(8), parameter :: csbar = 0.2
		real(8), parameter :: kappa = 0.0015
		real(8), parameter :: rw = 1.0
		real(8), parameter :: rs = 2.65
		real(8), parameter :: slope_20=0.36397d0 ! 20 degree 
		real(8), parameter :: csta = 0.6d0 
        integer k, kk,kkk,i, j,ii,jj, m, t,l, m1,n,slope_end_idx,nn,ll,n_count,riv_k
		integer iii, jjj, down_xk, down_yk,alfa !added 20240101

		S_angle=34./180.d0*3.14159d0

		!write (*,'(a,f,a,f)') 'conduct the slope sediment computation at time:', time - slo_sedi_cal_duration , 'dt for slope sediment computation =', dt_slo_sed 
!write (*,*) time, dt_slo_sed,t
		!if(slope_ero_switch==1)then ! gully/rill erosion on slope area (sediment path); added by Qin 20220730
		!ddt = 0.05d0
		n_count = 0
		cal_time=0.d0
do 
	n_count=n_count+1
	!write(*,*) 'n_count=', n_count
	ddt_slo = dt_slo_sed
	if(n_count*dt_slo_sed .ge. slo_sedi_cal_duration)then
	 	ddt_slo =slo_sedi_cal_duration-  (n_count-1)*dt_slo_sed
		if(ddt_slo .le. 0.01) then 
		ddt_slo = 0.01
		go to 8
		endif
		if(ddt_slo.gt.dt_slo_sed) ddt_slo = dt_slo_sed
	endif 
!pause '1'	
!$omp parallel do private(i,j,m,soildepth_fp_idx,gammaa_k,dgammaa,D_inf,dh,ns,V,h_dash,usta,thetr,emb_gully,alfa,fp,bbeta,Csb,z_star,Cae1,Cae,ristar,Ce,soildepth_fp_before,grad_H)
        do k = 1, slo_count
            i = slo_idx2i(k)
            j = slo_idx2j(k)
			if( domain(i,j) == 0 ) cycle
!modified 20240115			
			if(debris_switch.ne.0)then		
			!if(slo_grad(k).ge.0.25.and.zb(i,j)>100.d0 .and. riv(i,j).ne.1)cycle !modified for gofukuya river 20231101;20240121	
			!if(slo_grad(k).ge.0.25) cycle !20240126
			if(hki_g(k)>0.) cycle 	! no suspended sediment transportation taking place in debris flow grid cell; modified for gofukuya river 20240109	
			endif
!	write(*,*) "check1205 01"				
            !if(riv(i,j) == 1) cycle !modified 20230711 for inundation				
	    	ns = ns_slo_idx(k)
			do m = 1, Np
			 slo_qsisum(k,m) = 0.d0
			 slo_Dsi(k,m)  = 0.d0 
			 slo_Esi(k,m)  = 0.d0 
	        enddo
!	write(*,*) "check1205 02"				
			 gammaa_k=gammaa_idx(k)
			 dgammaa=1./(1.-gammaa_k)
			! h_surf(i,j) = hs_idx(k) - da_idx(k) ! modified 20231204
			! if(h_surf(i,j).le.0.) h_surf(i,j) = 0.d0
			! dh = h_surf(i,j)!(h-dap„ÅåË°®ÊµÅÊ∞¥„ÅÆÊ∞¥Ê∑±)
			!modified 20250618
			if(infil_w_depth(k)>da_idx(k))then
				D_inf=da_idx(k)
			else
				D_inf=infil_w_depth(k)
			endif	
			dh = hs_idx(k)-D_inf	
			if(dh.le.0.d0) dh=0.d0

            !modified20250503
			!if(riv(i,j)==1)then
			!water_v_cell(k) = hs_idx(k) *area
			!if(water_v_cell(k)<0.d0) water_v_cell(k)=0.d0
			!Qg(k) = sqrt(qs_ave_idx(1,k) ** 2.d0 + qs_ave_idx(2,k) ** 2.d0) * area
			!if(Qg(k).le. 0.d0) Qg(k) =0.d0
			!do m = 1, Np
			!slo_qsi(k,m)=Qg(k)*slo_ssi(k,m) 
			!enddo
			! all sediment in river cell supply into unit channel (without sediment inundation) 20250513
			if(riv(i,j)==1) then
			 water_v_cell(k) =0.d0
			 Qg(k)=0.d0
			 do m=1,Np
				slo_ssi(k,m) =0.d0
				slo_qsi(k,m) =0.d0
			 enddo
			else			
			!water_v_cell(k) = dh*area
			water_v_cell(k) = hs_idx(k)*area !20250514
			if(water_v_cell(k)<0.d0) water_v_cell(k)=0.d0
			Qg(k) = sqrt(qs_ave_idx(1,k) ** 2.d0 + qs_ave_idx(2,k) ** 2.d0) * area !20250514
			if(Qg(k).le. 0.d0) Qg(k) =0.d0
			if(dh.le.surflowdepth )then 
				!Qg(k)=0.d0
				do m = 1, Np  
				slo_Esi(k,m) = 0.
				slo_Dsi(k,m)=slo_ssi(k,m)*w0(m) !20250505
				slo_qsi(k,m) = slo_ssi(k,m)*Qg(k)
				if(slo_Dsi(k,m).le.0.d0) slo_Dsi(k,m) = 0.d0
					if (isnan(slo_ssi(k,m)).or.isnan(slo_Dsi(k,m))) then 
						write(*,*) k,m,slo_Dsi(k,m),slo_qsi(k,m) ,Qg(k),usta, dh,slo_grad(k),v
						stop "slope ssi or dsi is Nan"
					endif	
				enddo
			else !for dh>surflowdepth	
	!write(*,*) "check1205 05"
				!Qg(k) = sqrt(qsur_ave_temp_idx(1,k) ** 2.d0 + qsur_ave_temp_idx(2,k) ** 2.d0) * area
				if(slope_ero_switch==1)then !20230924 
					if(D_chan(k).le.0.d0 .or. B_chan(k).le.0.d0 .or. slo_grad(k).ge.slope_20) then !to prevent the erosion in urban area and the slope> 20 deg (deposition only) 20250402
						h_chan(k) = 0.d0
						q_chan(k)=0.d0
						do m = 1, Np 
							slo_qsi(k,m) = slo_ssi(k,m)*Qg(k)  
							slo_Dsi(k,m)  = slo_ssi(k,m) * w0(m) 		
							if(slo_Dsi(k,m).le.0.d0) slo_Dsi(k,m) = 0.d0
							slo_Esi(k,m) =0.d0	
						enddo								
					else	
	!modified 20240126
						if(len_slo_idx(1, k) > 0.d0 .and. len_slo_idx(2, k) > 0.d0 )then ! right-end side and down-end side boundaries of the basin do not have the width for x-flux and y-flux, respectively 
							q_chan(k) =  abs((qsur_ave_temp_idx(1,k) * area / len_slo_idx(1, k) * abs(slo_ang_cos(k)) + qsur_ave_temp_idx(2,k) * area / len_slo_idx(2, k) * abs(slo_ang_sin(k)))) / B_chan(k)
							!V = sqrt((qsur_ave_temp_idx(1,k) * area / (len_slo_idx(1, k) * dh)) ** 2.d0 + (qsur_ave_temp_idx(2,k) * area / (len_slo_idx(2, k) *dh)) ** 2.d0)
						elseif(len_slo_idx(1, k) == 0.d0 .and. len_slo_idx(2, k) > 0.d0 )then
							q_chan(k) =  abs(qsur_ave_temp_idx(2,k) * area / len_slo_idx(2, k) * slo_ang_sin(k)) / B_chan(k) 
							!V = abs((qsur_ave_temp_idx(2,k) * area / (len_slo_idx(2, k) *dh)))
						elseif(len_slo_idx(1, k) > 0.d0 .and. len_slo_idx(2, k) == 0.d0 )then
							q_chan(k) =  abs(qsur_ave_temp_idx(1,k) * area / len_slo_idx(1, k) * slo_ang_cos(k)) / B_chan(k)  
							!V = abs((qsur_ave_temp_idx(1,k) * area / (len_slo_idx(1, k) *dh)))
						else
							q_chan(k) = 0.d0 
							!V = 0.d0
						endif
						!modified 20250406
						h_chan(k) = (ns_MS(k) * q_chan(k) / sqrt(slo_grad(k) )) ** 0.6
						if(h_chan(k)>D_chan(k)) then
						h_chan(k) = D_chan(k)
						!elseif(h_chan(k)<dh)then
						! h_chan(k) = dh
						endif
						V = 1./ns_MS(k)*sqrt(slo_grad(k))*h_chan(k)**(2./3.)
						q_chan(k)=V*h_chan(k)
						usta = sqrt(grav * h_chan(k) *slo_grad(k))

						if(h_chan(k).lt.0.d0 .or. isnan(h_chan(k)))then
							write(*,*) k,ns_MS(k),q_chan(k),h_chan(k),B_chan(k),D_chan(k), len_slo_1d_idx(K),V,dh, qsur_ave_temp_idx(1,k),qsur_ave_temp_idx(2,k),len_slo_idx(2, k) , slo_ang_sin(k),len_slo_idx(1, k), slo_ang_cos(k)
							stop "water depth in sediment path is not correct "
						endif
						h_dash= 0.d0
						if(h_chan(k).gt.0.d0 .and. slo_grad(k).gt.0.d0)then
								thetr = atan(slo_grad(k))
								emb_gully = usta ** 2.d0 /((rs-1.) * csta / 2. * grav * cos(thetr) * (tan(S_angle)-tan(thetr)))!Exchange layer thickness
								h_dash = emb_gully/h_chan(k)
							if(h_dash.ge.0.9)then
								h_dash = 0.9d0
								emb_gully= h_dash*h_chan(k)
								Usta = sqrt(emb_gully*((rs-1.)*csta/2.*grav*cos(thetr)*(tan(S_angle)-tan(thetr))))
								grad_H= usta**2./grav/h_chan(k)
								V=1./ns_MS(k)*sqrt(grad_H)*h_chan(k)**(2./3.)
							endif
						!if(emb_gully.le.Emc) Usta=0.d0
						endif
				!------modified 20230420
				    endif !endif for urban or not
				!pause'2'

					if(h_chan(k).gt.0.d0.and.h_dash.gt.0.d0)then ! modified 20231226
						!added 20240101	  
						alfa=1.
						alfa = (slo_sur_zb(k)-zb_slo_idx(k))/emb_gully
						if(alfa.ge.1.) alfa=1.
						if(alfa <= 0.d0 .or.emb_gully.le.Emc) alfa=0. !turned off 20250407
						!if(alfa <= 0.d0) alfa=0.
						do m = 1, Np
							slo_qsi(k,m) = slo_ssi(k,m)*Qg(k)  
							slo_Dsi(k,m)  = slo_ssi(k,m) * w0(m) 		
							if(slo_Dsi(k,m).le.0.d0) slo_Dsi(k,m) = 0.d0
							slo_Esi(k,m) =0.d0	
							if(dsi(m).ge.h_chan(k)) cycle !added 20240131									
							if(fmslo(k,m)>0.d0 .and. alfa>0.d0)then 
								fp = fmslo(k,m)
								if (fp.lt.0.d0)then
									write(*,*) fp,k,m !, gully_sedi_dep(k,m)
									stop 'GSD of slope soil is negative'
								endif
								if (w0(m) < usta)then
									if(isuseq == 1)then  !Lane-Kalinske formula
										bbeta = 6.0*w0(m)/karmans/usta
										Csb = slo_ssi(k,m)*bbeta/(1.-exp(-1. * bbeta))
										slo_Dsi(k,m)  = Csb * w0(m)
										!if(slo_Dsi(k,m).le.1e-7) slo_Dsi(k,m) = 0.d0 !added 20240101	
										if(slo_Dsi(k,m).le.0.d0) slo_Dsi(k,m) = 0.d0								
										z_star = usta/w0(m)
										Cae1 = 5.55*(0.5*z_star*exp(-(1./z_star)**2.))** 1.61 * 1.0E-4 
										if(Cae1.ge.0.2d0) Cae1 = 0.2d0
										Cae = Cae1 *fp
										slo_Esi(k,m)  = Cae * w0(m)*alfa !modified 20240101 									
									elseif(isuseq == 2)then  !Density stratified flow formula
										ristar = (rs/rw - 1.) * csbar * grav * h_chan(k) / (V ** 2.) 
										Ce = kappa/ristar * V * csbar/w0(m)
										if(Ce.ge.0.3d0)then !modified 20231225
											Ce = 0.3d0
											ristar = kappa/Ce*V*csbar/w0(m)        !Once ristar is modified by this regulation for smallest size, the regulated ristar is applied to all classes.
										endif
											slo_Esi(k,m)  = fp * kappa / ristar * V * csbar * alfa !modified 20240101 										
									endif
								endif	
							endif

							slo_Esi(k,m) = slo_Esi(k,m)*area_chan(k)/area !modified 20250407
	!modified 20231226
							soildepth_fp_idx=(slo_sur_zb(k)-zb_slo_idx(k))*fmslo(k,m)
							!modified 20240101
							if(soildepth_fp_idx.le.0.) then
								slo_Esi(k,m)=0.d0
							else	
								soildepth_fp_before = soildepth_fp_idx
								soildepth_fp_idx = soildepth_fp_idx + (slo_Dsi(k,m) -slo_Esi(k,m) ) * ddt_slo * dgammaa
								if (soildepth_fp_idx.le.0.d0) then			
								slo_Esi(k,m)  = (soildepth_fp_before + slo_Dsi(k,m) *ddt_slo*dgammaa)/ddt_slo/dgammaa
								endif			
								if(isnan(slo_Esi(k,m)))then
									write(*,*) slo_Esi(k,m) ,slo_Dsi(k,m), k, V, ristar, dh, B_chan(k),usta,slo_grad(k),fgully(k,m),slo_sur_zb(k),zb_slo_idx(k)
									stop "E is nan" 
								endif
								if(slo_Esi(k,m) .gt.1.0E+2) then
									!if(ust_idx(k).gt.0.01) then
									write(*,'(a)') 'k      dsi       fp       w0      ssi      Esi       Dsi    v    I_chan    usta'
									write(*,*)k,  dsi(m), fp,w0(m),slo_ssi(k,m), slo_Esi(k,m), slo_Dsi(k,m), V, grad_H , usta
									!pause'check erosion rate'
								endif
								if(slo_Esi(k,m)<0.d0)then
									write(*,*) slo_Esi(k,m),soildepth_fp_before,k,m,fgully(k,m),v, grad_H ,alfa,slo_Dsi(k,m)
									stop 'E in gully is negative'
								endif
									!slo_Esi(k,m) = slo_Esi(k,m)*area_chan(k)/area !revise  20231226	turned off 20250407
							endif						
						enddo !end of do m = 1, Np	
					endif	
	!write(*,*) "check1205 07"						
				else ! without slope erosion 20230924
					 !V= Qg(k)/(width_idx(k)*dh)
					if(slo_grad(k).ge.slope_20 )then !added 20240131
					 V = 0.d0
					 usta = 0.d0
					 grad_H = 0.d0
					 h_dash = 0.d0
					else
					 !modified 20240126
					 if(len_slo_idx(1, k) > 0.d0 .and. len_slo_idx(2, k) > 0.d0 )then ! right-end side and down-end side boundaries of the basin do not have the width for x-flux and y-flux, respectively 
					  V = sqrt((qsur_ave_temp_idx(1,k) * area / len_slo_idx(1, k) ) ** 2.d0 + (qsur_ave_temp_idx(2,k) * area / len_slo_idx(2, k) ) ** 2.d0) / dh
					 elseif(len_slo_idx(1, k) == 0.d0 .and. len_slo_idx(2, k) > 0.d0 )then
					 V = abs(qsur_ave_temp_idx(2,k)) * area / (len_slo_idx(2, k) * dh)
					 elseif(len_slo_idx(1, k) > 0.d0 .and. len_slo_idx(2, k) == 0.d0 )then
					 V = abs(qsur_ave_temp_idx(1,k)) * area / (len_slo_idx(1, k) * dh)
					 else
					 V = 0.d0
					 endif
					 grad_H = (ns * v / dh ** (2./3.) ) ** 2.
					 usta = sqrt(grav * dh * grad_H) !modified 20231226 
					 h_dash = 0.d0
					!20240311
					!write(*,*) k,V,grad_H, dh,usta,qsur_ave_temp_idx(1,k),qsur_ave_temp_idx(2,k),len_slo_idx(1, k),len_slo_idx(2, k),ns
					!endif 		
	!write(*,*) "check1205 08"				 
						thetr=atan(grad_H) !added 20231226
						emb_gully = usta**2./((rs-1.)*csta/2.*grav*cos(thetr)*(tan(S_angle)-tan(thetr)))!Exchange layer thickness
						h_dash = emb_gully/dh
!modified 20231205						
						if (h_dash.ge.0.9d0) then !modified 20250407
						h_dash = 0.9d0
						emb_gully= h_dash*dh
						Usta = sqrt(emb_gully*((rs-1.)*csta/2.*grav*cos(thetr)*(tan(S_angle)-tan(thetr))))
						grad_H= usta**2./grav/dh
						V=1./ns*sqrt(grad_H)*dh**(2./3.)
						endif
	!write(*,*) "check1205 09"	
						!if(emb_gully.le.Emc) Usta=0.d0
					 endif !endif slope>20 degree
!added 20240101
						alfa=1.
						alfa = (slo_sur_zb(k)-zb_slo_idx(k))/emb_gully
						if(alfa.ge.1.) alfa=1.
						if(alfa <= 0.d0 .or.emb_gully.le.Emc) alfa=0. !turned off 20250407
						do m = 1, Np
							!slo_ssi(k,m) = slo_qsi(k,m)/Qg(k)  
							slo_qsi(k,m) = slo_ssi(k,m)*Qg(k)  
							slo_Dsi(k,m)  = slo_ssi(k,m) * w0(m) 
							!if(slo_Dsi(k,m).le.1e-7) slo_Dsi(k,m) = 0.d0 !added 20240101		
							if(slo_Dsi(k,m).le.0.d0) slo_Dsi(k,m) = 0.d0					
							slo_Esi(k,m) =0.d0
							if(alfa==0.) cycle	
							if(fmslo(k,m)==0.) cycle
			!				if(dzslo_idx(k).ge.soildepth_idx(k)) cycle
							fp = fmslo(k,m)
							if(dsi(m).ge.dh) cycle !added 20240131
							if (w0(m) .ge. usta) cycle
							if(isuseq == 1)then  !Lane-Kalinske formula
								bbeta = 6.0*w0(m)/karmans/usta
								Csb = slo_ssi(k,m)*bbeta/(1.-exp(-1. * bbeta))
								slo_Dsi(k,m)  = Csb * w0(m)
								!if(slo_Dsi(k,m).le.1e-7) slo_Dsi(k,m) = 0.d0 !added 20240101		
								if(slo_Dsi(k,m).le.0.d0) slo_Dsi(k,m) = 0.d0							
								z_star = usta/w0(m)
								Cae1 = 5.55*(0.5*z_star*exp(-(1./z_star)**2.))**1.61 * 1.0E-4 
								if(Cae1.ge.0.2d0) Cae1 = 0.2d0
								Cae = Cae1 *fp
								!if (riv(i,j)==1)then
								slo_Esi(k,m) = Cae * w0(m)*alfa !modified 20240101  	
								!else
								!slo_Esi(k,m) = Cae * w0(m)*0.5*alfa !modified 20240101 !modified erosion rate on slope area 20231111
								!endif						
							elseif(isuseq == 2)then  !Density stratified flow formula
								ristar = (rs/rw-1.)*csbar*grav*dh/(V**2.) 
								Ce = kappa/ristar*V*csbar/w0(m)
								if(Ce.ge.0.3d0)then !modified 20231225
									Ce = 0.3d0
									ristar = kappa/Ce*V*csbar/w0(m)        !Once ristar is modified by this regulation for smallest size, the regulated ristar is applied to all classes.
								endif
								slo_Esi(k,m) = fp*kappa/ristar*V*csbar*alfa !modified 20240101 								
							endif
							soildepth_fp_idx=(slo_sur_zb(k)-zb_slo_idx(k))*fmslo(k,m)
							!modified 20240101
							if(soildepth_fp_idx.le.0.) then
								slo_Esi(k,m)=0.d0
							else
								soildepth_fp_before = soildepth_fp_idx
								soildepth_fp_idx = soildepth_fp_idx + (slo_Dsi(k,m) -slo_Esi(k,m) )*ddt_slo*dgammaa
								if (soildepth_fp_idx.le.0.d0) then			
								slo_Esi(k,m)  = (soildepth_fp_before + slo_Dsi(k,m) *ddt_slo*dgammaa)/ddt_slo/dgammaa
								endif			
								if(isnan(slo_Esi(k,m)))then
								write(*,*) slo_Esi(k,m) ,slo_Dsi(k,m), k,V, ristar, dh,ns,len_slo_idx(1, k),len_slo_idx(2, k),qsur_ave_temp_idx(1,k),qsur_ave_temp_idx(2,k),usta,slo_grad(k),fmslo(k,m),slo_sur_zb(k),zb_slo_idx(k)
								stop "E is nan" 
								endif
								if(slo_Esi(k,m) .gt.1.0E+2) then
									!if(ust_idx(k).gt.0.01) then
									write(*,'(a)') 'k      dsi       fp       w0      ssi      Esi       Dsi    v    slope    usta'
									write(*,*)k,  dsi(m), fp,w0(m),slo_ssi(k,m), slo_Esi(k,m), slo_Dsi(k,m), V, slo_grad(k), usta,fmslo(k,m)
									!pause'check erosion rate'
								endif
								if(slo_Esi(k,m)<0.d0)then
								write(*,*) slo_Esi(k,m),soildepth_fp_before,k,m,alfa,v, slo_Dsi(k,m)
								stop 'E of slope is negative'
								endif
							endif
						enddo !end of do m = 1, Np						
	!write(*,*) "check1205 11"						
				endif ! for slope_ero_switch
			endif !endif for dh>surflow
		endif !endif for riv ==1	added 20250503
		enddo !end of do k = 1, slope_count	
	!write(*,*) "check1205 11"	
!pause'after deposition'
!pause '2'
!$omp parallel do private(m,n,riv_k,kk,kkk,down_xk,down_yk,i,j,ii,jj,ll,x_outflow,y_outflow,x_up_inflow ,y_up_inflow)
!--------------------trasportation of suspended sediment
		do k = 1, slo_count		
		    i=slo_idx2i(k)
			j=slo_idx2j(k)
		if( domain(i,j) == 0 ) cycle
		slo_sur_zb_before(k) = slo_sur_zb(k)
		if(isnan(water_v_cell(k)))then
		write(*,*) water_v_cell(k), hs_idx(k),da_idx(k)
		stop 'water volume in cell is nan'
		endif  
		if(dif_slo_idx(k)== 0)then
			do m = 1, Np
!modified 20240115			
			if(debris_switch.ne.0)then		
			!if(slo_grad(k).ge.0.25.and.zb(i,j)>100.d0 .and. riv(i,j).ne.1) exit !modified for gofukuya river 20231101;20240121	
			if(hki_g(k)>0.) exit	! no suspended sediment transportation taking place in debris flow grid cell; modified for gofukuya river 20240109	
			endif		
				do n = 1, 8
				if(up_slo_idx(k,n) == 0) exit 
					kkk = up_slo_idx(k, n)
					if(domain_riv_idx(kkk)==0)then !20250513
					slo_qsisum(k,m) = slo_qsisum(k,m) - slo_qsi(kkk,m)  !m3/s   !Á¨¶Âè∑„Å´Ê≥®ÅEΩ?!!!  ffd ÅEΩ?-„ÅÆ„Å®„ÅçÔøΩ??Á©çÔøΩ?+„ÅÆ„Å®„Åç‰æµÅEΩ? !modified 20230924
					endif
				end do		
				!if(riv(i,j)==1.and. qrs(i,j).ge.0.d0)then !modified 20250503
				if(riv(i,j)==1)then !all sediment in river cell should be delivered to unit channel 20250513
				slo_qsisum(k,m)=slo_qsisum(k,m)
				else
				slo_qsisum(k,m) = slo_qsisum(k,m) + slo_qsi(k,m)	 !modified 20230924
				endif
				if (isnan(slo_qsisum(k,m))) then
				write(*,*) k, m, riv(i,j), land(i,j), domain(i,j),slo_qsisum(k,m),slo_ssi(k,m),slo_qsi(k,m),qsur_ave_temp_idx(1,k),qsur_ave_temp_idx(2,k),water_v_cell(k)
				stop "check Nan of suspended sediment buguet in kinematic wave area"
			endif
			end do
	!write(*,*) "check1205 12"	

!!!!!modified for diffusion wave; by Qin 2023/4/26		
		elseif(dif_slo_idx(k)== 1)then
		kk=0
		kkk=0
			ii=slo_idx2i(k)-1
			jj=slo_idx2j(k)-1
			if(ii>0) kkk = slo_ij2idx(ii,j) !up:l=2 y direction inflow flux
			if(jj>0) kk = slo_ij2idx(i,jj) ! left: l=1 x direction inflow flux
!added 20240101			
			down_xk = down_slo_idx(1,k)
			down_yk = down_slo_idx(2,k)
!modified 20231206			
			do m =1,Np
			if(debris_switch.ne.0)then		
			!if(slo_grad(k).ge.0.25.and.zb(i,j)>100.d0 .and. riv(i,j).ne.1) exit !modified for gofukuya river 20231101;20240121	
			!if(slo_grad(k).ge.0.25) exit !20240126
			if(hki_g(k)>0.) exit	! no suspended sediment transportation taking place in debris flow grid cell; modified for gofukuya river 20240109	
			endif	
			x_up_inflow=0.d0
			y_up_inflow =0.d0
			x_outflow =0.d0
			y_outflow=0.d0
			!backward 20240126
			if(kk>0)then
			x_up_inflow = slo_ssi(kk,m) * qs_ave_idx(1,kk) !20250505
			endif
			if(kkk> 0)then
			y_up_inflow = slo_ssi(kkk,m) * qs_ave_idx(2,kkk) !20250505
			endif	
!added 20240101	
!upwind		
			if(qs_ave_idx(1,k) < 0.d0)then!20250505
			if(down_xk > 0)then
			x_outflow= slo_ssi(down_xk,m)*qs_ave_idx(1,down_xk) !20250505
			endif
			else		
			x_outflow= slo_ssi(k,m)*qs_ave_idx(1,k)		!20250505
			endif
			if(qs_ave_idx(2,k) < 0.d0)then!20250505
			if(down_yk > 0)then
			y_outflow= slo_ssi(down_yk,m)*qs_ave_idx(2, down_yk) !20250505
			endif
			else
			y_outflow= slo_ssi(k,m) * qs_ave_idx(2,k) !20250505
			endif			
			 slo_qsisum(k,m) = slo_qsisum(k,m)+ (x_outflow + y_outflow -x_up_inflow -y_up_inflow) *area !modified 20240101
			if (isnan(slo_qsisum(k,m))) then
				write(*,*) k, m, riv(i,j), land(i,j), x_outflow ,y_outflow ,x_up_inflow ,y_up_inflow,water_v_cell(k),domain(i,j),domain(i,jj),domain(ii,j),slo_qsisum(k,m),slo_ssi(k,m),qsur_ave_temp_idx(1,k),qsur_ave_temp_idx(2,k),slo_ssi(kk,m) ,qsur_ave_temp_idx(1,kk),slo_ssi(kkk,m) ,qsur_ave_temp_idx(1,kkk),qsur_ave_temp_idx(1,down_xk),slo_ssi(down_xk,m),slo_ssi(down_yk,m),qsur_ave_temp_idx(2,down_yk)
				stop "check Nan of suspended sediment buguet in diffusion wave area"
			endif	
			enddo
	!write(*,*) "check1205 13"				
		else
		stop "Please change the setting for kinematic approximation or diffusion approximation "
		endif
		! turnded off sediment inundation 20250513
		!if (riv(i,j)==1.and. qrs(i,j)<0.d0)then
		!	riv_k=riv_ij2idx(i,j) 
			!ll= link_to_riv(riv_k)			
		!	if(ll == 0) cycle     !in case ll is not defined 
		!	do m  =1, Np				  
		!		slo_qsisum(k,m) =slo_qsisum(k,m)+lin_to_slo_idx_di(riv_k,m)
		!	enddo
		!endif 	

		enddo		
!pause '3'
!write(*,*) "check1205 14"	
!$omp parallel do private(i,j,riv_k,m,soildepth_fp_before,gammaa_k,dgammaa,soildepth_fp_idx,over_ss,sum_fmslo,nn,d80)	  
  	do k = 1, slo_count
			i=slo_idx2i(k)
			j=slo_idx2j(k)
		if( domain(i,j) == 0 ) cycle
		riv_k=riv_ij2idx(i,j)
!modified 20240115			
			if(debris_switch.ne.0)then		
			!if(slo_grad(k).ge.0.25.and.zb(i,j)>100.d0 .and. riv(i,j).ne.1) cycle !modified for gofukuya river 20231101;20240121	
			!if(slo_grad(k).ge.0.25) cycle !20240126
			!if(hki_g(k)>0.) cycle	! no suspended sediment transportation taking place in debris flow grid cell; modified for gofukuya river 20240109	
			endif		
		ss_slope(k) = 0.d0
		qss_slope(k) = 0.d0
		gammaa_k=gammaa_idx(k)
		dgammaa=1./(1.-gammaa_k)
			do m = 1, Np
				soildepth_fp_before = (slo_sur_zb(k)-zb_slo_idx(k))*fmslo(k,m) 
				soildepth_fp_idx = soildepth_fp_before
			!	if(water_v_cell(k).le.0.d0.or.hki_g(k)>0.)then !modified for gofukuya 20240109
				if(riv(i,j)==1 .or. water_v_cell(k).le.0.d0)then
						if(slo_Dsi(k,m)>1.)then
						write(*,*) k, m, slo_Dsi(k,m),dsi(m), w0(m), slo_ssi(k,m), qrs(i,j), slo_qsisum(k,m)
					    stop 'slope deposition rate is too large check1'
						endif
					!slo_Dsi(k,m) =slo_Dsi(k,m)+ (-1.*slo_qsisum(k,m))/area !added 20231226	
					slo_Dsi(k,m) = (-1.*slo_qsisum(k,m))/area !20250511	
						slo_ssi(k,m) = 0.d0
					if(slo_Dsi(k,m).le.0.d0) slo_Dsi(k,m) = 0.d0
					slo_Esi(k,m) = 0.d0 !modified 20231226											
					slo_qsisum(k,m) = 0.d0
					slo_qsi(k,m) =0.d0
				else
					if(slo_Esi(k,m)<0.d0)then
					write(*,*) slo_Esi(k,m),soildepth_fp_before,k,m,fmslo(k,m)
					stop 'erosion rate is negative'
					endif
					slo_ssi(k,m) = slo_ssi(k,m) + (-1.* slo_qsisum(k,m) + (slo_Esi(k,m)- slo_Dsi(k,m))*area)*ddt_slo/water_v_cell(k)!modified for inundation 20230711
					if(slo_ssi(k,m).lt.0.d0) then
						slo_ssi(k,m) = slo_ssi(k,m) - (-1.* slo_qsisum(k,m) + (slo_Esi(k,m)- slo_Dsi(k,m))*area)*ddt_slo/water_v_cell(k)!modified for inundation 20230711
						slo_Dsi(k,m)  = -1.* slo_qsisum(k,m)/area + slo_Esi(k,m)      !D„Çí‰øÆÊ≠£„Åó„Å¶  !modified for inundation 20230711
						if(slo_Dsi(k,m)  < 0.) slo_Dsi(k,m)  = 0.d0 !„ÅæÅEΩ?D„ÅåÔøΩ?„ÅÆÊôÇÔøΩ?ÅEΩ‰∏ä‰∏ãÊµÅEøΩ?ÅEΩ„Ç¢„É≥„Éê„É©„É≥„Çπ„Å´„Çà„Çã„ÅÆ„ÅßD„ÅØ0„Å®„Å™ÅEΩ?
						!if(slo_Dsi(k,m).le.1e-7) slo_Dsi(k,m) = 0.d0 !added 20240101						
						slo_ssi(k,m) = slo_ssi(k,m) + (-1.* slo_qsisum(k,m) + (slo_Esi(k,m)- slo_Dsi(k,m))*area)*ddt_slo/water_v_cell(k) !modified for inundation 20230711
						if(slo_ssi(k,m).lt.0.d0) slo_ssi(k,m) = 0.d0 !modified 20240101	
					endif		

					if (isnan(slo_ssi(k,m))) then
						write(*,*) k, m, slo_ssi(k,m),water_v_cell(k),Qg(k),slo_qsisum(k,m),slo_Dsi(k,m),slo_Esi(k,m), area,h_chan(k) 
						stop "slope suspended sediment concentration is Nan"
					endif
				endif

				!if(riv(i,j)==1.and. qrs(i,j).ge.0.d0) then !added 20250503 all the sediment should supply into river if no overflow
				!slo_to_lin_idx_di(riv_k,m)=slo_to_lin_idx_di(riv_k,m)+slo_ssi(k,m)*water_v_cell(k)+slo_Dsi(k,m)*area*ddt_slo
				if(riv(i,j)==1) then !all sediment in river cell should be delivered to unit channel 20250513
				slo_to_lin_idx_di(riv_k,m)=slo_to_lin_idx_di(riv_k,m)+slo_Dsi(k,m)*area*ddt_slo
				slo_Dsi(k,m) =0.d0
				slo_ssi(k,m)=0.d0
				slo_qsi(k,m)=0.d0
				qss_slope(k)=0.d0
				else

				soildepth_fp_idx = soildepth_fp_idx + (slo_Dsi(k,m) -slo_Esi(k,m) )*ddt_slo*dgammaa
				!added 20231226
				if(soildepth_fp_idx .lt. 0.d0) then	
				slo_ssi(k,m) = slo_ssi(k,m) +(slo_Dsi(k,m) -slo_Esi(k,m) )*ddt_slo*area/water_v_cell(k) 	
				slo_Esi(k,m) = soildepth_fp_before / ddt_slo / dgammaa + slo_Dsi(k,m)	
				slo_ssi(k,m) = slo_ssi(k,m) +(slo_Esi(k,m) - slo_Dsi(k,m) )*ddt_slo*area/water_v_cell(k) 
				if(slo_ssi(k,m).le.0.d0) slo_ssi(k,m) = 0.d0	
				soildepth_fp_idx = 0.d0
				endif			
				dzslo_fp_idx(k,m) = (slo_Dsi(k,m) -slo_Esi(k,m) )*ddt_slo*dgammaa + dzslo_fp_idx(k,m)	
				ss_slope(k) = ss_slope(k) + slo_ssi(k,m)
				if (ss_slope(k)>thresh_ss) then
					!write(*,'(a,i)') "k=", k
					!write(*,'(a,f)') "ssi=", slo_ssi(k,m)
					!write(*,'(a,f)') "Esi=", slo_Esi(k,m)
					!write(*,'(a,f)') "Dsi=", slo_Dsi(k,m)
					!write(*,'(a,f)') "Qg=", Qg(k)
					!write(*,*) "slope suspended sediment concentration is larger than",thresh_ss 	
					!stop"slope suspended sediment concentration is larger than thresh_ss"	
					ss_slope(k)=ss_slope(k)-slo_ssi(k,m)
					over_ss= slo_ssi(k,m)-(thresh_ss-ss_slope(k))
					if(over_ss .le. 0.d0) over_ss = 0.d0
!modified 20231110	

					if(over_ss>0.d0) slo_Dsi(k,m) = slo_Dsi(k,m)+over_ss* water_v_cell(k)/area/ddt_slo		!20240101
						if(slo_Dsi(k,m)>1.)then
						write(*,*) k, m, slo_Dsi(k,m),dsi(m), w0(m), slo_ssi(k,m), slo_qsisum(k,m), over_ss, slo_Esi(k,m)
					   ! stop 'slope deposition rate is too large check3'
						endif
					if(slo_Dsi(k,m).le.0.d0) then
					 slo_Dsi(k,m) = 0.d0 !added 20240101
					 over_ss = 0.d0
					endif			
					slo_ssi(k,m)=thresh_ss-ss_slope(k)
					ss_slope(k) = thresh_ss
					if (slo_ssi(k,m).le.0.d0) slo_ssi(k,m)= 0.d0
!modified 20240101					
					dzslo_fp_idx(k,m) = dzslo_fp_idx(k,m)+over_ss* water_v_cell(k)/area*dgammaa	
				endif
				
				if(ss_slope(k).ge.1.d0)then 
					write(*,'(a,i)') "k=", k
					write(*,'(a,i)') "m=", m
					write(*,'(a,f)') "ss_slope=", ss_slope(k) 
					write(*,'(a,f)') "ssi=", slo_ssi(k,m)
					write(*,'(a,f)') "slo_qsisum=", slo_qsisum(k,m) 
					write(*,'(a,f)') "Esi=", slo_Esi(k,m)
					write(*,'(a,f)') "Dsi=", slo_Dsi(k,m)
					write(*,'(a,f)') "Qg=", Qg(k)
					write(*,'(a,f)') "water_vol=",water_v_cell(k)
					write(*,'(a,f)') "dzslo=",dzslo_fp_idx(k,m)
					stop "slope suspended sediment concentration is larger than 1"
				endif
				slo_qsi(k,m)  = slo_ssi(k,m) * Qg(k) !m3/s modified 20231226
				qss_slope(k) = qss_slope(k)+slo_qsi(k,m) 						

!added 20240115				
					!hs_idx(k)=hs_idx(k)+(slo_Dsi(k,m) -slo_Esi(k,m) )*ddt_slo*dgammaa								
			endif !endif for riv==1 sediment supply 20250503
				slo_sur_zb(k)= slo_sur_zb(k) + (slo_Dsi(k,m) -slo_Esi(k,m) )*ddt_slo*dgammaa				
				Ero_slo_vol(k) =  Ero_slo_vol(k)+(slo_Esi(k,m)-slo_Dsi(k,m)) *area*ddt_slo*dgammaa !minus value: deposition; modified 20230711 for inundation
				if(isnan(Ero_slo_vol(k)))then
					write(*,*) k, m, slo_Esi(k,m) ,slo_Dsi(k,m) ,slo_ssi(k,m),slo_qsisum(k,m),area_chan(k), q_chan(k),  water_v_cell(k),Qg(k),B_chan(k), l_chan(k), I_chan(k),dzslo_fp_idx(k,m)
					stop "erosion is NAN on the slope area"
				endif
			enddo !end of do m = 1, Np

	!write(*,*) "check1205 15"			 
!	  if(riv(i,j).ne.1.or. qrs(i,j)<0.d0) then !added 20250503
	  if(riv(i,j).ne.1) then !modified 20250513
		sum_fmslo =0.d0
		if(slo_sur_zb(k).gt.zb_slo_idx(k)) then	
			do m = 1,Np		
				soildepth_fp_idx=(slo_sur_zb_before(k) -zb_slo_idx(k))*fmslo(k,m)
				fmslo(k,m) = (soildepth_fp_idx+(slo_Dsi(k,m) -slo_Esi(k,m) )*ddt_slo*dgammaa)/(slo_sur_zb(k)-zb_slo_idx(k))
				if(fmslo(k,m).lt.0.d0) fmslo(k,m) = 0.d0
				sum_fmslo =sum_fmslo + fmslo(k,m)
				if(isnan(fmslo(k,m)).or.fmslo(k,m)<0.d0)then
					write(*,*) k,m, fmslo(k,m),fp,E,soildepth_fp_idx,(slo_sur_zb(k)-zb_slo_idx(k))*fmslo(k,m),soildepth_idx(k), V, ristar,usta
					stop "fm of slope is incorrect" 
				endif
			end do 	
		else
		do m = 1, Np
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
					d80 = d80+ fmslo(k,m)
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
	!consider vegetation cover effect;modified 20230710 for inundation sediment flow computation
			!ns_MS(k)=2.*ns_MS(k) !add the coeffient for the vegetation cover effect later 20230924
		endif!20230924	
		endif !20250503
	end do ! end of do k = 1, slo_count	

!pause '4'
	!write(*,*) "check1205 16"	
!$omp parallel do private(i,j)		
		  do k = 1,slo_count
		  	i = slo_idx2i(k)
			j = slo_idx2j(k)
			if( domain(i,j) == 0 ) cycle
		   soildepth_idx(k)=slo_sur_zb(k)-zb_slo_idx(k)
		   if(soildepth_idx(k).le.0.)then
		    soildepth_idx(k)=0.
		    da_idx(k)= 0.
			dm_idx(k) = 0.
		   	slo_sur_zb(k) = zb_slo_idx(k)
		   else
			da_idx(k)= (slo_sur_zb(k)-zb_slo_idx(k))*gammaa_idx(k)
			dm_idx(k) =  (slo_sur_zb(k)-zb_slo_idx(k)) * gammam_idx(k)
		   endif			
		   end do		   
	   
	!write(*,*) "check1205 18"	
!pause 'Finished dt_slo_sed time slope erosion calculation'
!modified 20231204	
8 continue	
		cal_time = cal_time +ddt_slo
	!write(*,*) "check1205 19"			
        if(cal_time.ge.slo_sedi_cal_duration) then
		!time = t*dt	
		!write (*,'(a,f,a,f,i,a)') 'sediment transport comuptation in slope area from time = ', time, 'to time = ', time - slo_sedi_cal_duration,  n_count,'th iteration'		
		!write(*,*) "check1205 20"			
		exit ! finish for this timestep
		endif
!pause 'finished dt time step slope erosion computation'		
enddo 
	end subroutine slo_sedi_cal !modified 20230920 

!Sediment flow exchange between slope and river;Added 20240305
	subroutine sed_exchange (sed_lin, hs_idx)
		use globals
		use sediment_mod
		implicit none
		type(sed_struct) sed_lin(link_count)
		real(8) hs_idx(slo_count)
		integer k, i, j, riv_k, ll, m
		real(8):: sedi
		inflow_sedi(:) = 0.d0
		overflow_sedi(:) =0.d0
		overflow_sedi_di(:,:) = 0.d0
		lin_to_slo_idx_di(:,:)=0.d0
		! Keep this loop serial: multiple river cells can map to the same unit-channel link.
		do k = 1, slo_count		
		    i=slo_idx2i(k)
			j=slo_idx2j(k)
			if( domain(i,j) == 0 ) cycle
			if(riv(i,j) == 0) cycle
			riv_k=riv_ij2idx(i,j) 
			ll= link_to_riv(riv_k)			
			if(ll == 0) cycle     !in case ll is not defined 
			if (isnan(qrs(i,j))) then
				write(*,*) k, m, riv(i,j), land(i,j), water_v_cell(k),domain(i,j),slo_qsisum(k,m),slo_ssi(k,m),qrs(slo_idx2i(k),slo_idx2j(k)),Qg(k),q_chan(k),slo_qsisum(k,m),slo_Dsi(k,m),slo_Esi(k,m), area_chan(k),area,B_chan(k),l_chan(k),h_chan(k),ns_MS(k),I_chan(k) 
				stop "river and cell flow exchange value is Nan"
			endif  
			do m = 1, Np
				if ( cut_overdepo_switch > 0 .and. overdepo_sedi_di(ll,m)>0.d0)then
					overflow_sedi_di(ll,m) = overdepo_sedi_di(ll,m) * width_lin(ll) * dis_riv_idx(riv_k)
					lin_to_slo_sum_di(ll,m) = lin_to_slo_sum_di(ll,m) + overflow_sedi_di(ll,m)
					lin_to_slo_sed_sum(ll)= lin_to_slo_sed_sum(ll) + overflow_sedi_di(ll,m)
					lin_to_slo_sed_sum_idx(riv_k) = lin_to_slo_sed_sum_idx(riv_k) + overflow_sedi_di(ll,m)
					overflow_sedi(riv_k) = overflow_sedi(riv_k) + overflow_sedi_di(ll,m) / dt
					lin_to_slo_idx_di(riv_k,m)=lin_to_slo_idx_di(riv_k,m)-overflow_sedi_di(ll,m) / dt
					overdepo_sedi_di(ll,m) = 0.d0 
				endif
				!sedi =slo_to_lin_idx_di(riv_k,m) !modified 20250503
				slo_to_lin_sed(ll,m) = slo_to_lin_sed(ll,m) + slo_to_lin_idx_di(riv_k,m)  !the sediment supply to unit channel
				slo_to_lin_sum_di(ll,m)=slo_to_lin_sum_di(ll,m)+slo_to_lin_idx_di(riv_k,m) 
				slo_to_lin_sed_sum(ll) = slo_to_lin_sed_sum(ll) + slo_to_lin_idx_di(riv_k,m) 
				slo_to_lin_sed_sum_idx(riv_k) = slo_to_lin_sed_sum_idx(riv_k) + slo_to_lin_idx_di(riv_k,m) 	
				inflow_sedi(riv_k)= inflow_sedi(riv_k) + slo_to_lin_idx_di(riv_k,m) / dt
				slo_to_lin_idx_di(riv_k,m)=0.d0	
			enddo	

			 if( qrs(i,j).lt.0.d0)then! sediment overflow (Only sediment loss from the unit channel due is considered; the inundation process in the slope area is ignored 20250513)
			  do m  =1, Np  
			  	lin_to_slo_idx_di(riv_k,m) =lin_to_slo_idx_di(riv_k,m)+ qrs(i,j)*area*sed_lin(ll)%ssi(m) !negative value
			  	sedi= lin_to_slo_idx_di(riv_k,m)*dt
				overflow_sedi_di(ll,m) = overflow_sedi_di(ll,m)- sedi
				lin_to_slo_sum_di(ll,m) = lin_to_slo_sum_di(ll,m)-sedi
				lin_to_slo_sed_sum(ll)= lin_to_slo_sed_sum(ll)-sedi
				lin_to_slo_sed_sum_idx(riv_k) = lin_to_slo_sed_sum_idx(riv_k)-sedi
				overflow_sedi(riv_k) = overflow_sedi(riv_k)- sedi / dt
			  enddo
			 endif 
	!write(*,*) "check1205 13.5"			 	
		enddo		
	end subroutine sed_exchange

!-----Debris flow  (added 2021/11/12)
    subroutine cal_Landslide ( hs_idx ) 
        use globals
        use sediment_mod
        implicit none	

		real(8), parameter :: PI = acos(-1.0d0)
		real(8), parameter :: rho = 1.d0   !=1.0d0(kg/m3)  density of water
		real(8) hs_idx(slo_count)
		real(8) :: Csta,sint,cost,tant,tanp
		real(8) :: h, dap, ha, hsp, D
		real(8) :: gravi, resi
		integer k, i

		csta = 1.d0 - lambda
		tanp = tan(PI/180.0d0*phi)      !=tan(phi)
		!n_LS = 0

!----Land slide (provided by Yamazaki and Egashira)
!$omp parallel do private(sint,cost,tant,dap,h,ha,D,hsp,gravi,resi)		
        do k = 1, slo_count
			!if(slo_grad(k).le.0.25.and.ns_slo_idx(k).lt.0.5)cycle !modified for gofukuya river 20231101;20240121	
			!if(slo_grad(k) .le. 0.25)cycle !modified for gofukuya river 20231101;20240127								
            if(riv(slo_idx2i(k),slo_idx2j(k)) == 1) cycle
!added 20240717		
			if(past_LS_switch>0)then	
			if(time<0.1)then
			if(past_ls(slo_idx2i(k),slo_idx2j(k)) == 1) LS_idx(k) = 3
			endif
			if(past_ls(slo_idx2i(k),slo_idx2j(k)) == 1) cycle 
			endif
!20240717			
			sf(k) = 0.d0
			if(LS_idx(k) > 0.99) cycle
			if(vol(k) > 0.1 ) cycle

			sint=sin(atan(slo_grad(k)))
			cost=cos(atan(slo_grad(k)))
			tant=slo_grad(k)
			
			dap = da_idx(k)
			ha = hs_idx(k)
			if (time < 0.1) ha = ha + L_rain_ini/1000.
			D = soildepth_idx(k)
			!----calculate water surface----
			if( ha > dap ) then
				h = ha - dap
				hsp = D
				pw(k) = pwc
			elseif( ha > pwc*D )then
				h = 0.d0
				hsp = (ha - pwc*D)/(lambda-pwc)
				pw(k) = pwc
			elseif( ha >= 0.d0 )then
				h = 0.d0
				hsp = 0.d0
				pw(k) = ha/D
			else
				h = 0.d0
				hsp = 0.d0
				pw(k) = 0.d0
			end if

			gravi=rho*grav*D*sint*(s/rho*(1.d0-lambda)+(1.d0-hsp/D)*pw(k)+hsp/D*lambda+h/D)
 !   		resi=rho*grav*D*cost*(s/rho*(1.d0-lambda)+(1.d0-hsp/D)*pw(k)-hsp/D*(1.d0-lambda))*tanp+c_dash(k)
     		resi=rho*grav*D*cost*(s/rho*(1.d0-lambda)+(1.d0-hsp/D)*pw(k)-hsp/D*(1.d0-lambda))*tanp+cohe
            sf(k) = gravi/resi !G/R		
			if(sf(k).ge.1.d0)then
			if(time<0.1)then
				LS_idx(k) = 2. ! the initial unstable meshes 20240630
			else	
				LS_idx(k) = 1.
			end if
			endif
		end do

	!$omp parallel do private(k) reduction(+ : LS_num) 	
	do k = 1, slo_count
		if(sf(k).ge.1.d0)then
			LS_num = LS_num + 1
		end if
	end do
	if (time < 0.1) then
		write(*,*) 'initial rainfall (mm) to eliminate unstable mesh = ', L_rain_ini
		write(*,*) 'Number_of_initial_unstable_mesh=', LS_num
		!pause
	else
	    write(*,*) 'Number_of_unstable_mesh=', LS_num
	end if
	end subroutine cal_Landslide

!----Debris flow_ mass point system (provided by Yamazaki and Egashira)
	subroutine cal_mspnt ( hs_idx )
        use globals
        use sediment_mod
        implicit none
		
		real(8), parameter :: PI = acos(-1.0d0)
		real(8), parameter :: deg2rad = PI/180.0d0
		real(8), parameter :: rad2deg = 180.0d0/PI
		real(8), parameter :: rhow = 1.d0   !=1.0d0(kg/m3)  density of water
		real(8) hs_idx(slo_count)
		real(8) :: Csta,sld,fld,sat,cco,cfo,ao,vo,c2,a,cc,cf,cstad,pcd,pfd,v
		real(8) :: h, dap, d_mp, ha, hsp, pw_mp, c_mp, rho, bsi, bbo, dz_mp !d_mp resotored 20250219 
		real(8) :: tant, cost, tanp, tane, dis
		real(8) :: delt, ldelt, alp, hkiarea
		integer k, l, m, i, j, kk, kkk, ii, jj, n_LS2

		csta = 1.d0 - lambda
		tanp = tan(PI/180.0d0*phi)      !=tan(phi)
		alp=0.01d0
		!allocate(k_LS(n_LS2))    !modified 20240724
		allocate(k_LS(n_LS))
		n_LS2 = 0
		!OMP must be prohibited for this loop 
		!!$omp parallel do private(k) reduction(+ : n_LS2)
!$omp single 
		do k = 1, slo_count
			if(sf(k)>=1.0d0)then
				n_LS2 = n_LS2 + 1
				k_LS(n_LS2) = k
			end if
		end do
		do l = 1, n_LS2
			i = slo_idx2i(k_LS(l))
            j = slo_idx2j(k_LS(l))
			k = k_LS(l)
			d_mp = soildepth_idx_deb(k)   !20250219
			if(d_mp == 0.d0) cycle
			dap = da_idx(k)
			ha = hs_idx(k)
			call cal_wl(ha, dap, d_mp, h, hsp, pw_mp)
			dis = dis_slo_1d_idx(k)
			sld=d_mp*(pc*csta)
			fld=d_mp*(pf*csta + (1.d0-csta)*hsp/d_mp + (1.d0-hsp/d_mp)*pw(k)) + h
			sat=csta          + (1.d0-csta)*hsp/d_mp + (1.d0-hsp/d_mp)*pw(k)
			cco=sld/(sld+fld)
			if(csta<=cco)then
!				write(*,'(a,10f10.3)')'csta<=cco,hs,d',cco,hsp,d_mp
				cycle
			end if
			cfo=d_mp*pf*csta/fld

			ao=(d_mp*sat + h ) * dis
			vo=ao*b_mp
			hkiarea = dis*b_mp
			hki_g(k) = hki_g(k) + 1.

			kk = k
			kkk = down_slo_1d_idx(k)

			do 
				i = slo_idx2i(kk)
				j = slo_idx2j(kk)
				ii = slo_idx2i(kkk)
				jj = slo_idx2j(kkk)
				d_mp=soildepth_idx_deb(kk)
				if(hki_g(kk)<=0.0001) d_mp = d_mp_ini   !caution !turned off at 20250122  !modified 20250219
				if( domain(i,j) == 0 ) exit
				if( d_mp == 0.d0 )exit
				if( riv(i,j) == 1 )then
!					write(*,*) 'kk,ao,vo',kk,ao,vo
					vo_total(slo_riv_idx(kk)) = vo    !„ÅìÔøΩ?ÅEΩk„ÅØÊñúÈù¢„ÅÆk„Åß„ÅÇ„Å£„Å¶riv„ÅÆk„Åß„ÅØ„Å™ÅEΩ?
					debris_total = debris_total + vo
					hki_area(slo_riv_idx(kk)) = hkiarea
					hki_total = hki_total + hkiarea
					hkiarea = 0.
!					write(*,*) 'cal_mspnt'
					exit 
				end if

				dis = dis_slo_1d_idx(kk)
				ao=vo/b_mp
				tant = (zb(i, j) - zb(ii, jj))/dis
!				tant = (zb(i, j)+dzslo_mspnt_idx(kk) - zb(ii, jj)-dzslo_mspnt_idx(kkk))/dis !trial
				tant = max(0.d0,tant)
				cost = cos(atan(tant))

				dap = da_idx(kk)
				ha = hs_idx(kk)
				call cal_wl(ha, dap, d_mp, h, hsp, pw_mp)
				c_mp = c_dash(kk)
				rho=(s-rhow)*cfo + rhow

				if(0.d0<d_mp)then
					sat=csta + (1.d0-csta)*hsp/d_mp + (1.d0-hsp/d_mp)*pw_mp
				else
					sat=1.d0
			    end if

				c2=0.d0 !c_mp/(1.d0*rho*g*cost*tanp)
				bsi=((s/rho-1.d0)*cco+c2)*tanp
				bbo= (s/rho-1.d0)*cco+1.d0
				tane=bsi/bbo

				delt=(tant-tane) / (1.d0 + tant*tane)
				ldelt=(ao/alp)**0.5d0 * delt

				if(tane <= tant)then !erosion
					!write(*,*) 'Erosion'
					if(hki_g(kk)<=0.0001)then     !20230319 to avoid double count
						dz_mp=d_mp !/cost
					else
						dz_mp=0.0001
						d_mp =0.0001
					endif
					if(d_mp*cost < ldelt)then
						ldelt=d_mp*cost
					end if
					if(hki_g(kk)<=0.0001)then
						a=ao + dis*ldelt*sat +dis*h   !20230319 to avoid double count
					else
						a=ao
					end if
					cc=(ao*cco + pc*csta*dis*ldelt) / a
					if(csta<cc)then
						cc=csta
					end if
					
					cf=(ao*(1.d0-cco)*cfo + pf*csta*dis*ldelt) / (a*(1.d0-cc))
					if(csta<cf)then
						cf=csta
					else if(cf<0.d0)then
						cf=0.d0
					end if

				else ! deposition
					!write(*,*) 'Deposition'
					dz_mp=zb(ii, jj) + dis*tane - zb(i, j)
!					dz_mp=zb(ii, jj) + dis*tane - zb(i, j) +dzslo_mspnt_idx(kkk)-dzslo_mspnt_idx(kk)
					if(dz_mp<0.d0)then
!					  write(*,*)'dz<0',zb(ii, jj),dis*tane,zb(i, j)
					 ! read(*,*)
					end if
					ldelt=max(ldelt,-dz_mp*cost,-ao/dis,-ao*cco/(csta*dis))
					a=ao + dis*ldelt
					if(a<=0.d0)then
					  !write(*,*)'a<0@ac',ao,a
					  a=0.d0
					  cc=0.d0
					  cf=0.d0
					  !read(*,*)
					  !ldelt=-ao/dx
					else				  
					  cstad=csta+(1.d0-csta)*cfo
					  pcd=csta/cstad
					  pfd=(cstad-csta)/cstad
	
					  cc=(ao*cco + pcd*dis*ldelt) / a
					  if(csta<cc)then
						!write(*,*)'csta<cc @ mp',cc
						cc=csta
					  else if(cc<0.d0)then
						!write(*,*)'cc<0.d0 @ mp',cc
						!write(*,*)'ao, cco, a' , ao, cco, a
						cc=0.d0
					  end if
				  
					  cf=(ao*(1.d0-cco)*cfo + pfd*dis*ldelt) / (a*(1.d0-cc))
					  if(csta<cf)then
						!write(*,*)'csta<cf @ mp',cf
						cf=csta
					  else if(cf<0.d0)then
						!write(*,*)'cf<0.d0 @ mp',cf
						!write(*,*)'ao, cco, a' , ao, cco, a
						cf=0.d0
					  end if
				  
					end if
				end if
				!write(*,*)'a, cc, cf' , a, cc, cf
				dz_mp=-ldelt/cost
				dzslo_mspnt_idx(kk) = dz_mp
				if (isnan(dzslo_mspnt_idx(kk))) then
					write(*,*) kk, dz_mp, ldelt, cost
				stop "dzslo_mspnt_idx(kk) is Nan"
			endif  
				!dzslo_mspnt_idx(kk) = -1.
				soildepth_idx_deb(kk) = soildepth_idx_deb(kk) + dz_mp
				soildepth_idx_deb(kk) = max(soildepth_idx_deb(kk),0.0001)
!				write(*,*)'kk, dz_mp, soil_depth' , kk, dzslo_mspnt_idx(kk),soildepth_idx(kk)
!20240724 tentative disable to prevent output error					
!				if(0.d0<=dz_mp)then
!					hs_idx(kk)=hs_idx(kk)+dz_mp
!				else
!					hs_idx(kk)=min(hs_idx(kk),soildepth_idx(kk))					
!				end if

				v=a*b_mp
				vol(kk)=vol(kk)+v
				vcc(kk)=vcc(kk)+v*cc
				vcf(kk)=vcf(kk)+v*(1.d0-cc)*cf
				if(hki_g(kk)<=0.0001)then              !added 20230317
					hkiarea = hkiarea + dis*b_mp
				end if
!				write(*,*)'k,kk,hki_g(kk),a,hkiarea,soil_depth,hs',k, kk, hki_g(kk), a, hkiarea,soildepth_idx(kk),hs_idx(kk)
				
				hki_g(kk) = hki_g(kk) + 1.
				if(a==0.d0)then
					write(*,'(a,5i5,5f10.3)')'exit:a=0',i,j,ii,jj,0,a,cc
					exit
				end if

				kk=kkk
				kkk=down_slo_1d_idx(kk)
				ao=a
				vo=v
				cco=cc
				cfo=cf
			end do
		end do
!$omp end single
		deallocate(k_LS)

	end subroutine cal_mspnt


	subroutine cal_wl(ha, dap, D, h, hsp, pw_mp)
		use globals
        use sediment_mod
		implicit none

		real(8) :: ha, dap, D, h, hsp, pw_mp

		!----calculate water surface----
		if( ha > dap ) then
			h = ha - dap
			hsp = D
			pw_mp = pwc
		elseif( ha > pwc*D )then
			h = 0.d0
			hsp = (ha - pwc*D)/(lambda-pwc)
			pw_mp = pwc
		elseif( ha >= 0.d0 )then
			h = 0.d0
			hsp = 0.d0
			pw_mp = ha/D
		else
			h = 0.d0
			hsp = 0.d0
			pw_mp = 0.d0
		end if
	
	end subroutine cal_wl


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


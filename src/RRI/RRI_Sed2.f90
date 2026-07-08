! RRI_Sediment
! The geometry parameters of river channel for sediment calculation
   subroutine riv_set4sedi
   use globals
   use sediment_mod
   use dam_mod
   implicit none
   integer i, j, ii, jj, k, kk, kkk,l,ll,n
   real(8) distance
   integer n_count,nn
	integer lll,kkkk, down_k ,m,f !20240601
	real(8) slope, thet, zbtemp !20240601

    !move1 ; sperate to RSR_chan_para 20240422
    !-----------------------added by yorozuya at 2016/08/03
    allocate( up_riv_idx(riv_count, 8))
    !-----------------------added by harada at 2020/11/17
    allocate( node_ups(riv_count) )
    !-----------------------added by Qin at 2021/4/28
    !allocate( damflg(riv_count))
    !-----------------------added by Qin at 2021/6/23
    allocate(riv_0th_idx(riv_count))
    allocate(divi_cell(riv_count))

    !move1 ;sperate to RSR_chan_para 20240422

    !move2; 20240422
    !-----------------------------------------------------------
    !!!!!Dividing unit channel 20230924 
	up_riv_idx(:, :) = 0
	divi_cell(:)=0
	node_ups(:) = 0 !20240603

    if(link_divi_switch==1)then
	do l = 1,sele_l_num
   ! i = sele_up_loci(l)
   ! j = sele_up_locj(l)
   ! ii=sele_down_loci(l)
   ! jj= sele_down_locj(l)
   ! modified for IRIC GUI
	i = ny-sele_loc(2,l)+1!modified 20250329
	j = sele_loc(1,l) 
	ii =ny-sele_loc(4,l)+1!modified 20250329)
	jj =sele_loc(3,l)
    k = riv_ij2idx(i,j)
	kk = riv_ij2idx(ii,jj)
	if(k.lt.1) then
        write (*, *) "check the i,j of section", l,"I", sele_loc(1,l),"J",sele_loc(2,l)
        stop
	endif	
	if(kk.lt.1) then
        write (*, *) "check the i,j of section", l,"I", sele_loc(3,l),"J",sele_loc(4,l)
        stop
	endif
    divi_cell(k) = 1
    divi_cell(kk) = 1
	!write(*,*)l, i,j,ii,jj, k,kk
    do
        do n = 1, merg_cell_num
        kk = down_riv_idx(k)
        if (kk == riv_ij2idx(ii,jj)) goto 1111
        k = kk
        enddo
        divi_cell(k) = 1
    enddo
    1111 continue
	enddo
    endif

    !-----------------------------------------------------------
	!2222 continue
    ! search for upstream gridcell (up_idx) added by yorozuya at 2016/08/03
    riv_count = 0
    link_count = 0
    do i = 1, ny
    do j = 1, nx

    n_count = 0

    if(domain(i,j).eq.0 .or. riv(i,j).ne.1) cycle
    ! domain(i, j) = 1 or 2 and riv(i, j) = 1
    riv_count = riv_count + 1

	riv_0th_idx(riv_count) = riv_0th(i,j) !added by Qin 2021/6/23
    !write(*,'(a,3i)') 'j/i/riv_count=', j, i, riv_count

    ! right
    ii = i
    jj = j + 1

    !write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
    !pause'right'

    if( dir(ii,jj).eq.16 ) then
    if(domain(ii,jj).eq.0 .or. riv(ii,jj).ne.1) goto 10
    n_count = n_count + 1 
    up_riv_idx(riv_count,n_count) = riv_ij2idx(ii, jj)
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
    if(domain(ii,jj).eq.0 .or. riv(ii,jj).ne.1) goto 20
    n_count = n_count + 1 
    up_riv_idx(riv_count,n_count) = riv_ij2idx(ii, jj)
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
    up_riv_idx(riv_count,n_count) = riv_ij2idx(ii, jj)
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
    up_riv_idx(riv_count,n_count) = riv_ij2idx(ii, jj)
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
    up_riv_idx(riv_count,n_count) = riv_ij2idx(ii, jj)
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
    up_riv_idx(riv_count,n_count) = riv_ij2idx(ii, jj)
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
    up_riv_idx(riv_count,n_count) = riv_ij2idx(ii, jj)
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
    up_riv_idx(riv_count,n_count) = riv_ij2idx(ii, jj)
    !write(*,'(a,5i)') 'i/j/ii/jj/dir=', j, i, jj, ii, dir(ii,jj)
    !pause'right up'
    endif

    80 continue

    !   ii = i
    !   jj = j
    !  if(dir(ii,jj).eq.0 ) then
    !   write(*,*) "dir(i, j) is error (", i, j, ")", dir(i, j)
    !   stop
    !  endif

    !write(*,*) riv_idx2j(riv_count), riv_idx2i(riv_count), n_count, (up_riv_idx(riv_count,nn), nn = 1, 8)
    !pause
    if(dam_switch == 1) then ! Added by Qin at 2021/4/21
    !kk = up_riv_idx(riv_count,1)
    !kk = riv_count-1
    if (damflg(riv_count).gt.0.and.n_count.eq.1) then ! Dam cell is set as a node
        !n_count = n_count + 1
        !up_riv_idx(riv_count, 1) = riv_ij2idx(ii, jj)
        link_count = link_count + 1
        node_ups(riv_count) = 1
        write(*,*) dam_switch, dam_num,damflg(riv_count),dam_iy(damflg(riv_count)), dam_ix(damflg(riv_count)),riv_count, n_count, link_count+1
    !elseif(damflg(riv_count)>0 .and.n_count>1)then
    !write(*,*) 'Dam location:', 'i=',dam_iy(damflg(riv_count)), 'j=',dam_ix(damflg(riv_count)), 'Dam cell is an confluence point'
    !stop "Change Dam loaciton"
    endif
    do f = 1, dam_num
    if(riv_count==down_riv_idx(dam_loc(f))) then
    !if(damflg(kk).gt.0.and. n_count.eq.1) then ! the downstream cell of the Dam is set as a node
        if(n_count.eq.1) then 
        !n_count = n_count+1
        !up_riv_idx(riv_count, 1) = riv_ij2idx(ii, jj)
        	link_count = link_count + 1
        	node_ups(riv_count) = 1			
    	endif
    !write(*,*) damflg(riv_count),riv_count, n_count, link_count
    endif
    enddo
    endif

	!for santanmria river old yamy dam; diversion channel added 20250329
	if(div_switch==2)then
	 do f = 1, div_id_max
		if(div_org_idx(f)==riv_count .and. n_count==1)then
		link_count = link_count + 1
		node_ups(riv_count) = 1
        write(*,*) div_switch, div_org_idx(f),riv_count, n_count, link_count
		elseif(down_riv_idx(div_org_idx(f))==riv_count .and. n_count==1)then
			link_count = link_count+1
			node_ups(riv_count) =1
		endif
	 enddo	
	endif

    !20230924 20231016
    if(link_divi_switch>0)then
    if(divi_cell(riv_count)==1.and.n_count<2.and.node_ups(riv_count)==0)then
        !if(n_count>1)then
        !riv_idx2i(riv_count) = i
        ! riv_idx2j(riv_count) = j
        !write(*,*) 'i=', i, 'j=', j, 'is a confluence point' 
        !stop "Change the location of the upstream and downstream cells for defining the division section"
        !else
        link_count = link_count + 1
        node_ups(riv_count) = 1
        !endif 
    endif 
    endif

    if(n_count .ge. 2) then
    link_count = link_count + n_count
    node_ups(riv_count) = n_count - 1   ! node_ups is the number of the inflow cells from the upstream to each target cell-1; node_ups = 0 when the cell is inside the link /  node_ups ge. 1 when the cell is the node
    end if
    enddo
    enddo

    link_count = link_count + 1    ! to add downstream end

    !-----allocate node cell  (added by Harada) 
    if(sed_switch==2)then
    write(*,*)'Number of unit channel =',link_count

    allocate( link_idx(riv_count, 8), link_idx_k(link_count), link_ups_k(riv_count) )
    allocate( Link_len(link_count), area_lin(link_count),link_0th_order(link_count), link_to_riv(riv_count) )
    allocate( up_riv_lin(link_count, 8), link_ij(ny,nx) )
    allocate( width_lin(link_count) )
    allocate( link_cell_num(link_count)) !---added by Qin 2021/08/18
    !allocate(chan_capa4wid_expan(link_count), depth_uni_chan(link_count)) !added 20240422
    allocate(wl4wid_expan(riv_count), width_lin_0(link_count),wid_expa_rate(link_count)) !added 20240422
	allocate(zb_riv_slope0_lin(link_count)) !20240601

    link_count = 0
    link_idx(:, :) = 0
    link_idx_k(:) = 0
	wid_expa_rate(:) =1.d0
    do k = 1, riv_count
        if(node_ups(k)==0.and.domain_riv_idx(k) .ne.2 ) cycle
        if(domain_riv_idx(k) .ne.2)then
        do n = 1, 8
            !if(up_riv_idx(k,n)==0) cycle
        if(up_riv_idx(k,n)==0) exit
            link_count = link_count + 1
            link_idx(k,n) = link_count
            link_idx_k(link_count) = up_riv_idx(k,n)  !downstream end cell in the link is set as link_idx
            !write(*,*) k, n, link_idx(k,n), link_idx_k(link_count)
        enddo
        elseif(domain_riv_idx(k) .eq.2)then!then downstream end unit channel for whole river basin 
        link_count = link_count + 1
        link_idx(k,1) = link_count
        link_idx_k(link_count) = k
        end if
        !write(*,*) k, link_count
		!max_width_added 20240229
    	if(width_idx(k)>max_width) max_width = width_idx(k)
    enddo

    !-----connect upstream and downstream node (added by Harada) 
    link_ups_k(:) = 0
    link_len(:) = 0.
    area_lin(:)=0.
    width_lin(:) =0.
    width_lin_0(:) = 0.
    link_0th_order(:) = 0
    link_to_riv(:) = 0
    link_cell_num(:)= 0  !--added by Qin 2021/8/18
    !chan_capa4wid_expan(:) = 0.
    !depth_uni_chan(:) = 0. 
	wl4wid_expan(:) = 0. !added 20240422
	
    do l = 1, link_count
        kkk = link_idx_k(l)
    !   width_dws = width_idx(kkk)
		
        do
        Link_len(l) = Link_len(l) + dis_riv_idx(kkk)
        width_lin(l) = width_lin(l) + width_idx(kkk)
		!depth_uni_chan(l) = depth_uni_chan(l) + depth_idx(kkk) + height_idx(kkk)
        link_cell_num(l) = link_cell_num(l)+1
        link_to_riv(kkk) = l
        if(node_ups(kkk) .ge. 1 ) then !if node_ups lager than 1, then target cell is a node, this cell should be the upstream start of this link (l)
            link_ups_k(l) = kkk
            link_0th_order(l) = 1
    !      width_ups = width_idx(kkk)
            !write(*,*) link_idx_k(l), link_ups_k(l), Link_len(l), link_0th_order(l)
            exit
        endif
        if(up_riv_idx(kkk,1) .eq. 0 ) then ! when node_ups =0, and up_rive_idx = 0, there is no link in upstream side of the target cell, this cell should be the upstream start of this link (l)
            link_ups_k(l) = kkk
            link_0th_order(l) = 0
    !      width_ups = width_idx(kkk)
            !write(*,*) link_idx_k(l), link_ups_k(l), Link_len(l), link_0th_order(l)
            exit
        end if  
        kkk = up_riv_idx(kkk,1) ! kkk is a link cell (node_ups = 0 and up_riv_idx .ne. 0), then go into the loop, to search the upstream end of this link
        end do
			
        !set river width (average) 
        width_lin(l) = width_lin(l)/real(link_cell_num(l)) !---modified by Qin 2021/8/18
        width_lin_0(l) = width_lin(l) ! the initial width of unit channels; added 20240422		
        area_lin(l) = width_lin(l)*link_len(l)
				
    end do

    !-----connect upstream and downstream link
    up_riv_lin(:,:) = 0
    ! link_to_riv(0) = 0
    do l = 1, link_count
        kk = link_ups_k(l)
        if(link_0th_order(l)== 0)then
            up_riv_lin(l, :) = 0
        else
        do n = 1, 8
            kkk = up_riv_idx(kk,n) 
            if(kkk==0) exit
            !if(kkk==0) cycle
            up_riv_lin(l, n) = link_to_riv(kkk)      !up_riv_lin=0 when no link is in upstream
            !write(*,*) up_riv_lin(l, n) , l, n, up_riv_idx(kk,n), link_ups_k(l)
        end do  
        endif
    end do
	
!---Eliminating the upstream end unit channel which steeper than the maximum slope threshold which given by user;added by Qin 20240601
	kill_n = 0	
	if(link_count>3)then

	 	do l = 1, link_count
			k = link_ups_k(l)
			kk = link_idx_k(l)
			kkk = down_riv_idx(kk)
			lll= link_to_riv(kkk)
			!modified the link slope by Qin 2021/5/25
			if(domain_riv_idx(kkk).eq.0 .or. damflg(kkk).gt.0) then
				m = 0
				zbtemp = 0.d0
				do n = 1, 8
					if (up_riv_lin(l,n) == 0) cycle
					m = m + 1
					ll = up_riv_lin(l,n)
					kkkk = link_idx_k(ll)
					zbtemp = zbtemp + zb_riv_idx(kkkk)
				enddo
				if(m.lt.1) then
				slope =	(zb_riv_idx(k) - zb_riv_idx(kk)) / link_len(l)
				else
                slope = (zbtemp/real(m)	- zb_riv_idx(kk)) / (link_len(l) + dis_riv_idx(kkkk))
				endif
			else
				slope = (zb_riv_idx(k)-zb_riv_idx(kkk))/(link_len(l) + dis_riv_idx(kkk))
			endif
			if(slope .le. min_slope) slope = min_slope     ! in case still the slope value is too small
			thet = atan(slope) * 180./3.14159
			zb_riv_slope0_lin(l) = thet ! degree	
		 if(zb_riv_slope0_lin(l).ge. max_slope .and. link_0th_order(l)==0)then
				kill_n = kill_n + 1
				n=0
				i = riv_idx2i(k)
				j = riv_idx2j(k)
				riv(i,j) = 0
				down_k = down_riv_idx(k)
				do
				n=n+1
				if(link_to_riv(down_k) .ne. l) exit 
				i = riv_idx2i(down_k)
				j = riv_idx2j(down_k)
				riv(i,j) = 0
				width(i,j) = 0.
				depth(i,j) = 0.
				height(i,j) = 0.
				len_riv(i,j) = 0.
				riv_0th(i,j) = 0
				area_ratio(i,j) = 0.
				down_k  = down_riv_idx(down_k)
				enddo	
				if(n.ne.link_cell_num(l))then
				write(*,*) 'l=',l,'killed cell number=', n, 'number of river cell in the eliminated unit channel= ',link_cell_num(l) 
				stop 'unit channel elimination is failed'
				endif
		 endif	
		end do
	endif	

	write(*,*) kill_n," of upstream end unit channels"," have been eliminated" 
	if(kill_n>0 .and. link_count>3)then
		deallocate( up_riv_idx)
		deallocate( node_ups)
		deallocate(riv_0th_idx)
		deallocate(divi_cell)
		deallocate( link_idx, link_ups_k)
		deallocate(link_idx_k, stat = k)
		deallocate( Link_len, area_lin,link_0th_order, link_to_riv)
		deallocate( up_riv_lin, link_ij)
		deallocate( width_lin)
		deallocate( link_cell_num) !---added by Qin 2021/08/18
		deallocate(wl4wid_expan, width_lin_0,wid_expa_rate) !added 20240422
		deallocate(zb_riv_slope0_lin) !20240601
!---
		deallocate( riv_idx2i, riv_idx2j )
		deallocate( down_riv_idx, domain_riv_idx )
		deallocate( width_idx, depth_idx )
		deallocate( height_idx, area_ratio_idx )
		deallocate( zb_riv_idx, dis_riv_idx )
		deallocate( dif_riv_idx )
		deallocate( sec_map_idx, len_riv_idx ) ! add v1.4
		deallocate( damflg) !off for dam 
		write(*,*) "finish deallocating"
	else
		open(1234, file = 'upcurdown.txt', status = 'unknown' )
		write(1234,'(a)') '     n    2i,    2j,    up,   cur,  down,  link'
		do k = 1, riv_count
		do n = 1, 8
		!write(1234,'(8i6)') n, riv_idx2i(k), riv_idx2j(k),  up_riv_idx(k,n), k, down_riv_idx(k), link_to_riv(k)
		write(1234,'(8i6)') n, riv_idx2j(k), ny+1-riv_idx2i(k), up_riv_idx(k,n), k, down_riv_idx(k), link_to_riv(k) !convert to i,j in IRIC GUI 20250505
		enddo
		enddo
		close(1234)

		link_ij(:,:) = 0 
		open( 2456, file = 'linkascii.txt', status = 'unknown' )
		do i = 1, ny
			do j = 1, nx
			if (riv(i,j).eq.1) then 
				k = riv_ij2idx(i,j)
				link_ij(i,j) = link_to_riv(k)
			else
				link_ij(i,j) = -9999
			endif 
		enddo
		write(2456,'(1000000i6)') (link_ij(i,j), j =1, nx)
		enddo   
		close(2456)	
		kill_n = -99	
	endif

    endif
    !move2 ; 20240422
    end subroutine riv_set4sedi
! calculation of sediment

!------------------------------------------------------
! slope to river connction (this subroutine is only for sed_link model)(moved from RRI_sub.f90 20240724)
subroutine slo_to_riv
  use globals
  use sediment_mod !added 20240422
  implicit none
  real(8) distance, slope1, slope2, diagonal
  integer k, kk, kkk, i, j, ii, jj, l,m,n
  integer steep_dirc, slo_xdire, slo_ydire
  allocate( slo_grad(slo_count) )
  allocate(slograd(ny,nx))
  allocate(l_chan(slo_count),slo_ang_sin(slo_count),slo_ang_cos(slo_count))

    l_chan(:) = 0.
    diagonal = sqrt(dx**2.+ dy**2.)
  !This loop was separated (modified 2022/3/6)
  do k = 1, slo_count 
    !Steepest gradient (slope) of the cell (Steepest direction) (4-direction) 
    slope1 = -9999.    !steepest
    slope2 = -9999.    !2nd steepest

    i = slo_idx2i(k)
    j = slo_idx2j(k)

    do l = 1, 4    !steepest search 
      if( l.eq.1 ) then
        ii = i
        jj = j + 1
        distance = dx
      elseif( l.eq.2 ) then
        ii = i
        jj = j - 1
        distance = dx
      elseif( l.eq.3 ) then
        ii = i + 1
        jj = j
        distance = dy
      elseif( l.eq.4 )then
        ii = i - 1
        jj = j
        distance = dy
      end if

      if( ii .gt. ny ) cycle
      if( jj .gt. nx ) cycle
      if( ii .lt. 1 ) cycle
      if( jj .lt. 1 ) cycle
      if( domain(ii,jj).eq.0 ) cycle

      if( (zb(i, j) - zb(ii, jj))/distance .ge. slope1 ) then
        slope1 = (zb(i, j) - zb(ii, jj))/distance
        if(l==1.or.l==2) then
        steep_dirc = 1     !x-direction is steepest 
        if(jj.gt.j)then
        slo_xdire = 1         ! positive x-direction: right
        else
        slo_xdire = 2         ! negative x-direction : left
        endif 
        else
        steep_dirc = 2    !y-direction is steepest
        if(ii.gt.i)then
        slo_ydire = 1         ! positive y-direction : down
        else
        slo_ydire = 2         ! negative y-direction: up
        endif 
        endif
      !  if(l==3.or.l==4) steep_dirc = 2     !y-direction is steepest
      end if
    end do

    do l = 1, 2    !2nd steepest search
      if(steep_dirc==1)then   !x-direction is steepest/ search y-direction only
        if( l.eq.1 ) then
          ii = i + 1
          jj = j
        elseif( l.eq.2 )then
          ii = i - 1
          jj = j
        end if
        distance = dy
      else                    !steep_dirc = 2: y-direction is steepest/ search x-direction only
        if( l.eq.1 ) then
          ii = i
          jj = j + 1
        elseif( l.eq.2 ) then
          ii = i
          jj = j - 1
        end if
        distance = dx
      endif

      if( ii .gt. ny ) cycle
      if( jj .gt. nx ) cycle
      if( ii .lt. 1 ) cycle
      if( jj .lt. 1 ) cycle
      if( domain(ii,jj).eq.0 ) cycle

      if( (zb(i, j) - zb(ii, jj))/distance .ge. slope2 ) then
        slope2 = (zb(i, j) - zb(ii, jj))/distance
        if(steep_dirc ==1)then
        if(l==1)then
        slo_ydire =1
        elseif(l==2)then
        slo_ydire = 2
        endif
        else
        if(l==1)then
        slo_xdire =1
        elseif(l==2)then
        slo_xdire = 2
        endif
        endif
      end if
    end do

!modified by Qin 20220803 revised by Qin20240226
    if(slope2 .le.0.d0) then
     if(slope1.le.0.d0)then
      slo_grad(k) = 0.00001
     else
      slo_grad(k) = slope1
     endif
     if(steep_dirc==1)then
      slo_ang_sin(k) = 0.      
      if(slo_xdire==1) then
       slo_ang_cos(k) = 1.
      elseif(slo_xdire==2) then
       slo_ang_cos(k) = -1.
      endif
     elseif(steep_dirc==2)then
      !slo_dy(k) = sqrt((dy*slo_grad(k))**2.+dy**2.)
      !slo_dx(k) = dx
      slo_ang_cos(k) = 0.  
      if(slo_ydire==1) then
       slo_ang_sin(k) = 1.
      elseif(slo_ydire==2) then
       slo_ang_sin(k) = -1.
      endif
     endif
    else
     slo_grad(k) = sqrt(slope1**2. + slope2**2.)   !Steepest gradient (slope) of the cell
      if(steep_dirc==1)then    
       if(slo_xdire==1) then      
        slo_ang_cos(k) = slope1/slo_grad(k)   
       elseif(slo_xdire==2)then
        slo_ang_cos(k) = -slope1/slo_grad(k)    
       endif    
       if(slo_ydire==1)then
          slo_ang_sin(k) = slope2/slo_grad(k)
        elseif(slo_ydire==2)then
         slo_ang_sin(k) = - slope2/slo_grad(k)
       endif
       
      elseif(steep_dirc==2)then     
       if(slo_xdire==1) then      
        slo_ang_cos(k) = slope2/slo_grad(k)
       elseif(slo_xdire==2)then
        slo_ang_cos(k) = -slope2/slo_grad(k)
       endif
       if(slo_ydire==1)then
         slo_ang_sin(k) = slope1/slo_grad(k)
       elseif(slo_ydire==2)then
         slo_ang_sin(k) = -slope1/slo_grad(k)
       endif
      endif
    endif 

    if(slo_grad(k) .gt. 9998.) slo_grad(k) = 0.00001
    slograd(i,j) = slo_grad(k)
    !modified for gofukuya-river 20240121
   ! if(zb(i,j) .le. 50.d0) ns_slo_idx(k) = 0.05 !20240126
    !if(riv(i,j)==1) then
    !   ns_slo_idx(k) = 0.05
    !else   
     ! if (slo_grad(k).lt.0.25.and.domain(i,j).ne.0) then
     !   if(zb(i,j)>100.d0)then
     !    ns_slo_idx(k) = ns_slo_idx(k)
       ! else  
       !   ns_slo_idx(k) = 0.05
       ! endif
      !endif
    !endif
    ! the length of sediment path(gully) in each grid cell
    l_chan(k) = dx/abs(slo_ang_cos(k))
    if (l_chan(k).gt.diagonal) l_chan(k) = dy/abs(slo_ang_sin(k))
    !modified 20240220
    !if (l_chan(k).le.diagonal)then
    !l_chan_x(k) = dx * abs(slo_ang_cos(k))/(slo_ang_cos(k))
    !l_chan_y (k) = l_chan(k) * slo_ang_sin(k)
    !else
    !l_chan(k) = dy/abs(slo_ang_sin(k))
    !l_chan_x(k) = l_chan(k) * slo_ang_cos(k)
    !l_chan_y(k) = dy * abs(slo_ang_sin(k))/(slo_ang_sin(k))
    !endif
  end do

  open(5888, file = 'gradient.out' )
    do i = 1, ny
        write(5888,'(10000000f14.5)') (slograd(i, j), j = 1, nx)
    enddo    
  close(5888)
end subroutine slo_to_riv

!-----sediment size distribution setting 
subroutine sed_distr2 (sed)
	use globals
	use sediment_mod	    
	use RRI_iric
    use iric
	implicit none

	type (sed_struct) :: sed(ny,nx)

	character*256 fi
	integer :: i,j,m,n,ierror, isnumber,nn,mm,ier,k,kk,nnn
	INTEGER :: tmpint, j_mix_dis_dep, fid, np2, np3,j_set_slo_GSD!modified for slope erosion
	real(8) :: di_temp, fm_temp, fm1_temp, ft_temp, fms_temp, fms2_temp
	real(8) :: fm_all, ft_all, fm1_all, fms_all,f90
	real(8) :: A1, AA
	!character(50) :: mix_label, cm

    dsi(:) =0.d0
    w0(:) =0.d0

	nm_cell = 9  !tentative
	np = 0 
!---modified for slope erosion
    allocate ( sdt_c(ny,nx) )
	allocate (sdt_s(ny,nx))
	sdt_c(:,:) = 1
	sdt_s(:,:) = 1
    call iric_read_cell_attr_int('GSD_c', sdt_c)
	call iric_read_cell_attr_int('GSD_s', sdt_s)
!---	
	call cg_iric_read_integer(cgns_f, "j_mix_dis_dep", j_mix_dis_dep, ier)
	
	CALL CG_IRIC_READ_FUNCTIONALSIZE(cgns_f, 'mixfile_fr',tmpint,ier)
          IF(ier == 0) THEN
             np = tmpint
             mm = np
          ENDIF
			if(np.le.0)then
				write(*,*) 'Please set sediment size distribution for mixed layer'
				stop
			end if
          allocate(xtmp(tmpint),ytmp(tmpint))
          !CALL CG_IRIC_READ_FUNCTIONAL(cgns_f, 'mixfile_fr',xtmp,ytmp,ier)

		  write(*,*)'number of particles:',np
		  allocate(ddist_mm(np))
		  allocate( pdist_m_100(np,0:nm_cell), pmk0(np,0:nm_cell) )

          call cg_iric_read_functionalwithname(cgns_f, 'mixfile_fr', 'diameter_k', xtmp, ier)

          do k=1,np
             ddist_mm(k) = xtmp(k)
          end do
 
          do n=0,nm_cell
             write(cm,'(i1)') n
             mix_label = 'Region'//trim(cm)          
             call cg_iric_read_functionalwithname(cgns_f, 'mixfile_fr', mix_label, ytmp, ier)          
             do k=1,np
                pdist_m_100(k,n) = ytmp(k)
             end do                    
          end do

          DEALLOCATE(xtmp, STAT = ier)
          DEALLOCATE(ytmp, STAT = ier)

          if( j_mix_dis_dep==1 ) then
             
             CALL CG_IRIC_READ_FUNCTIONALSIZE(cgns_f, 'mixfile_fr_d',tmpint,ier)
             allocate(xtmp(tmpint),ytmp(tmpint))
			 allocate( pdist_d_100(np,0:nm_cell), pdk0(np,0:nm_cell) )
			 IF(ier == 0) THEN
             	np2 = tmpint
          	 ENDIF
			if( np2==0)then
				write(*,*) 'Please set sediment size distribution for the depostied layer'
				stop
			end if
             if( np2/=np ) then
                write(*,*) "The sediment size class in deposited layer is different from one in the mixed layer."
                write(*,*) "The sediment size class (number of class and each diameter) must be same in both layers."
                !cg_iric_close(cgns_f, ier)
                stop
             end if
       
             do n=0,nm_cell
                write(cm,'(i1)') n
                mix_label = 'Region'//trim(cm)          
                call cg_iric_read_functionalwithname(cgns_f, 'mixfile_fr_d', mix_label, ytmp, ier)          
                do k=1,np
                   pdist_d_100(k,n) = ytmp(k)
                end do                    
             end do
       
             DEALLOCATE(xtmp, STAT = ier)
             DEALLOCATE(ytmp, STAT = ier)
             
          end if

		 ! if(debris_switch==1.or.slo_sedi_cal_switch==1)then !modified for slope erosion
		if(debris_switch==1.or.slo_sedi_cal_switch==1)then   
		 call cg_iric_read_integer(cgns_f, "j_set_slo_GSD", j_set_slo_GSD, ier)
		 if(j_set_slo_GSD==0)then 
			write(*,*) 'Please set sediment size distribution for the slope area'
			stop
		 else
             CALL CG_IRIC_READ_FUNCTIONALSIZE(cgns_f, 'mixfile_fr_slopearea',tmpint,ier)
             allocate(xtmp(tmpint),ytmp(tmpint))
			 allocate( pdist_slosoil_100(np,0:nm_cell) )

			 IF(ier == 0) THEN
             	np3 = tmpint
          	 ENDIF

			if( np3==0)then
				write(*,*) 'Please set sediment size distribution for the slope area'
				stop
			end if
             if( np3/=np ) then
                write(*,*) "The sediment size for the slope area is not set or the class is different from one in the mixed layer."
                write(*,*) "The sediment size class (number of class and each diameter) must be same with river sediment."
                !cg_iric_close(cgns_f, ier)
                stop
             end if
       
             do n=0,nm_cell
                write(cm,'(i1)') n
                mix_label = 'Region'//trim(cm)          
                call cg_iric_read_functionalwithname(cgns_f, 'mixfile_fr_slopearea', mix_label, ytmp, ier)          
                do k=1,np
                   pdist_slosoil_100(k,n) = ytmp(k)
                end do                    
             end do

			 DEALLOCATE(xtmp, STAT = ier)
             DEALLOCATE(ytmp, STAT = ier)

		  end if
		endif
    allocate (dsi(Np),w0(Np))!20220729
    allocate (ns_MS(slo_count))
	ns_MS(:) = 0.02
	nnn = 0

	do i = 1,ny
	  do j = 1,nx
		fm_all = 0.d0
		ft_all = 0.d0
		fm1_all = 0.d0
		fms_all = 0.d0
		sed(i,j)%dmean = 0.d0
		sed(i,j)%dmean_slo = 0.d0
		sed(i,j)%Bed_D60 = 0.d0 !added 20240422
		sed(i,j)%Bed_D90 =0.d0
		sed(i,j)%samp_dep = samp_dep_min !added 20240419
		sed(i,j)%stor_sed_sum = 0.d0 !20240424
		do m = 1, Np
			sed(i,j)%dsed(m) = ddist_mm(m)/1000.
			sed(i,j)%fm(m) = pdist_m_100(m,sdt_c(i,j)-1)
			if( j_mix_dis_dep==0 ) then
				sed(i,j)%fm1(m) = pdist_m_100(m,sdt_c(i,j)-1)
				sed(i,j)%ft(m) = pdist_m_100(m,sdt_c(i,j)-1)
			else
				sed(i,j)%fm1(m) = pdist_d_100(m,sdt_c(i,j)-1)
				sed(i,j)%ft(m) = pdist_d_100(m,sdt_c(i,j)-1)
			end if
			if(debris_switch==1.or.slo_sedi_cal_switch==1)then   !20240806
				sed(i,j)%fms(m) = pdist_slosoil_100(m,sdt_s(i,j)-1)
			else
				sed(i,j)%fms(m) = pdist_m_100(m,sdt_c(i,j)-1)
			end if
			sed(i,j)%qbi(m) = 0.d0
			sed(i,j)%qsi(m) = 0.d0
			sed(i,j)%qswi(m) = 0.d0
			sed(i,j)%qdi(m) = 0.d0
			sed(i,j)%ffd(m) = 0.d0
			sed(i,j)%qdsum(m) = 0.d0
		! for bedload, suspended load; GSD added by Qin 2021/6/25
			sed(i,j)%fqbi(m) = 0.d0
			sed(i,j)%fqsi(m) = 0.d0
			!sed(i,j)%fsur(m) =0.d0 !2021/11/12	
			sed(i,j)%Esisum(m) = 0.d0 !2021/11/26
			sed(i,j)%Dsisum(m) = 0.d0		
			sed(i,j)%stor_sed_sum_i(m) = 0.d0 ! added 20240424   						 
			fm_all = fm_all + sed(i,j)%fm(m) !ここでは%表示。例えば10
			fm1_all = fm1_all + sed(i,j)%fm1(m)
			ft_all = ft_all + sed(i,j)%ft(m)
			fms_all = fms_all + sed(i,j)%fms(m)
		end do
		sed(i,j)%Bed_D60 = sed(i,j)%dsed(Np) !added 20240422
		sed(i,j)%Bed_D90=sed(i,j)%dsed(Np)
		f90=0.d0
		nn = 0
		do m = 1, Np			
					   sed(i,j)%fm(m) = sed(i,j)%fm(m)/fm_all !ここでは%/100表示。例えば0.1
					   sed(i,j)%ft(m) = sed(i,j)%ft(m)/ft_all
					   sed(i,j)%fm1(m) = sed(i,j)%fm1(m)/fm1_all
					   sed(i,j)%fms(m) = sed(i,j)%fms(m)/fms_all
					    if(sed(i,j)%fms(m)>1.d0) then
!	                    write(*,*) 'input fmslo >1 ', i,j, m, sed(i,j)%dsed(m), fms_all,sed(i,j)%fms(m), sed(i,j)%fms(1),sed(i,j)%fms(2),sed(i,j)%fms(10)!modified 20240509
	                    write(*,*) 'input fmslo >1 ', j, ny+1-i, m, sed(i,j)%dsed(m), fms_all,sed(i,j)%fms(m), sed(i,j)%fms(1),sed(i,j)%fms(2),sed(i,j)%fms(10)!convert i,j in IRIC GUI 20250505
		                stop
	                    endif
						sed(i,j)%fsur(m)=sed(i,j)%fm(m)  !20240509
					   sed(i,j)%dmean = sed(i,j)%dsed(m) * sed(i,j)%fm(m) + sed(i,j)%dmean
					   sed(i,j)%dmean_slo = sed(i,j)%dsed(m) * sed(i,j)%fms(m) + sed(i,j)%dmean_slo
					!define D90 for the bed surface GSD estimation; added by Qin 20230421	
					   if (nn.ge.0)then !modified 20240510
						nn=nn+1
						f90= f90 + sed(i,j)%fm(m)
						!added 20240422
						if(nn>1) then
							if(f90>0.6 .and. (f90-sed(i,j)%fm(m)) .le. 0.6 )sed(i,j)%Bed_D60= sed(i,j)%dsed(m-1)+ (0.6d0-f90+sed(i,j)%fm(m))*(sed(i,j)%dsed(m)-sed(i,j)%dsed(m-1))/sed(i,j)%fm(m)
						else
							if(f90 > 0.6) sed(i,j)%Bed_D60 = sed(i,j)%dsed(m)
						endif	
						if(f90>0.9d0) then
						if(nn.le.1)then
						sed(i,j)%Bed_D90= sed(i,j)%dsed(m)
						else
						sed(i,j)%Bed_D90= sed(i,j)%dsed(m-1)+ (0.9d0-f90+sed(i,j)%fm(m))*(sed(i,j)%dsed(m)-sed(i,j)%dsed(m-1))/sed(i,j)%fm(m)
						endif
						nn=-9999
						f90=-9999
						endif	
					   endif				  		 
						do n = 1, Nl       !modified from insumber to Nl 20240724
					  		sed(i,j)%fd(m,n) = sed(i,j)%fm(m)
						enddo
		
					!	do n = isnumber + 1, Nl
					!  		sed(i,j)%fd(m,n) = sed(i,j)%fm1(m)
					!	enddo
		end do		
		
		! For uniformed sediment
				if(sed_type_switch .eq. 1) then
					 sed(i,j)%dmean = ds_river
					 sed(i,j)%dmean_slo = ds_river   !tentative
					 !sed(i,j)%dzbpr(m) = 0.d0
						
					 do m = 1, Np
						  sed(i,j)%fm(m) = 1.d0
						 sed(i,j)%ft(m) = 1.d0
						 sed(i,j)%fms(m) = 1.d0
						 
						 sed(i,j)%qbi(m) = 0.d0
						 sed(i,j)%qsi(m) = 0.d0
						 sed(i,j)%qswi(m) = 0.d0
					 end do			 	 
				end if
	  !-----------------Estimating the roughness of slope area based on Manning-Strickler; added by Qin 2022/7/29 
                        kk = slo_ij2idx(i, j) !modified 20240509
                        ns_MS(kk) = (2.d0*sed(i,j)%dmean)**(1./6.)/(7.66*grav**0.5)
                        if(ns_MS(kk)<0.02)then 
                            ns_MS(kk)=0.02 
                        elseif(ns_Ms(kk)>0.03)then
                            ns_MS(kk)=0.03
                        endif
                        !consider vegetation cover effect;modified 20230710 for inundation sediment flow computation
                        ns_MS(kk)=2.*ns_MS(kk)
	!---fall velocity of each sediment particle;added by Qin 20220802
                if(nnn==0)then
                	if(riv(i,j)==1)then
                 		nnn=nnn+1
			    		do m = 1, Np
                 			dsi(m) = sed(i,j)%dsed(m)
				 			A1 = 6.*kin_visc/dsi(m)
				 			AA = (2./3.)*(s-1.)*grav*dsi(m)
				 			w0(m) = -1.0*A1+sqrt(A1**2.+AA)	
			    		enddo 
                	endif
                endif
	  end do
	end do

	print*,"done: reading sediment distribution files"
deallocate(sdt_c)

end subroutine sed_distr2


!-----Slope sediment distribution setting
	  subroutine sed_distr ( i,j, sed, fi)
		use globals
		use sediment_mod
		implicit none
			  
		type (sed_struct) :: sed(ny,nx)
		!allocate (Bed_D90(ny,nx))
			  
		character*256 fi
		integer :: i,j,m,n,ierror, isnumber,nn
		real(8) :: di_temp, fm_temp, fm1_temp, ft_temp, fms_temp, fms2_temp
		real(8) :: fm_all, ft_all, fm1_all, fms_all,f90
		!write(*,'(2i,a)') i, j, fi
		!pause
		
		fm_all = 0.d0
		ft_all = 0.d0
		fm1_all = 0.d0
		fms_all = 0.d0
		sed(i,j)%dmean = 0.d0
		sed(i,j)%dmean_slo = 0.d0
		sed(i,j)%Bed_D60 = 0.d0 !added 20240422
		sed(i,j)%Bed_D90 =0.d0
		sed(i,j)%samp_dep = samp_dep_min !added 20240419
		sed(i,j)%stor_sed_sum = 0.d0 !20240424
		
		! For non-uniform sized sediment
				 if(sed_type_switch .eq. 2) then
				   open (11, file = fi, status = "old")
					read (11,*)
				   	if(num_of_landuse.ge.3)then					 
					do m = 1, Np
					 read(11,*) di_temp, fm_temp, fm1_temp, fms_temp,fms2_temp
					 sed(i,j)%dsed(m) = di_temp
					 sed(i,j)%fm(m) = fm_temp
					 sed(i,j)%fm1(m) = fm1_temp
					 sed(i,j)%ft(m) = fm_temp
        			 sed(i,j)%fms(m) = fms_temp			!size distribution for plowed and paddy field 
					 if(land(i,j)==3) sed(i,j)%fms(m) = fms2_temp		!size distribution for mountainous area 20211226	
					 !sed(i,j)%fsur(m)=fm_temp		!20240509 
					 sed(i,j)%dzbpr(m) = 0.d0
						
					 sed(i,j)%qbi(m) = 0.d0
					 sed(i,j)%qsi(m) = 0.d0
					 sed(i,j)%qswi(m) = 0.d0
					 sed(i,j)%qdi(m) = 0.d0
					 sed(i,j)%ffd(m) = 0.d0
					 sed(i,j)%qdsum(m) = 0.d0
		! for bedload, suspended load; GSD added by Qin 2021/6/25
					 sed(i,j)%fqbi(m) = 0.d0
					 sed(i,j)%fqsi(m) = 0.d0
					! sed(i,j)%fsur(m) =0.d0 !2021/11/12	
					 sed(i,j)%Esisum(m) = 0.d0 !2021/11/26
					 sed(i,j)%Dsisum(m) = 0.d0								 
					 fm_all = fm_all + sed(i,j)%fm(m) !ここでは%表示。例えば10
					 fm1_all = fm1_all + sed(i,j)%fm1(m)
					 ft_all = ft_all + sed(i,j)%ft(m)
					 fms_all = fms_all + sed(i,j)%fms(m)
					enddo

					else
					do m = 1, Np
					 read(11,*) di_temp, fm_temp, fm1_temp, fms_temp
					 sed(i,j)%dsed(m) = di_temp
					 sed(i,j)%fm(m) = fm_temp
					 sed(i,j)%fm1(m) = fm1_temp
					 sed(i,j)%ft(m) = fm_temp
        			 sed(i,j)%fms(m) = fms_temp			
					 sed(i,j)%fsur(m)=fm_temp	 
					 sed(i,j)%dzbpr(m) = 0.d0
						
					 sed(i,j)%qbi(m) = 0.d0
					 sed(i,j)%qsi(m) = 0.d0
					 sed(i,j)%qswi(m) = 0.d0
					 sed(i,j)%qdi(m) = 0.d0
					 sed(i,j)%ffd(m) = 0.d0
					 sed(i,j)%qdsum(m) = 0.d0
		! for bedload, suspended load; GSD added by Qin 2021/6/25
					 sed(i,j)%fqbi(m) = 0.d0
					 sed(i,j)%fqsi(m) = 0.d0
					 !sed(i,j)%fsur(m) =0.d0 !2021/11/12	
					 sed(i,j)%Esisum(m) = 0.d0 !2021/11/26
					 sed(i,j)%Dsisum(m) = 0.d0		
					 sed(i,j)%stor_sed_sum_i(m) = 0.d0 ! added 20240424   						 
					 fm_all = fm_all + sed(i,j)%fm(m) !ここでは%表示。例えば10
					 fm1_all = fm1_all + sed(i,j)%fm1(m)
					 ft_all = ft_all + sed(i,j)%ft(m)
					 fms_all = fms_all + sed(i,j)%fms(m)
					enddo
					endif
		
					read(11,*)
					read(11,*) isnumber
		
				  close(11)
				  	sed(i,j)%Bed_D60 = sed(i,j)%dsed(Np) !added 20240422
					sed(i,j)%Bed_D90=sed(i,j)%dsed(Np)
					f90=0.d0
					nn = 0
					 do m = 1, Np
					   sed(i,j)%fm(m) = sed(i,j)%fm(m)/fm_all !ここでは%/100表示。例えば0.1
					   sed(i,j)%ft(m) = sed(i,j)%ft(m)/ft_all
					   sed(i,j)%fm1(m) = sed(i,j)%fm1(m)/fm1_all
					   sed(i,j)%fms(m) = sed(i,j)%fms(m)/fms_all
					    if(sed(i,j)%fms(m)>1.d0) then
	                    !write(*,*) 'input fmslo >1 ', i,j, m, sed(i,j)%dsed(m), fms_all,sed(i,j)%fms(m), sed(i,j)%fms(1),sed(i,j)%fms(2),sed(i,j)%fms(10)!modified 20240509
	                    write(*,*) 'input fmslo >1 ', j, ny+1-i, m, sed(i,j)%dsed(m), fms_all,sed(i,j)%fms(m), sed(i,j)%fms(1),sed(i,j)%fms(2),sed(i,j)%fms(10)!convert i,j in IRIC GUI 20250505		                
						stop
	                    endif
						sed(i,j)%fsur(m)=sed(i,j)%fm(m)  !20240509
					   sed(i,j)%dmean = sed(i,j)%dsed(m) * sed(i,j)%fm(m) + sed(i,j)%dmean
					   sed(i,j)%dmean_slo = sed(i,j)%dsed(m) * sed(i,j)%fms(m) + sed(i,j)%dmean_slo
					!define D90 for the bed surface GSD estimation; added by Qin 20230421	
					   if (nn.ge.0)then !modified 20240510
						nn=nn+1
						f90= f90 + sed(i,j)%fm(m)
						!added 20240422
						if(nn>1) then
							if(f90>0.6 .and. (f90-sed(i,j)%fm(m)) .le. 0.6 )sed(i,j)%Bed_D60= sed(i,j)%dsed(m-1)+ (0.6d0-f90+sed(i,j)%fm(m))*(sed(i,j)%dsed(m)-sed(i,j)%dsed(m-1))/sed(i,j)%fm(m)
						else
							if(f90 > 0.6) sed(i,j)%Bed_D60 = sed(i,j)%dsed(m)
						endif	
						if(f90>0.9d0) then
						if(nn.le.1)then
						sed(i,j)%Bed_D90= sed(i,j)%dsed(m)
						else
						sed(i,j)%Bed_D90= sed(i,j)%dsed(m-1)+ (0.9d0-f90+sed(i,j)%fm(m))*(sed(i,j)%dsed(m)-sed(i,j)%dsed(m-1))/sed(i,j)%fm(m)
						endif
						nn=-9999
						f90=-9999
						endif	
					   endif						 
						do n = 1, isnumber
					  		sed(i,j)%fd(m,n) = sed(i,j)%fm(m)
						enddo
		
						do n = isnumber + 1, Nl
					  		sed(i,j)%fd(m,n) = sed(i,j)%fm1(m)
						enddo
					 end do		
					 	!write(*,*) fm_all,fm1_all,ft_all, fms_all !20240509
					   ! write(*,*) (sed(i,j)%fm(m), m = 1, Np)					
					   ! write(*,*) (sed(i,j)%fms(m), m = 1, Np)
 					   ! pause 
					 
		!----------------------------------------------------------------
		
		! For uniformed sediment
				elseif(sed_type_switch .eq. 1) then
					 sed(i,j)%dmean = ds_river
					 sed(i,j)%dmean_slo = ds_river   !tentative
					 !sed(i,j)%dzbpr(m) = 0.d0
						
					 do m = 1, Np
						  sed(i,j)%fm(m) = 1.d0
						 sed(i,j)%ft(m) = 1.d0
						 sed(i,j)%fms(m) = 1.d0
						 
						 sed(i,j)%qbi(m) = 0.d0
						 sed(i,j)%qsi(m) = 0.d0
						 sed(i,j)%qswi(m) = 0.d0
					 end do		
				else	 	 
				end if	
	end subroutine sed_distr

!----------------------------------------------------------------------------
! calculation of surface slope
subroutine cal_zs_slope
use globals
use sediment_mod !added 20240422
implicit none
integer i, j, ii, jj
real(8) distance11
real(8) slope11, thet11, thetr11

do i = 1, ny
 do j = 1, nx

distance11 = 0.0

!write(*,'(4i7,2f)') i, j, domain(i,j), dir(i,j), dx, dy

  if(domain(i,j).eq.0) cycle

  ! right
  if( dir(i,j).eq.1 ) then
   ii = i
   jj = j + 1
   distance11 = dx
  ! right down
  elseif( dir(i,j).eq.2 ) then
   ii = i + 1
   jj = j + 1
   distance11 = sqrt(dx*dx+dy*dy)
  ! down
  elseif( dir(i,j).eq.4 ) then
   ii = i + 1
   jj = j
   distance11 = dy
  ! left down
  elseif( dir(i,j).eq.8 ) then
   ii = i + 1
   jj = j - 1
   distance11 = sqrt(dx*dx+dy*dy)
  ! left
  elseif( dir(i,j).eq.16 ) then
   ii = i
   jj = j - 1
   distance11 = dx
  ! left up
  elseif( dir(i,j).eq.32 ) then
   ii = i - 1
   jj = j - 1
   distance11 = sqrt(dx*dx+dy*dy)
  ! up
  elseif( dir(i,j).eq.64 ) then
   ii = i - 1
   jj = j
   distance11 = dy
  ! right up
  elseif( dir(i,j).eq.128 ) then
   ii = i - 1
   jj = j + 1
   distance11 = sqrt(dx*dx+dy*dy)
  elseif( dir(i,j).eq.0 ) then
   ii = i
   jj = j
   distance11 = dx
  else
   write(*,*) "dir(i, j) is error (", i, j, ")", dir(i, j)
   stop
  endif

  if( ii .gt. ny ) cycle
  if( jj .gt. nx ) cycle
  if( ii .lt. 1 ) cycle
  if( jj .lt. 1 ) cycle
!  if( domain(ii,jj).eq.0 ) cycle

  slope11 = (zs(i,j) - zs(ii,jj))/distance11
  thetr11 = atan(slope11)
  thet11 = atan(slope11) * 180./3.14159265
  zs_slope(i,j) = thet11

!   if(thet11.gt.1.0) then
!    write(*,'(2i5,5f12.5)') i, j, thet11, slope11, zs(i,j), zs(ii,jj), distance11
!    pause
!   endif

 enddo
enddo

end subroutine cal_zs_slope

!----------------------------------------------------------------
 ! 2D -> 1D 'sed_struct'
      subroutine sub_sed_ij2idx ( a, a_idx )
	     use globals
	     use sediment_mod
	     implicit none
		 
	     type (sed_struct) a(ny,nx), a_idx(riv_count)
	     integer k, m, n
		 
	     do k = 1, riv_count
	         do m = 1, Np
		    	 a_idx(k)%dsed(m) = a (riv_idx2i(k), riv_idx2j(k))%dsed(m)
			     a_idx(k)%fm(m) = a (riv_idx2i(k), riv_idx2j(k))%fm(m)
			     a_idx(k)%ft(m) = a (riv_idx2i(k), riv_idx2j(k))%ft(m)
			     a_idx(k)%dzbpr(m) = a (riv_idx2i(k), riv_idx2j(k))%dzbpr(m)
				 a_idx(k)%fm1(m) = a (riv_idx2i(k), riv_idx2j(k))%fm1(m) ! 2021/10/04
				 
			     a_idx(k)%qbi(m) = a (riv_idx2i(k), riv_idx2j(k))%qbi(m)
			     a_idx(k)%qsi(m) = a (riv_idx2i(k), riv_idx2j(k))%qsi(m)
			     a_idx(k)%qswi(m) = a (riv_idx2i(k), riv_idx2j(k))%qswi(m)
		         a_idx(k)%qdi(m) = a (riv_idx2i(k), riv_idx2j(k))%qdi(m)
		         a_idx(k)%ffd(m) = a (riv_idx2i(k), riv_idx2j(k))%ffd(m)
		         a_idx(k)%qdsum(m) = a (riv_idx2i(k), riv_idx2j(k))%qdsum(m)

! for bedload, suspended load; GSD added by Qin 2021/6/25
				 a_idx(k)%fqbi(m) = a (riv_idx2i(k), riv_idx2j(k))%fqbi(m)
				 a_idx(k)%fqsi(m) = a (riv_idx2i(k), riv_idx2j(k))%fqsi(m)
				 a_idx(k)%fsur(m) = a(riv_idx2i(k),riv_idx2j(k))%fsur(m) !2021/11/12
				 a_idx(k)%Esisum(m) = a(riv_idx2i(k),riv_idx2j(k))%Esisum(m)!2021/11/26
				 a_idx(k)%Dsisum(m) = a(riv_idx2i(k),riv_idx2j(k))%Dsisum(m)	
				 a_idx(k)%stor_sed_sum_i(m) = a(riv_idx2i(k),riv_idx2j(k))%stor_sed_sum_i(m)!20240424			 
			     do n = 1,Nl
			      a_idx(k)%fd(m,n) = a (riv_idx2i(k), riv_idx2j(k))%fd(m,n)
		         end do				  
		  	 end do
	            a_idx(k)%dmean = a(riv_idx2i(k), riv_idx2j(k))%dmean
				a_idx(k)%Bed_D60= a(riv_idx2i(k),riv_idx2j(k))%Bed_D60 !added 20240422
				a_idx(k)%Bed_D90= a(riv_idx2i(k),riv_idx2j(k))%Bed_D90
				a_idx(k)%samp_dep = a(riv_idx2i(k),riv_idx2j(k))%samp_dep !20240419
				a_idx(k)%stor_sed_sum = a(riv_idx2i(k),riv_idx2j(k))%stor_sed_sum !20240424
	     end do
      end subroutine sub_sed_ij2idx
	  
!---------------------------------------------------------------------------------------------------------

! 1D -> 2D 'sed_struct'
      subroutine sub_sed_idx2ij ( a_idx, a )
         use globals
	     use sediment_mod
	     implicit none
		 
	     type (sed_struct) a(ny,nx), a_idx(riv_count)
	     integer k, m, n
		 
	     do k = 1, riv_count
		     do m = 1, Np
	             a (riv_idx2i(k), riv_idx2j(k))%dsed(m) = a_idx(k)%dsed(m)
			     a (riv_idx2i(k), riv_idx2j(k))%fm(m) = a_idx(k)%fm(m) 
			     a (riv_idx2i(k), riv_idx2j(k))%ft(m) = a_idx(k)%ft(m)
			     a (riv_idx2i(k), riv_idx2j(k))%dzbpr(m) = a_idx(k)%dzbpr(m)
			      
			     a (riv_idx2i(k), riv_idx2j(k))%qbi(m) = a_idx(k)%qbi(m)
			     a (riv_idx2i(k), riv_idx2j(k))%qsi(m) = a_idx(k)%qsi(m)
			     a (riv_idx2i(k), riv_idx2j(k))%qswi(m) = a_idx(k)%qswi(m)
			     a (riv_idx2i(k), riv_idx2j(k))%qdi(m) = a_idx(k)%qdi(m)
			     a (riv_idx2i(k), riv_idx2j(k))%ffd(m) = a_idx(k)%ffd(m)
			     a (riv_idx2i(k), riv_idx2j(k))%qdsum(m) = a_idx(k)%qdsum(m)
			     a (riv_idx2i(k), riv_idx2j(k))%ssi(m) = a_idx(k)%ssi(m) !20240619
! for bedload, suspended load; GSD added by Qin 2021/6/25
				  a (riv_idx2i(k), riv_idx2j(k))%fqbi(m) = a_idx(k)%fqbi(m) 
				 a (riv_idx2i(k), riv_idx2j(k))%fqsi(m) =  a_idx(k)%fqsi(m) 				 
				 a (riv_idx2i(k), riv_idx2j(k))%fsur(m) =  a_idx(k)%fsur(m) !2021/11/12
				 a (riv_idx2i(k), riv_idx2j(k))%Esisum(m) =  a_idx(k)%Esisum(m) !2021/11/26
				 a (riv_idx2i(k), riv_idx2j(k))%Dsisum(m) =  a_idx(k)%Dsisum(m) 
				 a (riv_idx2i(k), riv_idx2j(k))%stor_sed_sum_i(m) =  a_idx(k)%stor_sed_sum_i(m) !20240424 
			     do n = 1,Nl
				     a (riv_idx2i(k), riv_idx2j(k))%fd(m,n) = a_idx(k)%fd(m,n)
		         end do				 
	         end do
             a (riv_idx2i(k), riv_idx2j(k))%dmean = a_idx(k)%dmean	
			 a(riv_idx2i(k), riv_idx2j(k))%Bed_D60 = a_idx(k)%Bed_D60 !added 20240422				 
			 a(riv_idx2i(k), riv_idx2j(k))%Bed_D90 = a_idx(k)%Bed_D90	
			 a(riv_idx2i(k), riv_idx2j(k))%samp_dep = a_idx(k)%samp_dep !20240419	
			 a(riv_idx2i(k), riv_idx2j(k))%stor_sed_sum = a_idx(k)%stor_sed_sum !20240424	  
	     end do		 
	  end subroutine sub_sed_idx2ij	  
	  
!---------------------------------------------------------------------------------------------------------
 ! 1D -> link 'sed_struct'  added by harada   ! use downstream end information
      subroutine sub_sed_idxlin ( a_idx, a_lin )
		use globals
		use sediment_mod
		implicit none
		
		type (sed_struct) a_idx(riv_count)
		type (sed_struct2) a_lin(link_count)
		integer l, m, n
		
		do l = 1, link_count
			do m = 1, Np
				a_lin(l)%dsed(m) = a_idx (link_idx_k(l))%dsed(m)
				a_lin(l)%fm(m) = a_idx (link_idx_k(l))%fm(m)
				a_lin(l)%ft(m) = a_idx (link_idx_k(l))%ft(m)
				a_lin(l)%dzbpr(m) = a_idx (link_idx_k(l))%dzbpr(m)
				
				a_lin(l)%qbi(m) = a_idx (link_idx_k(l))%qbi(m)
				a_lin(l)%qsi(m) = a_idx (link_idx_k(l))%qsi(m)
				a_lin(l)%qswi(m) = a_idx (link_idx_k(l))%qswi(m)
				a_lin(l)%qdi(m) = a_idx (link_idx_k(l))%qdi(m)
				a_lin(l)%ffd(m) = a_idx (link_idx_k(l))%ffd(m)
				a_lin(l)%qdsum(m) = a_idx (link_idx_k(l))%qdsum(m)

! for bedload, suspended load; GSD added by Qin 2021/6/25
				a_lin(l)%fqbi(m) = a_idx(link_idx_k(l))%fqbi(m)
				a_lin(l)%fqsi(m) = a_idx(link_idx_k(l))%fqsi(m)							
				a_lin(l)%fsur(m) = a_idx(link_idx_k(l))%fsur(m) !2021/11/12
				a_lin(l)%Esisum(m) = a_idx(link_idx_k(l))%Esisum(m) !2021/11/26
				a_lin(l)%Dsisum(m) = a_idx(link_idx_k(l))%Dsisum(m)
				a_lin(l)%stor_sed_sum_i(m) = a_idx(link_idx_k(l))%stor_sed_sum_i(m) !20240424				
				do n = 1,Nl
					a_lin(l)%fd(m,n) = a_idx (link_idx_k(l))%fd(m,n)
				end do				  
		 	end do
		 	a_lin(l)%dmean = a_idx (link_idx_k(l))%dmean
			a_lin(l)%Bed_D60=a_idx(link_idx_k(l))%Bed_D60 !added 20240422			
			a_lin(l)%Bed_D90=a_idx(link_idx_k(l))%Bed_D90
			a_lin(l)%samp_dep=a_idx(link_idx_k(l))%samp_dep !20240419
			a_lin(l)%stor_sed_sum=a_idx(link_idx_k(l))%stor_sed_sum !20240424
		end do
	 end subroutine sub_sed_idxlin

!---------------------------------------------------------------------------------------------------------
! link -> 1D 'sed_struct'  added by harada  
	 subroutine sub_sed_linidx ( a_lin, a_idx )
		use globals
		use sediment_mod
		implicit none
		
		type (sed_struct) a_lin(link_count)
		type (sed_struct) a_idx(riv_count)
		integer k, m, n, l
		
		do k = 1, riv_count
			do m = 1, Np
				a_idx(k)%dsed(m) = a_lin(link_to_riv(k))%dsed(m)
				a_idx(k)%fm(m) = a_lin(link_to_riv(k))%fm(m) 
				a_idx(k)%ft(m) = a_lin(link_to_riv(k))%ft(m)
				a_idx(k)%dzbpr(m) = a_lin(link_to_riv(k))%dzbpr(m)
				 
				a_idx(k)%qbi(m) = a_lin(link_to_riv(k))%qbi(m)
				a_idx(k)%qsi(m) = a_lin(link_to_riv(k))%qsi(m)
				a_idx(k)%qswi(m) = a_lin(link_to_riv(k))%qswi(m)
				a_idx(k)%qdi(m) = a_lin(link_to_riv(k))%qdi(m)
				a_idx(k)%ffd(m) = a_lin(link_to_riv(k))%ffd(m)
				a_idx(k)%qdsum(m) = a_lin(link_to_riv(k))%qdsum(m)
				a_idx(k)%ssi(m) = a_lin(link_to_riv(k))%ssi(m)!---added by Qin

! for bedload, suspended load; GSD added by Qin 2021/6/25
				a_idx(k)%fqbi(m) = a_lin(link_to_riv(k))%fqbi(m)
				a_idx(k)%fqsi(m) = a_lin(link_to_riv(k))%fqsi(m)					
				a_idx(k)%fsur(m) = a_lin(link_to_riv(k))%fsur(m) !2021/11/12
				a_idx(k)%Esisum(m) = a_lin(link_to_riv(k))%Esisum(m)!2021/11/26
				a_idx(k)%Dsisum(m) = a_lin(link_to_riv(k))%Dsisum(m)
				a_idx(k)%stor_sed_sum_i(m) = a_lin(link_to_riv(k))%stor_sed_sum_i(m) ! 20240424
				do n = 1,Nl
					a_idx(k)%fd(m,n) = a_lin(link_to_riv(k))%fd(m,n)
				end do				 
			end do
			a_idx(k)%dmean = a_lin(link_to_riv(k))%dmean
			a_idx(k)%Bed_D60=a_lin(link_to_riv(k))%Bed_D60 ! added 20240422				
			a_idx(k)%Bed_D90=a_lin(link_to_riv(k))%Bed_D90	
			a_idx(k)%samp_dep=a_lin(link_to_riv(k))%samp_dep	!20240419			
			a_idx(k)%stor_sed_sum=a_lin(link_to_riv(k))%stor_sed_sum !20240424  
		end do		 
	 end subroutine sub_sed_linidx	 

!---------------------------------------------------------------------------------------------------------
!-------------------------------------
! 1D -> link    added by harada
subroutine sub_riv_idxlin( a_idx, a_lin )
  use globals
  use sediment_mod
  implicit none
  
  real(8) a_idx(riv_count), a_lin(link_count)
  integer l
  
  !a_idx(:) = 0.d0
  do l = 1, link_count
   a_lin(l) = a_idx(link_idx_k(l))   !downstream cell  (not average)
  enddo
  
  end subroutine sub_riv_idxlin

!-------------------------------------
! link -> 1D    added by harada
  subroutine sub_riv_linidx( a_lin, a_idx )
    use globals
    use sediment_mod !added 20240422
    implicit none
    
    real(8) a_idx(riv_count), a_lin(link_count)
    integer k
    
    !a_idx(:) = 0.d0
    do k = 1, riv_count
      a_idx(k) = a_lin(link_to_riv(k))
    enddo
    
    end subroutine sub_riv_linidx  

!--------Debris setting-------------------------------
subroutine debris_setting
		use globals
		use sediment_mod
		implicit none

		real(8), parameter :: PI = acos(-1.0d0)
    real(8), parameter :: rho = 1.d0   !=1.0d0(kg/m3)  density of water
		real(8) tanp, hsc0
    integer k

	tanp = tan(PI/180.0d0*phi)      !=tan(phi)

	do k = 1, slo_count
      !write(*,*) da_idx(k)
      !write(*,*) slo_grad(k), tanp
		  c_dash(k) = cohe/(rho*grav*da_idx(k)*cos(atan(slo_grad(k)))*tanp)
!      hsc0 = ((1.d0-tan(atan(slo_grad(k)))/tanp)*((1.d0-lambda)*s/rho+pw(i,j))+&
!      c2-h(i,j)/d(i,j)*tan(atan(grads(i,j)))/tanp) &
!      /((1.d0-tan(atan(grads(i,j)))/tanp)*((1.d0-lamda)+pw(i,j))+tan(atan(grads(i,j)))/tanp)
		soildepth_idx_deb(k) = d_mp_ini
    end do

    pc = 1.0d0 - pf

end subroutine debris_setting
!------------------------------------------------------


 ! Determination of river bed slope
      subroutine det_rivebedslope
		use globals
		use sediment_mod
		use dam_mod, only: dam_switch,damflg !added by Qin 2021/6/13
		implicit none

		real(8) slope, slope1,thet, thetr,thet1
		integer k, kk,kkk,m,n
!added parallel 20231028		
!!$omp parallel do private(slope, thet, thetr,kk)		   
		do k = 1, riv_count-1
		   kk = down_riv_idx(k)
		   slope = (zb_riv_idx(k) - zb_riv_idx(kk))/dis_riv_idx(k)
		   thetr = atan(slope)
		   thet = atan(slope) * 180./3.14159
		   zb_riv_slope_idx(k) = thet
		end do
   
		   zb_riv_slope_idx(riv_count) = zb_riv_slope_idx(riv_count-1)

!added by Qin 2021/6/9
		if(ave_slope_switch == 1)then
		   m = mod(riv_count,ave_cell_num)
!added parallel 20231028		   
!!$omp parallel do private(slope, slope1,thet, thetr,thet1,kk,n)			   
		do k = 1, riv_count-m
		   slope = 0.d0
		   kk = down_riv_idx(k)
		   n = 0
		   do 
			   slope1 = tan(zb_riv_slope_idx(kk)*3.14159/180.)
			   if(slope1 .le. min_slope) slope1 = min_slope
			   slope = slope + slope1
			   kk = down_riv_idx(kk)
			   n = n+1
			   if(n==ave_cell_num-1) exit
		   enddo	
			   slope1 = tan(zb_riv_slope_idx(k)*3.14159/180.)
			   if(slope1 .le. min_slope) slope1 = min_slope
			   slope = (slope1+slope)/real(ave_cell_num)
			   thetr = atan(slope)
			   thet = atan(slope) * 180./3.14159
			   zb_riv_slope_idx(k)= thet
	   enddo	
	   slope = 0.d0
	   !do k = riv_count-m-1, riv_count-1	
		   kk  = down_riv_idx(riv_count-m)
		   n = 0
		   do 
			   slope1 = tan(zb_riv_slope_idx(kk)*3.14159/180.)
			   if(slope1 .le. min_slope) slope1 = min_slope
			   slope = slope + slope1
			   kk = down_riv_idx(kk)
			   n = n+1
			   if(n==m) exit
		   enddo	
	   !enddo	
		   slope = (tan(zb_riv_slope_idx(riv_count-m)*3.14159/180.)+slope)/real(m+1)
		   thetr = atan(slope)
		   thet = atan(slope) * 180./3.14159
	   do k = riv_count-m, riv_count-1
			   kk = down_riv_idx(k)
			   zb_riv_slope_idx(kk) = thet
	   enddo	
!		zb_riv_slope_idx(riv_count) = zb_riv_slope_idx(riv_count-1)
	   endif	
!----added by Qin 2021/6/13
	   if(dam_switch==1)then
		do k = 1, riv_count
			m = 0
			thet = atan(min_slope)* 180./3.14159 !modified 20240427
			if(damflg(k).gt.0)then
				do n = 1,8
					kkk = up_riv_idx(k,n)
					if (domain_riv_idx(kkk).eq.0) cycle
					m = m+1
					thet = thet+ zb_riv_slope_idx(kkk)
				enddo
				zb_riv_slope_idx(k)= thet/real(m)
			endif
		enddo			
		endif

	 end subroutine det_rivebedslope	  

 ! Determination of river bed slope/ added by harada
	  subroutine det_rivebedslope2(hr_idx, hr_lin,water_v_lin, hr_idxa, qr_ave_idx)
		use globals
		use sediment_mod
		use dam_mod!, only: dam_switch,damflg ! added by Qin
		use omp_lib
		implicit none

		real(8) hr_idx(riv_count),hr_idxa(riv_count),qr_ave_idx(riv_count) !added 20240422
		real(8) hr_lin(link_count), water_v_lin(link_count)
		real(8) slope, thet, thetr,zbtemp, hrtemp,thet1, water_v, h_dam, slope1, slope11
		integer l,ll, lll,k, kk, kkk,kkkk,n, m, f

! water volume(m3) of each link modified 20230924	!moved to here 20240426
!!$omp parallel do private(k,water_v) schedule(static)
	do l = 1, link_count
		k = link_idx_k(l)
		water_v = 0.
		if (damflg(k).gt.0)then ! added by Qin
			water_v_lin(l) = dam_w_vol(damflg(k))
			area_lin(l) = dam_reserv_area(damflg(k))
		else	
		do
			if(hr_idx(k).le.0.d0) hr_idx(k) =0.
			water_v = water_v + dis_riv_idx(k) * width_idx(k) * hr_idx(k) !the actual water volume in the unit channel
			if(node_ups(k) .ge. 1 ) then
			water_v_lin(l)=water_v 
				exit
			end if
			if(up_riv_idx(k,1) .eq. 0 ) then
			water_v_lin(l)=water_v 
				exit
			end if
			k = up_riv_idx(k,1)
		end do
		if(isnan(water_v_lin(l))) then !check 2021/6/3
			write(*,*) water_v_lin(l), hr_idx(k),hr_idxa(k),hr_lin(l),k,l,link_0th_order(l),zb_riv_slope0_lin(l),zb_riv_idx(k)-zb_riv0_idx(k)
			stop "Water volume in unit channel is nan"
		endif
		endif
		if(water_v_lin(l).le.0.) then
		 water_v_lin(l)= 0.
		 hr_lin(l) = 0.
		endif
		hr_lin(l) = water_v_lin(l)/area_lin(l) ! water depth of unit channel taken by averaged value ;added 20240419
	end do	
!----set the water level of dam reservoir same as the average level of its neighbouring upstream cells; added by Qin;20240424		
	if(dam_switch.gt.0)then 
	   do f = 1, dam_num
	   	  k = dam_loc(f)
		  l = link_to_riv(k)
			h_dam = 0.d0
			m = 0
			do n = 1,8
				if (up_riv_lin(l,n)==0) cycle
				if (up_riv_lin(l,n)==0) exit
				m = m+1
				ll = up_riv_lin(l,n)
				kk = link_idx_k(ll)
				h_dam = h_dam+ hr_lin(ll)+zb_riv_idx(kk)
			enddo
			hr_lin(l) = h_dam/real(m)-zb_riv_idx(k)
			if(hr_lin(l).le.0.d0) hr_lin(l) = 0.d0 !modified by Qin 2021/11/06
		enddo	
	endif

!added parallel 20231028		
!!$omp parallel do private(slope, thet, thetr,zbtemp,hrtemp,ll,lll, k, kk, kkk,kkkk,n, m) schedule(static)
		do l = 1, link_count
			k = link_ups_k(l)
			kk = link_idx_k(l)
			kkk = down_riv_idx(kk)
			lll= link_to_riv(kkk)
			!modified the link slope by Qin 2021/5/25
			if(domain_riv_idx(kkk).eq.0 .or. damflg(kkk).gt.0) then
				m = 0
				zbtemp = 0.d0
				hrtemp = 0.d0
				do n = 1, 8
					if (up_riv_lin(l,n) == 0) cycle
					m = m + 1
					ll = up_riv_lin(l,n)
					kkkk = link_idx_k(ll)
					zbtemp = zbtemp + zb_riv_idx(kkkk)
					hrtemp = hrtemp +hr_idx(kkkk)
				enddo
				if(m.lt.1) then
				slope =	(zb_riv_idx(k) - zb_riv_idx(kk)) / link_len(l)
				else
                slope = (zbtemp/real(m)	- zb_riv_idx(kk)) / (link_len(l) + dis_riv_idx(kkkk))
				endif
			else
				slope = (zb_riv_idx(k)-zb_riv_idx(kkk))/(link_len(l) + dis_riv_idx(kkk))
			endif
			!if(slope .le. 1e-3) slope = (zb_riv_idx(kk) - zb_riv_ idx(kkk))/dis_riv_idx(kk)
			if(slope .le. min_slope) slope = min_slope     ! in case still the slope value is too small
			thetr = atan(slope)
			thet = atan(slope) * 180./3.14159
			zb_riv_slope_lin(l) = thet ! degree			
		end do
!----slope smoothing for the short unit channel---added by Qin 2021/8/18--
!OMP must be prohibited for this loop 
!!$omp parallel do private( thet, thet1,ll,lll, k, kk, n, m) 
		do l = 1, link_count
			k = link_idx_k(l)
			kk = down_riv_idx(k)
			if(link_cell_num(l) > min_num_cell_link) cycle 
			if(link_0th_order(l)== 0) cycle!---Do not conduct the slope smoothing for 0th  link 
			if( damflg(k) < 1)then  
				m = 0
				do n =1,8	
					if (up_riv_lin(l,n) == 0) cycle
					ll = up_riv_lin(l,n)
					thet1 = zb_riv_slope_lin(ll)
					m = m + 1
					if(m.eq.1)then
					 thet = thet1
					endif 
					if(m.gt.1.and.thet1.lt.thet)then 
					thet = thet1
					endif
				enddo	
	!			if(m.gt.1)then
	!			zb_riv_slope_lin(l)= thet
	!			else
				lll= link_to_riv(kk)
				thet1 = zb_riv_slope_lin(lll)
	!			slope =	(tan(thet*3.14159/180.)+ tan(thet1*3.14159/180.))/2.
	!			if(slope .le. 1e-3) slope = 1e-3 
				zb_riv_slope_lin(l) = (thet+thet1)/2.
				if(damflg(kk) > 0) then
				 zb_riv_slope_lin(l) = thet
				endif	
			else ! calculation of the slope of dam link
				m = 0
				thet = atan(min_slope)* 180./3.14159 !modified 20240427
				do n =1,8	
					if (up_riv_lin(l,n) == 0) cycle
					ll = up_riv_lin(l,n)
					thet1 = zb_riv_slope_lin(ll)
					m = m + 1
					if(m.eq.1) then
					thet = thet1
					endif
					if(m.gt.1.and.thet1.lt.thet) then 
					thet = thet1
					endif
				enddo	
				zb_riv_slope_lin(l)= thet
			endif
		enddo	
!		enddo
			
		!endif
	!OMP must be prohibited for this loop 
	!!$omp parallel do private(k,slope11,h_dam,m,n,ll,kk)
	!turned off at 20240419
	do l = 1, link_count
		k = link_idx_k(l)
		slope11 = tan(zb_riv_slope_lin(l)*3.14159/180.)!use the slope of link; modified by Qin
		if(qr_ave_idx(k) .le. 0.d0) qr_ave_idx(k) = 0.d0	!modified by Qin 2021/6/3
		if(slope11.le.min_slope)then !modified 20231025 !turned off at 20231102
		slope11 = min_slope
		endif
		hr_lin(l) = ( ns_river * qr_ave_idx(k)/ width_lin(l) / sqrt(slope11) )**(3./5.)	!use the width of link; modified by Qin
		if(hr_lin(l).le.0.d0) hr_lin(l) = 0.d0!modified by Qin 20230924 20231024
	!----set the water level of dam reservoir same as the average level of its neighbouring upstream cells; added by Qin	
		if(damflg(k).gt.0)then 
			h_dam = 0.d0
			m = 0
			do n = 1,8
				if (up_riv_lin(l,n)==0) cycle
				if (up_riv_lin(l,n)==0) exit
				m = m+1
				ll = up_riv_lin(l,n)
				kk = link_idx_k(ll)
				h_dam = h_dam+ hr_lin(ll)+zb_riv_idx(kk)
			enddo
			hr_lin(l) = h_dam/real(m)-zb_riv_idx(k)
			if(hr_lin(l).le.0.d0) hr_lin(l) = 0.d0 !modified by Qin 2021/11/06
		endif		
	end do

	!!$omp parallel do private(l) schedule(static)
	do k = 1, riv_count
		l = link_to_riv(k)
		hr_idxa(k) = hr_lin(l)
	enddo

	  end subroutine det_rivebedslope2
!---------------------------------------------------------------------------------------------------------
! Determining new dmean
 !     subroutine mean_diameter(dzb_temp, sed_idx, hr_idx, qr_ave_idx, ust_idx, width_idx, len_riv_idx, Nb_idx)
 !     subroutine mean_diameter(dzb_temp, sed_idx, hr_idx, qr_ave_idx, ust_idx,qsb_idx,qss_idx,Et_idx ) !revised by Qin 20210924
		subroutine mean_diameter(dzb_temp, sed_idx, qsb_idx,qss_idx)	
	     use globals
	     use sediment_mod
		 
	     implicit none
		 
	     type(sed_struct) sed_idx(riv_count)
	     real(8) dzb_temp(riv_count)
		 real(8) qsb_idx(riv_count), qss_idx(riv_count)! 2021/6/25;added by Qin
 !        integer Nb_idx(riv_count)
!	     real(8) width_idx(riv_count), len_riv_idx(riv_count)
!		 real(8) Et_idx(riv_count) !revised by Qin 202120924
	     !real(8) Embnew 
		 real(8) Etnew, Fmnew, Ftnew, dzb
	     real(8) FM, FT, DZBPR1, Fmall, Ftall, Fdall, Fsur_all,f90
		 !real(8) Fball, Fsall !
	     integer k, kk, m, n, Nbnew, Nbl,nn
	     !real(8) test
	     !real(8) fm_sum111(30)
	     !real(8) fm_sum222
	     !real(8) :: sfm1, sfm2, sfm3
	     !real(8) :: trans, trans1, vari, all
	    ! real(8) :: u_flax, Frn
         !integer k1d, k1u, k1
	     !real(8) :: x0, x1, y0, y1, slope
		 !integer inb
		real(8),parameter:: csta = 0.6d0
		real(8) De_depo, De_Emb,qdi, c_ave,Bed_depth
		integer DNl

!pause'in subroutine mean diameter'
		
        sed_idx(:)%dmean = 0.0
		!sed_idx(:)%Bed_D90=0.0
		c_ave  = csta/2.d0		

!added parallel 20231028 turned off 20240424
!!$omp parallel do private( Etnew, Fmnew, Ftnew, dzb,FM, FT, DZBPR1, Fmall, Ftall, Fdall, Fsur_all,f90, &
!!$omp kk, m, n, Nbnew, Nbl,nn,De_depo, De_Emb,qdi, Bed_depth, DNl )		
        do k = 1, riv_count
!		 inb = 0 !inb = 1 when it reached to riverbed
		 dzb = dzb_temp(k)
! Update the new deposition layer number		 
!if(Emb_idx(k).ne.0.d0) then

	     if(dzb.ge.0.) then !dzb>0 aggradation
			if (Et_idx(k)+dzb .le. Ed) then
			 Etnew = Et_idx(k) + dzb !eq(78,1)
			 Nbnew = Nb_idx(k)       !eq(80,1)
			else
			 Etnew = Et_idx(k) + dzb - Ed !eq(78,2)
			 Nbnew = Nb_idx(k) + 1        !eq(80,2)
			endif
		  else !dzb<0 degradation
			if(Et_idx(k)+dzb .le. 0.) then
			 Etnew = Ed + Et_idx(k) + dzb !eq(83,2)
			 Nbnew = Nb_idx(k) - 1        !eq(85,2)
			else
			 Etnew = Et_idx(k) + dzb !eq(83,1)
			 Nbnew = Nb_idx(k)       !eq(85,1)
			endif
		  endif
		  
		  if(Etnew .le. 0.) Etnew = dabs(dzb)*0.000001 !
		  !check here 2021/5/14
		  if(Nbnew .eq. 0) then
!modified by Qin 20211018			
			   Nbnew = Nl 
			   do m = 1, Np
				 do n = 1, Nl
					  sed_idx(k)%fd(m,n) = sed_idx(k)%fm1(m)
				 end do
				end do			  
		  elseif(Nbnew .eq. Nl+1) then
			  Nbnew = Nl
			do m = 1, Np
			 do n = 1, Nl
				 if (n.eq.Nl) then
				  sed_idx(k)%fd(m,n) = sed_idx(k)%fd(m,n)
				 else
				  sed_idx(k)%fd(m,n) = sed_idx(k)%fd(m,n+1)
				 endif
			 end do
			end do
		  endif

      do m = 1, Np
		    Nbl = Nb_idx(k)
		    FM = sed_idx(k)%fm(m)
		    FT = sed_idx(k)%ft(m)
		    DZBPR1 = sed_idx(k)%dzbpr(m)
			 
	    if (dzb .ge. 0.) then !dzb>0 aggradation
		    !Fmnew = (1.-dzb/Emb_idx(k))*FM + DZBPR1/Emb_idx(k) !eq(77)
			Fmnew = (1.-(1.-lambda)*dzb/(Emb_idx(k)*C_ave))*FM + DZBPR1*(1-lambda)/(Emb_idx(k)*c_ave) !modified by Qin 2021/7/16
			
		     if (Et_idx(k)+dzb .le. Ed) then
		       Ftnew = Et_idx(k)/Etnew*FT + dzb/Etnew*FM       !eq(79,1)
		       sed_idx(k)%fd(m,Nbnew) = sed_idx(k)%fd(m,Nbnew) !eq(81,1)
		     else
		       Ftnew = FM                                                    !eq(79,2)
		       sed_idx(k)%fd(m,Nbnew) = Et_idx(k)/Ed*FT+(1.-Et_idx(k)/Ed)*FM !eq(81,2)
		     endif
		else !dzb<0 degradation
		     if (Et_idx(k)+dzb .le. 0.) then
			 !  Fmnew = FM +  Et_idx(k)/Emb_idx(k)*FT - (Et_idx(k) + dzb)/Emb_idx(k)*sed_idx(k)%fd(m,Nbl) + DZBPR1/Emb_idx(k) !eq(82,2)
		   	   Fmnew = FM + ((1.-lambda)*Et_idx(k)/(Emb_idx(k)*c_ave))*FT -(1.-lambda)* (Et_idx(k) + dzb)/(Emb_idx(k)*C_ave)*sed_idx(k)%fd(m,Nbl)+ DZBPR1*(1-lambda)/(Emb_idx(k)*c_ave) !modified by Qin 2021/7/16
			   Ftnew = sed_idx(k)%fd(m,Nbl)                                                                                  !eq(84,2)
			 else
			   !Fmnew = FM - dzb/Emb_idx(k)*FT + DZBPR1/Emb_idx(k) !eq(82,1)
			   Fmnew = FM-((1.-lambda)*dzb/(Emb_idx(k)*c_ave))*FT + DZBPR1*(1-lambda)/(Emb_idx(k)*c_ave) !modified by Qin 2021/7/16
			   Ftnew = FT                                         !eq(84,1)
			 endif
		endif
		if(Fmnew.lt.0.0) then
	         sed_idx(k)%fm(m) = 0.0
		else
	         sed_idx(k)%fm(m) = Fmnew
		endif

		if(Ftnew.lt.0.0) then
	         sed_idx(k)%ft(m) = 0.0
		else
	         sed_idx(k)%ft(m) = Ftnew
		endif

	       Et_idx(k) = Etnew
	       Nb_idx(k) = Nbnew
      enddo

!      if(inb.eq.0) then 
	     Fmall = 0.
	     Ftall = 0.
	     Fdall = 0.
		Nbl = Nb_idx(k) 
	     do m = 1, Np
		   !  Nbl = Nb_idx(k)
		     Fmall = Fmall + sed_idx(k)%fm(m)
		     Ftall = Ftall + sed_idx(k)%ft(m)
		     Fdall = Fdall + sed_idx(k)%fd(m,Nbl)
	     end do
		! Nbl = Nb_idx(k)
!		 if(Ftall == 0.d0) then !modified by Qin 20211009
!			Nb_idx(k) = Nbl - 1
!			Et_idx(k) = Ed
!		 endif	
		 fsur_all=0.d0
		 nn = 0
		 f90=0.d0
		 De_Emb = Emb_idx(k)*c_ave*dlambda
		 Bed_depth = zb_riv_idx(k) - zb_roc_idx(k)+De_Emb
		 qdi = 0.d0

!check later 20231106	 
	     do m = 1, Np
		 if(Fdall ==0.d0)then

			if(Nb_idx(k)-1==0.d0) then
			sed_idx(k)%fd(m,Nbl) = sed_idx(k)%fm1(m)
			else
				Nb_idx(k)= Nb_idx(k)-1
				Nbl = Nb_idx(k)
			endif	
		 else
!		     Nbl = Nb_idx(k)
		     sed_idx(k)%fd(m,Nbl) = sed_idx(k)%fd(m,Nbl)/Fdall
		 endif	 
		     sed_idx(k)%ft(m) = sed_idx(k)%ft(m)/Ftall
!----revised by Qin 20210924
!			 if(Ftall==0.d0) sed_idx(k)%ft(m)=0.d0
			 if(Ftall==0.d0) then !modified by Qin 20211009
				sed_idx(k)%ft(m)=sed_idx(k)%fd(m,Nbl) 
			 endif	
	!		 if(Fdall==0.d0) sed_idx(k)%fd(m,Nbl)=sed_idx(k)%fd(m,Nbl-1)
			 if(fmall.ne.0.d0)then
		     sed_idx(k)%fm(m) = sed_idx(k)%fm(m)/Fmall
			 else
!			 sed_idx(k)%fm(m)=0.d0
 			 sed_idx(k)%fm(m)= sed_idx(k)%ft(m) !modified by Qin 20211009
			 endif	
	         sed_idx(k)%dmean = sed_idx(k)%dsed(m) * sed_idx(k)%fm(m) + sed_idx(k)%dmean!dmean is calculated by the  bedload layer
	
		! bedload, suspended load GSD; added by Qin 2021/6/25
			 sed_idx(k)%fqbi(m) = sed_idx(k)%qbi(m)/qsb_idx(k)
			 sed_idx(k)%fqsi(m) = sed_idx(k)%qsi(m)/qss_idx(k)
			if (qsb_idx(k)==0.d0) sed_idx(k)%fqbi(m) = 0.d0
			if (qss_idx(k)==0.d0) sed_idx(k)%fqsi(m) = 0.d0
		
!---added by Qin 20210924
			if(isnan(sed_idx(k)%fm(m)))then
				write(*,*) k, m, dzb,sed_idx(k)%fm(m),sed_idx(k)%ft(m),sed_idx(k)%fd(m,Nbl),sed_idx(k)%dzbpr(m),fdall,ftall,fmall
				write(*,*)
				write(*,*)sed_idx(k)%dzbpr(1),sed_idx(k)%dzbpr(2),sed_idx(k)%dzbpr(3),sed_idx(k)%dzbpr(4)
				write(*,*)
				write(*,*)ddt,dlambda,sed_idx(k)%ffd(m)
				stop"fm(m) is nan value"
				write(*,*)"fm(m) is nan value"
				endif
			if(isnan(sed_idx(k)%fm(m)).or.isnan(sed_idx(k)%ft(m)).or.isnan(sed_idx(k)%fd(m,Nbl)))then
			write(*,*) k, m, dzb,sed_idx(k)%fm(m),sed_idx(k)%ft(m),sed_idx(k)%fd(m,Nbl),sed_idx(k)%dzbpr(m),fdall,ftall,fmall, Et_idx(k)+dzb, Ed,Emb_idx(k),Nbl, Nbnew, Nb_idx(k), Nl
			stop"nan value"
			endif	

!surface GSD added by Qin 2021/11/12 !evaluate the fsur by taking the smaple of the surface bed material in depth of D90; D90 of bed surface at previous time step;modified by Qin 20230421 ; 
!set the minimum sampling depth as amp_dep_min given by users to avoid the sampling depth is too small; 20240419 
		  if(sed_idx(k)%Bed_D90 .gt. samp_dep_min) then
		  sed_idx(k)%samp_dep = sed_idx(k)%Bed_D90 
		  else
		  sed_idx(k)%samp_dep = samp_dep_min
		  endif
		  if(De_Emb+Et_idx(k).ge.sed_idx(k)%samp_dep)then
			if(De_Emb.ge.sed_idx(k)%samp_dep)then
		    sed_idx(k)%fsur(m) = sed_idx(k)%fm(m)
			else
			sed_idx(k)%fsur(m)= (Emb_idx(k)*c_ave*sed_idx(k)%fm(m)+sed_idx(k)%ft(m)*(sed_idx(k)%samp_dep-De_Emb)*(1.d0-lambda))/(sed_idx(k)%samp_dep*(1.d0-lambda))	
		  	endif
		  else

		   if(sed_idx(k)%samp_dep.ge. Bed_depth)then
		  	De_depo=Bed_depth- De_Emb -Et_idx(k)
		   else
			De_depo= sed_idx(k)%samp_dep- De_Emb -Et_idx(k)
		   endif	
		   if(de_depo.le.Ed)then
				DNl = 1
				if (de_depo.le.0.) DNl=0
		   else
			if(dmod(de_depo,Ed).le.0.001*Ed)then
				DNl = int(de_depo/Ed)
			else	
				DNl=int(De_depo/Ed)+1
			endif
		   endif

			if (DNL.gt.0)then
			do n = 1, DNl
			if(n==DNl)then	
			qdi = qdi+(1.d0-lambda)*(sed_idx(k)%fd(m,Nbl-n+1)*(1.d0-(dble(DNL)-dble(De_depo/Ed)))*Ed)
			else
			qdi = qdi + (1.d0 - lambda)*sed_idx(k)%fd(m,Nbl-n+1)*Ed
			endif
			enddo	
			endif
			sed_idx(k)%fsur(m) = (Emb_idx(k)*c_ave*sed_idx(k)%fm(m) + Et_idx(k)*sed_idx(k)%ft(m)*(1.d0-lambda)+ qdi)/(sed_idx(k)%samp_dep*(1.d0-lambda))
!			fsur_all = fsur_all + sed_idx(k)%fsur(m)
		  endif		
		    fsur_all = fsur_all + sed_idx(k)%fsur(m)						
	     end do
		 !surface GSD added by Qin 2021/11/12
		 if(fsur_all == 0.d0)then
		  do m = 1, Np
			sed_idx(k)%fsur(m)=sed_idx(k)%fm(m)
		  enddo	
		 else
		  do m = 1, Np
			sed_idx(k)%fsur(m)=sed_idx(k)%fsur(m)/fsur_all
		  enddo	
		 endif
!-------re-evaluating D90,D60; added by Qin 20230421
		 	nn= 0
			f90= 0.d0
		do m = 1,Np	
			nn=nn+1
		    f90= f90 + sed_idx(k)%fsur(m)
!added 20240422
			if(nn>1) then
				if(f90>0.6 .and. (f90 - sed_idx(k)%fsur(m)) .le. 0.6) sed_idx(k)%Bed_D60= sed_idx(k)%dsed(m-1)+ (0.6d0-f90+sed_idx(k)%fsur(m))*(sed_idx(k)%dsed(m)-sed_idx(k)%dsed(m-1))/sed_idx(k)%fsur(m)
			else
				if(f90 > 0.6) sed_idx(k)%Bed_D60 = sed_idx(k)%dsed(m)
			endif
			if(f90>0.9d0) then
			if(nn.le.1)then
			sed_idx(k)%Bed_D90= sed_idx(k)%dsed(m)
			else
			sed_idx(k)%Bed_D90= sed_idx(k)%dsed(m-1)+ (0.9d0-f90+sed_idx(k)%fsur(m))*(sed_idx(k)%dsed(m)-sed_idx(k)%dsed(m-1))/sed_idx(k)%fsur(m)
			endif
!			nn=-9999
!			f90=-9999
			exit
			endif	
		enddo	
!-----				 
      ! else
	   !  do m = 1, Np
		!     Nbl = 1
		!     sed_idx(k)%fd(m,Nbl) = 0.0
		!     sed_idx(k)%ft(m) = 0.0
		!     sed_idx(k)%fm(m) = 0.0
	    !     sed_idx(k)%dmean = 0.0
		!	 sed_idx(k)%fqbi(m) = 0.d0
		!	 sed_idx(k)%fqbi(m) = 0.d0
	    ! end do
	  ! endif
		 
if(abs(dzb).gt.0.1E+3) then   
write(*,'(a,i4,e12.3,i4,f)') 'k/dzb/Nb_idx(k)/Et_idx(k)=', k, dzb, Nb_idx(k), Et_idx(k) 
write(*,'(a)') '    m        dzbr        dsed          fm          ft         fd1         fd2         fd3'
do m = 1, Np
write(*,'(i5,e12.3,6f12.5)') m, sed_idx(k)%dzbpr(m), sed_idx(k)%dsed(m), sed_idx(k)%fm(m), sed_idx(k)%ft(m), sed_idx(k)%fd(m,Nbl), sed_idx(k)%fd(m,Nbl-1), sed_idx(k)%fd(m,Nbl-2)
enddo
!pause'in cal mean_diameter'
endif

!-------------determination of d50----------------------
!                fm_sum222 = 0.0
!                fm_sum111(:) = 0.0
!	     do m = 1, Np
!		fm_sum222 = fm_sum222 + sed_idx(k)%fm(m)
!		fm_sum111(m) = fm_sum222
!	     enddo
!
!	     do m = 1, Np
!	      if(fm_sum111(m).ge.0.5) then
!	        y1 = fm_sum111(m)
!	        y0 = fm_sum111(m-1)
!	        x1 = sed_idx(k)%dsed(m)
!	        x0 = sed_idx(k)%dsed(m-1)
!		slope = (y1-y0)/(x1-x0)
!		sed_idx(k)%dmean = x0 + (0.5 - y0)/slope
!		exit
!	      endif
!	     enddo
!
!
!endif !end if of if(Emb_idx(k).gt.Emc) then

	enddo !end do of do k = 1, riv_count
!pause'end of subroutine mean diameter'

      end subroutine mean_diameter
!------------------------------------------------------------------------------------------------------------------
		subroutine mean_diameter_lin(dzb_temp_lin, sed_lin, qsb_lin, qss_lin)
	     use globals
	     use sediment_mod

	     implicit none

	     type(sed_struct) sed_lin(link_count)
	     real(8) dzb_temp_lin(link_count)
		 real(8) qsb_lin(link_count), qss_lin(link_count)
		 real(8) Etnew, Fmnew, Ftnew, dzb
	     real(8) FM, FT, DZBPR1, Fmall, Ftall, Fdall, Fsur_all, f90
	     integer l, k, m, n, Nbnew, Nbl, nn
		real(8),parameter:: csta = 0.6d0
		real(8) De_depo, De_Emb, qdi, c_ave, Bed_depth
		integer DNl

		c_ave  = csta/2.d0

!$omp parallel do private(k,m,n,Nbnew,Nbl,nn,DNl,Etnew,Fmnew,Ftnew,dzb,FM,FT,DZBPR1, &
!$omp                     Fmall,Ftall,Fdall,Fsur_all,f90,De_depo,De_Emb,qdi,Bed_depth) schedule(static)
        do l = 1, link_count
          k = link_idx_k(l)
		  dzb = dzb_temp_lin(l)
          sed_lin(l)%dmean = 0.0

	      if(dzb.ge.0.) then
			if (Et_lin(l)+dzb .le. Ed) then
			  Etnew = Et_lin(l) + dzb
			  Nbnew = Nb_lin(l)
			else
			  Etnew = Et_lin(l) + dzb - Ed
			  Nbnew = Nb_lin(l) + 1
			endif
		  else
			if(Et_lin(l)+dzb .le. 0.) then
			  Etnew = Ed + Et_lin(l) + dzb
			  Nbnew = Nb_lin(l) - 1
			else
			  Etnew = Et_lin(l) + dzb
			  Nbnew = Nb_lin(l)
			endif
		  endif

		  if(Etnew .le. 0.) Etnew = dabs(dzb)*0.000001
		  if(Nbnew .eq. 0) then
		    Nbnew = Nl
			do m = 1, Np
			  do n = 1, Nl
			    sed_lin(l)%fd(m,n) = sed_lin(l)%fm1(m)
			  end do
			end do
		  elseif(Nbnew .eq. Nl+1) then
		    Nbnew = Nl
			do m = 1, Np
			  do n = 1, Nl
			    if (n.ne.Nl) sed_lin(l)%fd(m,n) = sed_lin(l)%fd(m,n+1)
			  end do
			end do
		  endif

          do m = 1, Np
		    Nbl = Nb_lin(l)
		    FM = sed_lin(l)%fm(m)
		    FT = sed_lin(l)%ft(m)
		    DZBPR1 = sed_lin(l)%dzbpr(m)

	        if (dzb .ge. 0.) then
			  Fmnew = (1.-(1.-lambda)*dzb/(Emb_lin(l)*C_ave))*FM + DZBPR1*(1-lambda)/(Emb_lin(l)*c_ave)
		      if (Et_lin(l)+dzb .le. Ed) then
		        Ftnew = Et_lin(l)/Etnew*FT + dzb/Etnew*FM
		        sed_lin(l)%fd(m,Nbnew) = sed_lin(l)%fd(m,Nbnew)
		      else
		        Ftnew = FM
		        sed_lin(l)%fd(m,Nbnew) = Et_lin(l)/Ed*FT+(1.-Et_lin(l)/Ed)*FM
		      endif
		    else
		      if (Et_lin(l)+dzb .le. 0.) then
		   	    Fmnew = FM + ((1.-lambda)*Et_lin(l)/(Emb_lin(l)*c_ave))*FT &
                         -(1.-lambda)*(Et_lin(l)+dzb)/(Emb_lin(l)*C_ave)*sed_lin(l)%fd(m,Nbl) &
                         + DZBPR1*(1-lambda)/(Emb_lin(l)*c_ave)
			    Ftnew = sed_lin(l)%fd(m,Nbl)
			  else
			    Fmnew = FM-((1.-lambda)*dzb/(Emb_lin(l)*c_ave))*FT + DZBPR1*(1-lambda)/(Emb_lin(l)*c_ave)
			    Ftnew = FT
			  endif
		    endif

		    if(Fmnew.lt.0.0) then
	          sed_lin(l)%fm(m) = 0.0
		    else
	          sed_lin(l)%fm(m) = Fmnew
		    endif

		    if(Ftnew.lt.0.0) then
	          sed_lin(l)%ft(m) = 0.0
		    else
	          sed_lin(l)%ft(m) = Ftnew
		    endif

	        Et_lin(l) = Etnew
	        Nb_lin(l) = Nbnew
          enddo

	      Fmall = 0.
	      Ftall = 0.
	      Fdall = 0.
		  Nbl = Nb_lin(l)
	      do m = 1, Np
		    Fmall = Fmall + sed_lin(l)%fm(m)
		    Ftall = Ftall + sed_lin(l)%ft(m)
		    Fdall = Fdall + sed_lin(l)%fd(m,Nbl)
	      end do

		  fsur_all=0.d0
		  nn = 0
		  f90=0.d0
		  De_Emb = Emb_lin(l)*c_ave*dlambda
		  Bed_depth = zb_riv_idx(k) - zb_roc_idx(k)+De_Emb
		  qdi = 0.d0

	      do m = 1, Np
		    if(Fdall ==0.d0)then
			  if(Nb_lin(l)-1==0.d0) then
			    sed_lin(l)%fd(m,Nbl) = sed_lin(l)%fm1(m)
			  else
			    Nb_lin(l)= Nb_lin(l)-1
			    Nbl = Nb_lin(l)
			  endif
		    else
		      sed_lin(l)%fd(m,Nbl) = sed_lin(l)%fd(m,Nbl)/Fdall
		    endif

		    sed_lin(l)%ft(m) = sed_lin(l)%ft(m)/Ftall
			if(Ftall==0.d0) sed_lin(l)%ft(m)=sed_lin(l)%fd(m,Nbl)
			if(fmall.ne.0.d0)then
		      sed_lin(l)%fm(m) = sed_lin(l)%fm(m)/Fmall
			else
  			  sed_lin(l)%fm(m)= sed_lin(l)%ft(m)
			endif
	        sed_lin(l)%dmean = sed_lin(l)%dsed(m) * sed_lin(l)%fm(m) + sed_lin(l)%dmean

			sed_lin(l)%fqbi(m) = sed_lin(l)%qbi(m)/qsb_lin(l)
			sed_lin(l)%fqsi(m) = sed_lin(l)%qsi(m)/qss_lin(l)
			if (qsb_lin(l)==0.d0) sed_lin(l)%fqbi(m) = 0.d0
			if (qss_lin(l)==0.d0) sed_lin(l)%fqsi(m) = 0.d0

			if(isnan(sed_lin(l)%fm(m)))then
			  write(*,*) k, m, dzb,sed_lin(l)%fm(m),sed_lin(l)%ft(m),sed_lin(l)%fd(m,Nbl),sed_lin(l)%dzbpr(m),fdall,ftall,fmall
			  write(*,*)
			  write(*,*)sed_lin(l)%dzbpr(1),sed_lin(l)%dzbpr(2),sed_lin(l)%dzbpr(3),sed_lin(l)%dzbpr(4)
			  write(*,*)
			  write(*,*)ddt,dlambda,sed_lin(l)%ffd(m)
			  stop"fm(m) is nan value"
			  write(*,*)"fm(m) is nan value"
			endif
			if(isnan(sed_lin(l)%fm(m)).or.isnan(sed_lin(l)%ft(m)).or.isnan(sed_lin(l)%fd(m,Nbl)))then
			  write(*,*) k, m, dzb,sed_lin(l)%fm(m),sed_lin(l)%ft(m),sed_lin(l)%fd(m,Nbl),sed_lin(l)%dzbpr(m),fdall,ftall,fmall, Et_lin(l)+dzb, Ed,Emb_lin(l),Nbl, Nbnew, Nb_lin(l), Nl
			  stop"nan value"
			endif

		    if(sed_lin(l)%Bed_D90 .gt. samp_dep_min) then
		      sed_lin(l)%samp_dep = sed_lin(l)%Bed_D90
		    else
		      sed_lin(l)%samp_dep = samp_dep_min
		    endif
		    if(De_Emb+Et_lin(l).ge.sed_lin(l)%samp_dep)then
			  if(De_Emb.ge.sed_lin(l)%samp_dep)then
		        sed_lin(l)%fsur(m) = sed_lin(l)%fm(m)
			  else
			    sed_lin(l)%fsur(m)= (Emb_lin(l)*c_ave*sed_lin(l)%fm(m)+sed_lin(l)%ft(m)*(sed_lin(l)%samp_dep-De_Emb)*(1.d0-lambda))/(sed_lin(l)%samp_dep*(1.d0-lambda))
		  	  endif
		    else
		      if(sed_lin(l)%samp_dep.ge. Bed_depth)then
		  	    De_depo=Bed_depth- De_Emb -Et_lin(l)
		      else
			    De_depo= sed_lin(l)%samp_dep- De_Emb -Et_lin(l)
		      endif
		      if(de_depo.le.Ed)then
			    DNl = 1
			    if (de_depo.le.0.) DNl=0
		      else
			    if(dmod(de_depo,Ed).le.0.001*Ed)then
				  DNl = int(de_depo/Ed)
			    else
				  DNl=int(De_depo/Ed)+1
			    endif
		      endif

			  if (DNL.gt.0)then
			    do n = 1, DNl
			      if(n==DNl)then
			        qdi = qdi+(1.d0-lambda)*(sed_lin(l)%fd(m,Nbl-n+1)*(1.d0-(dble(DNL)-dble(De_depo/Ed)))*Ed)
			      else
			        qdi = qdi + (1.d0 - lambda)*sed_lin(l)%fd(m,Nbl-n+1)*Ed
			      endif
			    enddo
			  endif
			  sed_lin(l)%fsur(m) = (Emb_lin(l)*c_ave*sed_lin(l)%fm(m) + Et_lin(l)*sed_lin(l)%ft(m)*(1.d0-lambda)+ qdi)/(sed_lin(l)%samp_dep*(1.d0-lambda))
		    endif
		    fsur_all = fsur_all + sed_lin(l)%fsur(m)
	      end do

		  if(fsur_all == 0.d0)then
		    do m = 1, Np
			  sed_lin(l)%fsur(m)=sed_lin(l)%fm(m)
		    enddo
		  else
		    do m = 1, Np
			  sed_lin(l)%fsur(m)=sed_lin(l)%fsur(m)/fsur_all
		    enddo
		  endif

		  nn= 0
		  f90= 0.d0
		  do m = 1,Np
			nn=nn+1
		    f90= f90 + sed_lin(l)%fsur(m)
			if(nn>1) then
			  if(f90>0.6 .and. (f90 - sed_lin(l)%fsur(m)) .le. 0.6) &
                sed_lin(l)%Bed_D60= sed_lin(l)%dsed(m-1)+ (0.6d0-f90+sed_lin(l)%fsur(m)) &
                *(sed_lin(l)%dsed(m)-sed_lin(l)%dsed(m-1))/sed_lin(l)%fsur(m)
			else
			  if(f90 > 0.6) sed_lin(l)%Bed_D60 = sed_lin(l)%dsed(m)
			endif
			if(f90>0.9d0) then
			  if(nn.le.1)then
			    sed_lin(l)%Bed_D90= sed_lin(l)%dsed(m)
			  else
			    sed_lin(l)%Bed_D90= sed_lin(l)%dsed(m-1)+ (0.9d0-f90+sed_lin(l)%fsur(m)) &
                  *(sed_lin(l)%dsed(m)-sed_lin(l)%dsed(m-1))/sed_lin(l)%fsur(m)
			  endif
			  exit
			endif
		  enddo

          Et_idx(k) = Et_lin(l)
          Nb_idx(k) = Nb_lin(l)

          if(abs(dzb).gt.0.1E+3) then
            write(*,'(a,i4,e12.3,i4,f)') 'k/dzb/Nb_idx(k)/Et_idx(k)=', k, dzb, Nb_lin(l), Et_lin(l)
            write(*,'(a)') '    m        dzbr        dsed          fm          ft         fd1         fd2         fd3'
            do m = 1, Np
              write(*,'(i5,e12.3,6f12.5)') m, sed_lin(l)%dzbpr(m), sed_lin(l)%dsed(m), sed_lin(l)%fm(m), sed_lin(l)%ft(m), sed_lin(l)%fd(m,Nbl), sed_lin(l)%fd(m,Nbl-1), sed_lin(l)%fd(m,Nbl-2)
            enddo
          endif
	    enddo
!$omp end parallel do

      end subroutine mean_diameter_lin
!------------------------------------------------------------------------------------------------------------------
! output sediment variation for specified stations
      subroutine sed_output (sed,sed_idx ,t_char)
         use globals
	     use sediment_mod
		 
	     implicit none
		 
	     type(sed_struct) sed(ny,nx)
		 type(sed_struct) sed_idx(riv_count)
         character*6 t_char
         integer i, ios, m, k,n,nn,nbl, maxloc
		 character*256 name(100)
! outfile
         character*256 outfile

! target point !modified 20260208
		 maxloc = 0
		 i = 1
		 do 
		 	if(itar(i)==0.and.jtar(i)==0)then
				exit
			else
				maxloc = i
				write(name(i), '(A,I0)') 'No.', maxloc
				i = i + 1
			end if
		 end do
         
! output
         do i = 1, maxloc
             outfile = trim(sdout_file) // trim(name(i)) // trim(t_char) // ".txt"
			 k = riv_ij2idx(itar(i), jtar(i))
             open(30, file = outfile)
			     write(30,*) "di        fm        ft    fsur      csi             f_qbi         f_qsi"

				 do m = 1, Np

			 write(30,'(1000es10.3)')sed(itar(i),jtar(i))%dsed(m),&
                            sed(itar(i),jtar(i))%fm(m),&
			    sed(itar(i),jtar(i))%ft(m),sed_idx(k)%fsur(m),sed_idx(k)%ssi(m),sed_idx(k)%fqbi(m),sed_idx(k)%fqsi(m)
		        end do
			 write(30,*) 
			 write(30,'(a, es15.7)') 'dmean=', sed(itar(i),jtar(i))%dmean
			 write(30,'(a,es15.7)') 'D90= ', sed(itar(i),jtar(i))%Bed_D90
			 write(30,'(a,es15.7)') 'D60=', sed(itar(i),jtar(i))%Bed_D60
			 write(30,'(a,es10.3)') 'The_sampling_depth_of_bed_surface_material=', sed(itar(i),jtar(i))%samp_dep
		     close(30)
    	 enddo
	  end subroutine sed_output

	  !------------added by Qin---
      subroutine sed_output2(sed,t_char, dzb_temp,sumdzb_idx,hr_idxa,qr_ave_idx,ust_idx,qsb_idx,qss_idx,sumqsb_idx,sumqss_idx)
		use globals
		use sediment_mod
		
		implicit none
		
		type(sed_struct) sed(ny,nx)
		real(8) dzb_temp(riv_count),sumdzb_idx(riv_count)!,zb_riv_slope_idx(riv_count)
		real(8) hr_idxa(riv_count), qr_ave_idx(riv_count)
		real(8) ust_idx(riv_count)
		real(8) qsb_idx(riv_count),qss_idx(riv_count), sumqsb_idx(riv_count),sumqss_idx(riv_count)
		character*6 t_char
		integer i, ios, k 
! outfile
		character*256 outfile

! target point
		character*20 name(100)
		integer target_i(100), target_j(100), maxloc
		
! read location file
		open( 5, file = sdt_location, status = "old")
			i = 1
			do
				read(5,*,iostat=ios) name(i), target_i(i), target_j(i)
				if(ios.ne.0) exit
				i = i + 1
			enddo
			maxloc = i-1
		close(5)
		
! output
		outfile = "./out/sedTar/target_sed_"// trim(t_char) // ".txt"
		open(301, file = outfile)
		write(301,'(a)') 'Target_site                   k        Q           Depth      slope    Zb       U*           Dzb      sumDzb   Dmean        Qsb         Qss         sumQsb      sumQss'
		do i = 1, maxloc
			k = riv_ij2idx(target_i(i), target_j(i))
		    write(301,'(a,i,2f12.3,f10.5,f12.3,f10.5,2f10.2,f10.5,4f12.3)') name(i), k, qr_ave_idx(k),hr_idxa(k),tan(zb_riv_slope_idx(k)*3.14159/180.),zb_riv_idx(k),ust_idx(k), dzb_temp(k) &
				, sumdzb_idx(k), sed(target_i(i),target_j(i))%dmean,qsb_idx(k),qss_idx(k), sumqsb_idx(k),sumqss_idx(k)	
		enddo
		close(301)
	 end subroutine sed_output2	  

	!modified 20240125 20260608
	subroutine sed_output3(sed,t_char, dzb_temp,sumdzb_lin,hr_idxa,qr_ave_idx,ust_idx,qsb_idx,qss_idx,sumqsb_idx,sumqss_idx)
		use globals
		use sediment_mod
		
		implicit none
		
		type(sed_struct) sed(ny,nx)
		real(8) dzb_temp(riv_count)
		real(8) hr_idxa(riv_count), qr_ave_idx(riv_count)
		real(8) ust_idx(riv_count)
		real(8) qsb_idx(riv_count),qss_idx(riv_count), sumqsb_idx(riv_count),sumqss_idx(riv_count)
		real(8) sumdzb_lin(link_count)!,zb_riv_slope_lin(link_count)
		character*6 t_char
        integer i, ios, m, k,n,nn,nbl, l, maxloc
		character*256 name(100)
		logical :: ex
		integer(kind=8) :: fsize

		inquire(file=outfile_sdout, exist=ex, size=fsize)

		open(301, file=outfile_sdout, status='unknown', action='write', position='append', iostat=ios)
		if (ios /= 0) stop 'Failed to open outfile_sdout'

		if (.not. ex .or. fsize == 0_8) then
    		if (j_drf == 0) then
        		write(301,'(a)') 'Target_site     time    i      j      k     l    qr     qb      qs     di        fm       ssi(m=1,Np)(Suspended Sediment concentration)'
    		else
        		write(301,'(a)') 'Target_site    time     i      j     k       l       qr      qb      qs      cw       di       fm      ssi(m=1,Np)(Suspended Sediment concentration)'
    		end if
		end if
		
! count number of target points !modified 20260208
		 maxloc = 0
		 i = 1
		 do 
		 	if(itar(i)==0.and.jtar(i)==0)then
				exit
			else
				maxloc = i
				write(name(i), '(A,I0)') 'No.', maxloc
				i = i + 1
			end if
		 end do
		
! output
		if(j_drf==0)then
			do i = 1, maxloc
				k = riv_ij2idx(itar(i), jtar(i))
				l = link_to_riv(k)
				write(301,'(2a10,4i6,500f15.8)')name(i), t_char, itar(i), jtar(i), k, l, qr_ave_idx(k),qsb_idx(k),qss_idx(k),sed(itar(i),jtar(i))%dmean, (sed(itar(i),jtar(i))%fm(m),m=1,Np), (sed(itar(i),jtar(i))%ssi(m),m=1,Np)
			enddo		
		else
			do i = 1, maxloc
				k = riv_ij2idx(itar(i), jtar(i))
				l = link_to_riv(k)
				write(301,'(2a10,4i6,100f15.8)')name(i),  t_char, itar(i), jtar(i), k, l, qr_ave_idx(k),qsb_idx(k),qss_idx(k),cw(l), sed(itar(i),jtar(i))%dmean, (sed(itar(i),jtar(i))%fm(m),m=1,Np), (sed(itar(i),jtar(i))%ssi(m),m=1,Np)
	
			!write(301,*)
			!write(301,'(a,2i,2f12.3,f10.5,3f12.3,f10.5,f15.8,f10.2,f10.5,4f12.3,f20.2, 20f15.8)') name(i), l,k,qr_ave_idx(k),hr_idxa(k), tan(zb_riv_slope_lin(l)*3.14159/180.),width_lin(l),Link_len(l),zb_riv_idx(k),ust_idx(k), dzb_temp(k) &
			!	, sumdzb_lin(l), sed(target_i(i),target_j(i))%dmean,qsb_idx(k),qss_idx(k),sumqsb_idx(k),sumqss_idx(k),area_lin(l), slo_to_lin_sed_sum(l),(slo_to_lin_sum_di(l,m), m = 1,Np)	
			enddo	
		endif

		close(301)
	 end subroutine sed_output3	 	 	 	  
	   	 	 	  
	  

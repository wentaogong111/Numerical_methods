
#include "paramesh_preprocessor.fh"


	subroutine advance_soln(mype,dt,qin,slope,snm)




!
! This file is a template describing how the solution can be
! initialized on the initial grid. Modify it for your own use.
!
!--------------------------------------------------------------
! include files for amr
		use paramesh_dimensions
		use physicaldata
		use tree

		integer :: mype
		integer :: im1, im2, jm1, jm2
		integer :: ixb1, ixb2, iyb1, iyb2
		real    :: dt
		double precision :: uvp, vup, hsup, hsvp
		double precision :: advection, pressure, roughness
		double precision :: dx, dy
		double precision :: xi, yi

		real, intent(in) :: qin, slope, snm

		Include 'mpif.h'

		integer :: nguard0
		
		integer, parameter :: iiia = 0
		
		double precision, parameter :: g		= 9.81d0
		double precision, parameter :: hmin		= 0.02d0
		double precision, parameter :: hmin2	= hmin*2.

		double precision :: h  (il_bnd:iu_bnd,jl_bnd:ju_bnd) 
		double precision :: wh (il_bnd:iu_bnd,jl_bnd:ju_bnd)
		double precision :: hn (il_bnd:iu_bnd,jl_bnd:ju_bnd) 
		double precision :: eta(il_bnd:iu_bnd,jl_bnd:ju_bnd)
		double precision :: u  (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		double precision :: v  (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		double precision :: wu (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		double precision :: wv (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		double precision :: qu (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		double precision :: qv (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)

!--------------------------------------------------------------

		nguard0 = nguard*npgs
		
		im1 = il_bnd+nguard
		im2 = iu_bnd-nguard
		jm1 = jl_bnd+nguard*k2d
		jm2 = ju_bnd-nguard*k2d

		call amr_timestep(dt,dtmin,dtmax,mype)

! loop over leaf grid blocks
		if( lnblocks.gt.0 ) then
			do lb=1,lnblocks
			
				if( nodetype(lb).eq.1 .or. advance_all_levels ) then

					if( neigh(1,1,lb)<-20 ) then
			!		if( neigh(1,1,lb)<  0 ) then
						ixb1 = 1
					else
						ixb1 = 0
					end if

					if( neigh(1,2,lb)<-20 ) then
			!		if( neigh(1,2,lb)<  0 ) then
						ixb2 = 0
					else
						ixb2 = 1
					end if

					if( neigh(1,3,lb)<-20 ) then
			!		if( neigh(1,3,lb)<  0 ) then
						iyb1 = 1
					else
						iyb1 = 0
					end if

					if( neigh(1,4,lb)<-20 ) then
			!		if( neigh(1,4,lb)<  0 ) then
						iyb2 = 0
					else
						iyb2 = 1
					end if
				
					if(nvar.gt.0) then
					
						k = 1
						
						do j=jl_bnd, ju_bnd
							do i=il_bnd, iu_bnd
								h(i,j)		= unk(1,i,j,k,lb)
								eta(i,j)	= unk(2,i,j,k,lb)	!+unk(3,i,j,k,lb)
								hn(i,j)		= h(i,j)+eta(i,j)
							end do
						end do
						
						do j=jl_bnd, ju_bnd
							do i=il_bnd, iu_bnd+1
								u(i,j) = facevarx(1,i,j,k,lb)
							end do
						end do
								
						do j=jl_bnd, ju_bnd+1
							do i=il_bnd, iu_bnd
								v(i,j) = facevary(1,i,j,k,lb)
							end do
						end do
																		
						dx = bsize(1,lb)/dble(nxb)
						dy = bsize(2,lb)/dble(nyb)
						
						do j=jm1,jm2
!							do i=im1+1,im2
!							do i=im1,im2+1
							do i=im1+ixb1,im2+ixb2
								hsup = (h(i-1,j)+h(i,j))*0.5d0
								
									vup = (v(i-1,j)+v(i,j)+v(i-1,j+1)+v(i,j+1))*0.25d0

!									wu(i,j) = u(i,j) - ( &
!												+g*(-h(i-1,j)+h(i,j))/dx+g*(-eta(i-1,j)+eta(i,j))/dx	&
!												+g*snm**2.d0/hsup*(4./3.)*u(i,j)*sqrt(u(i,j)**2.d0+vup**2.d0)			&
!												+( (u(i,j)+abs(u(i,j)))*(-u(i-1,j)+u(i,j))/dx			&
!												  +(u(i,j)-abs(u(i,j)))*(-u(i,j)+u(i+1,j))/dx			&
!												  +(vup+abs(vup))*(-u(i,j-1)+u(i,j))/dy				&
!												  +(vup-abs(vup))*(-u(i,j)+u(i,j+1))/dy )*0.5d0 )*dt

									advection = -( (u(i,j)+abs(u(i,j)))*(-u(i-1,j)+u(i,j))/dx			&
												  +(u(i,j)-abs(u(i,j)))*(-u(i,j)+u(i+1,j))/dx			&
												  +(vup+abs(vup))*(-u(i,j-1)+u(i,j))/dy					&
												  +(vup-abs(vup))*(-u(i,j)+u(i,j+1))/dy )*0.5d0*dt
									
									roughness = g*snm**2.d0/hsup*(4./3.)*sqrt(u(i,j)**2.d0+vup**2.d0)*dt
									
									if( (eta(i-1,j)>eta(i,j)												&
												.and. h(i-1,j)<=hmin2 .and. hn(i,j)<eta(i-1,j) )	&
										.or.																&
											(eta(i,j)>eta(i-1,j)											&
												.and. h(i,j)<=hmin2 .and. hn(i-1,j)<eta(i,j) ) ) then
				         	  			pressure = 0.d0
				         	  		else
			         	  				pressure = dt*(g*(-hn(i-1,j)+hn(i,j))/dx)
				         	  		end if

									wu(i,j) = (u(i,j)+advection-pressure)/(1.d0+roughness)
												  
									if( h(i-1,j)<=hmin2 .and. wu(i,j)>0.d0 ) wu(i,j) = 0.d0
									if( h(i  ,j)<=hmin2 .and. wu(i,j)<0.d0 ) wu(i,j) = 0.d0
												  
									qu(i,j) = wu(i,j)*hsup
	!								qu(i,j) = ((wu(i,j)+abs(wu(i,j)))*h(i-1,j)+(wu(i,j)-abs(wu(i,j)))*h(i,j))*0.5d0
								
							end do
						end do
						
!						do j=jm1,jm2
!							qu(im1  ,j) = qu(im1+1,j)
!							qu(im2+1,j) = qu(im2  ,j)
!							wu(im1  ,j) = wu(im1+1,j)
!							wu(im2+1,j) = wu(im2  ,j)
!						end do

!						do j=jm1+1,jm2
!						do j=jm1,jm2+1
						do j=jm1+iyb1,jm2+iyb2
							do i=im1,im2
								hsvp = (h(i,j-1)+h(i,j))*0.5d0
								
								
									uvp = (u(i,j-1)+u(i+1,j-1)+u(i,j)+u(i+1,j))*0.25d0
									
!									wv(i,j) = v(i,j) - ( &
!												+g*(-h(i,j-1)+h(i,j))/dy+g*(-eta(i,j-1)+eta(i,j))/dy	&
!												+g*snm**2.d0/hsvp**(4./3.)*v(i,j)*sqrt(v(i,j)**2.d0+uvp**2.d0)			&
!												+( (v(i,j)+abs(v(i,j)))*(-v(i,j-1)+v(i,j))/dy			&
!												  +(v(i,j)-abs(v(i,j)))*(-v(i,j)+v(i,j+1))/dy			&
!												  +(uvp+abs(uvp))*(-v(i-1,j)+v(i,j))/dx				&
!												  +(uvp-abs(uvp))*(-v(i,j)+v(i+1,j))/dx )*0.5d0 )*dt

									advection = -( (v(i,j)+abs(v(i,j)))*(-v(i,j-1)+v(i,j))/dy				&
												  +(v(i,j)-abs(v(i,j)))*(-v(i,j)+v(i,j+1))/dy			&
												  +(uvp+abs(uvp))*(-v(i-1,j)+v(i,j))/dx					&
												  +(uvp-abs(uvp))*(-v(i,j)+v(i+1,j))/dx )*0.5d0*dt
												  
									roughness = g*snm**2.d0/hsvp**(4./3.)*sqrt(v(i,j)**2.d0+uvp**2.d0)*dt
									
									if( (eta(i,j-1)>eta(i,j)								&
												.and. h(i,j-1)<=hmin2 .and. hn(i,j)<eta(i,j-1) )	&
										.or.													&
											(eta(i,j)>eta(i,j-1)								&
												.and. h(i,j)<=hmin2 .and. hn(i,j-1)<eta(i,j) ) ) then
				          	 			pressure = 0.d0
									else
										pressure = dt*(g*(-hn(i,j-1)+hn(i,j))/dy)
									end if

									wv(i,j) = (v(i,j)+advection-pressure)/(1.d0+roughness)
												  
									if( h(i,j-1)<=hmin2 .and. wv(i,j)>0.d0 ) wv(i,j) = 0.d0
									if( h(i,j  )<=hmin2 .and. wv(i,j)<0.d0 ) wv(i,j) = 0.d0
									
									qv(i,j) = wv(i,j)*hsvp
								!	qv(i,j) = ((wv(i,j)+abs(wv(i,j)))*h(i,j-1)+(wv(i,j)-abs(wv(i,j)))*h(i,j))*0.5d0
								
							end do
						end do

	!						do i=im1,im2
	!							qv(i,jm1  ) = qv(i,jm1+1)
	!							qv(i,jm2+1) = qv(i,jm2  )
	!							wv(i,jm1  ) = wv(i,jm1+1)
	!							wv(i,jm2+1) = wv(i,jm2  )
	!						end do

!						if( neigh(1,1,lb)<-20 ) then
						if( neigh(1,1,lb)<  0 ) then
							if( neigh(1,1,lb)==-21) then
								do j=jm1,jm2
									qu(im1,j) = 0.d0
									wu(im1,j) = 0.d0
								end do
							else if( neigh(1,1,lb)==-22 ) then
								do j=jm1,jm2
									yi =  bnd_box(1,2,lb) + dy*(real(j-nguard0)-.5)
							!		if( abs(yi)<5.d0 ) then
										qu(im1,j) = qin
										wu(im1,j) = qu(im1,j)/h(im1+1,j)
							!		else
							!			qu(im1,j) = 0.d0
							!			wu(im1,j) = 0.d0
							!		end if
								end do
							else				! ‚í‚©‚ç‚È‚¢‚Ì‚Å
				!				do j=jm1,jm2
				!					wu(im1,j) = u(im1,j)
				!					qu(im1,j) = ((wu(im1,j)+abs(wu(im1,j)))*h(im1-1,j)+(wu(im1,j)-abs(wu(im1,j)))*h(im1,j))*0.5d0
				!				end do
							end if
						end if

!						if( neigh(1,2,lb)<-20 ) then
						if( neigh(1,2,lb)<  0 ) then
							if( neigh(1,2,lb)==-21 ) then
								do j=jm1,jm2
									qu(im2+1,j) = 0.d0
									wu(im2+1,j) = 0.d0
								end do
							else if( neigh(1,2,lb)==-23 ) then
								do j=jm1,jm2
						!			yi =  bnd_box(1,2,lb) + dy*(real(j-nguard0)-.5)
						!			if( abs(yi)<5.d0 ) then
										qu(im2+1,j) = qu(im2,j)
										wu(im2+1,j) = qu(im2+1,j)/h(im2,j)
						!			else
						!				qu(im2+1,j) = 0.d0
						!				wu(im2+1,j) = 0.d0
						!			end if
								end do
							else
				!				do j=jm1,jm2
				!					wu(im2+1,j) = u(im2+1,j)
				!					qu(im2+1,j) = ((wu(im2+1,j)+abs(wu(im2+1,j)))*h(im2,j)+(wu(im2+1,j)-abs(wu(im2+1,j)))*h(im2+1,j))*0.5d0
				!				end do
							end if
						end if
						
!						if( neigh(1,3,lb)<-20 ) then
						if( neigh(1,3,lb)<  0 ) then
							if( neigh(1,3,lb)==-21 ) then
								do i=im1,im2
									qv(i,jm1) = 0.d0
									wv(i,jm1) = 0.d0
								end do
							else
				!				do i=im1,im2
				!					wv(i,jm1) = v(i,jm1)
				!					qv(i,jm1) = ((wv(i,jm1)+abs(wv(i,jm1)))*h(i,jm1-1)+(wv(i,jm1)-abs(wv(i,jm1)))*h(i,jm1))*0.5d0
				!				end do
							end if
						end if

!						if( neigh(1,4,lb)<-20 ) then
						if( neigh(1,4,lb)<  0 ) then
							if( neigh(1,4,lb)==-21 ) then
								do i=im1,im2
									qv(i,jm2+1) = 0.d0
									wv(i,jm2+1) = 0.d0
								end do
							else
				!				do i=im1,im2
				!					wv(i,jm2+1) = v(i,jm2+1)
				!					qv(i,jm2+1) = ((wv(i,jm2+1)+abs(wv(i,jm2+1)))*h(i,jm2)+(wv(i,jm2+1)-abs(wv(i,jm2+1)))*h(i,jm2+1))*0.5d0
				!				end do
							end if
						end if
						
						do j=jm1,jm2
							do i=im1,im2
								wh(i,j) = h(i,j)-((-qu(i,j)+qu(i+1,j))/dx+(-qv(i,j)+qv(i,j+1))/dy)*dt
!								wh(i,j) = h(i,j) + dt/(dx*dx)* ( h(i+1,j) + h(i-1,j) + &
!     &                   						h(i,j+1) + h(i,j-1) - h(i,j)*4.0 )

								if( wh(i,j)<=hmin ) then
									wh(i,j) = hmin
									if( wu(i  ,j)<0.d0 ) wu(i  ,j) = 0.d0
									if( wu(i+1,j)>0.d0 ) wu(i+1,j) = 0.d0
									if( wv(i,j  )<0.d0 ) wv(i,j  ) = 0.d0
									if( wv(i,j+1)>0.d0 ) wv(i,j+1) = 0.d0
								end if
							end do
						end do

						if( neigh(1,2,lb)==-23 ) then
							do j=jm1,jm2
								wh(im2  ,j) = wh(im2-1,j)
								wu(im2+1,j) = wu(im2  ,j)
							end do
						end if
						
						do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 

							do j=jm1,jm2
								do i=im1,im2
									unk(1,i,j,k,lb) = wh(i,j)
								end do
							end do

							do j=jm1,jm2
								do i=im1,im2+1
									facevarx(1,i,j,k,lb) = wu(i,j)
								end do
							end do

							do j=jm1,jm2+1
								do i=im1,im2
									facevary(1,i,j,k,lb) = wv(i,j)
								end do
							end do
							
						end do
						
					endif
					
				end if

			enddo ! end loop over grid blocks
		endif
				
	end subroutine advance_soln


	subroutine momentum(mype,dt,qin,slope,snm)

		use paramesh_dimensions
		use physicaldata
		use tree
		implicit none

		integer :: mype
		integer :: i, j, k, lb
		integer :: im1, im2, jm1, jm2
		integer :: ixb1, ixb2, iyb1, iyb2
		real    :: dt, dtmax, dtmin
		double precision :: uvp, vup, hsup, hsvp
		double precision :: advection, pressure, roughness
		double precision :: dx, dy
		double precision :: xi, yi

		real, intent(in) :: qin, slope, snm

		Include 'mpif.h'

		integer :: nguard0
		
		double precision, parameter :: g		= 9.81d0
		double precision, parameter :: hmin		= 0.02d0
		double precision, parameter :: hmin2	= hmin*2.

		double precision :: h  (il_bnd:iu_bnd,jl_bnd:ju_bnd) 
		double precision :: wh (il_bnd:iu_bnd,jl_bnd:ju_bnd)
		double precision :: hn (il_bnd:iu_bnd,jl_bnd:ju_bnd) 
		double precision :: eta(il_bnd:iu_bnd,jl_bnd:ju_bnd)
		double precision :: u  (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		double precision :: v  (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		double precision :: wu (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		double precision :: wv (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		double precision :: qu (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		double precision :: qv (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)

!--------------------------------------------------------------

		nguard0 = nguard*npgs
		
		im1 = il_bnd+nguard
		im2 = iu_bnd-nguard
		jm1 = jl_bnd+nguard*k2d
		jm2 = ju_bnd-nguard*k2d

		call amr_timestep(dt,dtmin,dtmax,mype)

! loop over leaf grid blocks
		if( lnblocks>0 ) then
!$omp parallel
!$omp do private(lb,i,j,k,ixb1,ixb2,iyb1,iyb2,h,hn,wh,u,v,qu,qv,wu,wv,eta,dx,dy,uvp,vup,hsup,hsvp,advection,pressure,roughness,xi,yi)
			do lb=1,lnblocks
			
				if( nodetype(lb).eq.1 .or. advance_all_levels ) then

					if( neigh(1,1,lb)<-20 ) then
			!		if( neigh(1,1,lb)<  0 ) then
						ixb1 = 1
					else
						ixb1 = 0
					end if

					if( neigh(1,2,lb)<-20 ) then
			!		if( neigh(1,2,lb)<  0 ) then
						ixb2 = 0
					else
						ixb2 = 1
					end if

					if( neigh(1,3,lb)<-20 ) then
			!		if( neigh(1,3,lb)<  0 ) then
						iyb1 = 1
					else
						iyb1 = 0
					end if

					if( neigh(1,4,lb)<-20 ) then
			!		if( neigh(1,4,lb)<  0 ) then
						iyb2 = 0
					else
						iyb2 = 1
					end if
				
					if( nvar>0 ) then
					
						k = 1
						
						do j=jl_bnd, ju_bnd
							do i=il_bnd, iu_bnd
								h(i,j)		= unk(1,i,j,k,lb)
								eta(i,j)	= unk(2,i,j,k,lb)	!+unk(3,i,j,k,lb)
								hn(i,j)		= h(i,j)+eta(i,j)
							end do
						end do
						
						do j=jl_bnd, ju_bnd
							do i=il_bnd, iu_bnd+1
								u(i,j) = facevarx(1,i,j,k,lb)
							end do
						end do
								
						do j=jl_bnd, ju_bnd+1
							do i=il_bnd, iu_bnd
								v(i,j) = facevary(1,i,j,k,lb)
							end do
						end do
																		
						dx = bsize(1,lb)/dble(nxb)
						dy = bsize(2,lb)/dble(nyb)
						
						do j=jm1,jm2
							do i=im1+ixb1,im2+ixb2
								hsup = (h(i-1,j)+h(i,j))*0.5d0
								
									vup = (v(i-1,j)+v(i,j)+v(i-1,j+1)+v(i,j+1))*0.25d0

!									wu(i,j) = u(i,j) - ( &
!												+g*(-h(i-1,j)+h(i,j))/dx+g*(-eta(i-1,j)+eta(i,j))/dx	&
!												+g*snm**2.d0/hsup*(4./3.)*u(i,j)*sqrt(u(i,j)**2.d0+vup**2.d0)			&
!												+( (u(i,j)+abs(u(i,j)))*(-u(i-1,j)+u(i,j))/dx			&
!												  +(u(i,j)-abs(u(i,j)))*(-u(i,j)+u(i+1,j))/dx			&
!												  +(vup+abs(vup))*(-u(i,j-1)+u(i,j))/dy				&
!												  +(vup-abs(vup))*(-u(i,j)+u(i,j+1))/dy )*0.5d0 )*dt

									advection = -( (u(i,j)+abs(u(i,j)))*(-u(i-1,j)+u(i,j))/dx			&
												  +(u(i,j)-abs(u(i,j)))*(-u(i,j)+u(i+1,j))/dx			&
												  +(vup+abs(vup))*(-u(i,j-1)+u(i,j))/dy					&
												  +(vup-abs(vup))*(-u(i,j)+u(i,j+1))/dy )*0.5d0*dt
									
									roughness = g*snm**2.d0/hsup*(4./3.)*sqrt(u(i,j)**2.d0+vup**2.d0)*dt
									
									if( (eta(i-1,j)>eta(i,j)												&
												.and. h(i-1,j)<=hmin2 .and. hn(i,j)<eta(i-1,j) )	&
										.or.																&
											(eta(i,j)>eta(i-1,j)											&
												.and. h(i,j)<=hmin2 .and. hn(i-1,j)<eta(i,j) ) ) then
				         	  			pressure = 0.d0
				         	  		else
			         	  				pressure = dt*(g*(-hn(i-1,j)+hn(i,j))/dx)
				         	  		end if

									wu(i,j) = (u(i,j)+advection-pressure)/(1.d0+roughness)
												  
									
									qu(i,j) = wu(i,j)*hsup
	!								qu(i,j) = ((wu(i,j)+abs(wu(i,j)))*h(i-1,j)+(wu(i,j)-abs(wu(i,j)))*h(i,j))*0.5d0
							
							end do
						end do
						
						do j=jm1+iyb1,jm2+iyb2
							do i=im1,im2
								hsvp = (h(i,j-1)+h(i,j))*0.5d0
								
									uvp = (u(i,j-1)+u(i+1,j-1)+u(i,j)+u(i+1,j))*0.25d0
									
!									wv(i,j) = v(i,j) - ( &
!												+g*(-h(i,j-1)+h(i,j))/dy+g*(-eta(i,j-1)+eta(i,j))/dy	&
!												+g*snm**2.d0/hsvp**(4./3.)*v(i,j)*sqrt(v(i,j)**2.d0+uvp**2.d0)			&
!												+( (v(i,j)+abs(v(i,j)))*(-v(i,j-1)+v(i,j))/dy			&
!												  +(v(i,j)-abs(v(i,j)))*(-v(i,j)+v(i,j+1))/dy			&
!												  +(uvp+abs(uvp))*(-v(i-1,j)+v(i,j))/dx				&
!												  +(uvp-abs(uvp))*(-v(i,j)+v(i+1,j))/dx )*0.5d0 )*dt

									advection = -( (v(i,j)+abs(v(i,j)))*(-v(i,j-1)+v(i,j))/dy				&
												  +(v(i,j)-abs(v(i,j)))*(-v(i,j)+v(i,j+1))/dy			&
												  +(uvp+abs(uvp))*(-v(i-1,j)+v(i,j))/dx					&
												  +(uvp-abs(uvp))*(-v(i,j)+v(i+1,j))/dx )*0.5d0*dt
												  
									roughness = g*snm**2.d0/hsvp**(4./3.)*sqrt(v(i,j)**2.d0+uvp**2.d0)*dt
									
									if( (eta(i,j-1)>eta(i,j)								&
												.and. h(i,j-1)<=hmin2 .and. hn(i,j)<eta(i,j-1) )	&
										.or.													&
											(eta(i,j)>eta(i,j-1)								&
												.and. h(i,j)<=hmin2 .and. hn(i,j-1)<eta(i,j) ) ) then
				          	 			pressure = 0.d0
									else
										pressure = dt*(g*(-hn(i,j-1)+hn(i,j))/dy)
									end if

									wv(i,j) = (v(i,j)+advection-pressure)/(1.d0+roughness)
												  
									
									
									qv(i,j) = wv(i,j)*hsvp
								!	qv(i,j) = ((wv(i,j)+abs(wv(i,j)))*h(i,j-1)+(wv(i,j)-abs(wv(i,j)))*h(i,j))*0.5d0
							
							end do
						end do

!						if( neigh(1,1,lb)<-20 ) then
						if( neigh(1,1,lb)<  0 ) then
							if( neigh(1,1,lb)==-21) then
								do j=jm1,jm2
									qu(im1,j) = 0.d0
									wu(im1,j) = 0.d0
								end do
							else if( neigh(1,1,lb)==-22 ) then
								do j=jm1,jm2
									yi =  bnd_box(1,2,lb) + dy*(real(j-nguard0)-.5)
							!		if( abs(yi)<5.d0 ) then
										qu(im1,j) = qin
										wu(im1,j) = qu(im1,j)/h(im1+1,j)
							!		else
							!			qu(im1,j) = 0.d0
							!			wu(im1,j) = 0.d0
							!		end if
								end do
							else				! ‚í‚©‚ç‚È‚¢‚Ì‚Å
								do j=jm1,jm2
				!					wu(im1,j) = u(im1,j)
				!					qu(im1,j) = ((wu(im1,j)+abs(wu(im1,j)))*h(im1-1,j)+(wu(im1,j)-abs(wu(im1,j)))*h(im1,j))*0.5d0
								end do
							end if
						end if

!						if( neigh(1,2,lb)<-20 ) then
						if( neigh(1,2,lb)<  0 ) then
							if( neigh(1,2,lb)==-21 ) then
								do j=jm1,jm2
									qu(im2+1,j) = 0.d0
									wu(im2+1,j) = 0.d0
								end do
							else if( neigh(1,2,lb)==-23 ) then
								do j=jm1,jm2
						!			yi =  bnd_box(1,2,lb) + dy*(real(j-nguard0)-.5)
						!			if( abs(yi)<5.d0 ) then
										qu(im2+1,j) = qu(im2,j)
										wu(im2+1,j) = qu(im2+1,j)/h(im2,j)
						!			else
						!				qu(im2+1,j) = 0.d0
						!				wu(im2+1,j) = 0.d0
						!			end if
								end do
							else
								do j=jm1,jm2
				!					wu(im2+1,j) = u(im2+1,j)
				!					qu(im2+1,j) = ((wu(im2+1,j)+abs(wu(im2+1,j)))*h(im2,j)+(wu(im2+1,j)-abs(wu(im2+1,j)))*h(im2+1,j))*0.5d0
								end do
							end if
						end if
						
!						if( neigh(1,3,lb)<-20 ) then
						if( neigh(1,3,lb)<  0 ) then
							if( neigh(1,3,lb)==-21 ) then
								do i=im1,im2
									qv(i,jm1) = 0.d0
									wv(i,jm1) = 0.d0
								end do
							else
								do i=im1,im2
				!					wv(i,jm1) = v(i,jm1)
				!					qv(i,jm1) = ((wv(i,jm1)+abs(wv(i,jm1)))*h(i,jm1-1)+(wv(i,jm1)-abs(wv(i,jm1)))*h(i,jm1))*0.5d0
								end do
							end if
						end if

!						if( neigh(1,4,lb)<-20 ) then
						if( neigh(1,4,lb)<  0 ) then
							if( neigh(1,4,lb)==-21 ) then
								do i=im1,im2
									qv(i,jm2+1) = 0.d0
									wv(i,jm2+1) = 0.d0
								end do
							else
								do i=im1,im2
				!					wv(i,jm2+1) = v(i,jm2+1)
				!					qv(i,jm2+1) = ((wv(i,jm2+1)+abs(wv(i,jm2+1)))*h(i,jm2)+(wv(i,jm2+1)-abs(wv(i,jm2+1)))*h(i,jm2+1))*0.5d0
								end do
							end if
						end if
						
						do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 

							do j=jm1,jm2
								do i=im1,im2+1
									facevarx(1,i,j,k,lb) = wu(i,j)
									facevarx(2,i,j,k,lb) = qu(i,j)
								end do
							end do

							do j=jm1,jm2+1
								do i=im1,im2
									facevary(1,i,j,k,lb) = wv(i,j)
									facevary(2,i,j,k,lb) = qv(i,j)
								end do
							end do
							
						end do
						
					endif
					
				end if

			enddo ! end loop over grid blocks
!$omp end do
!$omp end parallel
		endif
				
	end subroutine momentum


	subroutine continuity(mype,dt,qin,slope,snm)

		use paramesh_dimensions
		use physicaldata
		use tree
		implicit none

		integer :: mype
		integer :: i, j, k, lb
		integer :: im1, im2, jm1, jm2
		integer :: ixb1, ixb2, iyb1, iyb2
		real    :: dt, dtmax, dtmin
		double precision :: uvp, vup, hsup, hsvp
		double precision :: dx, dy
		double precision :: xi, yi

		real, intent(in) :: qin, slope, snm

		Include 'mpif.h'

		integer :: nguard0
		
		double precision, parameter :: g		= 9.81d0
		double precision, parameter :: hmin		= 0.02d0
		double precision, parameter :: hmin2	= hmin*2.

		double precision :: h  (il_bnd:iu_bnd,jl_bnd:ju_bnd) 
		double precision :: wh (il_bnd:iu_bnd,jl_bnd:ju_bnd)
		double precision :: u  (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		double precision :: v  (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		double precision :: qu (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		double precision :: qv (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)

!--------------------------------------------------------------

		nguard0 = nguard*npgs
		
		im1 = il_bnd+nguard
		im2 = iu_bnd-nguard
		jm1 = jl_bnd+nguard*k2d
		jm2 = ju_bnd-nguard*k2d

! loop over leaf grid blocks
		if( lnblocks>0 ) then
!$omp parallel
!$omp do private(lb,i,j,k,h,wh,u,v,qu,qv,dx,dy)
			do lb=1,lnblocks
			
				if( nodetype(lb).eq.1 .or. advance_all_levels ) then
				
					if( nvar>0 ) then
					
						k = 1
						
						do j=jl_bnd, ju_bnd
							do i=il_bnd, iu_bnd
								h(i,j)		= unk(1,i,j,k,lb)
							end do
						end do
						
						do j=jl_bnd, ju_bnd
							do i=il_bnd, iu_bnd+1
								u (i,j) = facevarx(1,i,j,k,lb)
								qu(i,j) = facevarx(2,i,j,k,lb)
							end do
						end do
								
						do j=jl_bnd, ju_bnd+1
							do i=il_bnd, iu_bnd
								v(i,j)  = facevary(1,i,j,k,lb)
								qv(i,j) = facevary(2,i,j,k,lb)
							end do
						end do
																		
						dx = bsize(1,lb)/dble(nxb)
						dy = bsize(2,lb)/dble(nyb)
						
						do j=jm1,jm2
							do i=im1,im2
								wh(i,j) = h(i,j)-((-qu(i,j)+qu(i+1,j))/dx+(-qv(i,j)+qv(i,j+1))/dy)*dt

								if( wh(i,j)<=hmin ) then
									wh(i,j) = hmin
									if( u(i  ,j)<0.d0 ) u(i  ,j) = 0.d0
									if( u(i+1,j)>0.d0 ) u(i+1,j) = 0.d0
									if( v(i,j  )<0.d0 ) v(i,j  ) = 0.d0
									if( v(i,j+1)>0.d0 ) v(i,j+1) = 0.d0
								end if
							end do
						end do

						if( neigh(1,2,lb)==-23 ) then
							do j=jm1,jm2
								wh(im2  ,j) = wh(im2-1,j)
								u (im2+1,j) =  u(im2  ,j)
							end do
						end if
						
						do k=kl_bnd+nguard*k3d,ku_bnd-nguard*k3d 

							do j=jm1,jm2
								do i=im1,im2
									unk(1,i,j,k,lb) = wh(i,j)
								end do
							end do

							do j=jm1,jm2
								do i=im1,im2+1
									facevarx(1,i,j,k,lb) = u(i,j)
								end do
							end do

							do j=jm1,jm2+1
								do i=im1,im2
									facevary(1,i,j,k,lb) = v(i,j)
								end do
							end do
							
						end do
						
					endif
					
				end if

			enddo ! end loop over grid blocks
!$omp end do
!$omp end parallel
		endif
				
	end subroutine continuity

	subroutine bed_evolution(mype,dt,diam,spec,snm)

		use paramesh_dimensions
		use physicaldata
		use tree

		integer :: mype
		integer :: im1, im2, jm1, jm2
		integer :: ixb1, ixb2, iyb1, iyb2
		real    :: dt
		real	:: up, vp
		real	:: dzdx, dzdy, dzdn, cosc, sinc, qbsc, tc, gam, qbnc, dex

		real, intent(in) :: diam, spec, snm

		Include 'mpif.h'

		integer :: nguard0
				
		real, parameter :: g		= 9.81d0
		real, parameter :: hmin		= 0.01d0
		real, parameter :: hmin2	= hmin*2.
		
		real, parameter :: tsc		= 0.05
		real, parameter :: muk		= 0.54
		real, parameter :: poro		= 0.4

		real :: h  (il_bnd:iu_bnd,jl_bnd:ju_bnd)
		real :: eta(il_bnd:iu_bnd,jl_bnd:ju_bnd)
		real :: u  (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		real :: v  (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)

		real :: vv  (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		real :: ts  (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		real :: qb  (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		real :: qbx (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		real :: qby (il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		real :: ssin(il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)
		real :: scos(il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1)

!--------------------------------------------------------------

		nguard0 = nguard*npgs
		
		im1 = il_bnd+nguard
		im2 = iu_bnd-nguard
		jm1 = jl_bnd+nguard*k2d
		jm2 = ju_bnd-nguard*k2d

! loop over leaf grid blocks
		if( lnblocks>0 ) then
			do lb=1,lnblocks
			
				if( nodetype(lb)==1 .or. advance_all_levels ) then
				
					if( neigh(1,1,lb)<-20 ) then
						ixb1 = 1
					else
						ixb1 = 0
					end if

					if( neigh(1,2,lb)<-20 ) then
						ixb2 = 0
					else
						ixb2 = 1
					end if

					if( neigh(1,3,lb)<-20 ) then
						iyb1 = 1
					else
						iyb1 = 0
					end if

					if( neigh(1,4,lb)<-20 ) then
						iyb2 = 0
					else
						iyb2 = 1
					end if
				
					if( nvar>0 ) then
					
						k = 1
						
						do j=jl_bnd, ju_bnd
							do i=il_bnd, iu_bnd
								h(i,j)		= unk(1,i,j,k,lb)
								eta(i,j)	= unk(2,i,j,k,lb)+unk(3,i,j,k,lb)
							end do
						end do
						
						do j=jl_bnd, ju_bnd
							do i=il_bnd, iu_bnd+1
								u(i,j) = facevarx(1,i,j,k,lb)
							end do
						end do
								
						do j=jl_bnd, ju_bnd+1
							do i=il_bnd, iu_bnd
								v(i,j) = facevary(1,i,j,k,lb)
							end do
						end do
																		
						dx = bsize(1,lb)/real(nxb)
						dy = bsize(2,lb)/real(nyb)
						
						do j=jl_bnd, ju_bnd
							do i=il_bnd, iu_bnd
								up = (u(i,j)+u(i+1,j))*0.5d0
								vp = (v(i,j)+v(i,j+1))*0.5d0
								
								if( h(i,j)<=hmin2 ) then
									vv(i,j)		= 0.d0
									ts(i,j)		= 0.d0
									ssin(i,j)	= 0.d0
									scos(i,j)	= 0.d0
									qb(i,j)		= 0.d0
								else
									vv(i,j)		= sqrt(up**2.d0+vp**2.d0)
									ts(i,j)		= (snm*vv(i,j))**2.d0/(spec*diam*h(i,j)**(1./3.))
									ssin(i,j)	= vp/vv(i,j)
									scos(i,j)	= up/vv(i,j)
									
									if( ts(i,j)>tsc ) then
										qb(i,j) = 8.d0*(ts(i,j)-tsc)**1.5*sqrt(spec*g*diam**3.)
									else
										qb(i,j) = 0.d0
									end if
								end if
							end do
						end do
						
						do j=jm1,jm2
							do i=im1+ixb1,im2+ixb2
								dzdx = (-eta(i-1,j)+eta(i,j))/dx
								dzdy = (-(eta(i-1,j-1)+eta(i,j-1))	&
										+(eta(i-1,j+1)+eta(i,j+1)))/dy*0.25d0
										
								cosc	= (scos(i-1,j)+scos(i,j))*0.5
								sinc	= (ssin(i-1,j)+ssin(i,j))*0.5
								dzdn	= -sinc*dzdx+cosc*dzdy
								qbsc	= (qb(i-1,j)+qb(i,j))*0.5d0
								tc		= (ts(i-1,j)+ts(i,j))*0.5d0
								
								if( tc>tsc ) then
									qbnc = -qbsc*sqrt(tsc/tc)/muk*dzdn
									
									qbx(i,j) = cosc*qbsc-sinc*qbnc
								else
									qbx(i,j) = 0.d0
								end if
								
								if( h(i-1,j)<=hmin2 .and. qbx(i,j)<0. ) qbx(i,j) = 0.d0
								if( h(i  ,j)<=hmin2 .and. qbx(i,j)>0. ) qbx(i,j) = 0.d0
								
							end do
						end do

						do j=jm1+iyb1,jm2+iyb2
							do i=im1,im2
								dzdx = (-(eta(i-1,j-1)+eta(i-1,j))			&
										+(eta(i+1,j-1)+eta(i+1,j)))/dx*0.25d0
								dzdy = (-eta(i,j-1)+eta(i,j))/dy
										
								cosc	= (scos(i,j-1)+scos(i,j))*0.5
								sinc	= (ssin(i,j-1)+ssin(i,j))*0.5
								dzdn	= -sinc*dzdx+cosc*dzdy
								qbsc	= (qb(i,j-1)+qb(i,j))*0.5d0
								tc		= (ts(i,j-1)+ts(i,j))*0.5d0
								
								if( tc>tsc ) then
									qbnc = -qbsc*sqrt(tsc/tc)/muk*dzdn
									
									qby(i,j) = sinc*qbsc+cosc*qbnc
								else
									qby(i,j) = 0.d0
								end if
								
								if( h(i,j-1)<=hmin2 .and. qby(i,j)<0. ) qby(i,j) = 0.d0
								if( h(i,j  )<=hmin2 .and. qby(i,j)>0. ) qby(i,j) = 0.d0
							
							end do
						end do

						if( neigh(1,1,lb)<-20 ) then
							if( neigh(1,1,lb)==-21) then
								do j=jm1,jm2
									qbx(im1,j) = 0.d0
								end do
							else if( neigh(1,1,lb)==-22 ) then
								do j=jm1,jm2
						!			yi =  bnd_box(1,2,lb) + dy*(real(j-nguard0)-.5)
						!			if( abs(yi)<5.d0 ) then
										qbx(im1,j) = qbx(im1+1,j)
						!			else
						!				qbx(im1,j) = 0.d0
						!			end if
								end do
							end if
						end if

						if( neigh(1,2,lb)<-20 ) then
							if( neigh(1,2,lb)==-21 ) then
								do j=jm1,jm2
									qbx(im2+1,j) = 0.d0
								end do
							else if( neigh(1,2,lb)==-23 ) then
								do j=jm1,jm2
						!			yi =  bnd_box(1,2,lb) + dy*(real(j-nguard0)-.5)
						!			if( abs(yi)<5.d0 ) then
										qbx(im2+1,j) = qbx(im2,j)
						!			else
						!				qu(im2+1,j) = 0.d0
						!				wu(im2+1,j) = 0.d0
						!			end if
								end do
							end if
						end if

						if( neigh(1,3,lb)<-20 ) then
							do i=im1,im2
								qby(i,jm1) = 0.d0
							end do
						end if

						if( neigh(1,4,lb)<-20 ) then
							do i=im1,im2
								qby(i,jm2+1) = 0.d0
							end do
						end if
						
						do j=jm1,jm2
							do i=im1,im2
								dex = -((-qbx(i,j)+qbx(i+1,j))/dx+(-qby(i,j)+qby(i,j+1))/dy)*dt/(1.d0-poro)
								
								eta(i,j) = eta(i,j)+dex
								
								unk(3,i,j,1,lb) = unk(3,i,j,1,lb)+dex
							end do
						end do

						if( neigh(1,2,lb)==-23 ) then
							do j=jm1,jm2
								unk(3,im2,j,1,lb) = unk(3,im2-1,j,1,lb)
							end do
						end if
						
					endif
					
				end if

			enddo ! end loop over grid blocks
		endif
				
	end subroutine bed_evolution
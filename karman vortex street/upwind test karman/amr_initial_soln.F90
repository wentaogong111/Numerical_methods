
#include "paramesh_preprocessor.fh"



	subroutine amr_initial_soln(qin,slope,snm)





!
! This file is a template describing how the solution can be
! initialized on the initial grid. Modify it for your own use.
!
!--------------------------------------------------------------
! include files for amr
		use paramesh_dimensions
		use physicaldata
		use tree

		Include 'mpif.h'

		integer :: nguard0

		real, intent(in)	:: qin, slope, snm
		real				:: hs0

		if( slope<=0.00000001 ) then
			hs0 = 1.
		else
			hs0 = (snm*qin/slope**0.5)**0.6
		end if

!--------------------------------------------------------------


		nguard0 = nguard*npgs

! loop over leaf grid blocks
		if(lnblocks.gt.0) then
		
			do lb=1,lnblocks
      
				if(nodetype(lb).eq.1 .or. advance_all_levels) then

					if(nvar.gt.0) then

! set values for unk
						dx = bsize(1,lb)/real(nxb)
						dy = bsize(2,lb)/real(nyb)

!		k = 1

						do k=kl_bnd+nguard0*k3d,ku_bnd-nguard0*k3d
							do j=jl_bnd+nguard0*k2d,ju_bnd-nguard0*k2d
								do i=il_bnd+nguard0,iu_bnd-nguard0
			!						unk(1,i,j,k,lb) = 5.0
							
									xi =  bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5) !+25.
									yi =  bnd_box(1,2,lb) + dy*(real(j-nguard0)-.5) !+10.
			              
									unk(2,i,j,k,lb) = slope*(g_xmax-xi)

									if( abs(xi)<2. .and. abs(yi)<2. ) then
										unk(2,i,j,k,lb) = unk(2,i,j,k,lb)+hs0*2.d0
										unk(1,i,j,k,lb) = 0.001
									else
										unk(1,i,j,k,lb) = hs0
									endif
			              
							!		unk(1,i,j,k,lb) = hs0	!- unk(2,i,j,k,lb)
			
									unk(3,i,j,k,lb) = 0.d0
		              
								enddo
							enddo
						enddo

					endif

					if(nvarcorn.gt.0) then

! set values for unk_n
						do k=kl_bnd+nguard0*k3d,ku_bnd+(nguard0+1)*k3d
							do j=jl_bnd+nguard0*k2d,ju_bnd+(nguard0+1)*k2d
								do i=il_bnd+nguard0,iu_bnd+(nguard0+1)
			!						unk_n(1,i,j,k,lb) = ???
			!						unk_n(2,i,j,k,lb) = ???
								enddo
							enddo
						enddo

					endif


					if(nfacevar.gt.0) then

!						k = 1

! set values for facevarx
						do k=kl_bnd+nguard0*k3d,ku_bnd+nguard0*k3d
							do j=jl_bnd+nguard0*k2d,ju_bnd+nguard0*k2d
								do i=il_bnd+nguard0,iu_bnd+nguard0+1
									hsup = (unk(1,i-1,j,k,lb)+unk(1,i,j,k,lb))*0.5
									if( hsup<=0.01 ) then
										facevarx(1,i,j,k,lb) = 0.d0
									else
										facevarx(1,i,j,k,lb) = hs0**(2./3.)*slope**0.5/snm
									end if
									facevarx(2,i,j,k,lb) = 0.d0
								enddo
							enddo
						enddo

! set values for facevary
						do k=kl_bnd+nguard0*k3d,ku_bnd+nguard0*k3d
							do j=jl_bnd+nguard0*k2d,ju_bnd+(nguard0+1)*k2d
								do i=il_bnd+nguard0,iu_bnd+nguard0
									facevary(1,i,j,k,lb) = 0.d0
									facevary(2,i,j,k,lb) = 0.d0
								enddo
							enddo
						enddo

! set values for facevarz
						do k=kl_bnd+nguard0*k3d,ku_bnd+(nguard0+1)*k3d
							do j=jl_bnd+nguard0*k2d,ju_bnd+nguard0*k2d
								do i=il_bnd+nguard0,iu_bnd+nguard0
									facevarz(1,i,j,k,lb) = 0.d0
								enddo
							enddo
						enddo

					endif

					if(nvaredge.gt.0) then

! set values for unk_e_x
						do k=kl_bnd+nguard0*k3d,ku_bnd+(nguard0+1)*k3d
							do j=jl_bnd+nguard0*k2d,ju_bnd+(nguard0+1)*k2d
								do i=il_bnd+nguard0,iu_bnd+nguard0
				! 					unk_e_x(1,i,j,k,lb) = ???
				!					unk_e_x(2,i,j,k,lb) = ???
								enddo
							enddo
						enddo

! set values for unk_e_y
						do k=kl_bnd+nguard0*k3d,ku_bnd+(nguard0+1)*k3d
							do j=jl_bnd+nguard0*k2d,ju_bnd+nguard0*k2d
								do i=il_bnd+nguard0,iu_bnd+(nguard0+1)
				!					unk_e_y(1,i,j,k,lb) = ???
				!					unk_e_y(2,i,j,k,lb) = ???
								enddo
							enddo
						enddo

! set values for unk_e_z
						do k=kl_bnd+nguard0*k3d,ku_bnd+nguard0*k3d
							do j=jl_bnd+nguard0*k2d,ju_bnd+(nguard0+1)*k2d
								do i=il_bnd+nguard0,iu_bnd+(nguard0+1)
				!					unk_e_z(1,i,j,k,lb) = ???
				!					unk_e_z(2,i,j,k,lb) = ???
								enddo
							enddo
						enddo

					endif

				endif

			enddo ! end loop over grid blocks
		endif

      return
      end subroutine amr_initial_soln

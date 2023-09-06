!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!
! This is a template main program for the multidimensional PARAMESH amr 
! package. It is suitable for use when permanent guardcell storage is 
! to be used for every grid block.
!


!----------------------------------------------------------------
!
!
! BLOCK 1



#include "paramesh_preprocessor.fh"

! include file to define physical qualities of the model and mesh
        use paramesh_dimensions
        use physicaldata

! include file defining the tree
        use tree


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE: IT IS EXTREMELY IMPORTANT TO EXPOSE ALL THE INTERFACES OF
! PARAMESH SUBROUTINES USED.
! PARAMESH3.0 WILL NOT WORK CORRECTLY WITHOUT DOING THIS.
! mpi-support, subroutine interfaces used by PARAMESH can be found in
! paramesh_mpi_interfaces, but none are needed for this example.

! identify interface blocks for routines called in this program
        use paramesh_interfaces, only : comm_start, &
     &                                  amr_initialize, &
     &                                  amr_refine_derefine, &
     &                                  amr_guardcell, &
     &                                  amr_prolong, &
     &                                  amr_restrict, &
     &                                  amr_test_refinement, &
     &                                  amr_close, &
     &                                  amr_checkpoint_wr

        Include 'mpif.h'

! local amr variables
        integer           :: nprocs, mype, ierr
        integer	:: three, four, five, six
        integer :: no_of_blocks, loop_count_max
        
        integer :: nguard0
        
        real :: tuk, optime
        real :: g_xmin, g_xmax, g_ymin, g_ymax, g_zmin, g_zmax

		real :: qin, slope, snm
		
		real :: diam, spec
		
		integer :: remote_block, remote_pe, iblk, cnodetype
        
        tuk = 0.5d0

		snm		= 0.03
		slope	= 0.001
		qin		= 1.d0
		
		diam	= 0.0001
		spec	= 1.65

!---------------------------------------------------------------
!
! BLOCK 2


! make initialization call for amr package

        call amr_initialize

	Call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
        Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

        if(mype.eq.0) write(*,*) 'Running on ',nprocs,' processors'

!---------------------------------------------------------------
!
! BLOCK 3


! set a limit on the refinement level
        lrefine_min = 1
!        lrefine_max = 7
        lrefine_max = 8
        
        no_of_blocks = 2**(lrefine_min-1)


! set up initial grid state.

! set up a single block covering the whole square domain
        lnblocks = 0
        if(mype.eq.0.) then
                lnblocks = 1
                coord(:,1) = 0.
				bnd_box(1,:,1) = -50.
				bnd_box(2,:,1) =  50.
                bsize(:,1) = bnd_box(2,:,1) - bnd_box(1,:,1)
!				bnd_box(1,1,1) = -100.
!				bnd_box(2,1,1) = 100.
!				bnd_box(1,2,1) = -200.
!				bnd_box(2,2,1) = 200.
!                bsize(:,1) = bnd_box(2,:,1) - bnd_box(1,:,1)
                nodetype(1) = 1
                lrefine(1) = 1

!                neigh(1,:,1) = 1                  ! initial block is its own
!                neigh(2,:,1) = 0                  ! neighbor ie periodic bc's
                neigh(1,:,1) = -21
				neigh(1,1,1) = -22
				neigh(1,2,1) = -23
                neigh(2,:,1) = 0		! Š„‚è“–‚½‚Á‚Ä‚éƒvƒƒZƒX”Ô†??

                refine(1)=.true.
        endif
        
		g_xmin = bnd_box(1,1,1)
		g_xmax = g_xmin+bsize(1,1)
		g_ymin = bnd_box(1,2,1)
		g_ymax = g_ymin+bsize(2,1)
		g_zmin = bnd_box(1,3,1)
		g_zmax = g_ymin+bsize(3,1)

		grid_xmin = g_xmin
		grid_xmax = g_xmax
		grid_ymin = g_ymin
		grid_ymax = g_ymax
		grid_zmin = g_zmin
		grid_zmax = g_zmax

		write(*,*) g_xmin, g_xmax, g_ymin, g_ymax, g_zmin, g_zmax
		write(*,*) grid_xmin, grid_xmax, grid_ymin, grid_ymax, grid_zmin, grid_zmax

		boundary_index = -21
!		boundary_index = 1

! x boundaries
		boundary_box(1,2:3,1:2) = -1.e10
		boundary_box(2,2:3,1:2) =  1.e10
		boundary_box(1,1,1) = -1.e10
		boundary_box(2,1,1) = g_xmin
		boundary_box(1,1,2) = g_xmax
		boundary_box(2,1,2) = 1.e10

! y boundaries
		if(ndim>=2) then
			three = (2*k2d) + 1
			four  = three + k2d
			boundary_box(1,1,three:four) = -1.e10
			boundary_box(2,1,three:four) =  1.e10
			boundary_box(1,3,three:four) = -1.e10
			boundary_box(2,3,three:four) =  1.e10
			boundary_box(1,2,three) = -1.e10
			boundary_box(2,2,three) = g_ymin
			boundary_box(1,2,four) = g_ymax
			boundary_box(2,2,four) = 1.e10
		endif

! z boundaries
		if(ndim.eq.3) then
			five = (4*k3d) + 1
			six  = five + k3d
			boundary_box(1,1:2,five:six) = -1.e10
			boundary_box(2,1:2,five:six) =  1.e10
			boundary_box(1,3,five) = -1.e10
			boundary_box(2,3,five) = g_zmin
			boundary_box(1,3,six) = g_zmax
			boundary_box(2,3,six) = 1.e10
		endif

		Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! Now cycle over blocks, refining all existing leaf blocks.

		loop_count_max = int(log(real(no_of_blocks))/log(2.)+1.)

        do loop_count=1,4	!loop_count_max

                refine(1:lnblocks) = .true.

! refine grid and apply morton reordering to grid blocks if necessary
                call amr_refine_derefine

        enddo


! Output details of the current grid
	if(mype.eq.0) write(*,*) 'pe / blk / blk-coords / blk-sizes'
        do l=1,lnblocks
        	if( nodetype(l)==1 ) then
		        write(*,100) mype,l,(coord(i,l),i=1,ndim),(bsize(j,l),j=1,ndim),neigh(1,1,l),neigh(1,2,l),neigh(1,3,l),neigh(1,4,l)
		    end if
        enddo
100     format(i2,1x,i4,4(2x,f13.7),4i5)

!		write(*,*) lperiodicx, lperiodicy

!---------------------------------------------------------------
!
! BLOCK 4

! Establish solution on the initial grid.

	call amr_initial_soln(qin,slope,snm)

!---------------------------------------------------------------
!
! BLOCK 5


! exchange guardcell information
        iopt = 1

        nlayers = nguard
        call amr_guardcell(mype,iopt,nlayers)


!        do lb=1,lnblocks
!          if(coord(1,lb).eq.5.0.and.coord(2,lb).eq.5.0) then
!            do j=1,nyb+2*nguard
!            write(*,50) j,(unk(1,i,j,1,lb),i=1,nxb+2*nguard)
!50          format(1x,i3,6(2x,f11.8))
!            enddo
!          endif
!        enddo
        

!----------------------------------------------------------
!
! BLOCK 6

        minstp = 1
        maxstp = 250
!        maxstp = 750
        maxstp = 30000
        

! loop over time steps.

!		istep = 0
!        call paramesh2paraview( istep )

!        do istep = minstp, maxstp

		istep	= 0
        time	= 0.d0
        optime	= 0.d0
        
        nguard0 = nguard*npgs
        
		call amr_timestep(dt,dtmin,dtmax,mype)
		
		write(*,*) nfluxes, nfaces
		write(*,*) il_bnd, iu_bnd, jl_bnd, ju_bnd, kl_bnd, ku_bnd
		write(*,*) il_bndi, iu_bndi, jl_bndi, ju_bndi, kl_bndi, ku_bndi
		
		do

! advance the solution using a User provided routine
!			call advance_soln(mype,dt,qin,slope,snm)

			call momentum(mype,dt,qin,slope,snm)
			
			call continuity(mype,dt,qin,slope,snm)
            	call momentum(mype,dt,qin,slope,snm)
			
			call continuity(mype,dt,qin,slope,snm)
            	call momentum(mype,dt,qin,slope,snm)
			
			call continuity(mype,dt,qin,slope,snm)
			
!			call bed_evolution(mype,dt,diam,spec,snm)

        	if (.not.advance_all_levels) then
! A valid solution will be required on the parents of leaf blocks
! when refinement testing is done. See the comment before the call
! to amr_test_refinement.
        		iempty = 0 
        		iopt = 1
        		call amr_restrict(mype,iopt,iempty)
        	endif
!
! test to see if additional refinement or derefinement is necessary
! note - if the solution is only being advanced on leaf nodes then
! a call to amr_restrict must come before this call to ensure that
! the refinement test can be done on parents of leafs also. This avoids
! a potential refinement/derefinement flip-flop happenning on successive
! timesteps.

        	call amr_test_refinement(mype,lrefine_min,lrefine_max)

! refine grid and apply morton reordering to grid blocks if necessary
        	call amr_refine_derefine

! prolong solution to any new leaf blocks if necessary
        	call amr_prolong(mype,iopt,nlayers)

! exchange guardcell information
        	call amr_guardcell(mype,iopt,nlayers)

!        	if(mod(istep,10).eq.0) then
			if( optime>tuk .or. istep==0 ) then
        		call paramesh2paraview( istep )
        		write(*,*) 'Time= ', time, dt
        		optime = optime-tuk
        		istep = istep+1
        	end if

!        if(mype.eq.0) write(*,*) 'iteration ',istep, &
!     &                           ' no of blocks = ',lnblocks

			if( istep>200 ) exit
			
			time	= time+dt
			optime	= optime+dt

        end do
! END OF BLOCK 6
!---------------------------------------------------------------
 
        if(mype.eq.0) write(*,*) 'pe / blk / blk-coords / blk-sizes'
        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        do l=1,lnblocks
        	if( nodetype(l)==1 ) then
		        write(*,100) mype,l,(coord(i,l),i=1,ndim),(bsize(j,l),j=1,ndim),neigh(1,1,l),neigh(1,2,l),neigh(1,3,l),neigh(1,4,l)
		    end if
        enddo
500     format(i2,1x,i4,4(2x,f13.7),4i5)

        call amr_checkpoint_wr(24,.false.)

        call amr_close

        stop
        end

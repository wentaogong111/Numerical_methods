!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_test_refinement(mype,llrefine_min,llrefine_max)



!--------------------------------------------------------------------------
! 
! This is a template to assist in constructing the routine AMR_TEST_REFINEMENT
! for use in your application. In this illustration we use the workspace
! array WORK to store the data which is used in computing the error measure
! at each grid point. This gives us the freedom to extend the testing
! beyond the normal bounds of individual blocks, since WORK is declared
! with NGUARD_WORK guard cells at each boundary, which can be set to a
! larger number than NGUARD.

! Arguments:
!   mype           integer      local processor number
!   llrefine_min   integer      minimum refinement level to be permitted
!   llrefine_max   integer      maximum refinement level to be permitted

!--------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      use paramesh_interfaces, only : amr_guardcell

      integer, intent(in)    ::  mype,llrefine_min,llrefine_max

      Include 'mpif.h'

      real :: error(ilw:iuw,jlw:juw,klw:kuw)
      
      integer :: ndel

      ndel = (nguard_work - nguard)*npgs

!-----------------------------------------------------------
!
! Re-initialize the refinement and derefinement flag arrays
      refine(:)   = .false.
      derefine(:) = .false.



!
! Set up the workspace array WORK to store the variable we wish to examine 
! in order to test the refinement level.


! Set up the workspace array to store the current solution.
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
      if(nodetype(lb).eq.1.or.nodetype(lb).eq.2) then
            do k=kl_bnd+ndel*k3d,ku_bnd+ndel*k3d
            do j=jl_bnd+ndel*k2d,ju_bnd+ndel*k2d
            do i=il_bnd+ndel,iu_bnd+ndel
              work(i,j,k,lb,1) = &
     &                 unk(1,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)+unk(2,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
!				work(i,j,k,lb,1) = dsqrt( &
!						((facevarx(1,i+1,j,k,lb)+facevarx(1,i,j,k,lb))*0.5d0)**2.d0	&
!					   +((facevary(1,i,j+1,k,lb)+facevary(1,i,j,k,lb))*0.5d0)**2.d0 )
            end do
            end do
            end do
      endif
      end do
      endif



!
! Fill the guard cell layers of the workspace array.
      iopt=2
      nlayers=1
      call amr_guardcell(mype,iopt,nlayers)


!
! Error limits which control the refinement and derefinement requests below.
!      ctore = 0.06
!      ctode = 0.015
!      ctore = 0.2
!      ctode = 0.01
!      ctore = 0.06
!      ctode = 0.01

	! vorticity threshold for Karman vortex street

      ctore = 0.5
      ctode = 0.02

	! bed angle

!      ctore = 0.05
!      ctode = 0.02


!
! Loop over all leaf blocks and all parents of leaf blocks
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
      if(nodetype(lb).eq.1.or.nodetype(lb).eq.2) then

	      dx = bsize(1,lb)/real(nxb)
	      dy = bsize(2,lb)/real(nyb)

!
! User provided routine which returns an array error, which has some error
! measure computed for each grid cell of the current grid block, based on 
! some computation on the input array WORK(:,:,:,lb,1).
!       call error_measure( error, work(1,1,1,lb,1) )

		k = 1
       error(:,:,:) = 0.
!       do k=klw,kuw
!       do j=jlw+1,juw-1
!       do i=ilw+1,iuw-1
		do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
		do i=il_bnd+nguard,iu_bnd-nguard
         error1 = abs(work(i+1,j,k,lb,1)-work(i,j,k,lb,1))/dx
         error2 = abs(work(i-1,j,k,lb,1)-work(i,j,k,lb,1))/dx
         error3 = abs(work(i,j+1,k,lb,1)-work(i,j,k,lb,1))/dy
         error4 = abs(work(i,j-1,k,lb,1)-work(i,j,k,lb,1))/dy
         error_num = max( error1,error2,error3,error4 )
!!         error_den = max( work(i,j,k,lb,1)  ,work(i+1,j,k,lb,1), &
!!     &                    work(i-1,j,k,lb,1),work(i,j+1,k,lb,1), &
!!     &                    work(i,j-1,k,lb,1), 1.0e-6 )

!         error(i,j,k) = error_num	!/error_den

			error1 = (-(facevarx(1,i,j-1,k,lb)+facevarx(1,i+1,j-1,k,lb))+(facevarx(1,i,j+1,k,lb)+facevarx(1,i+1,j+1,k,lb)))*0.25/dy
			error2 = (-(facevary(1,i-1,j,k,lb)+facevary(1,i-1,j+1,k,lb))+(facevary(1,i+1,j,k,lb)+facevary(1,i+1,j+1,k,lb)))*0.25/dx
			
			error(i,j,k) = abs(error2-error1)

       enddo
       enddo
!       enddo  
       error_max = maxval( error )

! Does the error measure on this block anywhere exceed the limit which 
! should trigger refinement?

      if( lrefine(lb).lt.llrefine_max ) then
        if ( error_max .ge. ctore) refine(lb) = .true.
      endif


! Can we derefine this block?

      if( lrefine(lb).gt.llrefine_min .and. (.not.refine(lb)) ) then
        if ( error_max .lt. ctode) derefine(lb) = .true.
      endif



      endif
      end do                                   ! end of loop over blocks
      endif


!-----------------------------------------------------------

      return
      end subroutine amr_test_refinement

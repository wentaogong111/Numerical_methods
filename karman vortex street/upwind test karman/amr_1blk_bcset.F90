!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

	subroutine amr_1blk_bcset(mype,ibc,lb,pe, &
		&    idest,iopt,iface,jface,kface,surrblks)


!------------------------------------------------------------------------
!
! This is a template which is designed to show how to set up the
! boundary conditions which you desire. PARAMESH calls this routine
! when setting guardcell data on external boundaries. You must 
! modify this routine to set up your boundary conditions.
!
!
! This routine sets guardcell values at external boundaries in the case
! where a single block is having its guardcells filled.
!
! It can be assumed in writing this routine, that all guardcells for this
! block which are not across an external boundary have already been
! properly filled.
!
!
!------------------------------------------------------------------------
!
!
! Arguments:
!      mype             local processor
!      ibc              the integer specifying the particular boundary
!                        condition to be imposed
!      lb               block number of selected block
!      pe               processor on which block lb is located
!      idest            selects the storage space in data_1blk.fh which is to
!                        be used in this call. If the leaf node is having its
!                        guardcells filled then set this to 1, if its parent
!                        is being filled set it to 2.
!      iface            a selector setting designating whether the guarcells
!                        to be set are at the left, center or right sections
!                        of the i index range, eg
!                           iface = -1      left end
!                                 =  0      middle
!                                 = +1      right. For example, if iface=-1,
!                        the i index applied when filling unk will run
!                        from 1:nguard, if iface=0 from 1+nguard:nxb+nguard,
!                        and if iface=+1 from nxb+nguard+1:nxb+2*nguard.
!      jface            a selector setting designating whether the guarcells
!                        to be set are at the left, center or right sections
!                        of the j index range.
!      kface            a selector setting designating whether the guarcells
!                        to be set are at the left, center or right sections
!                        of the k index range.
!
!
! Written :     Peter MacNeice          August 1998
! Modified:     Peter MacNeice          January 2001
!------------------------------------------------------------------------

		use paramesh_dimensions
		use physicaldata
		use tree
		use workspace

		Include 'mpif.h'

		integer, intent(in) :: mype,ibc,lb,pe
		integer, intent(in) :: idest,iopt,iface,jface,kface
		integer, intent(in) :: surrblks(:,:,:,:)




		integer :: i,j,k
		integer :: is, js

!---------------------------------------------------------------------------
! Section to be modified by user

! Which boundary condition has been specified? ibc is the value
! of NEIGH for the current block face. If ibc is less than or equal
! to -20 you are on an external boundary. This if test should 
! conditionally execute the appropriate code for each of the different 
! boundary conditions you wish to impose.
! The example here serves to indicate the index ranges which need to be set.

!-------------------------
		if( ibc==-21 ) then
!-------------------------
! Boundary condition 1
!
!
!--------------------
			if( iopt==1 ) then
!--------------------


! Do cell centered data
          if(nvar.gt.0) then

!  if iface = -1 then limits are     i1 = 1           , i2 = nguard
!  if iface =  0 then limits are     i1 = nguard+1    , i2 = nxb+nguard
!  if iface = +1 then limits are     i1 = nxb+nguard+1, i2 = nxb+2*nguard
!  if jface = -1 then limits are     j1 = 1           , j2 = nguard
!  if jface =  0 then limits are     j1 = nguard+1    , j2 = nyb+nguard
!  if jface = +1 then limits are     j1 = nyb+nguard+1, j2 = nyb+2*nguard
!  if kface = -1 then limits are     k1 = 1           , k2 = nguard
!  if kface =  0 then limits are     k1 = nguard+1    , k2 = nzb+nguard
!  if kface = +1 then limits are     k1 = nzb+nguard+1, k2 = nzb+2*nguard
!  if ndim < 2, j1 = j2 = 1
!  if ndim < 3, k1 = k2 = 1
				
				if( iface==-1 ) then
					is = nguard
					i1 = 1
					i2 = nguard
				else if( iface==1 ) then
					is = -1
					i1 = nxb+nguard+1
					i2 = nxb+2*nguard
				else
					is = 0
					i1 = nguard+1
					i2 = nxb+nguard
				end if
				
				if( jface==-1 ) then
					js = nguard
					j1 = 1
					j2 = nguard
				else if( jface==1 ) then
					js = -1
					j1 = nyb+nguard+1
					j2 = nyb+2*nguard
				else
					js = 0
					j1 = nguard+1
					j2 = nyb+nguard
				end if
				
				k1 = 1
				k2 = 1
				
		!		write(*,'(a3,5i5)') 'lb= ', lb, neigh(1,1,lb), neigh(1,2,lb), neigh(1,3,lb), neigh(1,4,lb)
		!		write(*,'(a6,i3,a6,i3)') 'iface= ', iface, 'jface= ', jface
		!		write(*,'(4i5)') i1,i2,j1,j2
				
				do k = k1,k2
					do j = j1,j2
						do i = i1,i2
			!				unk1(:,i,j,k,idest) =  5.               !<<<<< USER EDIT
			!				unk1(1,i,j,k,idest) =  5.               !<<<<< USER EDIT
			!				unk1(2,i,j,k,idest) =  0.               !<<<<< USER EDIT
							unk1(1,i,j,k,idest) =  unk(1,i+is,j+js,k,lb)
							unk1(2,i,j,k,idest) =  unk(2,i+is,j+js,k,lb)
						enddo
					enddo
				enddo

          endif                          ! end of nvar if test



! Do cell corner data
          if(nvarcorn.gt.0) then

!  if iface = -1 then limits are     i1 = 1           , i2 = nguard
!  if iface =  0 then limits are     i1 = nguard+1    , i2 = nxb+nguard+1
!  if iface = +1 then limits are     i1 = nxb+nguard+2, i2 = nxb+2*nguard+1
!  if jface = -1 then limits are     j1 = 1           , j2 = nguard
!  if jface =  0 then limits are     j1 = nguard+1    , j2 = nyb+nguard+1
!  if jface = +1 then limits are     j1 = nyb+nguard+2, j2 = nyb+2*nguard+1
!  if kface = -1 then limits are     k1 = 1           , k2 = nguard
!  if kface =  0 then limits are     k1 = nguard+1    , k2 = nzb+nguard+1
!  if kface = +1 then limits are     k1 = nzb+nguard+2, k2 = nzb+2*nguard+1
!  if ndim < 2, j1 = j2 = 1
!  if ndim < 3, k1 = k2 = 1
           do k = k1,k2
           do j = j1,j2
           do i = i1,i2
!             unk_n1(:,i,j,k,idest) =  ????               !<<<<< USER EDIT
           enddo
           enddo
           enddo

          endif                          ! end of nvarcorn if test

! Do cell face centered  data
          if(nfacevar.gt.0) then

! facevarx1
!  if iface = -1 then limits are     i1 = 1           , i2 = nguard
!  if iface =  0 then limits are     i1 = nguard+1    , i2 = nxb+nguard+1
!  if iface = +1 then limits are     i1 = nxb+nguard+2, i2 = nxb+2*nguard+1
!  if jface = -1 then limits are     j1 = 1           , j2 = nguard
!  if jface =  0 then limits are     j1 = nguard+1    , j2 = nyb+nguard
!  if jface = +1 then limits are     j1 = nyb+nguard+1, j2 = nyb+2*nguard
!  if kface = -1 then limits are     k1 = 1           , k2 = nguard
!  if kface =  0 then limits are     k1 = nguard+1    , k2 = nzb+nguard
!  if kface = +1 then limits are     k1 = nzb+nguard+1, k2 = nzb+2*nguard+1
!  if ndim < 2, j1 = j2 = 1
!  if ndim < 3, k1 = k2 = 1

				if( iface==-1 ) then
					is = nguard+2
					i1 = 1
					i2 = nguard
				else if( iface==1 ) then
					is = -2
					i1 = nxb+nguard+2
					i2 = nxb+2*nguard+1
				else
					is = 0
					i1 = nguard+1
					i2 = nxb+nguard+1
				end if
				
				if( jface==-1 ) then
					js = nguard
					j1 = 1
					j2 = nguard
				else if( jface==1 ) then
					js = -1
					j1 = nyb+nguard+1
					j2 = nyb+2*nguard
				else
					js = 0
					j1 = nguard+1
					j2 = nyb+nguard
				end if
				
				k1 = 1
				k2 = 1

           do k = k1,k2
           do j = j1,j2
           do i = i1,i2
!             facevarx(1,i,j,k,lb) =  0.               !<<<<< USER EDIT
             facevarx1(:,i,j,k,idest) =  facevarx(1,i+is,j+js,k,lb)               !<<<<< USER EDIT
           enddo
           enddo
           enddo


! facevary1
!  if iface = -1 then limits are     i1 = 1           , i2 = nguard
!  if iface =  0 then limits are     i1 = nguard+1    , i2 = nxb+nguard
!  if iface = +1 then limits are     i1 = nxb+nguard+1, i2 = nxb+2*nguard
!  if jface = -1 then limits are     j1 = 1           , j2 = nguard
!  if jface =  0 then limits are     j1 = nguard+1    , j2 = nyb+nguard+1
!  if jface = +1 then limits are     j1 = nyb+nguard+2, j2 = nyb+2*nguard+1
!  if kface = -1 then limits are     k1 = 1           , k2 = nguard
!  if kface =  0 then limits are     k1 = nguard+1    , k2 = nzb+nguard
!  if kface = +1 then limits are     k1 = nzb+nguard+1, k2 = nzb+2*nguard
!  if ndim < 3, k1 = k2 = 1
           if(ndim.ge.2) then


				if( iface==-1 ) then
					is = nguard
					i1 = 1
					i2 = nguard
				else if( iface==1 ) then
					is = -1
					i1 = nxb+nguard+1
					i2 = nxb+2*nguard
				else
					is = 0
					i1 = nguard+1
					i2 = nxb+nguard
				end if
				
				if( jface==-1 ) then
					js = nguard+2
					j1 = 1
					j2 = nguard
				else if( jface==1 ) then
					js = -2
					j1 = nyb+nguard+2
					j2 = nyb+2*nguard+1
				else
					js = 0
					j1 = nguard+1
					j2 = nyb+nguard+1
				end if

           do k = k1,k2
           do j = j1,j2
           do i = i1,i2
!             facevary(1,i,j,k,lb) =  0.               !<<<<< USER EDIT
             facevary1(:,i,j,k,idest) =  -facevary(1,i+is,j+js,k,lb)               !<<<<< USER EDIT
           enddo
           enddo
           enddo
           endif

! facevarz1
!  if iface = -1 then limits are     i1 = 1           , i2 = nguard
!  if iface =  0 then limits are     i1 = nguard+1    , i2 = nxb+nguard
!  if iface = +1 then limits are     i1 = nxb+nguard+1, i2 = nxb+2*nguard
!  if jface = -1 then limits are     j1 = 1           , j2 = nguard
!  if jface =  0 then limits are     j1 = nguard+1    , j2 = nyb+nguard
!  if jface = +1 then limits are     j1 = nyb+nguard+1, j2 = nyb+2*nguard
!  if kface = -1 then limits are     k1 = 1           , k2 = nguard
!  if kface =  0 then limits are     k1 = nguard+1    , k2 = nzb+nguard+1
!  if kface = +1 then limits are     k1 = nzb+nguard+2, k2 = nzb+2*nguard+1
           do k = k1,k2
           do j = j1,j2
           do i = i1,i2
!             facevarz1(:,i,j,k,idest) =  ????               !<<<<< USER EDIT
             facevarz1(1,i,j,k,idest) =  0.d0               !<<<<< USER EDIT
           enddo
           enddo
           enddo


          endif                          ! end of nfacevar if test


! Do cell edge centered data
          if(nvaredge.gt.0) then

! unk_e_x1
!  if iface = -1 then limits are     i1 = 1           , i2 = nguard
!  if iface =  0 then limits are     i1 = nguard+1    , i2 = nxb+nguard
!  if iface = +1 then limits are     i1 = nxb+nguard+1, i2 = nxb+2*nguard
!  if jface = -1 then limits are     j1 = 1           , j2 = nguard
!  if jface =  0 then limits are     j1 = nguard+1    , j2 = nyb+nguard+1
!  if jface = +1 then limits are     j1 = nyb+nguard+2, j2 = nyb+2*nguard+1
!  if kface = -1 then limits are     k1 = 1           , k2 = nguard
!  if kface =  0 then limits are     k1 = nguard+1    , k2 = nzb+nguard+1
!  if kface = +1 then limits are     k1 = nzb+nguard+2, k2 = nzb+2*nguard+1
!  if ndim < 2, j1 = j2 = 1
!  if ndim < 3, k1 = k2 = 1
           do k = k1,k2
           do j = j1,j2
           do i = i1,i2
!             unk_e_x1(:,i,j,k,idest) =  ????               !<<<<< USER EDIT
           enddo
           enddo
           enddo

! finally unk_e_y1
!
!  if iface = -1 then limits are     i1 = 1           , i2 = nguard
!  if iface =  0 then limits are     i1 = nguard+1    , i2 = nxb+nguard+1
!  if iface = +1 then limits are     i1 = nxb+nguard+2, i2 = nxb+2*nguard+1
!  if jface = -1 then limits are     j1 = 1           , j2 = nguard
!  if jface =  0 then limits are     j1 = nguard+1    , j2 = nyb+nguard
!  if jface = +1 then limits are     j1 = nyb+nguard+1, j2 = nyb+2*nguard
!  if kface = -1 then limits are     k1 = 1           , k2 = nguard
!  if kface =  0 then limits are     k1 = nguard+1    , k2 = nzb+nguard+1
!  if kface = +1 then limits are     k1 = nzb+nguard+2, k2 = nzb+2*nguard+1
!  if ndim < 2, j1 = j2 = 1
!  if ndim < 3, k1 = k2 = 1
           do k = k1,k2
           do j = j1,j2
           do i = i1,i2
!             unk_e_y1(:,i,j,k,idest) =  ????               !<<<<< USER EDIT
           enddo
           enddo
           enddo

! finally unk_e_z1
!
!  if iface = -1 then limits are     i1 = 1           , i2 = nguard
!  if iface =  0 then limits are     i1 = nguard+1    , i2 = nxb+nguard+1
!  if iface = +1 then limits are     i1 = nxb+nguard+2, i2 = nxb+2*nguard+1
!  if jface = -1 then limits are     j1 = 1           , j2 = nguard
!  if jface =  0 then limits are     j1 = nguard+1    , j2 = nyb+nguard+1
!  if jface = +1 then limits are     j1 = nyb+nguard+2, j2 = nyb+2*nguard+1
!  if kface = -1 then limits are     k1 = 1           , k2 = nguard
!  if kface =  0 then limits are     k1 = nguard+1    , k2 = nzb+nguard
!  if kface = +1 then limits are     k1 = nzb+nguard+1, k2 = nzb+2*nguard
!  if ndim < 2, j1 = j2 = 1
!  if ndim < 3, k1 = k2 = 1
           do k = k1,k2
           do j = j1,j2
           do i = i1,i2
!             unk_e_z1(:,i,j,k,idest) =  ????               !<<<<< USER EDIT
           enddo
           enddo
           enddo

          endif                          ! end of nvaredge if test



!--------------------
        elseif(iopt.ge.2) then
!--------------------


! Operate on work array

! if iface = -1 then limits are   i1 = 1                , i2 = nguard_work
! if iface =  0 then limits are   i1 = nguard_work+1    , i2 = nxb+nguard_work
! if iface = +1 then limits are   i1 = nxb+nguard_work+1, i2 = nxb+2*nguard_work
! if jface = -1 then limits are   j1 = 1                , j2 = nguard_work
! if jface =  0 then limits are   j1 = nguard_work+1    , j2 = nyb+nguard_work+1
! if jface = +1 then limits are   j1 = nyb+nguard_work+1, j2 = nyb+2*nguard_work
! if kface = -1 then limits are   k1 = 1                , k2 = nguard_work
! if kface =  0 then limits are   k1 = nguard_work+1    , k2 = nzb+nguard_work
! if kface = +1 then limits are   k1 = nzb+nguard_work+1, k2 = nzb+2*nguard_work
! if ndim < 2, j1 = j2 = 1
! if ndim < 3, k1 = k2 = 1
           
				if( iface==-1 ) then
					is = nguard_work
					i1 = 1
					i2 = nguard_work
				else if( iface==1 ) then
					is = -1
					i1 = nxb+nguard_work+1
					i2 = nxb+2*nguard_work
				else
					is = 0
					i1 = nguard_work+1
					i2 = nxb+nguard_work
				end if
				
				if( jface==-1 ) then
					js = nguard_work
					j1 = 1
					j2 = nguard_work
				else if( jface==1 ) then
					js = -1
					j1 = nyb+nguard_work+1
					j2 = nyb+2*nguard_work
				else
					js = 0
					j1 = nguard_work+1
					j2 = nyb+nguard_work
				end if
				
				k1 = 1
				k2 = 1
				
				do k = k1,k2
					do j = j1,j2
						do i = i1,i2
			!				unk1(:,i,j,k,idest) =  5.               !<<<<< USER EDIT
			!				unk1(1,i,j,k,idest) =  5.               !<<<<< USER EDIT
			!				unk1(2,i,j,k,idest) =  0.               !<<<<< USER EDIT
							work1(i,j,k,idest) = work(i+is,j+js,k,lb,iopt-1)
						enddo
					enddo
				enddo

!--------------------
        endif
!--------------------


!-------------------------
!       elseif(ibc.eq. ???? ) then
!-------------------------
! Boundary condition 2
!

!-------------------------
       endif                            ! end of test of bc flag
!-------------------------

!	call amr_1blk_to_perm(lcc,lfc,lec,lnc,lb,iopt,idest)


! End of Section to be modified by user
!---------------------------------------------------------------------------

      return
      end subroutine amr_1blk_bcset

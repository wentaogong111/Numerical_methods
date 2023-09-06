
	subroutine paramesh2paraview( itime ) 
		use paramesh_dimensions
		use physicaldata
		use tree
      
		integer, intent(in) :: itime

		integer :: l, ll, lvis, mm
		integer :: nxy
		integer :: im1, im2, jm1, jm2
		integer :: i, j, ij
		real	:: dx, dy, xi, yi
		character(300) :: fpout, fpout2
		character (len=6)  :: fnum_string
		real :: fave
		
		
		write (fnum_string, '(i6.6)') itime

		im1 = il_bnd+nguard
		im2 = iu_bnd-nguard
		jm1 = jl_bnd+nguard*k2d
		jm2 = ju_bnd-nguard*k2d
		
		nxy = nxb*nyb
		
		lvis = 0
		do l=1,lnblocks
			if( nodetype(l)==1 ) then
				lvis = lvis+1
			end if
		end do
		
		fpout  = 'paramesh_vtk_' // fnum_string // '.vtk'
		fpout2 = 'paramesh_vtk_detail_' // fnum_string // '.vtk'

		open( 55, file=fpout, status='unknown' )
		
			write(55,'(a26)') '# vtk DataFile Version 3.0'
			write(55,'(a5)') 'Test '
			write(55,'(a5)') 'ASCII'
			write(55,'(a25)') 'DATASET UNSTRUCTURED_GRID'
			
!			write(55,'(a6,i7,a10)') 'POINTS',lnblocks*4,'FLOAT'
			write(55,'(a6,i7,a10)') 'POINTS',lvis*4,'FLOAT'
			
			do l=1,lnblocks
				if( nodetype(l)==1 ) then
					write(55,*) (coord(1,l)-bsize(1,l)*0.5d0), (coord(2,l)-bsize(2,l)*0.5d0), 0.d0
					write(55,*) (coord(1,l)+bsize(1,l)*0.5d0), (coord(2,l)-bsize(2,l)*0.5d0), 0.d0
					write(55,*) (coord(1,l)+bsize(1,l)*0.5d0), (coord(2,l)+bsize(2,l)*0.5d0), 0.d0
					write(55,*) (coord(1,l)-bsize(1,l)*0.5d0), (coord(2,l)+bsize(2,l)*0.5d0), 0.d0
				end if
			end do
			
!			write(55,'(a5,2i8)') 'CELLS', lnblocks, lnblocks*5
			write(55,'(a5,2i8)') 'CELLS', lvis, lvis*5
			
!			do l=0,lnblocks-1
!				ll = l*4
!				write(55,'(i3,4i8)') 4, ll, ll+1, ll+2, ll+3
!			end do

			mm = 0
			do l=1,lnblocks
				if( nodetype(l)==1 ) then
					ll = mm*4
					write(55,'(i3,4i8)') 4, ll, ll+1, ll+2, ll+3
					mm = mm+1
				end if
			end do
			
			write(55,'(a10,i8)') 'CELL_TYPES', lvis	!lnblocks
			
!			do l=1,lnblocks
			do l=1,lvis
				write(55,'(i3)') 7
			end do
			
			write(55,'(a10,i8)') 'POINT_DATA', lvis*4	!lnblocks*4
			write(55,'(a23,i5)') 'SCALARS ELEVATION FLOAT',1
			write(55,'(a20)') 'LOOKUP_TABLE DEFAULT'
			
			do l=1,lnblocks
				if( nodetype(l)==1 ) then
	!				write(55,*) 0.d0
	!				write(55,*) 0.d0
	!				write(55,*) 0.d0
	!				write(55,*) 0.d0
					write(55,*) unk(2,im1,jm1,1,l)
					write(55,*) unk(2,im2,jm1,1,l)
					write(55,*) unk(2,im2,jm2,1,l)
					write(55,*) unk(2,im1,jm2,1,l)
				end if
			end do
			
			write(55,'(a23)') 'VECTORS velocity double'
			
			do l=1,lnblocks
				if( nodetype(l)==1 ) then
					write(55,*) (facevarx(1,im1  ,jm1,1,l)+facevarx(1,im1  ,jm1-1,1,l))*0.5,	&
								(facevary(1,im1,jm1  ,1,l)+facevary(1,im1-1,jm1  ,1,l))*0.5, 0.
					write(55,*) (facevarx(1,im2+1,jm1,1,l)+facevarx(1,im2+1,jm1-1,1,l))*0.5,	&
								(facevary(1,im2,jm1  ,1,l)+facevary(1,im2+1,jm1  ,1,l))*0.5, 0.
					write(55,*) (facevarx(1,im2+1,jm2,1,l)+facevarx(1,im2+1,jm2+1,1,l))*0.5,	&
								(facevary(1,im2,jm2+1,1,l)+facevary(1,im2+1,jm2+1,1,l))*0.5, 0.
					write(55,*) (facevarx(1,im1  ,jm2,1,l)+facevarx(1,im1  ,jm2+1,1,l))*0.5,	&
								(facevary(1,im1,jm2+1,1,l)+facevary(1,im1-1,jm2+1,1,l))*0.5, 0.
				end if
			end do
			
			write(55,'(a9,i8)') 'CELL_DATA', lvis	!lnblocks
			write(55,'(a7,2x,a12,2x,a7,i4)') 'SCALARS', 'DD', 'double', 1
			write(55,'(a20)') 'LOOKUP_TABLE DEFAULT'
			
			do l=1,lnblocks
				if( nodetype(l)==1 ) then
					
					fave = 0.d0
					
					do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
						do i=il_bnd+nguard,iu_bnd-nguard
!					do j=jl_bnd,ju_bnd
!						do i=il_bnd,iu_bnd
							fave = fave + unk(1,i,j,1,l)
						end do
					end do
					
					fave = fave/dble(nxb*nyb)
						
					write(55,'(f15.9)') fave
				end if
			end do

			write(55,'(a7,2x,a12,2x,a7,i4)') 'SCALARS', 'HH', 'double', 1
			write(55,'(a20)') 'LOOKUP_TABLE DEFAULT'
			
			do l=1,lnblocks
				if( nodetype(l)==1 ) then
					
					fave = 0.d0
					
					do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
						do i=il_bnd+nguard,iu_bnd-nguard
!					do j=jl_bnd,ju_bnd
!						do i=il_bnd,iu_bnd
							fave = fave + unk(1,i,j,1,l)+unk(2,i,j,1,l)
						end do
					end do
					
					fave = fave/dble(nxb*nyb)
						
					write(55,'(f15.9)') fave
				end if
			end do

!			write(55,'(a7,2x,a12,2x,a7,i4)') 'SCALARS', 'BC', 'double', 1
!			write(55,'(a20)') 'LOOKUP_TABLE DEFAULT'

!			do l=1,lnblocks
!				if( nodetype(l)==1 ) then
				
!!					fave = 0.d0
					
!!					do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
!!						do i=il_bnd+nguard,iu_bnd-nguard
!!							fave = fave + sqrt(			&
!!												(facevarx(1,i,j,1,l)+facevarx(1,i+1,j,1,l))**2.d0		&
!!											   +(facevary(1,i,j,1,l)+facevary(1,i,j+1,1,l))**2.d0			)
!!						end do
!!					end do
					
!!					fave = fave/dble(nxb*nyb)

!					if( neigh(1,1,l)<-20 .or. neigh(1,2,l)<-20	&
!						.or. neigh(1,3,l)<-20 .or. neigh(1,4,l)<-20 ) then
!!					if( ( neigh(1,1,l)>-20.and.neigh(1,1,l)<0 ) .or. 	&
!!						( neigh(1,2,l)>-20.and.neigh(1,2,l)<0 ) .or.	&
!!						( neigh(1,3,l)>-20.and.neigh(1,3,l)<0 ) .or.	&
!!						( neigh(1,4,l)>-20.and.neigh(1,4,l)<0 ) ) then
!!						fave = 1.d0
!						fave = min(neigh(1,1,l),neigh(1,2,l))
!						fave = min(fave,float(neigh(1,3,l)))
!						fave = min(fave,float(neigh(1,4,l)))
!					else
!						fave = 0.d0
!					end if
!						
!					write(55,'(f15.9)') fave
!				end if
!			end do

			write(55,'(a7,2x,a12,2x,a7,i4)') 'SCALARS', 'DZ', 'double', 1
			write(55,'(a20)') 'LOOKUP_TABLE DEFAULT'

			do l=1,lnblocks
				if( nodetype(l)==1 ) then
					
					fave = 0.d0
					
					do j=jl_bnd+nguard*k2d,ju_bnd-nguard*k2d 
						do i=il_bnd+nguard,iu_bnd-nguard
							fave = fave + unk(3,i,j,1,l)
						end do
					end do
					
					fave = fave/dble(nxb*nyb)
					write(55,'(f15.9)') fave
				end if
			end do
			
		close(55)

		open( 56, file=fpout2, status='unknown' )
		
			write(56,'(a26)') '# vtk DataFile Version 3.0'
			write(56,'(a5)') 'Test '
			write(56,'(a5)') 'ASCII'
			write(56,'(a25)') 'DATASET UNSTRUCTURED_GRID'
			
			write(56,'(a6,i7,a10)') 'POINTS',lvis*nxy*4,'FLOAT'
			
			do l=1,lnblocks
				if( nodetype(l)==1 ) then
					
					dx = bsize(1,l)/real(nxb)
					dy = bsize(2,l)/real(nyb)
					
					do j=1,nyb
						do i=1,nxb
							xi = dx*float(i-1)+(coord(1,l)-bsize(1,l)*0.5d0)
							yi = dy*float(j-1)+(coord(2,l)-bsize(2,l)*0.5d0)
							
							write(56,*) xi   , yi   , 0.d0
							write(56,*) xi+dx, yi   , 0.d0
							write(56,*) xi+dx, yi+dy, 0.d0
							write(56,*) xi   , yi+dy, 0.d0
						end do
					end do
				end if
			end do
			
			write(56,'(a5,2i8)') 'CELLS', lvis*nxy, lvis*nxy*5
			
			mm = 0
			do l=1,lnblocks
				if( nodetype(l)==1 ) then
					ll = mm*nxy*4
					do j=1,nyb
						do i=1,nxb
							ij = (nxb*(j-1)+i-1)*4
							write(56,'(i3,4i8)') 4, ll+ij, ll+ij+1, ll+ij+2, ll+ij+3
						end do
					end do
					mm = mm+1
				end if
			end do
			
			write(56,'(a10,i8)') 'CELL_TYPES', lvis*nxy
			
			do l=1,lvis*nxy
				write(56,'(i3)') 7
			end do
			
			write(56,'(a10,i8)') 'POINT_DATA', lvis*nxy*4
			write(56,'(a23,i5)') 'SCALARS ELEVATION FLOAT',1
			write(56,'(a20)') 'LOOKUP_TABLE DEFAULT'
			
			do l=1,lnblocks
				if( nodetype(l)==1 ) then
					do j=jm1,jm2 
						do i=im1,im2
							write(56,*) unk(2,i,j,1,l)
							write(56,*) unk(2,i,j,1,l)
							write(56,*) unk(2,i,j,1,l)
							write(56,*) unk(2,i,j,1,l)
						end do
					end do
				end if
			end do
			
			write(56,'(a23)') 'VECTORS velocity double'
			
			do l=1,lnblocks
				if( nodetype(l)==1 ) then
					do j=jm1,jm2
						do i=im1,im2
							write(56,*) (facevarx(1,i  ,j,1,l)+facevarx(1,i  ,j-1,1,l))*0.5,	&
										(facevary(1,i,j  ,1,l)+facevary(1,i-1,j  ,1,l))*0.5, 0.
							write(56,*) (facevarx(1,i+1,j,1,l)+facevarx(1,i+1,j-1,1,l))*0.5,	&
										(facevary(1,i,j  ,1,l)+facevary(1,i+1,j  ,1,l))*0.5, 0.
							write(56,*) (facevarx(1,i+1,j,1,l)+facevarx(1,i+1,j+1,1,l))*0.5,	&
										(facevary(1,i,j+1,1,l)+facevary(1,i+1,j+1,1,l))*0.5, 0.
							write(56,*) (facevarx(1,i  ,j,1,l)+facevarx(1,i  ,j+1,1,l))*0.5,	&
										(facevary(1,i,j+1,1,l)+facevary(1,i-1,j+1,1,l))*0.5, 0.
						end do
					end do
				end if
			end do
			
			write(56,'(a9,i8)') 'CELL_DATA', lvis*nxy
			write(56,'(a7,2x,a12,2x,a7,i4)') 'SCALARS', 'DD', 'double', 1
			write(56,'(a20)') 'LOOKUP_TABLE DEFAULT'
			
			do l=1,lnblocks
				if( nodetype(l)==1 ) then
					
					do j=jm1,jm2 
						do i=im1,im2
							write(56,'(f15.9)') unk(1,i,j,1,l)
						end do
					end do

				end if
			end do

			write(56,'(a7,2x,a12,2x,a7,i4)') 'SCALARS', 'DZ', 'double', 1
			write(56,'(a20)') 'LOOKUP_TABLE DEFAULT'
			
			do l=1,lnblocks
				if( nodetype(l)==1 ) then
					
					do j=jm1,jm2 
						do i=im1,im2
							write(56,'(f15.9)') unk(2,i,j,1,l)
						end do
					end do

				end if
			end do
			
		close(56)

	end subroutine paramesh2paraview

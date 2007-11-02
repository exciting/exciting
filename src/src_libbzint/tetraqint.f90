!BOP
!
! !ROUTINE: tetraqint 
!
! !INTERFACE:
      subroutine tetraqint(nik,nt,nb,ebd,tetc,wtet,wqk,v,efer,oc,opm,&
     &                     intopm)
!      
! !DESCRIPTION:
!
!   This subroutine calculates integrals of the form
!
!\begin{equation}
!\mathcal{O}(\vec{k})=\int\limits_{BZ}{\langle \Phi_n(\vec{k})|X(\vec{q})|%
!\Phi_{n'}(\vec{k}-\vec{q})\rangle\Theta[\pm(\epsilon_F-\epsilon_{n'
!,\vec{k}-\vec{q}})]d^3q
!\end{equation}
!
!using the improved tetrahedron method

! !USES:
 
      use tetra_internal
       
      implicit none      
       
! !INPUT PARAMETERS:
 
      integer(4), intent(in) :: nik        ! Number of irreducible 
!                                             k-points
       
      integer(4), intent(in) :: nt         ! Number of tetrahedra
       
      integer(4), intent(in) :: nb         ! Number of bands
       
      real(8), target, intent(in) :: ebd(nb,nik)  ! Band energies
       
      integer(4), target, intent(in) :: tetc(4,*)  ! id. numbers of 
!                                                    the corners
!                                                    of the tetrahedra
   
      integer(4), target, intent(in) :: wtet(*)    ! weight of each 
!                                                    tetrahedron
       
      integer(4), target, intent(in) :: wqk(nik,nik)    ! weight of 
!                                                         each q
!                                                         point
       
      real(8), intent(in)    :: v         ! the volume of the tetrahedra
 
      real(8), intent(in)    :: efer       ! fermi energy
       
      logical, intent(in)    :: oc         ! .true. if the integration
!                                             runs over occupied states, 
!                                             .false. over excited ones
       
      real(8), intent(in)    :: opm(nik,nik,*) ! the values of the operator
!                                                 that has to be integrated
       
! !OUTPUT PARAMETERS:
       
      real(8), intent(out)   :: intopm    ! the value of the integral

!  
! !LOCAL VARIABLES:
 
      integer(4) :: ik,ib,jk,normfac
       
      real(8) :: fermi
       
      real(8), allocatable :: kw(:,:)    ! Weight of each k-point

      real(8), allocatable :: mop(:,:)   ! Weighted mean value of the
!                                           operator opm, summed over q
       
! !SYSTEM ROUTINES:
       
      intrinsic size
      
      external intw
       
! !REVISION HISTORY:
!
!   Created: 5th. March 2004 by RGA
!
!EOP
!BOC
 
      nirkp = nik
      ntet  = nt
      nband = nb
      vt = v
      tetcorn => tetc(1:4,1:ntet)
      qweig => wqk(1:nik,1:nik)
      tetweig => wtet(1:ntet)

      allocate(kw(nirkp,nband))
      allocate(mop(nirkp,nband))
      
      do ik=1,nirkp
        do ib=1,nband
          normfac=0
          mop(ik,nband)=0.0d0
          do jk=1,nirkp
            normfac=normfac+qweig(ik,jk)
            mop(ik,nband)=mop(ik,nband)+dble(qweig(ik,jk))*opm(ik,jk,ib)
          enddo
          mop(ik,nband)=mop(ik,nband)/dble(normfac)
        enddo
      enddo
          
      if(oc)then
        eband   => ebd(1:nband,1:nik)
        fermi = efer
      else
        allocate(eband(nband,nik))
        do ik=1,nik
          do ib=1,nband
            eband(ib,ik)=-1.0d0*ebd(ib,ik)
          enddo
        enddo
        fermi=-1.0d0*efer
      endif
      
      call intw(fermi,kw)
      
      intopm = 0
      
      do ik=1,nik
        do ib=1,nband
          intopm=intopm+kw(ik,ib)*mop(ik,ib)
        enddo
      enddo
      
      end subroutine tetraqint
      
!EOC

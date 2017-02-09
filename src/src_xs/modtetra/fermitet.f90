#ifdef TETRA
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: fermitet
!
! !INTERFACE:
      subroutine fermitet(nik,nbd,eband,ntet,tetc,wtet,vt,nel,sp,efer,efden,lprt)
!
! !DESCRIPTION:
!  This subroutine calculated the Fermi energy with the tetrahedron method
!
! !USES:
      
!<sag>
      use control, only: tetradbglv
!</sag>
      implicit none     
      
! !INPUT PARAMETERS:

      integer(4), intent(in) :: nik ! Number of irreducible k-points
      
      integer(4), intent(in) :: nbd ! Maximum number of bands
      
      real(8), intent(in) :: eband(nbd,nik) ! Band energies
      
      integer(4), intent(in) :: ntet        ! Number of tetrahedra
      
      integer(4), intent(in) :: tetc(4,*)   ! id. numbers of the corners
!                                             of the tetrahedra
  
      integer(4), intent(in) :: wtet(*)    ! weight of each tetrahedron
      
      real(8), intent(in)    :: vt         ! the volume of the tetrahedra

      real(8), intent(in)    :: nel        ! number of electrons
      
      logical, intent(in)    :: sp         ! .true. for spin polarized
!                                             case
      logical, intent(in) :: lprt
      
! !OUTPUT PARAMETERS:      
      
      real(8), intent(out)   :: efer       ! the fermi energy
      real(8), intent(out)   :: efden      ! the density of states at efer

! !REVISION HISTORY:
!
! Created 10th. March 2004 by RGA
!
! !LOCAL VARIABLES:

      integer(4) :: ik                 ! (Counter) Runs over k-points
      integer(4) :: ib                 ! (Counter) Runs over energies
      integer(4) :: it

      real(8) :: fact   ! Number of electrons per state.
      real(8) :: emin   ! Lower limit of the energy interval around ef.
      real(8) :: emax   ! Higher limit of the energy interval around ef.
      real(8) :: eint   ! New approximation to ef.
      real(8) :: ocmin  ! Number of electrons up to emin.
      real(8) :: ocmax  ! Number of electrons up to emax.
      real(8) :: ocint  ! Number of electrons up to eint.
      real(8) :: df     ! DOS at ef.
      
      real(8), external :: idos
      real(8), external :: dostet
      integer, parameter :: nitmax=200
      real(8), parameter :: epsb=1.d-4  ! Tolerance for bisection method
      real(8), parameter :: epss=1.d-8  ! Tolerance for secant method

!EOP
!BOC
      fact=1.0d0
      if(.not.sp)fact = fact + 1.0d0
!
!     Set emin to the lowest band energy, and emax to the highest.
!
      emin=100.00
      emax=-100.00
      do ik=1,nik
        do ib=1,nbd
          if(eband(ib,ik).lt.emin) emin=eband(ib,ik)
          if((eband(ib,ik).gt.emax).and.(eband(ib,ik).lt.5.0d+2))       &
     &                emax=eband(ib,ik)
        enddo
      enddo
!
!     Calculate the corresponding occupations.
!      
      ocmin=fact*idos(nik,nbd,eband,ntet,tetc,wtet,vt,emin)
      ocmax=fact*idos(nik,nbd,eband,ntet,tetc,wtet,vt,emax)
!     
!     If the maximum occupation is lower than the number of electrons
!     exit with the corresponding error message. 
!          
      if(ocmax.le.nel) stop 'not enough bands'
      
      ocint=0.0d0
      it=0

      if(lprt)then
        write(6,*)'fermi: bisection method'
        write(6,1) '#it',"emin","emax","eint",                  &
     &      "ocmin","ocmax","ocint"
      endif
!
! Use bisection method to determine solver the equation N( Efermi ) = Nel
!
      do it=1,nitmax
        eint=0.50d0*(emax+emin)
        ocint = fact*idos(nik,nbd,eband,ntet,tetc,wtet,vt,eint)

!        if(lprt) write(6,2) it, emin, emax, eint, ocmin, ocmax, ocint
          if(lprt) write(6,2) it, emin, emax, eint, ocmin-nel, ocmax-nel, ocint-nel
        
        if((abs(ocint-nel).le. epsb)) exit 
        if(ocint.ge.nel)then
          emax=eint
          ocmax=ocint
        else
          emin=eint
          ocmin=ocint
        endif
      enddo
      
      if(it.ge.nitmax) stop "fermi: fail to converge"

      df=dostet(nik,nbd,eband,ntet,tetc,wtet,vt,eint)
      
      if(df.gt.1.0d-9)then
        if(ocint.ge.nel)then
          emax=eint
          ocmax=ocint
        else
          emin=eint
          ocmin=ocint
        endif
        if(abs(emax-emin).gt.epsb)then
          eint=0.50d0*(emax+emin)
          ocint = fact*idos(nik,nbd,eband,ntet,tetc,wtet,vt,eint)
          if(ocint.lt.nel)then
            emin=eint
            ocmin=ocint
            eint=0.50d0*(emax+emin)
            ocint = fact*idos(nik,nbd,eband,ntet,tetc,wtet,vt,eint)
          endif  
        endif  
      
        df=dostet(nik,nbd,eband,ntet,tetc,wtet,vt,eint)
      endif

      if(lprt) write(6,3)'Efermi,DOS(Efermi) = ',eint, df
!
! For insulator (including semiconductor, set fermi energy as the middle of gap)
!
      if(df.lt.1.0d-9)then
        emin=maxval(eband(1:nbd,1:nik),mask=eband(1:nbd,1:nik).lt. &
     &              eint)
        emax=minval(eband(1:nbd,1:nik),mask=eband(1:nbd,1:nik).gt. &
     &              eint)
        eint=0.5d0*(emin+emax)
!<sag>
        if (tetradbglv > 0) then
!</sag>
           write(6,*)emin,emax,eint
!<sag>
        end if
!</sag>
      else
!
! For metals, Aply the secant method for highest precision in the Fermi energy             
!
        if(lprt)then
          write(6,*)
          write(6,*)'fermi: secant method'
          write(6,1) '#it',"emin","emax","eint",                  &
     &        "ocmin","ocmax","ocint"
        endif
        do it=1,nitmax
          eint=emin+(emax-emin)/(ocmax-ocmin)*(nel-ocmin)
          ocint = fact*idos(nik,nbd,eband,ntet,tetc,wtet,vt,eint)

          if(lprt) write(6,2) it, emin, emax, eint, ocmin-nel, ocmax-nel, ocint-nel
        
          if(ocint.ge.nel)then
            emax=eint
            ocmax=ocint
          else
            emin=eint
            ocmin=ocint
          endif
          if((abs(ocint-nel).le. epss).and.(ocint.ge.nel)) exit
        enddo
        df=dostet(nik,nbd,eband,ntet,tetc,wtet,vt,eint)
        if(it.ge.nitmax)then
          if(lprt) write(6,3)'Efermi,DOS(Efermi) = ',eint, df
          stop "fermi: secant method fail to converge"
        endif  
      endif        
      efer=eint
      efden=df
      if(lprt) write(6,3)'Efermi,DOS(Efermi) = ',efer, df

  1   format(A5,6(A12,8x))
  2   format(I5,6g20.12)
  3   format(A,2g16.8) 
      end subroutine fermitet
!EOC        
#endif

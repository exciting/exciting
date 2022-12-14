!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: occupy
! !INTERFACE:
!
!
Subroutine occupy
! !USES:
      Use modinput
      Use modmain
      use mod_opt_tetra
      use mod_kpointset

! !DESCRIPTION:
!   Finds the Fermi energy and sets the occupation numbers for the
!   second-variational states using the routine {\tt fermi}.
!
! !REVISION HISTORY:
!   Created February 2004 (JKD)
!   Modifiactions for tetrahedron method, November 2007 (RGA alias
!     Ricardo Gomez-Abal)
!   Modifications for tetrahedron method, 2007-2010 (Sagmeister)
!   Modifications for tetrahedron method, 2011 (DIN)
!   Simplicistic method for systems with gap added, 2013 (STK)
!EOP
!BOC
      implicit none
! local variables
      integer, parameter :: maxit = 1000
      real(8), parameter :: de0=1.d0
      integer :: ik, ist, it, nvm
      real(8) :: e0, e1, chg, x, t1
! external functions
      real(8) :: sdelta, stheta
      real(8) :: egap
      real(8) :: dfde(nstsv,nkpt)
      logical :: lspin

      character(1024) :: message

      external sdelta, stheta
      real(8), external :: dostet_exciting, dos_at_energy

      type( k_set) :: kset
      type( t_set) :: tetra
      
      if ( input%groundstate%stypenumber .ge. 0 ) then
         t1 = 1.d0 / input%groundstate%swidth

!     next lines taken in part from libbzint (STK)
!
!!    nvm is the number of bands for an insulating system
!!    since for a system with gap, the procedure to determine the
!!    band gap can be unstable, just try first whether it is an
!!    insulating system, but such a simplicistic way to determine the Fermi energy
!!    is valid only for no spin polarized cases
!
         if (.not.associated(input%groundstate%spin)) then
           nvm  = nint(chgval/occmax)
           e0 = maxval(evalsv(nvm,:))
           e1 = minval(evalsv(nvm+1,:))
           efermi = 0.5*(e0 + e1)

           fermidos = 0.d0
           chg = 0.d0
           Do ik = 1, nkpt
              Do ist = 1, nstsv
                 x = (evalsv(ist, ik)-efermi) * t1
                 fermidos = fermidos + wkpt (ik) * sdelta (input%groundstate%stypenumber, x) * t1
                 occsv (ist, ik) = occmax * stheta (input%groundstate%stypenumber, -x)
                 chg = chg + wkpt (ik) * occsv (ist, ik)
              End Do
           End Do
           fermidos = fermidos * occmax
           if ((e1 .ge. e0) .and. (abs(chg - chgval) .lt. input%groundstate%epsocc)) then
!            Write (*, '("Info(occupy): System has gap, simplicistic method used in determining efermi and occupation")')
             goto 10
           endif
         end if

! find minimum and maximum eigenvalues
         e0 = evalsv (1, 1)
         e1 = e0
         Do ik = 1, nkpt
            Do ist = 1, nstsv
               e0 = Min (e0, evalsv(ist, ik))
               e1 = Max (e1, evalsv(ist, ik))
            End Do
         End Do
!
! determine the Fermi energy using the bisection method
!
         Do it = 1, maxit
            efermi = 0.5d0 * (e0+e1)
            chg = 0.d0
            Do ik = 1, nkpt
               Do ist = 1, nstsv
                  x = (efermi-evalsv(ist, ik)) * t1
                  occsv (ist, ik) = occmax * stheta &
                   & (input%groundstate%stypenumber, x)
                  chg = chg + wkpt (ik) * occsv (ist, ik)
               End Do
            End Do
            If (chg .Lt. chgval) Then
               e0 = efermi
            Else
               e1 = efermi
            End If
            If ((e1-e0) .Lt. input%groundstate%epsocc) Go To 10
         End Do
         Write (*,*)
         Write (*, '("Error(occupy): could not find Fermi energy")')
         Write (*,*)
         Stop
10       Continue
! find the density of states at the Fermi surface in units of
! states/Hartree/spin/unit cell
         fermidos = dos_at_energy( nkpt, wkpt, nstsv, evalsv, occmax, efermi, &
                      input%groundstate%stypenumber, input%groundstate%swidth)
         Do ik = 1, nkpt
            If (occsv(nstsv, ik) .Gt. input%groundstate%epsocc) Then
               call warning('Warning(occupy):')
               Write (message, '(" Not enough empty states for k-point ", I6)') ik
               call warning(message)
            End If
         End Do

      else  if (input%groundstate%stypenumber==-1) then
         !------------------------------------------------------------
         ! Use the tetrahedron integration method (LIBBZINT library)
         !------------------------------------------------------------
         ! Calculate the Fermi energy
         lspin = associated(input%groundstate%spin)
         call fermi_exciting(lspin, &
                             chgval, &
                             nstsv, nkpt, evalsv, &
                             ntet, tnodes, wtet, tvol, &
                             efermi, egap, fermidos)
         ! write(*,*) 'occupy: ', efermi, egap, fermidos
         ! Calculate state occupation numbers
         call tetiw(nkpt, ntet, nstsv, evalsv, tnodes, wtet, tvol, efermi, occsv)
         do ik = 1, nkpt
           do ist = 1, nstsv
             occsv(ist,ik) = dble(occmax)/wkpt(ik)*occsv(ist,ik)
           end do
         end do

      else  if (input%groundstate%stypenumber==-2) then
         call generate_k_vectors( kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, input%groundstate%reducek, uselibzint=.false.)
         call opt_tetra_init( tetra, kset, 2, reduce=.true.)
         !--------------------------------------
         ! Use the improved tetrahedron method
         !--------------------------------------
	 nvm  = nint(chgval/occmax)
         efermi = 0.5d0*(maxval( evalsv( nvm, :)) + minval( evalsv( nvm+1, :)))
         call opt_tetra_efermi( tetra, chgval/dble(occmax), nkpt, nstsv, evalsv, efermi, occsv, ef0=efermi, df0=efermi)
         do ik = 1, nkpt
           occsv(:,ik) = dble(occmax)/wkpt(ik)*occsv(:,ik)
         end do
         !write(*,*) 'occsv=', occsv(:,1)

         call opt_tetra_wgt_delta( tetra, nkpt, nstsv, evalsv, 1, (/efermi/), dfde)
         fermidos = sum(dfde)
         call opt_tetra_destroy( tetra)
         call delete_k_vectors( kset)
         !write(*,*) 'dos at Ef=', fermidos

      End If ! modified tetrahedron integration method

      Return
End Subroutine
!EOC

!> The total charge \(\rho\) in the system is given by
!> \[ \rho = \sum_{n,{\bf k}} w_{\bf k} \, f_{n{\bf k}} \;, \]
!> where \(f_{n{\bf k}} = \theta(\epsilon_{\rm F} - \epsilon_{n{\bf k}})\) is the occupation of state \(\Psi_{n{\bf k}}\)
!> and \(\epsilon_{\rm F}\) is the Fermi energy.
!> The Fermi energy is iteratively refined using the bisection method until the resulting occupation numbers yield
!> the desired charge within a given tolerance.
!> \(\theta(x)\) is approximated with a smooth function according to `stype` and `swidth`.
subroutine find_fermi( nkpt, wkpt, nst, eval, chg, maxocc, stype, swidth, epschg, efermi, occ)
  use precision, only: dp
  use modmpi, only: terminate_if_false
  !> number of k-points
  integer, intent(in) :: nkpt
  !> k-point weights
  real(dp), intent(in) :: wkpt(nkpt)
  !> number of states
  integer, intent(in) :: nst
  !> eigenenergies
  real(dp), intent(in) :: eval(nst,nkpt)
  !> charge to fill
  real(dp), intent(in) :: chg
  !> maximum occupation per state
  real(dp), intent(in) :: maxocc
  !> smearing type
  integer, intent(in) :: stype
  !> smearing width
  real(dp), intent(in) :: swidth
  !> tolerance for occupied charged
  real(dp), intent(in) :: epschg
  !> Fermi energy
  real(dp), intent(out) :: efermi
  !> occupation numbers
  real(dp), intent(out) :: occ(nst,nkpt)

  integer, parameter :: maxiter = 100

  integer :: nvm, it, ip, ist
  real(dp) :: e0, e1, c, x, drift
  logical :: smaller, greater

  real(dp), external :: stheta, sdelta

  ! find Fermi level and occupations
  drift = 0.1_dp
  smaller = .false.; greater = .false.

  ! initial guess for energy interval from estimate of VBM and CBm
  nvm = nint( chg/maxocc)
  e0 = maxval( eval(nvm,:))
  e1 = minval( eval(nvm+1,:))
  do it = 1, maxiter
    efermi = 0.5*(e0 + e1)
    c = 0.d0
    do ip = 1, nkpt
      do ist = 1, nst
        x = (efermi - eval(ist,ip))/swidth 
        occ(ist,ip) = maxocc*stheta( stype, x)
        c = c + wkpt(ip)*occ(ist,ip)
      end do
    end do
    if( abs(c - chg) < epschg) exit
    if( c < chg) then
      smaller = .true.
      e0 = efermi
      e1 = e1 + (chg-c)*drift
    else
      greater = .true.
      e1 = efermi
      e0 = e0 + (chg-c)*drift
    end if
    if( smaller .and. greater) drift = 0.d0
  end do
  call terminate_if_false( it <= maxiter, &
    '(find_fermi) Fermi level search did not converge.')
end subroutine find_fermi

!> Find the Fermi level and occupation number response upon any perturbation for a given set of eigenenergies and k-points
!> using smearing and the bisection method.
!> The response of the total charge \(\rho\) in the system is given by
!> \[ \delta\rho = \sum_{n,{\bf k}} w_{\bf k} \, \delta f_{n{\bf k}} \;, \]
!> where \(\delta f_{n{\bf k}} = (\delta \epsilon_{\rm F} - \delta \epsilon_{n{\bf k}}) \delta(\epsilon_{\rm F} - \epsilon_{n{\bf k}})\)
!> is the occupation response of state \(\Psi_{n{\bf k}}\) and \(\delta\epsilon_{\rm F}\) is the Fermi energy response.
!> The Fermi energy response is iteratively refined using the bisection method until the resulting occupation number responses yield
!> the desired charge response within a given tolerance.
!> \(\delta(x)\) is approximated with a smooth function according to `stype` and `swidth`.
subroutine find_dfermi( nkpt, wkpt, nst, eval, deval, dchg, maxocc, efermi, stype, swidth, epschg, defermi, docc)
  use precision, only: dp
  use modmpi, only: terminate_if_false
  !> number of k-points
  integer, intent(in) :: nkpt
  !> k-point weights
  real(dp), intent(in) :: wkpt(nkpt)
  !> number of states
  integer, intent(in) :: nst
  !> eigenenergies
  real(dp), intent(in) :: eval(nst,nkpt)
  !> eigenenergy responses
  real(dp), intent(in) :: deval(nst,nkpt)
  !> charge response
  real(dp), intent(in) :: dchg
  !> maximum occupation per state
  real(dp), intent(in) :: maxocc
  !> Fermi energy
  real(dp), intent(in) :: efermi
  !> smearing type
  integer, intent(in) :: stype
  !> smearing width
  real(dp), intent(in) :: swidth
  !> tolerance for occupied charged
  real(dp), intent(in) :: epschg
  !> Fermi energy response
  real(dp), intent(out) :: defermi
  !> occupation number responses
  real(dp), intent(out) :: docc(nst,nkpt)

  integer, parameter :: maxiter = 100

  integer :: nvm, it, ip, ist
  real(dp) :: e0, e1, dc, x, dx, drift
  logical :: smaller, greater

  real(dp), external :: stheta, sdelta

  ! find Fermi level and occupations
  drift = 0.1d0
  smaller = .false.; greater = .false.

  ! initial guess of energy interval
  e0 = minval( deval)
  e1 = maxval( deval)
  do it = 1, maxiter
    defermi = 0.5*(e0 + e1)
    if( it == 1) defermi = 0.d0
    dc = 0.d0
    do ip = 1, nkpt
      do ist = 1, nst
        x = (efermi - eval(ist,ip))/swidth 
        dx = (defermi - deval(ist,ip))/swidth 
        docc(ist,ip) = maxocc*sdelta( stype, x)*dx
        dc = dc + wkpt(ip)*docc(ist,ip)
      end do
    end do
    if( abs(dc - dchg) < epschg) exit
    if( dc < dchg) then
      smaller = .true.
      e0 = defermi
      e1 = e1 + (dchg-dc)*drift
    else
      greater = .true.
      e1 = defermi
      e0 = e0 + (dchg-dc)*drift
    end if
    if( smaller .and. greater) drift = 0.d0
  end do
  call terminate_if_false( it <= maxiter, &
    '(find_dfermi) Fermi level response search did not converge.')
end subroutine find_dfermi

!> Evaluate the density of states at a given energy level using smearing given by
!> \[ D(\epsilon) = \sum_{n,{\bf k}} w_{\bf k}\, \delta(\epsilon_{n{\bf k}} - \epsilon) \;. \]
!> \(\delta(x)\) is approximated with a smooth function according to `stype` and `swidth`.
function dos_at_energy( nkpt, wkpt, nst, eval, maxocc, energy, stype, swidth) result(dos)
  use precision, only: dp
  !> number of k-points
  integer, intent(in) :: nkpt
  !> k-point weights
  real(dp), intent(in) :: wkpt(nkpt)
  !> number of states
  integer, intent(in) :: nst
  !> eigenenergies
  real(dp), intent(in) :: eval(nst,nkpt)
  !> maximum occupation per state
  real(dp), intent(in) :: maxocc
  !> energy at which to evaluate the DOS
  real(dp), intent(in) :: energy
  !> smearing type
  integer, intent(in) :: stype
  !> smearing width
  real(dp), intent(in) :: swidth
  !> DOS at given energy
  real(dp) :: dos

  integer :: ip, ist
  real(dp) :: x

  real(dp), external :: sdelta

  dos = 0._dp
  do ip = 1, nkpt
    do ist = 1, nst
      x = (eval(ist,ip) - energy)/swidth
      dos = dos + wkpt(ip)*sdelta( stype, x)/swidth
    end do
  end do
  dos = dos*maxocc
end function dos_at_energy

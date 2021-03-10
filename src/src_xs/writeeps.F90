! Copyright(C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_writeeps

  implicit none

  contains

    subroutine writeeps(iq, iop1, iop2, w, eps, fn)
      use mod_lattice, only: omega
      use mod_misc, only: version
      use modinput
      use modmpi
      use modbse, only: nk_bse
      use modxs, only: escale, unitout, unit1
      use modxs, only: ivgmt, vqlmt, vgcmt, vqcmt
      use modxs, only: sptclg, ivgigq
      use fox_wxml
      use m_getunit

      implicit none

      ! Arguments
      integer, intent(in) :: iq, iop1, iop2
      real(8), intent(in) :: w(:)
      complex(8), intent(in) :: eps(:)
      character(*), intent(in) :: fn

      ! Local variables
      type(xmlf_t), save :: xf
      character(256) :: buffer
      character(*), parameter :: thisnam = 'writeeps'
      integer :: n, iw, igqmt
      real(8), allocatable :: imeps(:), kkeps(:)
      character(128) :: fhdf5, group, ciq
      if(any(shape(w) .ne. shape(eps))) then
        write(unitout, '(a)') 'Error(' // thisnam // '): input&
          & arrays have diffenrent shape'
        call terminate
      end if

      n = size(w)

      allocate(imeps(n), kkeps(n))

      ! Kramers-Kronig transform imaginary part
      imeps(:) = aimag(eps(:))

      !write(*,*) "iop1,iop2,n", iop1, iop2, n
      !write(*,*) "w"
      !write(*,'(g10.3)') w
      !write(*,*) "imeps"
      !write(*,'(g10.3)') imeps
      Call kramkron(iop1, iop2, 1.d-8, n, w, imeps, kkeps)

      igqmt = ivgigq(ivgmt(1,iq),ivgmt(2,iq),ivgmt(3,iq),iq)

      Call getunit(unit1)
      Open(unit1, File=trim(fn), Action='write')
      write(unit1, '("# Macroscopic dielectric function epsm")')
      write(unit1, '("#")')
      write(unit1, '("# Momentum transfer Q=G+q in lattice cooridnates")')
      write(unit1, '("# G:",3i4)') ivgmt(1:3,iq) 
      write(unit1, '("# q:",3f12.7)') vqlmt(1:3,iq) 
      write(unit1, '("# Momentum transfer Q=G+q in Cartesian cooridnates")')
      write(unit1, '("# G:",3f12.7)') vgcmt(1:3,iq) 
      write(unit1, '("# q:",3f12.7)') vqcmt(1:3,iq) 
      write(unit1, '("# Norm2(G+q)",f12.7)') norm2(vgcmt(:,iq)+vqcmt(:,iq))
      write(unit1, '("#")')
      write(unit1, '("# Energy scale=", f12.6)') escale
      write(unit1, '("# Used broadening in scaled energy units:", f12.6)')&
        & escale*input%xs%broad
      write(unit1, '("#")')
      write(unit1, '("# Number of k-points=", i8)') nk_bse
      write(unit1, '("# Unit cell volume [au]=", f12.6)') omega
      write(unit1, '("# Coulomb potential v(Q) [au]=", f18.6)') sptclg(igqmt,iq)**2
      write(unit1, '("#")')
      write(unit1, '("#",a22,1x,a23,1x,a23,1x,a23)')&
        & "Frequency", "Re(epsm)", "Im(epsm)", "Re(epsm) form KKT"
      write(unit1, '(SP,E23.16,1x,E23.16,1x,E23.16,1x,E23.16)')&
        & (w(iw)*escale, eps(iw), kkeps(iw), iw=1, n)
      Close(unit1)

      deallocate(imeps, kkeps)
   end subroutine writeeps

end module m_writeeps

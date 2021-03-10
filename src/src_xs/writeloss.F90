! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_writeloss

  implicit none

  contains

    subroutine writeloss(iq, w, loss, fn)
      use modmpi
      Use FoX_wxml
      Use modmain, Only: version
      Use mod_lattice
      use constants, only: pi
      Use mod_charge_and_moment
      Use modxs
      Use m_getunit
      !Use m_writevars

      implicit none

      ! Arguments
      integer(4), intent(in) :: iq
      real(8), intent(in) :: w(:)
      real(8), intent(in) :: loss(:)
      character(*), intent(in) :: fn

      ! Local variables
      real(8) :: nval 
      character(*), parameter :: thisnam = 'writeloss'
      type(xmlf_t), save :: xf
      character(256) :: buffer
      integer(4) :: n, iw, igqmt, un

      if(any(shape(w) .ne. shape(loss))) then
         write(unitout, '(a)') 'Error(' // thisnam // '):&
           & Input arrays have different shapes'
         call terminate
      end if

      n = size(w,1)

      Call getunit(un)
      Open(un, File=trim(fn), Action='write')

      ! Include dynamical structure factor
      ! Dynamical structure factor; expression taken from Weissker, PRL 2006
      ! Units of dynamical structure factor are Hartree^-1
      ! Rescaling of units if electron volts are selected

      ! Get Gmt+qmt index
      igqmt = ivgigq(ivgmt(1, iq), ivgmt(2, iq), ivgmt(3, iq), iq)

      ! Compute avarage valence density = total valence charge / unti cell volume
      nval = chgval/omega

      write(un, '("#",1x,"Loss function L(Q,w) and dynamical structure factor S(Q,w)")')
      write(un, '("#")')
      write(un, '("#",1x,"L(Q,w)= -Im(1/eps_m(Q,w) [1]")')
      write(un, '("#",1x,"S(Q,w)= L(Q,w)*|Q|^2*/(4 pi^2 n_val) [Energy^-1]")')
      write(un, '("#")')
      write(un, '("# Momentum transfer Q=G+q in lattice coordinates")')
      write(un, '("# G:",3i4)') ivgmt(1:3,iq) 
      write(un, '("# q:",3f12.7)') vqlmt(1:3,iq) 
      write(un, '("# Momentum transfer Q=G+q in Cartesian coordinates")')
      write(un, '("# G:",3f12.7)') vgcmt(1:3,iq) 
      write(un, '("# q:",3f12.7)') vqcmt(1:3,iq) 
      write(un, '("# Norm2(G+q)",f12.7)') norm2(vgcmt(:,iq)+vqcmt(:,iq))
      write(un, '("#")')
      write(un, '("#",a22,1x,a23,1x,a23)')&
        & "Frequency/(eV/hbar)", "L(Q,w)", "S(Q,w)*1000"
      write(un, '(SP,E23.16,1x,E23.16,1x,E23.16)')&
        &(w(iw)*escale, loss(iw), &
        & 1000.0d0*loss(iw)/escale*(gqc(igqmt, iq)**2/(4.d0*pi**2*nval)), iw=1,  n)

      ! write relevant parameters to file
      !Call writevars(un, iq, iq)
      Close(un)
    End Subroutine writeloss

End Module m_writeloss

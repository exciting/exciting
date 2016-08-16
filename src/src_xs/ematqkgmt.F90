! Copyright(C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
subroutine ematqkgmt(iq, ik, igq, integrals)
  use modinput, only: input
  use mod_constants, only: zzero, zone, fourpi
  use mod_atoms, only: natmtot, nspecies, natoms, idxas
  use modxs, only: apwmaxsize, lomaxsize, istu2, istl2,&
                 & xiou, sfacgq, nst2, nst1,&
                 & apwsize, losize, cmtfun, cmtfun0,&
                 & cpumtaa
#ifdef USEOMP
  use omp_lib
#endif

  implicit none

  ! Arguments
  integer, intent(in) :: iq, ik, igq
  complex(8) :: integrals(apwmaxsize+lomaxsize,apwmaxsize+lomaxsize,natmtot)

  ! Local variables
  character(*), parameter :: thisnam = 'ematqkgmt'
  integer :: is, ia, ias
  integer :: lmax1, lmax3, ikt, zmsize, whichthread
  complex(8), allocatable :: zm(:,:)
  complex(8) :: prefactor
  real(8) :: cmt0, cmt1

#ifdef USEOMP
  whichthread=omp_get_thread_num()
#else
  whichthread=0
#endif

  ikt = ik
  lmax1 = input%xs%lmaxapwwf
  lmax3 = lmax1
  zmsize=apwmaxsize+lomaxsize

  allocate(zm(1:istu2-istl2+1,zmsize))

  xiou(:, :, igq) = zzero

  ! Loop over species and atoms
  do is = 1, nspecies
    do ia = 1, natoms(is)

      ias = idxas(ia, is)
      call timesec(cmt0)
      !---------------------------!
      !     apw-apw contribution  !
      !---------------------------!
      prefactor=fourpi*conjg(sfacgq(igq, ias, iq))

      call zgemm('n', 'n', nst2, apwsize(is)+losize(is), apwsize(is)+losize(is),&
        & zone, cmtfun(1,1,ias), nst2, integrals(1,1,ias), apwmaxsize+lomaxsize,&
        & zzero, zm, nst2)
      call zgemm('n', 't', nst1, nst2, apwsize(is)+losize(is),&
        & prefactor, cmtfun0(1,1,ias), nst1, zm, nst2, zone, xiou(1,1,igq), nst1)

      call timesec(cmt1)
      if(whichthread.eq.0) then
        cpumtaa = cpumtaa + cmt1 - cmt0
      endif

    ! End loop over species and atoms
    end do ! ia
  end do ! is

  deallocate(zm)
end subroutine ematqkgmt

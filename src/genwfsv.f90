
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genwfsv
! !INTERFACE:
subroutine genwfsv(tocc,ngp,igpig,evalsvp,apwalm,evecfv,evecsv,wfmt,wfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tocc    : .true. if only occupied wavefunctions are required (in,logical)
!   ngp     : number of G+p-vectors (in,integer)
!   igpig   : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   evalsvp : second-variational eigenvalue for every state (in,real(nstsv))
!   apwalm  : APW matching coefficients
!             (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evecfv  : first-variational eigenvectors (in,complex(nmatmax,nstfv))
!   evecsv  : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   wfmt    : muffin-tin part of the wavefunctions for every state in spherical
!             coordinates (out,complex(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
!   wfir    : interstitial part of the wavefunctions for every state
!             (out,complex(ngrtot,nspinor,nstsv))
! !DESCRIPTION:
!   Calculates the second-variational spinor wavefunctions in both the
!   muffin-tin and interstitial regions for every state of a particular
!   $k$-point. The wavefunctions in both regions are stored on a real-space
!   grid. A coarse radial mesh is assumed in the muffin-tins with with angular
!   momentum cut-off of {\tt lmaxvr}. If {\tt tocc} is {\tt .true.}, then only
!   the occupied states (those below the Fermi energy) are calculated.
!
! !REVISION HISTORY:
!   Created November 2004 (Sharma)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tocc
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: evalsvp(nstsv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(out) :: wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(out) :: wfir(ngrtot,nspinor,nstsv)
! local variables
integer ispn,is,ia,ias
integer i,j,n,ist,igp,ifg
real(8) t1
complex(8) zt1
! allocatable arrays
logical, allocatable :: done(:)
complex(8), allocatable :: wfmt1(:,:)
complex(8), allocatable :: wfmt2(:,:,:)
allocate(done(nstfv))
allocate(wfmt1(lmmaxvr,nrcmtmax))
allocate(wfmt2(lmmaxvr,nrcmtmax,nstfv))
!--------------------------------!
!     muffin-tin wavefunction    !
!--------------------------------!
do is=1,nspecies
  n=lmmaxvr*nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    done(:)=.false.
    do j=1,nstsv
      if ((.not.tocc).or.((tocc).and.(evalsvp(j).lt.efermi))) then
        if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
          wfmt(:,:,ias,:,j)=0.d0
          i=0
          do ispn=1,nspinor
            do ist=1,nstfv
              i=i+1
              zt1=evecsv(i,j)
              if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
                if (.not.done(ist)) then
                  call wavefmt(lradstp,lmaxvr,is,ia,ngp,apwalm,evecfv(1,ist), &
                   lmmaxvr,wfmt1)
! convert from spherical harmonics to spherical coordinates
                  call zgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,zone,zbshtapw, &
                   lmmaxapw,wfmt1,lmmaxvr,zzero,wfmt2(1,1,ist),lmmaxvr)
                  done(ist)=.true.
                end if
! add to spinor wavefunction
                call zaxpy(n,zt1,wfmt2(1,1,ist),1,wfmt(1,1,ias,ispn,j),1)
              end if
! end loop over first-variational states
            end do
! end loop over spin
          end do
        else
! spin-unpolarised wavefunction
          call wavefmt(lradstp,lmaxvr,is,ia,ngp,apwalm,evecfv(1,j),lmmaxvr, &
           wfmt1)
! convert from spherical harmonics to spherical coordinates
          call zgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,zone,zbshtapw,lmmaxapw, &
           wfmt1,lmmaxvr,zzero,wfmt(1,1,ias,1,j),lmmaxvr)
        end if
      end if
! end loop over second-variational states
    end do
! end loops over atoms and species
  end do
end do
!-----------------------------------!
!     interstitial wavefunction     !
!-----------------------------------!
t1=1.d0/sqrt(omega)
do j=1,nstsv
  wfir(:,:,j)=0.d0
  if ((.not.tocc).or.((tocc).and.(evalsvp(j).lt.efermi))) then
    if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
      i=0
      do ispn=1,nspinor
        do ist=1,nstfv
          i=i+1
          zt1=evecsv(i,j)
          if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
            zt1=t1*zt1
            do igp=1,ngp
              ifg=igfft(igpig(igp))
              wfir(ifg,ispn,j)=wfir(ifg,ispn,j)+zt1*evecfv(igp,ist)
            end do
          end if
        end do
      end do
    else
! spin-unpolarised wavefunction
      do igp=1,ngp
        ifg=igfft(igpig(igp))
        wfir(ifg,1,j)=t1*evecfv(igp,j)
      end do
    end if
! Fourier transform wavefunction to real-space
    do ispn=1,nspinor
      call zfftifc(3,ngrid,1,wfir(1,ispn,j))
    end do
  end if
end do
deallocate(done,wfmt1,wfmt2)
return
end subroutine
!EOC


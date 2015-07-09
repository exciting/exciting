!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
!BOP
! !ROUTINE: calcpmatgw
! !INTERFACE:


Subroutine pmatradcor
! !USES:
      Use modmain
      Use modinput
      Use modxs
      Use modgw, only: ncmax,nclm,ncore,ucore,ripacor,ripcora

! !DESCRIPTION:
!
!
! !REVISION HISTORY:
!
! Created based on src_xs/pmatrad.f90 (S.Sagmeister) July 2011 (DIN)
!
! !LOCAL VARIABLES:
      Implicit None
!     local variables
      Integer :: is, ia, ias, nr, ir
      Integer :: l1, m1, lm1, l3, m3, lm3
      Integer :: icor, io, j
!     automatic arrays
      Real(8) :: r2(nrmtmax), fr(nrmtmax), gr(nrmtmax), cf(3,nrmtmax)
!     allocatable arrays
      Real(8), Allocatable :: fapw(:), dapwfr(:,:,:,:,:)
      Real(8), Allocatable :: fcor(:), dcorfr(:,:,:,:,:)
!EOP
!BOC

!     allocate local arrays for radial derivatives
      Allocate(fapw(nrmtmax))
      Allocate(dapwfr(lmmaxapw,nrmtmax,3,apwordmax,lmmaxapw))
      fapw(:) = 0.d0
      dapwfr(:,:,:,:,:) = 0.d0

      Allocate(fcor(nrmtmax))
      Allocate(dcorfr(lmmaxapw,nrmtmax,3,ncmax,nclm))
      fcor(:) = 0.d0
      dcorfr(:,:,:,:,:) = 0.d0   
      
      ripacor(:,:,:,:,:,:) = 0.d0
      ripcora(:,:,:,:,:,:) = 0.d0

!     begin loop over species
      Do is = 1, nspecies
         nr = nrmt(is)
         Do ir = 1, nr
            ! calculate r^2
            r2(ir) = spr(ir,is)**2
         End Do
        ! begin loop over atoms
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            !--------------------!
            !     derivatives    !
            !--------------------!
            ! APW functions
            do l1 = 0, input%groundstate%lmaxapw
               do m1 = -l1, l1
                  lm1 = idxlm(l1,m1)
                  do io = 1, apword(l1,is)
                     fapw(1:nr) = apwfr(1:nr,1,io,l1,ias)
                     call gradzfmtr(input%groundstate%lmaxapw, nr, &
                    & spr(1,is), l1, m1, lmmaxapw, nrmtmax, fapw,  &
                    & dapwfr(1,1,1,io,lm1))
                  end do
               end do
            end do
            ! core functions
            do icor = 1, ncore(is)
               l1 = spl(icor,is)
               do m1 = -l1, l1
                  lm1 = idxlm(l1,m1)
                  fcor(1:nr) = ucore(1:nr,1,icor,ias)
                  call gradzfmtr(input%groundstate%lmaxapw, nr,  &
                 &  spr(1,is), l1, m1, lmmaxapw, nrmtmax, fcor,  &
                 &  dcorfr(1,1,1,icor,lm1))
               end do
            end do
            !----------------------------!
            !     APW-core-orbital       !
            !----------------------------!
            do l1 = 0, input%groundstate%lmaxapw
               do m1 = -l1, l1
                  lm1 = idxlm(l1,m1)
                  do io = 1, apword (l1, is)
                     do icor = 1, ncore(is)
                        l3 = spl(icor,is)
                        do m3 = -l3, l3
                           lm3 = idxlm(l3,m3)
                           do j = 1, 3
                              fr(1:nr) = r2(1:nr)*apwfr(1:nr,1,io,l1,ias)* &
                             &           dcorfr(lm1,1:nr,j,icor,lm3)
                              call fderiv(-1, nr, spr(1,is), fr, gr, cf)
                              ripacor(io,lm1,icor,lm3,ias,j) = gr(nr)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
           !----------------------------!
           !     core-orbital-APW      !
           !----------------------------!
               Do icor = 1, ncore(is)
                  l1 = spl(icor,is)
                  Do m1 = -l1, l1
                     lm1 = idxlm(l1,m1)
                     Do l3 = 0, input%groundstate%lmaxapw
                        Do m3 = -l3, l3
                           lm3 = idxlm(l3,m3)
                           Do io = 1, apword(l3,is)
                              Do j = 1, 3
                                 fr(1:nr) = r2(1:nr)*ucore(1:nr,1,icor,ias)* &
                                &           dapwfr(lm1,1:nr,j,io,lm3)
                                 Call fderiv(-1, nr, spr(1,is), fr, gr, cf)
                                 ripcora(icor,lm1,io,lm3,ias,j) = gr(nr)
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
        ! end loops over atoms and species
         End Do
      End Do
  ! deallocate
      deallocate (fapw, dapwfr)
      deallocate (fcor, dcorfr)
End Subroutine pmatradcor

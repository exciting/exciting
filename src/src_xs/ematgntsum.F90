! Copyright(C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine ematgntsum(iq, igq, integrals)
  use modinput, only: input
  use mod_misc, only: filext
  use mod_constants, only: zzero, zil
  use mod_muffin_tin, only: idxlm
  use mod_atoms, only: natmtot, nspecies, natoms, idxas
  use mod_APW_LO, only: lolmax, apwordmax, nlomax, apword,&
                      & nlorb, lorbl
  use modxs, only: apwmaxsize, apwsize, lomaxsize, ylmgq,&
                 & xsgnt, rilolo, riaa, riloa
  use m_findgntn0, only: l1map, l2map, l3map,&
                       & m1map, m2map, m3map,&
                       & l1shape, l2shape, l3shape,&
                       & m1shape, m2shape, m3shape
  use m_getunit

  implicit none

  ! Arguments
  integer, intent(in) :: iq, igq
  complex(8) :: integrals(apwmaxsize+lomaxsize, apwmaxsize+lomaxsize, natmtot)

  ! Local variables
  integer :: is, ia, ias, iaug1, iaug2
  integer :: l1, l2, l3, m2, lm2
  integer :: ilo, ilo1, ilo2, io, io1, io2
  integer :: lmax1, lmax2, lmax3, lmmax1, lmmax2, lmmax3
  integer :: u1, u2, u3, u4
  integer :: m1, m3, lm1, lm3, cl1, cm1, cl2, cm2, cl3, cm3
  complex(8), dimension(:,:,:,:), allocatable :: intrgaa, intrgloa, intrglolo, intrgalo

  ! Set lm related local variables
  lmax1 = max(input%xs%lmaxapwwf, lolmax)
  lmax2 = input%xs%lmaxemat
  ! lmax1 and lmax3 should be the same
  lmax3 = lmax1
  lmmax1 = (lmax1+1) ** 2
  lmmax2 = (lmax2+1) ** 2
  lmmax3 = (lmax3+1) ** 2

  ! Allocate arrays for radial integrals and Bessel functions
  allocate(intrgaa(lmmax1, apwordmax, lmmax3, apwordmax))
  allocate(intrgloa(-lolmax:lolmax, nlomax, lmmax3, apwordmax))
  allocate(intrgalo(-lolmax:lolmax, nlomax, lmmax3, apwordmax))
  allocate(intrglolo(-lolmax:lolmax, nlomax, -lolmax:lolmax, nlomax))

  integrals=zzero

  ! Debug output
  if(input%xs%dbglev .gt. 2) then
    ! Apw-apw
    call getunit(u1)
    open(unit=u1, file='IRADGAUNTaa'//filext, form='formatted', action='write', status='replace')
    write(u1, '(a)') 'igq, ias, lm1, io1, lm3, io2,   intrgaa'
    write(u1, '(a)') '------------------------------------------------------'
    ! Apw-lo
    call getunit(u2)
    open(unit=u2, file='IRADGAUNTalo'//filext, form='formatted', action='write', status='replace')
    write(u2, '(a)') 'igq, ias, m3, ilo, lm1, io,     intrgalo'
    write(u2, '(a)') '------------------------------------------------------'
    ! Lo-apw
    call getunit(u3)
    open(unit=u3, file='IRADGAUNTloa'//filext, form='formatted', action='write', status='replace')
    write(u3, '(a)') 'igq, ias, m1, ilo, lm3, io,     intrgloa'
    write(u3, '(a)') '------------------------------------------------------'
    ! Lo-lo
    call getunit(u4)
    open(unit=u4, file='IRADGAUNTlolo'//filext, form='formatted', action='write', status='replace')
    write(u4, '(a)') 'igq, ias, m1, ilo1, m3, ilo2,   intrglolo'
    write(u4, '(a)') '------------------------------------------------------'
  end if

  ! Begin loop over species
  do is = 1, nspecies
    ! Begin loop over atoms
    do ia = 1, natoms(is)

      intrgaa(:, :, :, :) = zzero
      intrgloa(:, :, :, :) = zzero
      intrgalo(:, :, :, :) = zzero
      intrglolo(:, :, :, :) = zzero

      ias = idxas(ia, is)
      !---------------------------!
      !     apw-apw integrals     !
      !---------------------------!
      do cl1 = 1, l1shape
        l1 = l1map(cl1)
        do cm1 = 1, m1shape(l1)
          m1 = m1map(l1, cm1)
          lm1 = idxlm(l1, m1)
          do io1 = 1, apword(l1, is)
            do cl2 = 1, l2shape(l1, m1)
              l3 = l2map(l1, m1, cl2)
              do cm2 = 1, m2shape(l1, m1, l3)
                m3 = m2map(l1, m1, l3, cm2)
                lm3 = idxlm(l3, m3)
                do io2 = 1, apword(l3, is)
                  do cl3 = 1, l3shape(l1, m1, l3, m3)
                    l2 = l3map(l1, m1, l3, m3, cl3)
                    do cm3 = 1, m3shape(l1, m1, l3, m3, l2)
                      m2 = m3map(l1, m1, l3, m3, l2, cm3)
                      lm2 = idxlm(l2, m2)
                      intrgaa(lm1, io1, lm3, io2) =&
                        & intrgaa(lm1, io1, lm3, io2)&
                        &+ conjg(zil(l2)) * riaa(l1, io1, l3, io2, l2, ias, igq)&
                        &* conjg(ylmgq(lm2, igq, iq)) * xsgnt(lm1, lm2, lm3)
                    end do
                  end do
                  ! Debug output
                  if(input%xs%dbglev .gt. 2) then
                     write(u1, '(6i5, 2g18.10)') igq, ias, &
                       & lm1, io1, lm3, io2, intrgaa(lm1, io1, lm3, io2)
                  end if
                end do
              end do
            end do
          end do
        end do
      end do

      iaug2=0
      do l3=0, input%xs%lmaxapwwf
        do m3=-l3, l3
          do io2=1, apword(l3, is)
            iaug2=iaug2+1
            lm3=idxlm(l3, m3)
            iaug1=0
            do l1=0, input%xs%lmaxapwwf
              do m1=-l1, l1
                  do io1=1, apword(l1, is)
                    iaug1=iaug1+1
                    lm1=idxlm(l1, m1)
                    integrals(iaug2, iaug1, ias)=intrgaa(lm1, io1, lm3, io2)
                  end do
                end do
              end do
          end do
        end do
      end do

      !-------------------------------------!
      !     apw-local-orbital integrals     !
      !-------------------------------------!
      do cl1 = 1, l1shape
        l1 = l1map(cl1)
        do cm1 = 1, m1shape(l1)
          m1 = m1map(l1, cm1)
          lm1 = idxlm(l1, m1)
          do io = 1, apword(l1, is)
            do ilo = 1, nlorb(is)
              l3 = lorbl(ilo, is)
              do cm2 = 1, m2shape(l1, m1, l3)
                m3 = m2map(l1, m1, l3, cm2)
                lm3 = idxlm(l3, m3)
                do cl3 = 1, l3shape(l1, m1, l3, m3)
                  l2 = l3map(l1, m1, l3, m3, cl3)
                  do cm3 = 1, m3shape(l1, m1, l3, m3, l2)
                    m2 = m3map(l1, m1, l3, m3, l2, cm3)
                    lm2 = idxlm(l2, m2)
                    intrgalo(m3, ilo, lm1, io) = intrgalo(m3, ilo, lm1, io)&
                      &+ conjg(zil(l2)) * riloa(ilo, l1, io, l2, ias, igq)&
                      &* conjg(ylmgq(lm2, igq, iq)) * xsgnt(lm1, lm2, lm3)
                  end do
                end do
                ! Debug output
                if(input%xs%dbglev .gt. 2) then
                   write(u2, '(6i5, 2g18.10)') igq, ias, &
                     & m3, ilo, lm1, io, intrgalo(m3, ilo, lm1, io)
                end if
              end do
            end do
          end do
        end do
      end do

      iaug2=0
      do ilo = 1, nlorb(is)
        l3 = lorbl(ilo, is)
        do m3=-l3, l3
          iaug2=iaug2+1
          lm3=idxlm(l3, m3)
         
          iaug1=0
          do l1=0, input%xs%lmaxapwwf
            do m1=-l1, l1
              do io=1, apword(l1, is)
                iaug1=iaug1+1
                lm1=idxlm(l1, m1)
                integrals(apwsize(is)+iaug2, iaug1, ias)=intrgalo(m3, ilo, lm1, io)
              end do
            end do
           end do

        end do
      end do

      !-------------------------------------!
      !     local-orbital-apw integrals     !
      !-------------------------------------!
      do ilo = 1, nlorb(is)
        l1 = lorbl(ilo, is)
        do cm1 = 1, m1shape(l1)
          m1 = m1map(l1, cm1)
          lm1 = idxlm(l1, m1)
          do cl2 = 1, l2shape(l1, m1)
            l3 = l2map(l1, m1, cl2)
            do cm2 = 1, m2shape(l1, m1, l3)
              m3 = m2map(l1, m1, l3, cm2)
              lm3 = idxlm(l3, m3)
              do io = 1, apword(l3, is)
                do cl3 = 1, l3shape(l1, m1, l3, m3)
                  l2 = l3map(l1, m1, l3, m3, cl3)
                  do cm3 = 1, m3shape(l1, m1, l3, m3, l2)
                    m2 = m3map(l1, m1, l3, m3, l2, cm3)
                    lm2 = idxlm(l2, m2)
                    intrgloa(m1, ilo, lm3, io) = intrgloa(m1, ilo, lm3, io)&
                      &+ conjg(zil(l2)) * riloa(ilo, l3, io, l2, ias, igq)&
                      &* conjg(ylmgq(lm2, igq, iq)) * xsgnt(lm1, lm2, lm3)
                  end do
                end do
                ! Debug output
                if(input%xs%dbglev .gt. 2) then
                   write(u3, '(6i5, 2g18.10)') igq, ias,&
                     & m1, ilo, lm3, io, intrgloa(m1, ilo, lm3, io)
                end if
              end do
            end do
          end do
        end do
      end do

      iaug2=0
      do ilo = 1, nlorb(is)
        l3 = lorbl(ilo, is)
        do m3=-l3, l3
          iaug2=iaug2+1
          lm3=idxlm(l3, m3)

          iaug1=0
          do l1=0, input%xs%lmaxapwwf
            do m1=-l1, l1
              do io=1, apword(l1, is)
                iaug1=iaug1+1
                lm1=idxlm(l1, m1)
                integrals(iaug1, apwsize(is)+iaug2, ias)=intrgloa(m3, ilo, lm1, io)
              end do
            end do
           end do

        end do
      end do

      !-----------------------------------------------!
      !     local-orbital-local-orbital integrals     !
      !-----------------------------------------------!
      do ilo1 = 1, nlorb(is)
        l1 = lorbl(ilo1, is)
        do cm1 = 1, m1shape(l1)
          m1 = m1map(l1, cm1)
          lm1 = idxlm(l1, m1)
          do ilo2 = 1, nlorb(is)
            l3 = lorbl(ilo2, is)
            do cm2 = 1, m2shape(l1, m1, l3)
              m3 = m2map(l1, m1, l3, cm2)
              lm3 = idxlm(l3, m3)
              do cl3 = 1, l3shape(l1, m1, l3, m3)
                l2 = l3map(l1, m1, l3, m3, cl3)
                do cm3 = 1, m3shape(l1, m1, l3, m3, l2)
                  m2 = m3map(l1, m1, l3, m3, l2, cm3)
                  lm2 = idxlm(l2, m2)
                  intrglolo(m1, ilo1, m3, ilo2) = intrglolo(m1, ilo1, m3, ilo2)&
                    &+ conjg(zil(l2)) * rilolo(ilo1, ilo2, l2, ias, igq)&
                    &* conjg(ylmgq(lm2, igq, iq)) * xsgnt(lm1, lm2, lm3)
                end do
              end do
              ! Debug output
              if(input%xs%dbglev .gt. 2) then
                 write(u4, '(6i5, 2g18.10)') igq, ias, m1, &
                   & ilo1, m3, ilo2, intrglolo(m1, ilo1, m3, ilo2)
              end if
            end do
          end do
        end do
      end do

      iaug2=0
      do ilo2 = 1, nlorb(is)
        l3 = lorbl(ilo2, is)
        do m3=-l3, l3
          iaug2=iaug2+1

          iaug1=0
          do ilo1 = 1, nlorb(is)
            l1 = lorbl(ilo1, is)
            do m1=-l1, l1
              iaug1=iaug1+1
              integrals(apwsize(is)+iaug2, apwsize(is)+iaug1, ias)=&
                & intrglolo(m1, ilo1, m3, ilo2)
            end do
          end do

        end do
      end do

    ! End loops over atoms and species
    end do
  end do

  ! Deallocate
  if(allocated(intrgaa)) deallocate(intrgaa)
  if(allocated(intrgloa)) deallocate(intrgloa)
  if(allocated(intrgalo)) deallocate(intrgalo)
  if(allocated(intrglolo)) deallocate(intrglolo)

  ! Debug code
  if(input%xs%dbglev .gt. 2) then
    ! Close files
    close(u1)
    close(u2)
    close(u3)
    close(u4)
  end if
      
end subroutine ematgntsum

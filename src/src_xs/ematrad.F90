! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: ematrad
! !INTERFACE:
subroutine ematrad(iq)
! !USES:
  use mod_misc, only: filext
  use mod_APW_LO, only: lolmax, apwordmax, nlomax, apword,&
                      & apwfr, nlorb, lorbl, lofr
  use mod_atoms, only: natmtot, nspecies, spr, natoms, idxas
  use mod_muffin_tin, only: nrmtmax, nrmt
  use modinput, only: input
  use modxs, only: riaa, riloa, rilolo, ngq, gqc
  use m_getunit

! !DESCRIPTION:
! This routine is used in the construction of the plane wave matrix elements.
! It calculates the involved radial integrals for the APW-APW, APW-LO and LO-LO
! combinations. (Compare the R quantities section 8.1 of Sagmeisters thesis.)
!
! !REVISION HISTORY:
! Added to documentation scheme. (Aurich)
!EOP
!BOC

  implicit none

  ! Arguments
  integer, intent(in) :: iq

  ! Local variables
  integer :: is, ia, ias, nr, ir, igq
  integer :: l1, l2, l3,lio,liomax
  integer :: ilo, ilo1, ilo2, io, io1, io2
  real(8) :: t1
  integer :: lmax1, lmax2, lmax3
  integer :: u11, u22, u33
  ! Automatic arrays
  real(8) :: r2(nrmtmax), fr(nrmtmax), gr(nrmtmax), cf(3, nrmtmax)
  ! Allocatable arrays
  real(8), allocatable :: jl(:, :), jhelp(:)
  integer, allocatable :: lio2l(:), lio2io(:)

  lmax1 = max(input%xs%lmaxapwwf, lolmax)
  lmax2 = input%xs%lmaxemat

  ! lmax1 and lmax3 should be the same!
  lmax3 = lmax1

  ! Allocate arrays for radial integrals and Bessel functions
  if(allocated(riaa)) deallocate(riaa)
  if(allocated(riloa)) deallocate(riloa)
  if(allocated(rilolo)) deallocate(rilolo)
  allocate(riaa(0:lmax1, apwordmax, 0:lmax3, apwordmax, 0:lmax2, natmtot, ngq(iq)))
  allocate(riloa(nlomax, 0:lmax3, apwordmax, 0:lmax2, natmtot, ngq(iq)))
  allocate(rilolo(nlomax, nlomax, 0:lmax2, natmtot, ngq(iq)))

  ! Allocate temporary arrays
  allocate(jl(nrmtmax,0:lmax2))
  allocate(jhelp(0:lmax2))
  allocate(lio2l((lmax1+1)*apwordmax))
  allocate(lio2io((lmax1+1)*apwordmax))

  jl(:, :) = 0.d0
  jhelp(:) = 0.d0

  ! Zero arrays for radial integrals
  riaa(:, :, :, :, :, :, :) = 0.d0
  riloa(:, :, :, :, :, :) = 0.d0
  rilolo(:, :, :, :, :) = 0.d0

  if(input%xs%dbglev .gt. 1) then

    ! Apw-apw
    call getunit(u11)
    open(unit=u11, file='IRADaa'//filext, form='formatted', action='write',&
      & status='replace')
    write(u11, '(a)') 'igq, ias, l1, io1, l3, io2, l2	 iraa'
    write(u11, '(a)') '-----------------------------------------------------'
    ! Lo-apw
    call getunit(u22)
    open(unit=u22, file='IRADalo'//filext, form='formatted', action='write',&
      & status='replace')
         write(u22, '(a)') 'igq, ias, ilo, l1, l3, io, l2,	 iralo'
         write(u22, '(a)') '-----------------------------------------------------'
    ! Lo-lo
    call getunit(u33)
    open(unit=u33, file='IRADlolo'//filext, form='formatted', action='write',&
      & status='replace')
    write(u33, '(a)') 'igq, ias, ilo1, l1, ilo2, l3, l2,   irlolo'
    write(u33, '(a)') '-----------------------------------------------------'

  end if

  ! Begin loop over G+q vectors
  do igq = 1, ngq(iq)

    ! Begin loop over species
    do is = 1, nspecies

      nr = nrmt(is)
      do ir = 1, nr
        ! Calculate r^2
        r2(ir) = spr(ir, is) ** 2
        ! Calculate spherical Bessel functions of first kind j_l(|g+q|r_a)
        call sbessel(lmax2, gqc(igq, iq)*spr(ir, is), jhelp)
        jl(ir,:) = jhelp(:)
      end do

      lio=0
      do l1 = 0, lmax1
        do io1 = 1, apword(l1, is)
          lio=lio+1
          lio2l(lio)=l1
          lio2io(lio)=io1
        end do
      end do
      liomax=lio

      ! Begin loop over atoms
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        !----------------!
        !     apw-apw    !
        !----------------!
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(lio,l1,io1,l3,io2,ir,t1,fr,cf,gr)
!$OMP DO
#endif
        do lio=1,liomax
          l1=lio2l(lio)
          io1=lio2io(lio)
          do l3 = 0, lmax3
            do io2 = 1, apword(l3, is)
              do l2 = 0, lmax2
                do ir = 1, nr
                  t1 = apwfr(ir, 1, io1, l1, ias) * apwfr(ir, 1, io2, l3, ias) * r2(ir)
                  fr(ir) = t1 * jl(ir,l2)
                end do
                call fderiv(-1, nr, spr(1, is), fr, gr, cf)
                riaa(l3, io2, l1, io1, l2, ias, igq) = gr(nr)
              end do ! l2
            end do ! io2
          end do ! l3
        end do ! l1
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

        !----------------------------!
        !     local-orbital-apw      !
        !----------------------------!
        do ilo = 1, nlorb(is)
          l1 = lorbl(ilo, is)
          do l3 = 0, lmax3
            do io = 1, apword(l3, is)
              do l2 = 0, lmax2
                do ir = 1, nr
                  t1 = lofr(ir, 1, ilo, ias) * apwfr(ir, 1, io, l3, ias) * r2(ir)
                  fr(ir) = t1 * jl(ir,l2)
                end do
                call fderiv(-1, nr, spr(1, is), fr, gr, cf)
                riloa(ilo, l3, io, l2, ias, igq) = gr(nr)
              end do ! l2
            end do ! io
          end do ! l3
        end do ! ilo
        !------------------------------------!
        !     local-orbital-local-orbital    !
        !------------------------------------!
        do ilo1 = 1, nlorb(is)
          l1 = lorbl(ilo1, is)
          do ilo2 = 1, nlorb(is)
            l3 = lorbl(ilo2, is)
            do l2 = 0, lmax2
              do ir = 1, nr
                t1 = lofr(ir, 1, ilo1, ias) * lofr(ir, 1, ilo2, ias) * r2(ir)
                fr(ir) = t1 * jl(ir,l2)
              end do
              call fderiv(-1, nr, spr(1, is), fr, gr, cf)
              rilolo(ilo1, ilo2, l2, ias, igq) = gr(nr)
            end do ! l2
          end do ! ilo2
        end do ! ilo1
        
        !****************************************
        ! Debugging segment with output to files
        !****************************************
        if(input%xs%dbglev .gt. 1) then
          !----------------!
          !     apw-apw    !
          !----------------!
          do l1 = 0, lmax1
            do io1 = 1, apword(l1, is)
              do l3 = 0, lmax3
                do io2 = 1, apword(l3, is)
                  do l2 = 0, lmax2
                    write(u11, '(7i5, g18.10)') igq, ias, l1, io1, l3, io2, l2,&
                      & riaa(l1, io1, l3, io2, l2, ias, igq)
                  end do
                end do ! io2
              end do ! l3
            end do ! io1
          end do ! l1
          !----------------------------!
          !     local-orbital-apw      !
          !----------------------------!
          do ilo = 1, nlorb(is)
            l1 = lorbl(ilo, is)
            do l3 = 0, lmax3
              do io = 1, apword(l3, is)
                do l2 = 0, lmax2
                  write(u22, '(7i5, g18.10)') igq, ias, ilo, l1, l3, io, l2,&
                    & riloa(ilo, l3, io, l2, ias, igq) 
                end do ! l2
              end do ! io
            end do ! l3
          end do ! ilo
          !------------------------------------!
          !     local-orbital-local-orbital    !
          !------------------------------------!
          do ilo1 = 1, nlorb(is)
            l1 = lorbl(ilo1, is)
            do ilo2 = 1, nlorb(is)
              l3 = lorbl(ilo2, is)
              do l2 = 0, lmax2
                write(u33, '(7i5, g18.10)') igq, ias, ilo1, l1, ilo2, l3, l2,&
                  & rilolo(ilo1, ilo2, l2, ias, igq)
              end do ! l2
            end do ! ilo2
          end do ! ilo1
        !**************************
        ! End of debugging segment
        !**************************
        endif
      ! End loops over atoms and species
      end do
    end do
  ! End loop over G+q vectors
  end do

  ! Deallocate
  deallocate(jl, jhelp)
  deallocate(lio2l,lio2io)

  if(input%xs%dbglev .gt. 1) then
    ! Close files
    close(u11)
    close(u22)
    close(u33)
  end if
end subroutine ematrad
!EOC

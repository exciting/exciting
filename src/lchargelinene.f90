!BOP
! !ROUTINE: lchargelinene
! !INTERFACE:
subroutine lchargelinene
! !USES:
      use modinput
      use modmain
      use modmpi, only: rank
! !DESCRIPTION:
!
! Calculate the "Optimal Energy Parameters" by S.Bluegel
!
! !REVISION HISTORY:
!   
!EOP
!BOC
  implicit none
  ! local variables
  integer :: ik,l,m,lm
  integer :: lmax,lmmax
  integer :: is,ia,ias,ist
  integer :: ispn
  integer :: ilo, io1, io2, ja, jas
  real(8) :: t1,t2,t3,t4,t5

  Real (8), Allocatable :: el(:,:)
  Real (8), Allocatable :: lcharge (:, :, :, :, :)
  Real (8), Allocatable :: lstate (:, :, :)
  Complex (8), Allocatable :: evecfv(:,:,:)
  Complex (8), Allocatable :: evecsv(:,:)
  Complex (8), Allocatable :: dmat (:, :, :, :, :)
  Complex (8), Allocatable :: apwalm (:, :, :, :, :)
  Logical :: done (natmmax)

  ! set lmax for maximum accuracy
  lmax=input%groundstate%lmaxapw
  lmmax = (lmax+1) ** 2
  
  Allocate (evecfv(nmatmax, nstfv, nspnfv))
  Allocate (evecsv(nstsv, nstsv))
  Allocate (dmat(lmmax, lmmax, nspinor, nspinor, nstsv))
  Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
  Allocate (lcharge(0:lmax, nspinor, natmtot, nstsv, nkpt))
  Allocate (el(0:input%groundstate%lmaxapw,natmtot))
  Allocate (lstate(0:input%groundstate%lmaxapw,nstsv,natmtot))
  
! loop over k-points
  Do ik = 1, nkpt
    
    Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
    Call getevecsv (vkl(:, ik), evecsv)
  
    ! find the matching coefficients
    Do ispn = 1, nspnfv
      Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :,ispn, ik), &
     &  sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, ispn))
    End Do ! ispn

    Do is = 1, nspecies
      done (:) = .False.
      Do ia = 1, natoms (is)
        If ( .Not. done(ia)) Then
          ias = idxas (ia, is)
          
! generate the density matrix
          Call gendmat (.True., .True., 0, lmax, is, ia, ngk(:, ik), &
         &  apwalm, evecfv, evecsv, lmmax, dmat)

! determine "l-like" charges
          Do ist = 1, nstsv
            Do ispn = 1, nspinor
              Do l = 0, lmax
                t1 = 0.d0
                Do m = - l, l
                  lm = idxlm (l, m)
                  t1 = t1 + dble(dmat(lm, lm, ispn, ispn, ist))
                End Do
                lcharge (l, ispn, ias, ist, ik) = t1
              End Do
            End Do ! ispn
          End Do ! ist
          done (ia) = .True.

! copy to equivalent atoms
          Do ja = 1, natoms (is)
            If ((.not.done(ja)) .and.(eqatoms(ia,ja,is))) Then
              jas = idxas(ja, is)
              do ist = 1, nstsv
                  lcharge(:,:,jas,ist,ik) = lcharge(:,:,ias,ist,ik)
              end do
              done(ja) = .True.
            End If
          End Do
          
        End If ! done
      End Do ! ia
    End Do ! is
  End Do ! ik
  Deallocate (dmat, apwalm)
  Deallocate (evecfv, evecsv)

  if (rank==0) then
    open(32,file='LCHARGE.OUT',action='write')
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        write(32,*) 'ATOM ', ias
        write(32,*) '# state    l=0     l=1     l=2     l=3     l=4'
        do ist = 1, nstsv
          do l = 0, input%groundstate%lmaxapw
            t1 = 0.d0
            do ik = 1, nkpt
              t1 = t1 + wkpt(ik)*lcharge(l,1,ias,ist,ik)
            end do
            lstate(l,ist,ias) = t1
          end do
          write(32,'(i8,5f12.6)') ist, (lstate(l,ist,ias), l = 0, 4)
        end do
       write(32,*)
     end do
    end do
    close(32)
  end if

!-----------------------
! Calculate l-energies
!-----------------------

! for each atom
  do is = 1, nspecies
    done (:) = .False.
    do ia = 1, natoms (is)
      if ( .Not. done(ia)) Then
        ias = idxas (ia, is)

! for each angular momentum
        do l = 0, lmax
! sum over k-points
          t4 = 0.d0
          t5 = 0.d0
          do ik = 1, nkpt
            t1 = 0.d0
            t2 = 0.d0
! sum over states
            do ist = 1, nstsv
              if (abs(occsv(ist,ik)) .gt. input%groundstate%epsocc) Then
                ! get rid of nspinor index
                t3 = sum(lcharge(l,:,ias,ist,ik))
                t1 = t1+evalsv(ist,ik)*t3*occsv(ist,ik)
                t2 = t2+t3*occsv(ist,ik)
              end if
            end do ! ist
            t4 = t4 + t1*wkpt(ik)
            t5 = t5 + t2*wkpt(ik)
          end do ! ik
          el(l,ias) = t4/t5
! check if previous radial functions have same default energies
          do io1 = 1, apword (l, is)
            if (apwve(io1, l, is)) then
! for too low and too high values of el we keep the default values
              if ( (el(l,ias) > -0.5d0) .and. &
             &     (el(l,ias) <  0.5d0)) then
                apwe(io1,l,ias) = el(l,ias)
              else
                apwe(io1,l,ias) = apwe0(io1,l,is)
              end if
            end if
          end do ! io1
       end do ! l

!---------------------------------!
!     local-orbital functions     !
!---------------------------------!
        do ilo = 1, nlorb (is)
          do io1 = 1, lorbord (ilo, is)
            if (lorbve(io1, ilo, is)) then
              l = lorbl (ilo, is)

! if lo radial functions have same default energies as apw
              do io2 = 1, apword (l, is)

                if (apwve(io2, l, is)) then

                  if ( abs(lorbe0(io1,ilo,is)-apwe0(io2,l,is) ) < 1.d-4) then
                    lorbe(io1,ilo,ias) = apwe(io2,l,ias)
                  else
                    lorbe(io1,ilo,ias) = el(l,ias)
                  end if

                end if

              end do ! io2

            end if
          end do ! io1
        end do ! ilo
        done (ia) = .True.

!       copy to equivalent atoms
        do ja = 1, natoms (is)
           if (( .not. done(ja)) .and. (eqatoms(ia, ja, is))) then
              jas = idxas (ja, is)
              do l = 0, input%groundstate%lmaxapw
                 do io1 = 1, apword (l, is)
                    apwe (io1, l, jas) = apwe (io1, l, ias)
                 end do
              end do
              do ilo = 1, nlorb (is)
                 do io1 = 1, lorbord (ilo, is)
                    lorbe (io1, ilo, jas) = lorbe (io1, ilo, ias)
                 end do
              end do
              done (ja) = .true.
           end if
        end do

      end if ! done
    end do ! ia
  end do ! is

  deallocate(lcharge,el)
  return
end subroutine
!EOC

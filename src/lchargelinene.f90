!BOP
! !ROUTINE: lchargelinene
! !INTERFACE:
subroutine lchargelinene
! !USES:
      use modinput
      use modmain
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
  integer :: ispn, jspn
  integer :: ilo, io1, io2, ja, jas
  real(8) :: el,t1,t2,t3,t4,t5

  Real (8), Allocatable :: lcharge (:, :, :, :, :)
  Real (8), Allocatable :: elm (:, :)
  Complex (8), Allocatable :: ulm (:, :, :)
  Complex (8), Allocatable :: a (:, :)
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
                  t1 = t1 + dble (dmat(lm, lm, ispn, ispn, ist))
                End Do
                lcharge (l, ispn, ias, ist, ik) = t1
              End Do
            End Do ! ispn
          End Do ! ist
          done (ia) = .True.

! copy to equivalent atoms
          Do ja = 1, natoms (is)
            If (( .Not. done(ja)) .And. (eqatoms(ia, ja, is))) Then
              jas = idxas (ja, is)
              lcharge (:, ispn, jas, ist, ik) = lcharge (:, ispn, ias, ist, ik)
              done (ja) = .True.
            End If
          End Do
          
        End If ! done
      End Do ! ia
    End Do ! is
  End Do ! ik
  Deallocate (dmat, apwalm)
  Deallocate (evecfv, evecsv)

  open(33,file='LCHARGELINENE.OUT',action='write')
  ! for each atom
  do is = 1, nspecies
    done (:) = .False.
    do ia = 1, natoms (is)
      if ( .Not. done(ia)) Then
        ias = idxas (ia, is)
       
        write(33,*)
        write(33,*)'#############################'
        write(33,*)'# atom ', ias
        write(33,*)'#############################'
        write(33,*)

!-----------------------!
!     APW functions     !
!-----------------------!

!       for each angular momentum
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
                t3 = 0.d0
                do ispn = 1, nspinor
                  t3 = t3+lcharge(l,ispn,ias,ist,ik)
                end do
                if (evalsv(ist,ik)>-1.d0) then
                  t1 = t1+evalsv(ist,ik)*t3*occsv(ist,ik)
                  t2 = t2+t3*occsv(ist,ik)
                end if
              end if
            end do ! ist
            t4 = t4 + t1*wkpt(ik)
            t5 = t5 + t2*wkpt(ik)
          end do ! ik

          el = t4/t5
          elcharge(l,ias) = el
          write(33,'("    l=",i4,"    partial charge=",f12.6,"    energy",f12.6)') l, t5, el

! check if previous radial functions have same default energies
          do io1 = 1, apword (l, is)
            if (apwve(io1, l, is)) then
              apwe (io1, l, ias) = el
            end if
          end do ! io1

        end do ! l
        write(33,*)

!---------------------------------!
!     local-orbital functions     !
!---------------------------------!
        do ilo = 1, nlorb (is)
          do io1 = 1, lorbord (ilo, is)
            if (lorbve(io1, ilo, is)) then
              l = lorbl (ilo, is)
! check if previous radial functions have same default energies
              do io2 = 1, apword (l, is)
                if (apwve(io2, l, is)) then
                  if (abs(lorbe0(io1, ilo, is)-apwe0(io2, l, is)) .lt. 1.d-4) then
                    lorbe (io1, ilo, ias) = apwe (io2, l, ias)
                    exit
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
  close(33)

  deallocate(lcharge)
  return
end subroutine
!EOC

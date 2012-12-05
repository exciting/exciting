!BOP
! !ROUTINE: optlinene
! !INTERFACE:
subroutine optlinene
! !USES:
      use modinput
      use modmain
! !DESCRIPTION:
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
  real(8) :: t1,t2,t3,t4,t5

  Real (8), Allocatable :: lcharge (:, :, :, :, :)
  Real (8), Allocatable :: elk (:, :, :)
  Real (8), Allocatable :: elm (:, :)
  Real (8), Allocatable :: temp (:)
  Complex (8), Allocatable :: ulm (:, :, :)
  Complex (8), Allocatable :: a (:, :)
  Complex (8), Allocatable :: evecfv(:,:,:)
  Complex (8), Allocatable :: evecsv(:,:)
  Complex (8), Allocatable :: dmat (:, :, :, :, :)
  Complex (8), Allocatable :: apwalm (:, :, :, :, :)
  
  Logical :: lmirep=.false.
    
  ! set lmax for maximum accuracy
  lmax=input%groundstate%lmaxapw
  lmmax = (lmax+1) ** 2
  
  Allocate (evecfv(nmatmax, nstfv, nspnfv))
  Allocate (evecsv(nstsv, nstsv))
  Allocate (dmat(lmmax, lmmax, nspinor, nspinor, nstsv))
  Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
  Allocate (lcharge(0:lmax, nspinor, natmtot, nstsv, nkpt))

! representation basis of the symmetry group at each atomic site
  If (lmirep) Then
    Allocate (elm(lmmax, natmtot))
    Allocate (ulm(lmmax, lmmax, natmtot))
    Allocate (a(lmmax, lmmax))
    Call genlmirep (lmax, lmmax, elm, ulm)
  End If
  
! loop over k-points
  Do ik = 1, nkpt
    
    Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
    Call getevecsv (vkl(:, ik), evecsv)
  
    ! find the matching coefficients
    Do ispn = 1, nspnfv
      Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :,ispn, ik), &
     & sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, ispn))
    End Do ! ispn

    Do is = 1, nspecies
      Do ia = 1, natoms (is)
        ias = idxas (ia, is)

! generate the density matrix
        Call gendmat (.True., .True., 0, lmax, is, ia, ngk(:, ik), &
       &  apwalm, evecfv, evecsv, lmmax, dmat)

! convert (l,m) part to an irreducible representation if required
        If (lmirep) Then
           Do ist = 1, nstsv
              Do ispn = 1, nspinor
                 Do jspn = 1, nspinor
                    Call zgemm ('N', 'N', lmmax, lmmax, lmmax, &
                   & zone, ulm(:, :, ias), lmmax, dmat(:, :, &
                   & ispn, jspn, ist), lmmax, zzero, a, lmmax)
                    Call zgemm ('N', 'C', lmmax, lmmax, lmmax, &
                   & zone, a, lmmax, ulm(:, :, ias), lmmax, &
                   & zzero, dmat(:, :, ispn, jspn, ist), lmmax)
                 End Do
              End Do
           End Do
        End If

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

      End Do ! ia
    End Do ! is
  End Do ! ik

  Deallocate (dmat, apwalm)
  Deallocate (evecfv, evecsv)
  If (lmirep) Deallocate (elm, ulm, a)

!------------------------------------------------------------
! Calculate the "Optimal Energy Parameters" by S.Bluegel
!------------------------------------------------------------

  open(33,file='OPTLINENE.OUT',action='write')
!
! K-dependent linearization energies
!  
  allocate(elk(0:lmax,natmtot,nkpt))
  elk(:,:,:)=0.d0

  ! for each atom
  Do is = 1, nspecies
    Do ia = 1, natoms (is)
      ias = idxas (ia, is)
      write(33,*)'#############################'
      write(33,*)'# atom ', ias
      write(33,*)'#############################'
      ! for each angular momentum
      do l = 0, lmax
        write(33,*)'l = ', l
        t4 = 0.d0
        t5 = 0.d0
        ! sum over k-points
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
              t1 = t1+evalsv(ist,ik)*t3*occsv(ist,ik)
              t2 = t2+t3*occsv(ist,ik)
            end if
          end do ! ist
          elk(l,ias,ik) = t1/t2
          write(33,'("    ik=",i4,f12.6)') ik, elk(l,ias,ik)
          t4 = t4 + t1*wkpt(ik)
          t5 = t5 + t2*wkpt(ik)
        end do ! ik
        write(33,'("Integrated value:")')
        write(33,'("l=",i4,"    charge=",f12.6,"    energy",f12.6)') l, t5, t4/t5
        write(33,*)'-----------------'
        write(33,*)
      end do ! l
      write(33,*)
    end do ! ia
  end do ! is
  deallocate(elk)
  close(33)

!--------------------------------
!
!--------------------------------
  open(34,file='OPTLINENE_2.OUT',action='write')
  allocate(temp(0:lmax))
  temp(:)=0.d0
  Do is = 1, nspecies
    Do ia = 1, natoms (is)
      ias = idxas (ia, is)
      write(34,*)'# atom ', ias
      ! for each state
      do ist = 1, nstsv
        ! for each angular momentum
        do l = 0, lmax
          t1=0.d0
          ! sum over k-points
          do ik = 1, nkpt
            t2 = 0.d0
            do ispn = 1, nspinor
              t2 = t2+lcharge(l,ispn,ias,ist,ik)
            end do
            t1 = t1+t2*wkpt(ik)
          end do ! ik
          temp(l)=t1
        end do ! l
        t1 = 0.d0
        do ik = 1, nkpt    
          t1 = t1+evalsv(ist,ik)*wkpt(ik)
        end do ! ik
        write(34,*) t1, temp(0:lmax)
      end do ! ist
      write(34,*)
    end do ! ia
  end do ! is
  deallocate(temp)
  close(34)

  deallocate(lcharge)
  return
end subroutine
!EOC

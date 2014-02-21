!BOP
! !ROUTINE: lchargelinene
! !INTERFACE:
subroutine lchargelinene
! !USES:
      use modinput
      use modmain
      use modmpi, only: rank
! !DESCRIPTION:
!   Calculate the "Optimal Energy Parameters" by S.Bluegel
! !REVISION HISTORY:
!   Created October 2013 (DIN)
!EOP
!BOC
  implicit none
  ! local variables
  integer :: ik,l,m,lm
  integer :: lmax,lmmax
  integer :: is,ia,ias,ist,n,i
  integer :: ispn
  integer :: ilo, io1, io2, ja, jas
  real(8) :: t1,t2
  real(8) :: de, de_min
  logical :: newset

  integer, allocatable  :: nl(:)
  Real (8), Allocatable :: el(:,:,:)
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
  
  Allocate (evecfv(nmatmax,nstfv,nspnfv))
  Allocate (evecsv(nstsv,nstsv))
  Allocate (dmat(lmmax,lmmax,nspinor,nspinor,nstsv))
  Allocate (apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
  Allocate (lcharge(0:lmax,nspinor,natmtot,nstsv,nkpt))
  Allocate (el(nstsv,0:input%groundstate%lmaxapw,natmtot))
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
        write(32,*) '# state    l=0         l=1         l=2         l=3         l=4'
        do ist = 1, nstsv
          do l = 0, input%groundstate%lmaxapw
            t1 = 0.d0
            do ik = 1, nkpt
              t1 = t1 + wkpt(ik)*sum(lcharge(l,:,ias,ist,ik))
            end do
            lstate(l,ist,ias) = t1
          end do ! l
          write(32,'(i8,5f12.6)') ist, (lstate(l,ist,ias), l = 0, 4)
        end do ! ist
       write(32,*)
     end do ! ia
    end do ! is
    close(32)
  end if

  allocate(nl(0:lmax))
  nl(:) = 0

! for each atom
  do is = 1, nspecies
    done (:) = .False.
    do ia = 1, natoms (is)
      if ( .Not. done(ia)) Then
        ias = idxas (ia, is)

!-----------------------
! Calculate l-energies
!-----------------------
        do l = 0, lmax
! sum over states
          t1 = 0.d0
          t2 = 0.d0
          newset = .false.
          do ist = 1, nstsv
            if (lstate(l,ist,ias) > 0.1d0) then
              newset = .true.
! integrate over BZ
              do ik = 1, nkpt
                t1 = t1+evalsv(ist,ik)*sum(lcharge(l,:,ias,ist,ik))*wkpt(ik)
              end do
              t2 = t2+lstate(l,ist,ias)
            else
              if (newset) then
                nl(l) = nl(l)+1
                el(nl(l),l,ias) = t1/t2
                t1 = 0.d0
                t2 = 0.d0
                newset = .false.
              end if
            end if
          end do ! ist
        end do ! l

!---------------------------------!
!      augmented functions         !
!---------------------------------!
        do l = 0, lmax 
          do io1 = 1, apword (l, is)
            apwe(io1,l,ias) = apwe0(io1,l,is)
            if (apwve(io1,l,is)) then
                n = 1
                de_min = dabs(el(1,l,ias)-apwe0(io1,l,is))
                do i = 2, nl(l)
                  de = dabs(el(i,l,ias)-apwe0(io1,l,is))
                  if (de < de_min) then
                    n = i
                    de_min = de 
                  end if
                end do
                apwe(io1,l,ias) = el(n,l,ias)
            end if
          end do ! io1
        end do ! l

!---------------------------------!
!      local-orbital functions    !
!---------------------------------!
       do ilo = 1, nlorb (is)
         do io1 = 1, lorbord (ilo, is)
           lorbe(io1,ilo,ias) = lorbe0(io1,ilo,is)
           if (lorbve(io1, ilo, is)) then
             l = lorbl (ilo, is)
! if lo has the same default energies as apw
             do io2 = 1, apword (l, is)
               if (abs(lorbe0(io1,ilo,is)-apwe0(io2,l,is) )<1.d-4) then
                 if (apwve(io2,l,is)) then
                   lorbe(io1,ilo,ias) = apwe(io2,l,ias)
                 else
                   lorbe(io1,ilo,ias) = lorbe0(io1,ilo,is)
                 end if
! energies are different
               else
                 n = 1
                 de_min = dabs(el(1,l,ias)-lorbe0(io1,ilo,is))
                 do i = 2, nl(l)
                   de = dabs(el(i,l,ias)-lorbe0(io1,ilo,is))
                   if (de < de_min) then
                     n = i
                     de_min = de 
                   end if
                 end do
                 lorbe(io1,ilo,ias) = el(n,l,ias)
               end if
             end do ! io2
           end if
         end do ! io1
       end do ! ilo
       done (ia) = .True.

!---------------------------------!
!      copy to equivalent atoms  !
!---------------------------------!
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
  deallocate(lstate,lcharge,nl,el)

  return
end subroutine
!EOC

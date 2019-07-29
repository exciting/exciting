!
!
!
!
! Copyright (C) 2013 exciting team
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: hmlint
! !INTERFACE:
!
!
Subroutine hmlint(mt_h)
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Calculates the "muffin-tin" Hamiltonian.
!
!EOP
!BOC
      Implicit None
      Type (MTHamiltonianList), intent(InOut) :: mt_h
!      Type (MTHamiltonianList) :: mt_h
! local variables
      Integer :: is, ia, ias, nr, ir, if1,if3,inonz,ireset1,ireset3
      Integer :: l1, l2, l3, m2, lm2, m1, m3, lm1, lm3
      Integer :: ilo, ilo1, ilo2, io, io1, io2, nalo1, maxnlo
      Real (8) :: t1,t2,angular
      Real (8), allocatable :: haaintegrals(:,:,:,:,:),halointegrals(:,:,:,:),hlolointegrals(:,:,:)
      Real (8) :: rmtable(nrmtmax),r2inv(nrmtmax)
      complex(8) :: zsum
! automatic arrays
      Real (8) :: r2 (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, nrmtmax),a,rm
      integer, allocatable :: lfromlm(:),mfromlm(:)
! fine structure constant
      Real (8), Parameter :: alpha = 1.d0 / 137.03599911d0
! electron g factor
      Real (8), Parameter :: ge = 2.0023193043718d0  
      Real (8), Parameter :: ga4 = ge * alpha / 4.d0      
      
      Type (apw_lo_basis_type) :: mt_basis
 

      mt_basis%lofr=>lofr
      mt_basis%apwfr=>apwfr
         

        call MTRedirect(mt_h%main,mt_h%spinless)        
        call mt_kin(veffmt,mt_basis,mt_h)
        call mt_pot(veffmt,mt_basis,mt_h)

! now the magnetic field
        if (associated(input%groundstate%spin)) then

          call MTInit(mt_h%alpha,mt_h%maxaa,mt_h%maxnlo)        
          call MTInit(mt_h%beta,mt_h%maxaa,mt_h%maxnlo)

! is the problem noncollinear?
          if (ncmag) then 
            call MTInit(mt_h%ab,mt_h%maxaa,mt_h%maxnlo)            
            call MTRedirect(mt_h%main,mt_h%alpha)
  
! adding the z component of the external magnetic field to bxcmt
            do is=1,nspecies
              do ia=1,natoms(is)
                ias=idxas(ia,is)
                bxcmt (1, :, ias, 3) = bxcmt (1, :, ias, 3) + ga4 * (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(3)+input%groundstate%spin%bfieldc(3))/y00
              enddo
            enddo
  
            call mt_pot(bxcmt(:,:,:,3),mt_basis,mt_h)
  
! removing the z component of the external magnetic field to bxcmt / restoring bxcmt to the initial state
            do is=1,nspecies
              do ia=1,natoms(is)
                ias=idxas(ia,is)
                bxcmt (1, :, ias, 3) = bxcmt (1, :, ias, 3) - ga4 * (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(3)+input%groundstate%spin%bfieldc(3))/y00
              enddo
            enddo
  
            mt_h%beta%aa=-mt_h%alpha%aa
! do we have local orbitals?
            if (associated(mt_h%beta%lolo)) then
              mt_h%beta%alo=-mt_h%alpha%alo
              mt_h%beta%loa=-mt_h%alpha%loa
              mt_h%beta%lolo=-mt_h%alpha%lolo
            endif
  
  
! adding the y component of the external magnetic field to bxcmt
            call MTRedirect(mt_h%main,mt_h%ab)
            do is=1,nspecies
              do ia=1,natoms(is)
                ias=idxas(ia,is)
                bxcmt (1, :, ias, 2) = bxcmt (1, :, ias, 2) + ga4 * (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(2)+input%groundstate%spin%bfieldc(2))/y00
              enddo
            enddo
  
            call mt_pot(bxcmt(:,:,:,2),mt_basis,mt_h)
  
! removing the y component of the external magnetic field to bxcmt / restoring bxcmt to the initial state
            do is=1,nspecies
              do ia=1,natoms(is)
                ias=idxas(ia,is)
                bxcmt (1, :, ias, 2) = bxcmt (1, :, ias, 2) - ga4 * (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(2)+input%groundstate%spin%bfieldc(2))/y00
              enddo
            enddo
  
! scale the alpha-beta block of the Hamiltonian by i, because we are handling the y component of the magnetic field now
            mt_h%ab%aa=-zi*mt_h%ab%aa
! do we have local orbitals?
            if (associated(mt_h%beta%lolo)) then
              mt_h%ab%alo=-zi*mt_h%ab%alo
              mt_h%ab%loa=-zi*mt_h%ab%loa
              mt_h%ab%lolo=-zi*mt_h%ab%lolo
            endif
  
  
! adding the x component of the external magnetic field to bxcmt
            call MTRedirect(mt_h%main,mt_h%ab)
            do is=1,nspecies
              do ia=1,natoms(is)
                ias=idxas(ia,is)
                bxcmt (1, :, ias, 1) = bxcmt (1, :, ias, 1) + ga4 * (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(1)+input%groundstate%spin%bfieldc(1))/y00
              enddo
            enddo
  
            call mt_pot(bxcmt(:,:,:,1),mt_basis,mt_h)
  
! removing the x component of the external magnetic field to bxcmt / restoring bxcmt to the initial state
            do is=1,nspecies
              do ia=1,natoms(is)
                ias=idxas(ia,is)
                bxcmt (1, :, ias, 1) = bxcmt (1, :, ias, 1) - ga4 * (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(1)+input%groundstate%spin%bfieldc(1))/y00
              enddo
            enddo
  
! the collinear case
          else
            call MTRedirect(mt_h%main,mt_h%alpha)
            do is=1,nspecies
              do ia=1,natoms(is)
                ias=idxas(ia,is)
                bxcmt (1, :, ias, 1) = bxcmt (1, :, ias, 1) + ga4 * (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(3)+input%groundstate%spin%bfieldc(3))/y00
              enddo
            enddo
  
            call mt_pot(bxcmt(:,:,:,1),mt_basis,mt_h)
  
            do is=1,nspecies
              do ia=1,natoms(is)
                ias=idxas(ia,is)
                bxcmt (1, :, ias, 1) = bxcmt (1, :, ias, 1) - ga4 * (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(3)+input%groundstate%spin%bfieldc(3))/y00
              enddo
            enddo
  
            mt_h%beta%aa=-mt_h%alpha%aa
! do we have local orbitals?
            if (associated(mt_h%beta%lolo)) then
              mt_h%beta%alo=-mt_h%alpha%alo
              mt_h%beta%loa=-mt_h%alpha%loa
              mt_h%beta%lolo=-mt_h%alpha%lolo
            endif
          endif
  
! add the spin orbit interaction if requested
          if (isspinorb()) call mt_so(veffmt,mt_basis,mt_basis,mt_h,level_zora)
        endif

!deallocate(mt_h%losize)

      Return
End Subroutine
!EOC

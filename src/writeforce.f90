!
!-----------------------------------------------------------------------------80.
! Copyright (C) 2013- exciting team
! This file is distributed under the terms of the GNU General Public License.
! Created       on 15-10-2013 Pasquale Pavone (exciting team)
! Modified      on 15-11-2013 Pasquale Pavone (exciting team)
! Last modified on 14-06-2014 Pasquale Pavone (exciting team)
!-----------------------------------------------------------------------------80
!
subroutine writeforce(fnum,verbosity)
    use modmain
    use modinput
    implicit none
! arguments
    integer, intent(in) :: fnum
    integer, intent(in) :: verbosity
! local variables
    integer :: is, ia, ias, i, j
    logical :: totlock
    character*3 :: clabel(3)
    
    if (verbosity < 1) return

    totlock = .False.
    if ( associated(input%relax) ) then 
        do is = 1, nspecies
            do ia = 1, natoms (is)
                do i = 1, 3
                    if (input%structure%speciesarray(is)%species%atomarray(ia)%atom%lockxyz(i)) then
                        totlock = .True.
                        go to 10
                    end if
                end do
            end do
        end do
    end if

10  i=0
    write(fnum,*)
    if (input%groundstate%tfibs) then
        if (totlock) then
           write(fnum,'(" Total atomic forces including IBS (cartesian) + constraints (cartesian) :")')
        else
           write(fnum,'(" Total atomic forces including IBS (cartesian) :")')
        end if
    else
        if (totlock) then
           write(fnum,'(" Total atomic forces (cartesian) + constraints (cartesian) :")')
        else
           write(fnum,'(" Total atomic forces (cartesian) :")')
        end if
    end if 
    do is = 1, nspecies
        do ia = 1, natoms (is)
            ias = idxas (ia, is)
            i=i+1
            if (totlock) then
                do j = 1, 3
                    clabel(j) = "F"
                    if (input%structure%speciesarray(is)%species%atomarray(ia)%atom%lockxyz(j)) clabel(j) = "T"
                end do
                write(fnum,'(" atom ",I5,2x,A2,T18,": ",3F14.8,3x,3A3)') &
               &  i, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
               &  forcetot(:,ias), clabel(:)
            else
                write(fnum,'(" atom ",I5,2x,A2,T18,": ",3F14.8)') &
               &  i, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
               &  forcetot(:,ias)
            end if
        end do
    end do

   i=0
   if (verbosity > 1) then
        write(fnum,*)
        if (input%groundstate%tfibs) then
            write(fnum,'(" Atomic force components including IBS (cartesian) :")')
        else
            write(fnum,'(" Atomic force components (cartesian) :")')
        end if
        do is = 1, nspecies
            do ia = 1, natoms (is)
                ias = idxas (ia, is)
                i=i+1
                write(fnum,'(" atom ",I5,2x,A2,T18,": ",3F14.8,"   HF force")') &
               &  i, trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
               &  forcehf(:,ias)
                write(fnum,'(                  T18,": ",3F14.8,"   core correction")') &
               &  forcecr(:,ias)
                if (input%groundstate%tfibs) &
               &write(fnum,'(                  T18,": ",3F14.8,"   IBS correction")') &
               &  forceibs(:,ias)
                If ( input%groundstate%vdWcorrection .Ne. "none" ) Then
                   If ( input%groundstate%vdWcorrection .Eq. "DFTD2" ) Then
                      write(fnum,'(                  T18,": ",3F14.8,"   DFT-D2 force")') &
                           &  force_disp(:,ias)
                   Else If ( input%groundstate%vdWcorrection .Eq. "TSvdW" ) Then
                      write(fnum,'(                  T18,": ",3F14.8,"   TS-vdW force")') &
                           &  force_disp(:,ias)
                   End If
                End If

             end do
        end do
    end if
    
    return
end subroutine

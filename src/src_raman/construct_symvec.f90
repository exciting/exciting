module m_symvec 
!              
integer, allocatable :: elements(:,:)                                                          
integer :: no_perm
!              

contains

subroutine construct_symvec(ev)
!
use raman_symmetry
use mod_symmetry
use mod_atoms
!
implicit none
! arguments
Complex (8), intent(out) :: ev(3*natmtot, 3*natmtot)
! local variables
real(8), allocatable :: vec_vib(:), symvec_phon(:, :)
real(8) :: norm_factor
integer :: i, k, j, ivc, igr, ia, ja, ip, is, iat, m, n, ia_gr, iev
integer :: no_vec_orth
integer :: icomp(3*natmtot)
logical :: new, acousticx, acousticy, acousticz
logical :: acoustic(3*natmtot)
!
ev (:, :) = zzero
!
! Projection into phonon symmetry subspaces
! +/-(1,0,0)T, +/-(0,0,1)T and +/-(0,0,1)T vectors on one atom of every sym group are used as basis
!
if (.not. sym_out) then
   write(*,*)
   write(*,*) ' The construction of symmetry vectors seems to be not meaningful'
   write(*,*) ' Do you really want to do this?'
   write(*,*)
endif
!
! write out symmetry groups of atoms
write(*,*)
write(*,'(" There are ",i2," symmetry groups of atoms :")') gr
do i = 1,gr
   write(*,'(" Group ",i2," contains atoms :",20i4)') i,(gr_atoms(i,j),j=1,gr_atoms_no(i))
enddo
!
! compute all possibilities for combining gr times +1 and -1, excluding linear dependent combinations
allocate( elements(2**(gr-1), gr) )
elements = 1
no_perm = 1
call permute(gr,1)
write(*,'(/," There are ",i3," possibilities to distribute +/- vectors over the ",&
&         i3," symmetry groups of atoms...",/)') no_perm, gr
!do i = 1,no_perm
!   write(*,*) 'perm ',i
!   write(*,*) elements(i,:)
!enddo
!
allocate( vec_vib(3*natmtot) )
allocate( symvec_phon(3*natmtot, 3*no_vec_nonorth) )
iev = 0
!
do k = 1, cl
  if (vib_mode(k)) then
    n = 0
    ! loop over components of an R3 vector
    do ivc = 1, 3
       ! loop over gr elements of all possibilities for combining +1's and -1's
       do ip = 1, no_perm
          ! loop over inequivalent groups of atoms
          do ia = 1, natmtot
             vec_vib = 0.d0
             ! sum over all contributions from different groups of atoms
             do igr = 1, gr
                ia_gr = mod(ia,gr_atoms_no(igr)) + 1
                ! loop over SOPs
                do j = 1, nsymcrys
                  ja = atom_sop(gr_atoms(igr,ia_gr), j)
                ! this is the projection operator (involving only characters of IREPs), 
                ! working on a vector times +1 or -1 (all possibilities over groups of atoms)
                  vec_vib(3*ja-2:3*ja) = vec_vib(3*ja-2:3*ja) &
    &              + dble(charact(class(j),k))*sopmat(1:3, ivc, j) &
    &                *dble(elements(ip,igr))
                enddo
             enddo
             ! normalize vector
             norm_factor = 0.d0
             do i = 1, 3*natmtot
               norm_factor = norm_factor + vec_vib(i)**2
             enddo
             if (norm_factor .gt. eps) norm_factor = 1.d0 / sqrt(norm_factor)
             vec_vib = vec_vib*norm_factor
             ! check if we got something new ---------------------------------------------------
             if (any(abs(vec_vib) .gt. eps)) then                                              !
               new = .true.                                                                    !
               do m = 1, n                                                                     !
                 if (all(abs(vec_vib(:) - symvec_phon(:, m)) .lt. eps) .or.  &                 !
    &                all(abs(vec_vib(:) + symvec_phon(:, m)) .lt. eps))      &                 !
    &              new = .false.                                                               !
               enddo                                                                           !
               if (new) then                                                                   !
                  n = n + 1                                                                    !
                  symvec_phon(:, n) = vec_vib(:)                                               !
               endif                                                                           !
             endif                                                                             !
             ! ---------------------------------------------------------------------------------
          enddo ! end loop over atoms
       enddo ! end loop over +/-1 combinations
    enddo ! end loop over vector components
!
! orthonormalize the collection of vectors
!   call dgeqrf(3*natmtot,n,symvec_phon,3*natmtot,tau,work,lwork,info)
    call schmidt(symvec_phon(:,1:n), 3*natmtot, n, no_vec_orth)
!
    write(*,*)
    write(*,'(" Symmetry vectors of IREP ",i2," (",a2,")")') k, irep_ch(k)
    write(*,*) 'n = ', no_vec_orth
    do m = 1, no_vec_orth
      write(*,*) ' Vector ', m
      do ia = 1, natmtot
        write(*,'(" Atom ",i2," : ",3f10.4)') ia, symvec_phon(3*ia-2:3*ia, m)
      enddo
    enddo
!
    if (no_vec_orth .gt. vib_ireps(k)*nint(dble(charact(1, k)))) then
       write(*,'(/," The construction of symmetry vectors for this system failed.")')
       write(*,'(  " IREP : ",i1,"   (",a3,")"/," no. phonon modes : ",i3,/," no. symvecs found : ",i3,/)') &
          &              k, irep_ch(k), vib_ireps(k), no_vec_orth
       write(*,'(" Please use the other options of input%properties%raman%getphonon",/)')
       stop
    endif
!
! check for acoustic modes
    do m = 1, no_vec_orth
!      acoustic(iev+m) = .false.
!      do i = 1, 3*natmtot
!         if (abs(symvec_phon(i, m)) .le. eps) then
!            icomp(i) = 0
!         elseif (symvec_phon(i, m) .gt. eps) then
!            icomp(i) = 1
!         else
!            icomp(i) = -1
!         endif
!      enddo
!      acousticx = .false.; acousticy = .false.; acousticz = .false.
!      if (all(icomp(1:3*natmtot-2:3) .eq.  0) .or. &
!       &  all(icomp(1:3*natmtot-2:3) .eq.  1) .or. &
!       &  all(icomp(1:3*natmtot-2:3) .eq. -1) ) acousticx = .true.
!      if (all(icomp(2:3*natmtot-1:3) .eq.  0) .or. &
!       &  all(icomp(2:3*natmtot-1:3) .eq.  1) .or. &
!       &  all(icomp(2:3*natmtot-1:3) .eq. -1) ) acousticy = .true.
!      if (all(icomp(3:3*natmtot:3) .eq.  0) .or. &
!       &  all(icomp(3:3*natmtot:3) .eq.  1) .or. &
!       &  all(icomp(3:3*natmtot:3) .eq. -1) ) acousticz = .true.
!      if (acousticx .and. acousticy .and. acousticz) acoustic(iev+m) = .true.
       call check_acoustic(symvec_phon(:, m), acoustic(iev+m))
    enddo
!
! weight by sqrt of masses and renormalize
    do m = 1, no_vec_orth
       iat = 0
       norm_factor = 0.d0
       do is = 1, nspecies
          do ia = 1, natoms(is)
             iat = iat + 1
             symvec_phon(3*iat-2:3*iat, m) = symvec_phon(3*iat-2:3*iat, m) / sqrt(spmass(is))
             norm_factor = norm_factor + symvec_phon(3*iat-2, m)**2 + symvec_phon(3*iat-1, m)**2 + &
               &                         symvec_phon(3*iat, m)**2
          enddo
       enddo
       if (norm_factor .gt. 1.d-12) norm_factor = 1.d0 / sqrt(norm_factor)
       symvec_phon(:, m) = symvec_phon(:, m) * norm_factor
    enddo
! save on array ev(:, :)
    do m = 1, no_vec_orth
       ev(:, iev+m) = cmplx(symvec_phon(:, m), 0., 8)
    enddo
    iev = iev + no_vec_orth
!
  endif
enddo
!
! rearrange ev(:, :) in order to have the acoustic modes in ev(:, 1:3)
! for systems with more than 1 atom
if (natmtot .gt. 1) then
   do i = 1, 3*natmtot
      if (acoustic(i) .and. i .gt. 3) then
         do j = 1, 3
            if (.not. acoustic(j)) then
               ! swap i and j
               symvec_phon(:, 1) = dble(ev(:, i))
               ev(:, i) = ev(:, j)
               ev(:, j) = cmplx(symvec_phon(:, 1), 0.d0, 8)
               acoustic(i) = .false.
               acoustic(j) = .true.
            endif
         enddo
      endif
   enddo
endif
!
deallocate( symvec_phon, vec_vib )
!
!
end subroutine construct_symvec

recursive subroutine permute(no_el,i_el)
implicit none
integer, intent(in) :: no_el
integer, intent(in) :: i_el
integer :: ii,el_temp(no_el)
!
el_temp(:) = elements(no_perm,:)
do ii = i_el+1,no_el
   if (ii .le. no_el) then
      no_perm = no_perm + 1
      elements(no_perm,:) = el_temp(:)
      elements(no_perm,ii) = -1
      call permute(no_el,ii)
   endif
enddo
return
end subroutine permute

end module m_symvec

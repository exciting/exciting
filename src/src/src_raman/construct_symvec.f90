! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created July 2014 by Stefan Kontur
!
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
real(8), allocatable :: vec_vib(:), symvec_phon(:, :), work(:), tau(:)
real(8) :: norm_factor
complex(8) :: ztmp(3*natmtot)
integer :: i, k, j, ivc, igr, ia, ja, ip, is, iat, m, n, ia_gr, iev
integer :: no_vec_orth, info
integer :: icomp(3*natmtot)
integer, allocatable :: indep(:)
logical :: new, acousticx, acousticy, acousticz
logical :: acoustic(3*natmtot)
!
ev (:, :) = zzero
!
! open symmetry output file
open(unit=13,file='RAMAN_SYM.OUT',status='unknown',form='formatted',position='append')
!
! Projection into phonon symmetry subspaces
! +/-(1,0,0)T, +/-(0,0,1)T and +/-(0,0,1)T vectors on one atom of every sym group are used as basis
!
write(13, '(/," Construction of symmetry vectors ",/, &
       &      " -------------------------------- ")')
if (.not. sym_out) then
   write(13,'(/," The construction of symmetry vectors might be not meaningful due to the low symmetry.",/, &
   &            " Check any result below carefully!",/)')
endif
!
! write out symmetry groups of atoms
write(13,'(/," There are ",i2," symmetry groups of atoms :")') gr
do i = 1,gr
   write(13,'(" Group ",i2," contains atoms :",20i4)') i,(gr_atoms(i,j),j=1,gr_atoms_no(i))
enddo
if (gr .gt. 20) then
    write(13,'(/," The construction of symmetry vectors for this system is not feasible.")')
    write(13,'(" Please use the other options of input%properties%raman%getphonon",/)')
    stop
endif
!
! compute all possibilities for combining gr times +1 and -1, excluding linear dependent combinations
allocate( elements(2**(gr-1), gr) )
elements = 1
no_perm = 1
call permute(gr,1)
write(13,'(/," There are ",i3," possibilities to distribute +/- unit vectors over the ",&
&         i3," symmetry groups of atoms.",/)') no_perm, gr
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
                do j = 1, numsop
                  ja = atom_sop(gr_atoms(igr,ia_gr), j)
                ! this is the projection operator (involving only characters of IREPs), 
                ! working on a vector times +1 or -1 (all possibilities over groups of atoms)
                  vec_vib(3*ja-2:3*ja) = vec_vib(3*ja-2:3*ja) &
    &              + dble(charact(class(j),k))*sopmatc(1:3, ivc, j) &
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
    &                all(abs(vec_vib(:) + symvec_phon(:, m)) .lt. eps)) then                   !
                    new = .false.                                                              !
                    exit                                                                       !
                 endif                                                                         !
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
! eliminate linearly dependent vectors
    allocate( indep(n) )
    call findindep(3*natmtot, n, symvec_phon(:,1:n), no_vec_orth, indep)
    do m = 1, no_vec_orth
       symvec_phon(:, m) = symvec_phon(:, indep(m))
    enddo
    deallocate( indep )
! orthonormalize the remaining collection
    call schmidt(symvec_phon(:,1:no_vec_orth), 3*natmtot, no_vec_orth)
    
!
    write(13,'(/," Symmetry vectors of IREP ",i2," (",a4,")")') k, irep_ch(k)
    write(13,'("  n = ",i3)') no_vec_orth
    do m = 1, no_vec_orth
      write(13,'("  Vector ",i3)') m
      do ia = 1, natmtot
        write(13,'("  Atom ",i3," : ",3f10.4)') ia, symvec_phon(3*ia-2:3*ia, m)
      enddo
    enddo
!
    if (no_vec_orth .gt. vib_ireps(k)*nint(dble(charact(1, k)))) then
       write(13,'(/," The construction of symmetry vectors for this system failed.")')
       write(13,'(  " IREP : ",i1,"   (",a3,")"/," no. phonon modes : ",i3,/," no. symvecs found : ",i3,/)') &
          &              k, irep_ch(k), vib_ireps(k)*nint(dble(charact(1, k))), no_vec_orth
       write(13,'(" Please use the other options of input%properties%raman%getphonon",/)')
       stop
    endif
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
!
if (natmtot .gt. 1) then
   do i = 4, 3*natmtot
      call check_acoustic(ev(:, i), acoustic(i))
      if (acoustic(i)) then
         do j = 1, 3
            call check_acoustic(ev(:, j), acoustic(j))
            if (.not. acoustic(j)) then
               ! swap i and j
               ztmp = ev(:, i)
               ev(:, i) = ev(:, j)
               ev(:, j) = ztmp
               acoustic(i) = .false.
               acoustic(j) = .true.
            endif
         enddo
      endif
   enddo
endif
!
write(13, '(/," --- Construction of symmetry vectors done. --- ",//)')
close(13)
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


subroutine schmidt(mat,n,m)
! orthogonalizes the m vectors in matrix mat(n,m)
! using a simple Schmidt orthogonalization
! the resulting nvo orthonormal vectors are returned
! again in mat(n,1:nvo)
use raman_params
implicit none
integer, intent(in) :: n,m
real(8), intent(inout) :: mat(n,m)
real(8) :: tempmat(n,m)
integer :: i,j
!
tempmat = 0.d0
do i = 1,m
   tempmat(:,i) = mat(:,i)
   do j = 1,i-1
      tempmat(:,i) = tempmat(:,i) - dot_product(mat(:,i),tempmat(:,j))*tempmat(:,j) / &
                   &                dot_product(tempmat(:,j),tempmat(:,j))
   enddo
   mat(:,i) = fnorm(tempmat(:,i),n)*tempmat(:,i)
enddo
!
end subroutine schmidt


function fnorm(vec,n)
! computes the norm or a vector
implicit none
integer, intent(in) :: n
real(8), intent(in) :: vec(n)
real(8) :: fnorm
integer :: i
!
fnorm = 0.d0
do i = 1,n
   fnorm = fnorm + vec(i)*vec(i)
enddo
!write(*,*) 'norm ',norm
if (fnorm .gt. 1.d-8) fnorm = 1.d0 / sqrt(fnorm)
return
end function fnorm



subroutine findindep(m, n, mat, k, indep)
! compute row echelon form of a general (m,n)-matrix and
! determine its linearly independet columns
! IN: m,n,mat(m,n)
! OUT: number of linearly independent columns k
!      position of lin indep columns in original matrix: first k entries in indep(n)
use raman_params, only : eps
implicit none
integer, intent(in) :: m, n
real(8), intent(in) :: mat(m, n)
integer, intent(out) :: k
integer, intent(out) :: indep(n)
! local variables
real(8) :: tempmat(m, n)
integer, dimension(m) :: ipiv
integer :: i,j,n1
!
tempmat(:,:) = mat(:,:)
indep(:) = 0
!
! form row echelon form
do i = 1,m
   call pivot(m, n, tempmat,ipiv)
   if (ipiv(i) .le. m) then
      if (abs(tempmat(i,ipiv(i))) .gt. eps) then
         tempmat(i,:) = tempmat(i,:) / tempmat(i,ipiv(i))
      endif
   else
     cycle
   endif
   do j = i+1,m
      tempmat(j,:) = tempmat(j,:) - tempmat(j,ipiv(i))*tempmat(i,:)
   enddo
enddo
!
! find leading ones
k = 0
n1 = 1
do i = 1, m
   do j = n1, n
      if (abs(tempmat(i,j)-1.d0) .lt. eps) then
         k = k + 1
         indep(k) = j
         n1 = j
         exit
      endif
   enddo
enddo
!
!
end subroutine findindep



subroutine pivot(m, n, mat, ipiv)
! pivoting of a general (m,n)-matrix
use raman_params, only : eps
implicit none
integer, intent(in) :: m, n
real(8), dimension(m,n), intent(inout) :: mat
integer, dimension(m), intent(out) :: ipiv
real(8), dimension(n) :: temp
integer :: i,j,k
!
ipiv = m+1
do i = 1,m
   do j = 1,n
      if (abs(mat(i,j)) .gt. eps) then
         ipiv(i) = j
         exit
      endif
   enddo
enddo
do i = 1,m
   do j = i+1,m
      if (ipiv(i) .gt. ipiv(j)) then
         temp(:) = mat(j,:)
         mat(j,:) = mat(i,:)
         mat(i,:) = temp(:)
         k = ipiv(j)
         ipiv(j) = ipiv(i)
         ipiv(i) = k
      endif
   enddo
enddo
end subroutine pivot


end module m_symvec

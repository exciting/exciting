! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created July 2014 by Stefan Kontur
!
! -----------------------------------------------------------------------------
subroutine construct_chartabl                                                 !
! -----------------------------------------------------------------------------
!
! References:
! Computation of Character Tables:
! J. McKay in: J.Leech (ed.), Computational problems in abstract algebra,
! Oxford, 1970. (Oxford, 1967), Pergamon Press.
!
!
use modinput
use mod_symmetry
use raman_symmetry
use mod_atoms
use m_eigvec
use raman_params, only : eps
implicit none
real(8) :: rdum1,rdum2,rdum3
integer :: i,j,k,fi,numrep,tr,sop_1,sop_2,sop_3,lwork,info
integer :: i1,i2
integer :: maxrot,maxrot_cl,minrot,minrot_cl,inv_cl,ia,ja,is, num_m, sigmah_cl
real(8) :: transl(3,48)
real(8) :: conj(3,3), tmpmat(3, 3), invmat(3, 3)
real(8) :: E(3,3)
real(8) :: sopmat_12(3,3),const_c(48,48,48),dim_(48)
real(8) :: char_vec(48),char_raman(48),vec_one(3)
real(8) :: char_rot(48),char_equiv(48)
complex(8) :: freq_fac_vec(48),freq_fac_raman(48),freq_fac_rot(48),freq_fac_vib(48)
real(8) :: dirvec(3), atpos_sop(3)
real(8), allocatable :: rwork(:),indet_u(:),atpos(:,:)
complex(8), allocatable :: mat_phi(:,:),eigval(:),work(:)
character(175) :: vec_sym,rot_sym,raman_sym,vib_sym
logical :: flag(48),representative(48),grouped(100)
integer :: rottype(48) 
integer :: rottabl(-3:3,-1:1),rotsum(-6:6)
character(2), dimension(-6:6) :: rotchar = (/ 'S3','  ','S4','S6',' m',' i','  ',' E','C2','C3','C4','  ','C6' /)
character(2) :: cl_rotchar(48)
logical :: centric,cubic,mix
character(5) :: hmsymb,schsymb,chnum
character(3) :: sense
character(1) :: mul_index
! external functions
      Real (8) :: r3mdet
      External r3mdet
! one vector
vec_one(:) = (/ 1.,1.,1. /)
! one matrix
E(:,1) = (/ 1.,0.,0. /)
E(:,2) = (/ 0.,1.,0. /)
E(:,3) = (/ 0.,0.,1. /)
!
! open symmetry output file
open(unit=13,file='RAMAN_SYM.OUT',status='unknown',form='formatted')
!
! lookup table for rotation types
rottabl(:,1) = (/ 0,0,2,3,4,6,1 /)
rottabl(:,-1) = (/ -1,-6,-4,-3,-2,0,0 /)
!
call init0
!
allocate( atpos(natmtot,3) )
i = 1
do is = 1, nspecies
   Do ia = 1, natoms (is)
      atpos(i, 1) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(1)
      atpos(i, 2) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(2)
      atpos(i, 3) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(3)
      i = i + 1
   enddo
enddo
numsop = nsymcrys
!
class = 0
rotsum = 0
!
representative = .true.
fi = 1
do i = 1,numsop
   ! use the index of SOPs referring to crystal symmetry
   sopmat(:, :, i) = dble( symlat(:, :, lsplsymc(i)) )
   sopmatc(:, :, i) = symlatc(:, :, lsplsymc(i))
   transl(:, i) = vtlsymc(:, i)
! check for i
   if (all(sopmat(:,:,i) .eq. -E)) fi = 2
! check for occurrences of (R,r1) and (R,r2)
   do j = 1,i-1
      if (all(sopmat(:,:,i) .eq. sopmat(:,:,j))) then
         representative(i) = .false.
      endif
   enddo
enddo
!
! determine classes of the group
!
cl = 0; flag = .true.
do i = 1,numsop
   if (flag(i)) then
!     element i is first element of a new class
      cl = cl + 1
      class(i) = cl
      flag(i) = .false.
      do j = 1,numsop
         if ((j .ne. i) .and. flag(j)) then
!           testing element j...
            do k = 1,numsop
               ! compute B = XAX^-1
               call r3minv (sopmat(:, :, k), invmat(:, :))
               call r3mm (sopmat(:, :, i), invmat(:, :), tmpmat(:, :))
               call r3mm (sopmat(:, :, k), tmpmat(:, :), conj(:, :))
               ! is SOP j equal to B?
               if (all(conj(:, :) .eq. sopmat(:,:,j))) then
               !  add j to class cl
                  class(j) = cl
               ! do not test this element again
                  flag(j) = .false.
               ! done with j
                  exit
               endif
            enddo
         endif
      enddo
   endif
enddo
!
! output classes
! and analyze SOPs
!
elem_cl = 0; cl_rotchar = ''
rotsum = 0
maxrot = -6; maxrot_cl = 0
minrot = 6; minrot_cl = 0
inv_cl = 0
sigmah_cl = 0
write(13, '(" Symmetry analysis of given system and SOPs",/ &
 &         ," ------------------------------------------",/)') 
write(13, '(/," The ",i2," SOPs form ",i2," classes",/)') numsop,cl
write(13,*)
do i = 1,cl
   write(13, '(/," Class ",i2," contains: ",/)') i
   num_m = 0
   do j = 1,numsop
      if (class(j) .eq. i) then
         tr = nint(trace(sopmat(:,:,j)))
         rottype(j) = rottabl(tr,symlatd(lsplsymc(j)))
         rotsum(rottype(j)) = rotsum(rottype(j)) + 1
         if (rottype(j) .eq. -2) num_m = num_m + 1
         if (rottype(j) .gt. maxrot) then
            maxrot_cl = i
            maxrot = rottype(j)
         endif
         if (rottype(j) .lt. minrot) then
            minrot_cl = i
            minrot = rottype(j)
         endif
         if (rottype(j) .eq. -1) inv_cl = i
         elem_cl(i) = elem_cl(i) + 1
         cl_rotchar(i) = rotchar(rottype(j))
         ! output characteristics of SOPs, seperated in classes
         ! compute direction of rotation axis
         call eigvec(sopmat(:,:,j),dirvec,symlatd(lsplsymc(j)))
         ! for Cn and Sn with n > 2 compute pos/neg sense of rotation
         ! M. Boisen and G. Gibbs, Mathematical Crystallography, Rev. Mineral. 15 (1990); theorem TA3.9
         if (abs(rottype(j)) .gt. 2) then
            rdum1 = dirvec(1)*sopmat(3,2,j)
            rdum2 = dirvec(2)*sopmat(3,1,j)
            rdum3 = dirvec(3)*sopmat(2,1,j)
            sense = 'neg'
            if (symlatd(lsplsymc(j)) .gt. 0) then
               if ((rdum3 - rdum2) .gt. 0.d0 .or. rdum1 .gt. 0.d0) sense = 'pos'
            else
               if ((rdum2 - rdum3) .gt. 0.d0 .or. -rdum1 .gt. 0.d0) sense = 'pos'
            endif
            write(13, '(" SOP No. ",i2," :",/,3(" ( ",3f7.3," )",/)," Type : ",a2,/,&
              &         " Rotation axis : ",3f7.3,/," Sense : ",a3,/)') &
              &       j,((sopmat(i1,i2,j),i2=1,3),i1=1,3),cl_rotchar(i),dirvec,sense
         elseif (abs(rottype(j)) .eq. 2) then
            write(13, '(" SOP No. ",i2," :",/,3(" ( ",3f7.3," )",/)," Type : ",a2,/,&
              &         " Rotation axis : ",3f7.3,/)') &
              &       j,((sopmat(i1,i2,j),i2=1,3),i1=1,3),cl_rotchar(i),dirvec
         else
            write(13, '(" SOP No. ",i2," :",/,3(" ( ",3f7.3," )",/)," Type : ",a2,/)') &
              &       j,((sopmat(i1,i2,j),i2=1,3),i1=1,3),cl_rotchar(i)
         endif
      endif
   enddo
   if (num_m .eq. 1) sigmah_cl = i
enddo
!
! Compute Character Table for the Group:
!
! (a) determine structure constants c_ijk by class multiplication
!     C_i C_j = Sum_k c_ijk C_k
!
const_c = 0.d0
do i = 1,cl
 do j = 1,cl
  do sop_1 = 1,numsop
   if (class(sop_1) .eq. i) then
    do sop_2 = 1,numsop
     if (class(sop_2) .eq. j) then 
      call r3mm (sopmat(:,:,sop_1), sopmat(:,:,sop_2), sopmat_12)
!     sopmat_12 = matmat(sopmat(:,:,sop_1),sopmat(:,:,sop_2))
      do k = 1,cl
       do sop_3 = 1,numsop
        if (class(sop_3) .eq. k) then
         if (all(sopmat_12 .eq. sopmat(:,:,sop_3))) const_c(i,j,k) = const_c(i,j,k) + 1.d0
        endif
       enddo
      enddo
     endif
    enddo
   endif
  enddo
 enddo
enddo
do i = 1,cl
 do j = 1,cl
  do k = 1,cl
   const_c(i,j,k) = const_c(i,j,k) / dble(elem_cl(k))
  enddo
 enddo
enddo
!
! (b) compute matrix \Phi = Sum_i u_i M_i
!     with the indeterminants u_i being random numbers
!     [-> replaced by u_i = 3/4*(3/pi)**(i-1) to obtain reproducible results]
!     and M_i(j,k) = c_ijk
!
allocate( indet_u(cl) )
allocate( mat_phi(cl,cl) )
mat_phi = cmplx(0.0,0.0)
!call random_seed
!call random_number(indet_u)
indet_u(1) = 0.75
do i = 2,cl
   indet_u(i) = indet_u(i-1)/pi*3.
enddo
do i = 1,cl
 do j = 1,cl
  do k = 1,cl
   mat_phi(j,k) = mat_phi(j,k) + cmplx(const_c(i,j,k)*indet_u(i),0.d0)
  enddo
 enddo
enddo
!
! (c) compute right eigenvectors of \Phi
!
allocate( charact(cl,cl) )
allocate( eigval(cl) )
lwork = 4*cl
allocate( work(lwork) )
allocate( rwork(2*cl) )
call zgeev('N','V',cl,mat_phi(:,:),cl,eigval,charact,cl,charact,cl,work,lwork,rwork,info)
deallocate( eigval,work,rwork )
!
! normalize eigenvectors in order to have charact(1) = 1
!
do j = 1,cl
   charact(:,j) = charact(:,j) / charact(1,j)
enddo
!
! compute dimension of representation
!
dim_ = 0.d0
do j = 1,cl
  do i = 1,cl
   dim_(j) = dim_(j) + abs(charact(i,j))**2 / dble(elem_cl(i))
  enddo
  dim_(j) = sqrt( dble(numsop) / dim_(j) )
enddo
! characters
do i = 1,cl
   do j = 1,cl
      charact(i,j) = charact(i,j)*dim_(j)/dble(elem_cl(i))
   enddo
enddo
!
! derive (rough) notation of irreducible representations:
! (the arbitrary numbering is ignored)
!
! (i) check for cubic groups
cubic = .false.
if (rotsum(3) .eq. 8) cubic = .true.
irep_ch = ''
do j = 1,cl
! (ii) main notation
   ! dimension 1: A or B
   if (nint(dble(charact(1,j))) .eq. 1) then
      if (.not.cubic .and. maxrot .ge. 2) then
         if (-minrot .eq. 2*maxrot .and. sigmah_cl .eq. 0) then
         ! for Dnd,S2n group take from improper rotation
            if (dble(charact(minrot_cl,j)) .gt. 0.d0) then
               irep_ch(j) = trim(irep_ch(j))//'A'
            else
               irep_ch(j) = trim(irep_ch(j))//'B'
            endif
         elseif (maxrot .eq. 2 .and. rotsum(2) .eq. 3 .and. -minrot .ne. 2*maxrot) then
         ! for D2 and D2h consider all C2 SOPs (which are in seperate classes)
            k = 0
            do i = 1,numsop
               if (rottype(i) .eq. 2) k = k + nint(dble(charact(class(i),j)))
            enddo
            if (k .eq. 3) then 
               irep_ch(j) = trim(irep_ch(j))//'A'
            else
               irep_ch(j) = trim(irep_ch(j))//'B'
            endif
         else
         ! else from proper rotation
            if (dble(charact(maxrot_cl,j)) .gt. 0.d0) then
               irep_ch(j) = trim(irep_ch(j))//'A'
            else
               irep_ch(j) = trim(irep_ch(j))//'B'
            endif
         endif
      else
         ! cubic groups
         irep_ch(j) = trim(irep_ch(j))//'A'
      endif
   ! higher dimensions
   elseif (nint(dble(charact(1,j))) .eq. 2) then
      irep_ch(j) = trim(irep_ch(j))//'E'
   elseif (nint(dble(charact(1,j))) .eq. 3) then
      irep_ch(j) = trim(irep_ch(j))//'T'
   endif
! (iii) gerade or ungerade?
   if (fi .eq. 2) then
      if (nint(dble(charact(inv_cl,j))) .gt. 0) then
         irep_ch(j) = trim(irep_ch(j))//'g'
      else
         irep_ch(j) = trim(irep_ch(j))//'u'
      endif
   endif
! (iv) prime or double prime?
   ! C3h and D3h
   if (maxrot .eq. 3 .and. sigmah_cl .gt. 0) then
      if (nint(dble(charact(sigmah_cl,j))) .gt. 0) then
         irep_ch(j) = trim(irep_ch(j))//''''
      else
         irep_ch(j) = trim(irep_ch(j))//''''''
      endif
   endif
   ! Cs
   if (sigmah_cl .gt. 0 .and. cl .eq. 2) then
      if (nint(dble(charact(sigmah_cl,j))) .gt. 0) then
         irep_ch(j) = trim(irep_ch(j))//''''
      else
         irep_ch(j) = trim(irep_ch(j))//''''''
      endif
   endif
enddo
irep_ch(:) = adjustl(irep_ch(:))
! 
! account for 1-dim conjugate complex IREPs
do i = 1, cl
   if (all(abs(aimag(charact(:, i))) .lt. eps)) cycle
   do j = 1, i
      if (all(aimag(charact(:, i)) + aimag(charact(:, j)) .lt. eps)) then
         irep_ch(i)(1:1) = 'E'
         irep_ch(i) = trim(irep_ch(i))//'a'
         irep_ch(i) = adjustl(irep_ch(i))
         irep_ch(j)(1:1) = 'E'
         irep_ch(j) = trim(irep_ch(j))//'b'
         irep_ch(j) = adjustl(irep_ch(j))
         exit
      endif
   enddo
enddo
!
! add indices 1, 2 (, 3) in case there are ambiguous IREP names
if (maxrot .eq. 2 .and. rotsum(2) .eq. 3 .and. -minrot .ne. 2*maxrot) then
 ! D2 and D2h: arbitrary numbering of B IREPs
 i1 = 1; i2 = 1
 do i = 1, cl
  if (irep_ch(i)(1:2) .eq. 'B ' .or. irep_ch(i)(1:2) .eq. 'Bg') then
     write(mul_index, '(i1.1)') i1
     irep_ch(i) = irep_ch(i)(1:1)//mul_index//irep_ch(i)(2:3) 
     i1 = i1 + 1
  endif
  if (irep_ch(i)(1:2) .eq. 'Bu') then
     write(mul_index, '(i1.1)') i2
     irep_ch(i) = irep_ch(i)(1:1)//mul_index//irep_ch(i)(2:3) 
     i2 = i2 + 1
  endif
 enddo
else
 ! all other cases
 do i = 1, cl
  do j = i, cl
   if (irep_ch(j) .eq. irep_ch(i)) then
    ! same name, need for index according to C_n with descending n
    do i1 = 6, -6, -1
     do i2 = 1, numsop
      if (rottype(i2) .eq. i1) then
       if (abs(dble(charact(class(i2), j))-dble(charact(class(i2), i))) .lt. eps) then
        exit
       else
        if (dble(charact(class(i2), i)) .gt. 0.d0) then
         irep_ch(i) = irep_ch(i)(1:1)//'1'//irep_ch(i)(2:3)
         irep_ch(j) = irep_ch(j)(1:1)//'2'//irep_ch(j)(2:3)
         goto 13
        else
         irep_ch(i) = irep_ch(i)(1:1)//'2'//irep_ch(i)(2:3)
         irep_ch(j) = irep_ch(j)(1:1)//'1'//irep_ch(j)(2:3)
         goto 13
        endif
       endif
      endif
     enddo
    enddo
   endif
  enddo
  13 continue
 enddo
endif
         
!
! Determine Point Group from SOPs
!
! reduce group in case I is present
if (fi .eq. 2) then
   do i = 1,numsop-1
      if (representative(i)) then
         do j = i+1,numsop
            if (representative(j)) then
               if (all(sopmat(:,:,j) .eq. -sopmat(:,:,i))) then
                  representative(j) = .false.
               endif
            endif
         enddo
      endif
   enddo
endif 
!
! reduce to representative group in order to determine point group
! and sum over properties of rotation part in representative group; only integer tr and det occur
!
numrep = 0
rotsum = 0
do i = 1,numsop
   if (representative(i)) then
      numrep = numrep + 1
      tr = nint(trace(sopmat(:,:,i)))
      rottype(i) = rottabl(tr,symlatd(lsplsymc(i)))
      rotsum(rottype(i)) = rotsum(rottype(i)) + 1
   endif
enddo
!
!
! look-up table taken from Grosse-Kunstleve, Acta Cryst. (1999) A55, 383-395
hmsymb = 'NONE '; schsymb = 'NONE '
sym_out = .true.
!
if (fi .eq. 2) then
   centric = .true.
elseif (fi .eq. 1) then
   centric = .false.
endif
!
if (rotsum(3)+rotsum(-3) .eq. 8) then
   if (numrep .eq. 24) then
      if (.not.centric) then
         if (rotsum(4) .eq. 6) then
            hmsymb = '432  '
            schsymb = 'O    '
            goto 17
         elseif (rotsum(-4) .eq. 6) then
            hmsymb = '-43m '
            schsymb = 'Td   '
            goto 17
         endif
      elseif (centric) then
         hmsymb = 'm-3m '
         schsymb = 'Oh   '
         goto 17
      endif
   elseif (numrep .eq. 12) then
      if (.not.centric) then
         hmsymb = '23   '
         schsymb = 'T    '
         goto 17
      elseif (centric) then
         hmsymb = 'm-3  '
         schsymb = 'Th   '
         goto 17
      endif
   endif
endif
!
if (rotsum(6)+rotsum(-6) .eq. 2) then
   if (numrep .eq. 12) then
      if (.not.centric) then
         if (rotsum(6) .eq. 2) then
            if (rotsum(2) .eq. 7) then
               hmsymb = '622  '
               schsymb = 'D6   '
               goto 17
            elseif (rotsum(-2) .eq. 6) then
               hmsymb = '6mm  '
               schsymb = 'C6v  '
               goto 17
            endif
         elseif (rotsum(-6) .eq. 2) then
            hmsymb = '-6m2 '
            schsymb = 'D3h  '
            goto 17
         endif
      elseif (centric) then
         hmsymb = '6/mmm'
         schsymb = 'D6h  '
         goto 17
      endif
   elseif (numrep .eq. 6) then
     if (.not.centric) then
         if (rotsum(6) .eq. 2) then
            hmsymb = '6    '
            schsymb = 'C6   '
            goto 17
         elseif (rotsum(-6) .eq. 2) then
            hmsymb = '-6   '
            schsymb = 'C3h  '
            goto 17
         endif
      elseif (centric) then
         hmsymb = '6/m  '
         schsymb = 'C6h  '
         goto 17
      endif
   endif
endif
!
if (rotsum(3)+rotsum(-3) .eq. 2) then
   if (numrep .eq. 3) then
      if (.not.centric) then
         hmsymb = '3    '
         schsymb = 'C3   '
         goto 17
      elseif (centric) then
         hmsymb = '-3   '
         schsymb = 'S6   '
         goto 17
      endif
   elseif (numrep .eq. 6) then
      if (.not.centric) then
         if (rotsum(2) .eq. 3) then
            hmsymb = '32   '
            schsymb = 'D3   '
            goto 17
         elseif (rotsum(-2) .eq. 3) then
            hmsymb = '3m   '
            schsymb = 'C3v  '
            goto 17
         endif
      elseif (centric) then
         hmsymb = '-3m  '
         schsymb = 'D3d  '
         goto 17
      endif
   endif
endif
!
if (rotsum(4)+rotsum(-4) .eq. 2) then
   if (numrep .eq. 4) then
      if (.not.centric) then
         if (rotsum(4) .eq. 2) then
            hmsymb = '4    '
            schsymb = 'C4   '
            goto 17
         elseif (rotsum(-4) .eq. 2) then
            hmsymb = '-4   '
            schsymb = 'S4   '
            goto 17
         endif
      elseif (centric) then
         hmsymb = '4/m  '
         schsymb = 'C4h  '
         goto 17
      endif
   elseif (numrep .eq. 8) then
      if (.not.centric) then
         if (rotsum(4) .eq. 2) then
            if (rotsum(2) .eq. 5) then
               hmsymb = '422  '
               schsymb = 'D4   '
               goto 17
            elseif (rotsum(-2) .eq. 4) then
               hmsymb = '4mm  '
               schsymb = 'C4v  '
               goto 17
            endif
         elseif (rotsum(-4) .eq. 2) then
            hmsymb = '-4m2 '
            schsymb = 'D2d  '
            goto 17
         endif
      elseif (centric) then
         hmsymb = '4/mmm'
         schsymb = 'D4h  '
         goto 17
      endif
   endif
endif
!
if (rotsum(2)+rotsum(-2) .eq. 3) then
   if (.not.centric) then
      if (rotsum(2) .eq. 3) then
         hmsymb = '222  '
         schsymb = 'D2   '
         goto 17
      elseif (rotsum(-2) .eq. 2) then
         hmsymb = 'mm2  '
         schsymb = 'C2v  '
         goto 17
      endif
   elseif (centric) then
      hmsymb = 'mmm  '
      schsymb = 'D2h  '
      goto 17
   endif
elseif (rotsum(2)+rotsum(-2) .eq. 1) then
   if (.not.centric) then
      if (rotsum(2) .eq. 1) then
         hmsymb = '2    '
         schsymb = 'C2   '
         goto 17
      elseif (rotsum(-2) .eq. 1) then
         hmsymb = 'm    '
         schsymb = 'Cs   '
         goto 17
      endif
   elseif (centric) then
      hmsymb = '2/m  '
      schsymb = 'C2h  '
      goto 17
   endif
elseif (numsop .eq. 1) then
   sym_out = .false.
   if (.not.centric) then
      hmsymb = '1    '
      schsymb = 'C1   '
   elseif (centric) then
      hmsymb = '-1   '
      schsymb = 'Ci   '
   endif
endif
17 continue
!
write(13, '(/," The (isomorphic) point group is ",a5," (Hermann-Maugin) or ", &
  &                                               a5," (Schoenflies)")') hmsymb,schsymb
write(13,*)
write(13, '(" Character Table",/, &
      &     " ---------------",/)')
write(13, '("              ",20("         ",i2,a2),/)') (elem_cl(i),cl_rotchar(i),i=1,cl)
do j = 1,cl
   write(13,'(i5,a10,20(" ",f5.2,",",f5.2," "))') &
 &       j,irep_ch(j),(dble(charact(i,j)),aimag(charact(i,j)),i=1,cl)
enddo
write(13,*)
! --- end of character table code
!
! determine Raman active ireps
!
allocate( atom_sop(natmtot, numsop) )
!
! for convenience we continue working with a re-defined array of equivalent atoms,
! one index for all atoms
do j = 1,numsop
   i = 0
   k = 0
   Do is = 1, nspecies
      Do ia = 1, natoms (is)
         i = i + 1
         atom_sop(i, j) = ieqatom(ia, is, j) + k
      enddo
      k = k + natoms (is)
   enddo
enddo
!
! determine groups of atoms which are mixed by SOPs
allocate( gr_atoms_no(natmtot) )
allocate( gr_atoms(natmtot,natmtot) )
gr = 0; gr_atoms_no = 0; gr_atoms = 0; grouped = .false.
do ia = 1,natmtot
   if (.not. grouped(ia)) then
      ! start new group
      gr = gr + 1
      gr_atoms_no(gr) = 1
      gr_atoms(gr,gr_atoms_no(gr)) = ia
      ! check for other members of the group
      do ja = ia+1,natmtot
         mix = .false.
         do j = 1,numsop
            if (atom_sop(ia, j) .eq. ja) mix = .true.
         enddo
         if (mix) then
            ! add to group
            gr_atoms_no(gr) = gr_atoms_no(gr) + 1
            gr_atoms(gr,gr_atoms_no(gr)) = ja
            ! do not consider further
            grouped(ja) = .true.
         endif
      enddo
   endif
enddo
         
! compute characters for polar and axial vectors, and atom equivalence
char_equiv = 0.d0
do i = 1,cl
   do j = 1,numsop
      if (class(j) .eq. i) exit
   enddo
   ! radial vector
   char_vec(i) = trace(sopmat(:,:,j))
   ! axial vector = rotation = anti-symmetric components of a rank 2 tensor
   char_rot(i) = dble(symlatd(lsplsymc(j)))*char_vec(i)
   ! symmetric tensor componenets
   char_raman(i) = char_vec(i)*char_vec(i) - char_rot(i)
   ! conformational equivalence
   do ia = 1,natmtot
      if (atom_sop(ia, j) .eq. ia) char_equiv(i) = char_equiv(i) + 1.d0
   enddo
enddo
! reduction
freq_fac_vec = zzero
freq_fac_raman = zzero
freq_fac_rot = zzero
freq_fac_vib = zzero
vec_sym = ''; rot_sym = ''; raman_sym = ''; vib_sym = ''
vib_mode = .false.
raman_active = .false.
no_vec_nonorth = 0
do j = 1,cl
   do i = 1,cl
      freq_fac_vec(j) = freq_fac_vec(j) + conjg(charact(i,j))*char_vec(i)*dble(elem_cl(i))
      freq_fac_raman(j) = freq_fac_raman(j) + conjg(charact(i,j))*char_raman(i)*dble(elem_cl(i))
      freq_fac_rot(j) = freq_fac_rot(j) + conjg(charact(i,j))*char_rot(i)*dble(elem_cl(i))
      freq_fac_vib(j) = freq_fac_vib(j) + conjg(charact(i,j))*char_equiv(i)*char_vec(i)*dble(elem_cl(i))
   enddo
   if (nint(dble(freq_fac_vec(j))/numsop) .gt. 0) then
      write(chnum,'(i5)') nint(dble(freq_fac_vec(j))/numsop)
      !write(vec_sym,*) trim(vec_sym),trim(chnum),' '
      !write(chnum,'(i2)') j
      !write(vec_sym,*) trim(vec_sym),trim(irep_ch(j)),'   '
      vec_sym = trim(adjustl(vec_sym))//'   '//trim(adjustl(chnum))//trim(adjustl(irep_ch(j)))
      vec_sym = adjustl(vec_sym)
   endif
   if (nint(dble(freq_fac_rot(j))/numsop) .gt. 0) then
      write(chnum,'(i5)') nint(dble(freq_fac_rot(j))/numsop)
      !write(rot_sym,*) trim(rot_sym),trim(chnum),' '
      !write(chnum,'(i2)') j
      !write(rot_sym,*) trim(rot_sym),trim(irep_ch(j)),'   '
      rot_sym = trim(adjustl(rot_sym))//'   '//trim(adjustl(chnum))//trim(adjustl(irep_ch(j)))
      rot_sym = adjustl(rot_sym)
   endif
   if (nint(dble(freq_fac_raman(j))/numsop) .gt. 0) then
      raman_active(j) = .true.
      write(chnum,'(i5)') nint(dble(freq_fac_raman(j))/numsop)
      !write(raman_sym,*) trim(raman_sym),trim(chnum),' '
      !write(chnum,'(i2)') j
      !write(raman_sym,*) trim(raman_sym),trim(irep_ch(j)),'   '
      raman_sym = trim(adjustl(raman_sym))//'   '//trim(adjustl(chnum))//trim(adjustl(irep_ch(j)))
      raman_sym = adjustl(raman_sym)
   endif
   if (nint(dble(freq_fac_vib(j))/numsop) .gt. 0) then
      vib_mode(j) = .true.
      no_vec_nonorth = no_vec_nonorth + 2**(nint(dble(freq_fac_vib(j))/numsop)-1)
      vib_ireps(j) = nint(dble(freq_fac_vib(j))/numsop)
      if (vib_ireps(j) .gt. nint(dble(charact(1, j)))) sym_out = .false.
      write(chnum,'(i5)') nint(dble(freq_fac_vib(j))/numsop)
      !write(vib_sym,*) trim(vib_sym),trim(chnum),' '
      !write(chnum,'(i2)') j
      !write(vib_sym,*) trim(vib_sym),trim(irep_ch(j)),'   '
      vib_sym = trim(adjustl(vib_sym))//'   '//trim(adjustl(chnum))//trim(adjustl(irep_ch(j)))
      vib_sym = adjustl(vib_sym)
   endif
enddo
write(13, '(" Symmetries of phonons are : ",a)') trim(vib_sym)
write(13, '(" IR active                 : ",a)') trim(vec_sym)
write(13, '(" Raman active              : ",a)') trim(raman_sym)
write(13, '(" Rotations                 : ",a,/)') trim(rot_sym)
!
!
!
!
close(13)
!
!
! ------------------------------------------------------------------------------
contains
function trace(mat)
!
! compute trace of 3x3 matrix
!
implicit none
real(8) :: mat(3,3),trace
trace = mat(1,1) + mat(2,2) + mat(3,3)
return
end function trace
!
!
!
end subroutine construct_chartabl

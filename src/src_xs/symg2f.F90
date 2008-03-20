
! Copyright (C) 2007-2008 S. Sagmeister, J. K. Dewhurst, S. Sharma and 
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symg2f
! !INTERFACE:
subroutine symg2f(vpl,ngp,igpig,fg)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   vpl     :vpl
!   ngp     :ngp
!   igpig   :igpig
!   fg      :fg
! !DESCRIPTION:
!   Symmetrises a real scalar function given in $G$-space
!   $$ f_{\bf q}({\bf G},{\bf G'}) = \frac{1}{\Omega}\int {\rm d}^3r 
!   e^{-i{\bf (G+q)r}} f({\bf r},{\bf r'}) e^{i{\bf (G'}+{\bf q}){\bf r'}} $$
!   of which its real-space representation is invariant under application of
!   a symmetry  operation to both spatial variables, i.e.:
!   $$ f({\bf r},{\bf r'}) = f(\{\alpha|\tau\}{\bf r},\{\alpha|\tau\}{\bf r'}).
!   $$
!   The symmetry operations are restricted to the subset that leaves the
!   ${\bf q}$-vector unaltered
!   $$ {\bf q}=\alpha^{-1}{\bf q} + {\bf G}_{\alpha},  $$
!   building up the small (little) group of {\bf q}.
!   For a function obeying the latter symmetry property one can derive a
!   symmetry relation for the Fourier-coefficients as well:
!   $$ f_{\bf q}({\bf G},{\bf G'}) =  
!   e^{i({\bf G'}-{\bf G})\tau_{\alpha^{-1}}} 
!   f_{\bf q}(\alpha[{\bf G}+{\bf G}_{\alpha}],
!   \alpha[{\bf G'}+{\bf G}_{\alpha}]) =: f^{(\alpha)}_{\bf q}({\bf G},
!   {\bf G'}).$$
!   Since this property may not exactly be fulfilled numerically, an average
!   of the rotated expressions $f^{(\alpha)}$ is taken:
!   $$  \bar{f}_{\bf q}({\bf G},{\bf G'})= \frac{1}{N_{\bf q}}
!   \sum_{\alpha}^{\mathcal{G}({\bf q})}
!   f^{(\alpha)}_{\bf q}({\bf G},{\bf G'}) $$
!   where $\mathcal{G}({\bf q})$ denotes the small (little) group of ${\bf q}$
!   and $N_{\bf q}$ its number of elements.
!
! !REVISION HISTORY:
!   Created March 2008 (SAG)
!EOP
!BOC
  implicit none
  ! arguments
  real(8), intent(in) :: vpl(3)
  integer, intent(in) :: ngp
  integer, intent(in) :: igpig(ngp)
  complex(8), intent(inout) :: fg(ngp,ngp)
  ! local variables
  integer j,isym,lspl,ilspl,sym(3,3)
  integer iv1(3),iv2(3),igp1,igp2,jgp1,jgp2
  integer :: nscq,scq(nsymcrys),ivgscq(3,nsymcrys),ivgs(3)
  real(8) vtc(3),t1
  complex(8) zt1
  ! allocatable arrays
  complex(8), allocatable :: fg2(:,:)

  allocate(fg2(ngp,ngp))
  call findgroupq(vpl,epslat,symlat,nsymcrys,lsplsymc,nscq,scq,ivgscq)

  write(*,*) 'symg2f: group of q has elements:',nscq

  fg2(:,:)=zzero
  ! loop over crystal symmetries of group of q
  do j=1,nscq
     isym=scq(j)
     vtc=matmul(avec,vtlsymc(:,isym))
     lspl=lsplsymc(isym)
     ilspl=isymlat(lspl)
     sym(:,:)=symlat(:,:,ilspl)
     ivgs(:)=ivgscq(:,isym)
     do igp1=1,ngp
        ! (s^-1)^LT ( G + G_s ) or in Cartesian s ( G + G_s )
        iv1=matmul(transpose(sym),ivg(:,igpig(igp1))+ivgs)
        jgp1=ivgig(iv1(1),iv1(2),iv1(3))
        do igp2=1,ngp
           iv2=matmul(transpose(sym),ivg(:,igpig(igp2))+ivgs)
           jgp2=ivgig(iv2(1),iv2(2),iv2(3))
           ! complex phase factor for translation
           t1=-dot_product(vgc(:,igp1)-vgc(:,igp2),vtc(:))
           zt1=cmplx(cos(t1),sin(t1),8)
           fg2(jgp1,jgp2) = fg2(jgp1,jgp2) + zt1*fg(igp1,igp2)
        end do
     end do
  end do
  ! normalize
  fg(:,:)=fg2(:,:)/dble(nscq)
  deallocate(fg2)
end subroutine symg2f

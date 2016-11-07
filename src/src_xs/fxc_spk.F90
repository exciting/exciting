!
!
!
!Complex (8), allocatable :: chimi(:,:)
Module m_fxc_spk
      Implicit None
!
Contains
!
!BOP
! !ROUTINE: fxc_spk
! !INTERFACE:
!
!
      Subroutine fxc_spk (fxctype,msiz,ngtot,nw, chim, chim_w, fxc, fxc_w)
! !USES:
         use modmpi
         Use mod_constants, Only: fourpi, zzero, zone
         Use modxs, Only: unitout
         Use invert
         Use m_dyson
         Use modinput
! !INPUT/OUTPUT PARAMETERS:
!   msiz  : matrix size of local field effects (in,integer)
!   chim : model static density response function (e.g., scissors corrected) (in,complex(:,:))
!   fxc   : xc-kernel  (out,complex(:,:))
! !DESCRIPTION:
!   Static spk xc-kernel; PRL 107, 186401 (2011), PRL 114, 146402 (2015)
!   Calculates the symmetrized BO and RBO xc-kernels.
!
! !REVISION HISTORY:
!   Created October 2015 (SR)
!EOP
!BOC
         Implicit None
    ! arguments
         Integer, Intent (In) :: fxctype
         Integer, Intent (In) :: msiz
         Integer, Optional, Intent (In) :: ngtot
         Integer, Optional, Intent (In) :: nw
    ! true if all G-components of fxc are to be considered
         Complex (8), Optional, Intent (In) :: chim (:, :)
         Complex (8), Optional, Intent (In) :: chim_w (:, :, :)
         Complex (8), Optional, Intent (Out) :: fxc (:, :)
         Complex (8), Optional, Intent (Out) :: fxc_w (:, :, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'fxc_spk'
         Integer :: sh (2), ig, iw, n, ig1, ig2, fxctype2
         Complex (8), allocatable :: aux(:,:), chimi(:,:), sqrmat(:,:)
         Complex (8), allocatable :: mat(:,:), chim2(:,:),fxcm(:,:)
         Complex (8), allocatable :: ker(:,:), chi(:,:), idfrpa_mat(:,:)
         Complex (8) :: t1, t2, zdotu, hcor, fxcs,idfrpa,chimodrpa
         Logical :: fullmat

         allocate(aux(msiz,msiz), chimi(msiz,msiz),sqrmat(msiz,msiz))
         allocate(mat(msiz,msiz), chim2(msiz,msiz),ker(msiz,msiz))
         allocate(chi(msiz,msiz),idfrpa_mat(msiz,msiz))
         aux=zzero
         chimi=zzero
         sqrmat=zzero
         mat=zzero              
         chim2=zzero
         ker=zzero
         chi=zzero
         idfrpa_mat=zzero

         Select Case (Abs(fxctype))
         Case(9)
            ! fxc = e^(-1)/chi_KS_00
            fxc (:, :) = zzero
            t1 = 1.d0/chim(1,1)

             If ( msiz.eq.1 ) Then
                fxc (1, 1) =  t1 - 0.5d0 + sqrt(0.25d0 - t1)
             Else
                Call zinvert_lapack(chim, chimi)

                ! build aux = (1-chimi-1/chim00)/2
                aux(:,:) = -chimi(:,:)/2.d0
                Do ig = 1, msiz
                    aux(ig,ig) = aux(ig,ig) + (1.d0 - t1)/2.d0
                End Do

                ! calculate sqrt(aux**2-chimi/chim00)
                call zgemm('n','n',msiz,msiz,msiz,zone,aux,msiz,aux,&
                & msiz,zzero,mat,msiz)
                mat(:,:) = mat(:,:) - chimi(:,:)*t1
                Call sqrtmatrix(msiz,mat,sqrmat,1)

                fxc(:,:) = - aux(:,:) + sqrmat(:,:)
             End If
         Case(10)
            ! fxc =  1/(e_M * chi_KS_00)
            n = msiz
            ! Set up modified Coulomb kernel
            ker = zzero
            do ig = 2, n
                ker(ig,ig) = zone
            end do
            ! Calculate \bar(\chi)_RPA i.e., modified chi in RPA
            Call dyson(n,chim(1:n,1:n),ker(1:n,1:n),chi(1:n,1:n))
            chimodrpa = chi(1,1)

            t1 = zone/chimodrpa
            t2 = zone/chim(1,1)
            ! Calculate BO_SCALAR kernel
            fxc = zzero
            fxc(1,1) = 0.5*(t1+t2-1.0d0) + sqrt(0.25*(t1+t2-1.0d0)**2-t1*t2)
         Case(11)
            ! fxc =  1/(e_{M,RPA} * \bar\chi_{RPA,00})
            n = msiz
            ! Calculate 1/\epsilon_M(RPA)
            ker=zzero
            do ig = 1, n
                ker(ig,ig) = zone
            end do
            Call dyson(n,chim(1:n,1:n),ker(1:n,1:n),chi(1:n,1:n))
            idfrpa = zone + chi(1,1)
            
            ! Calculate \bar\chi_RPA
            chi = zzero
            ker(1,1) = zzero
            Call dyson(n,chim(1:n,1:n),ker(1:n,1:n),chi(1:n,1:n))
            chimodrpa = chi(1,1)

            ! Calculate RBO kernel
            fxc = zzero
            fxc(1,1) = idfrpa/chimodrpa
         Case Default
            Write (unitout,*)
            Write (unitout, '("Error(fxc_spk): spktype not implemented for analytical calculation")')
            Write (unitout,*)
            Call terminate

         End Select
         deallocate (aux,chimi,sqrmat,mat,chim2,ker,chi,idfrpa_mat)
       End Subroutine fxc_spk

End Module m_fxc_spk
!EOC









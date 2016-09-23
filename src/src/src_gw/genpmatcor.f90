!BOP
! !ROUTINE: genpmat
! !INTERFACE:
subroutine genpmatcor(vpl,ngp,apwalm,evecfv,evecsv,pmc)
! !USES:
    use modinput
    use modmain
    use modgw
! !INPUT/OUTPUT PARAMETERS:
!   vpl    : k-vector (in,real(3))
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evecfv : first-variational eigenvector (in,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   pmc    : momentum matrix elements (out,complex(3,nstsv,nstsv))
! !DESCRIPTION:
!   Calculates the momentum matrix elements
!   $$ p_{ij}=\langle\Psi_{i,{\bf k}}|-i\nabla|\Psi_{j,{\bf k}}\rangle. $$
!
! !REVISION HISTORY:
!   Created August 2006 (RGA)
!   Revisited June 2011 (DIN)
!EOP
!BOC
    implicit none
    ! arguments
    real(8), intent(in)     :: vpl(3)
    integer, intent(in)     :: ngp
    complex(8), intent(in)  :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
    complex(8), intent(in)  :: evecfv(nmatmax,nstfv)
    complex(8), intent(in)  :: evecsv(nstsv,nstsv)
    complex(8), intent(out) :: pmc(3,ncg,nstsv)
!   local variables
    integer :: ic, icg, ie2, ispn, i, j, n
    integer :: is, ia, ias, l, m, lm, ir, irc
    real(8) :: arg, t1, t2
    real(8), dimension(nrcmtmax) :: fr, fr1, fr2, gr
    real(8), dimension(3,nrcmtmax) :: cf
    complex(8) :: phs, zt1, zv(3)
!   allocatable arrays
    complex(8), allocatable :: wfmt(:,:)
    complex(8), allocatable :: gwfmt(:,:,:)
    complex(8), allocatable :: pm(:,:,:)

!-------------------------------------------------------------------------------    
    
    allocate(wfmt(lmmaxapw,nrcmtmax))
    allocate(gwfmt(lmmaxapw,nrcmtmax,3))
    allocate(pm(3,ncg,nstfv))
    pm(:,:,:) = zzero
    
    do ie2 = 1, nstsv
    
      icg = 0
      do is = 1, nspecies
        n = lmmaxvr*nrcmt (is)
        do ia = 1, natoms(is)
        
          ias = idxas(ia,is)
          arg = atposl(1,ia,is)*vpl(1)+ &
          &     atposl(2,ia,is)*vpl(2)+ &
          &     atposl(3,ia,is)*vpl(3)
          phs = cmplx(cos(2.0d0*pi*arg),-sin(2.0d0*pi*arg),8)
 
          ! calculate the wavefunction
          call wavefmt(input%groundstate%lradstep, &
          &  input%groundstate%lmaxapw,is,ia,ngp,apwalm, &
          &  evecfv(:,ie2),lmmaxapw,wfmt)
          
          ! calculate the gradient
          call gradzfmt(input%groundstate%lmaxapw,nrcmt(is), &
          &  rcmt(:,is),lmmaxapw,nrcmtmax,wfmt,gwfmt)
          
          do ic = 1, ncore(is)
            l = spl(ic,is)
            irc = 0
            do ir = 1, nrmt(is), input%groundstate%lradstep
              irc = irc+1
              fr(irc) = ucore(ir,1,ic,ias)*spr(ir,is)*spr(ir,is)
            end do ! ir
            do m = -l, l
              lm = idxlm(l,m)
              ! combined atom+lm index
              icg = icg+1
              do i = 1, 3
                do irc = 1, nrcmt(is)
                  fr1(irc) = fr(irc)*dble(gwfmt(lm,irc,i))
                  fr2(irc) = fr(irc)*aimag(gwfmt(lm,irc,i))
                enddo  
                call fderiv(-1,nrcmt(is),rcmt(1,is),fr1,gr,cf)
                t1 = gr(nrcmt(is))
                call fderiv(-1,nrcmt(is),rcmt(1,is),fr2,gr,cf)
                t2 = gr(nrcmt(is))
                pm(i,icg,ie2) = phs*cmplx(t1,t2,8)
              enddo ! i
            end do ! m
          end do ! ic
          
        end do ! ia
      end do ! is
    end do ! ie2
    
    ! multiply by -i
    pm(1:3,:,:) = -zi*pm(1:3,:,:)

! NEED TO BE CHECKED !!    
    ! compute the second-variational momentum matrix elements
    if (input%groundstate%tevecsv) then
      do icg = 1, ncg
        do j = 1, nstsv
          zv(:) = 0.d0
          do ispn = 1, nspinor
            l = (ispn-1)*nstfv
            do ie2 = 1, nstfv
              l = l+1
              zt1 = evecsv(l,j)
              zv(:) = zv(:) + zt1*pm(:,icg,ie2)
            end do
          end do
          pmc(:,icg,j) = zv(:)
        end do
      end do
    else
       pmc(:,:,:) = pm(:,:,:)
    end if
    deallocate(pm)

    return
end subroutine
!EOC

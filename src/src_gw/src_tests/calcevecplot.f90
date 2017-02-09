!BOP
!!ROUTINE: calcevecplot

!!INTERFACE:

    subroutine calcevecplot(ik,ib,evec)

!!DESCRIPTION:
!
!This subroutine calculates the real space representation of the product
!of two eigenvectors in the line joining the two given atoms for ploting. 
!
!!USES:
    use modinput
    use modmain
    use modgw
    use mod_rpath

!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: ik
    integer(4), intent(in) :: ib

!!OUTPUT PARAMETERS:
    complex(8), intent(out) :: evec(*)
      
!!LOCAL VARIABLES:
    integer(4) :: is1, is2, ia1, ia2
    integer(4) :: ir, jr, kr
    integer(4) :: i, j, igp
    real(8)    :: phs, phsat, ri(3)
     
    complex(8), allocatable :: apwalm(:,:,:,:)
    complex(8), allocatable :: evecmtlm(:,:)
    complex(8), allocatable :: evecmt(:)
    complex(8), allocatable :: evecsv(:,:,:)
    complex(8), allocatable :: zylm(:)
    
!EOP
!BOC
    is1 = rpath%atom1(1)
    ia1 = rpath%atom1(2)
    
    is2 = rpath%atom2(1)
    ia2 = rpath%atom2(2)
    
    evec(1:rpath%nptot) = zzero
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
    allocate(zylm(lmmaxapw))
    allocate(evecmtlm(lmmaxapw,nrmtmax))
    allocate(evecmt(nrmtmax))
    
    allocate(evecsv(nmatmax,nstfv,nspinor))
    call getevecsvgw_new('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
    
    call match(Gkset%ngk(1,ik), &
    &          Gkset%gkc(:,1,ik), &
    &          Gkset%tpgkc(:,:,1,ik), &
    &          Gkset%sfacgk(:,:,1,ik),&
    &          apwalm)

    call ylm(rpath%rvec,input%groundstate%lmaxapw,zylm)

    ! calculate the radial wavefunctions of atom 1
    call wavefmt(1,input%groundstate%lmaxapw,is1,ia1,Gkset%ngk(1,ik), &
    &            apwalm,evecsv(:,ib,1),lmmaxapw,evecmtlm)
    call zgemv('T',lmmaxapw,nrmt(is1),zone,evecmtlm,lmmaxapw,zylm,1, &
    &          zzero,evecmt,1)
      
    do ir = 1, nrmt(is1)
      evec(ir) = evecmt(ir)
    end do
      
    ! Calculate the phase of the plane waves due to the change of origin
    do igp = 1, Gkset%ngk(1,ik)
      phsat = 2.0d0*pi*(Gkset%vgkl(1,igp,1,ik)*atposl(1,ia1,is1)+ &
      &                 Gkset%vgkl(2,igp,1,ik)*atposl(2,ia1,is1)+ &
      &                 Gkset%vgkl(3,igp,1,ik)*atposl(3,ia1,is1))
      ! last point is excluded
      do i = nrmt(is1)+1, nrmt(is1)+rpath%nptir-1
        ri(1:3) = rpath%rvec(1:3)/rpath%rlen*rpath%r(i,1)
        phs = phsat+(Gkset%vgkc(1,igp,1,ik)*ri(1)+ &
        &            Gkset%vgkc(2,igp,1,ik)*ri(2)+ &
        &            Gkset%vgkc(3,igp,1,ik)*ri(3))
        evec(i) = evec(i)+ &
        &         evecsv(igp,ib,1)/sqrt(omega)*cmplx(dcos(phs),dsin(phs),8)
      end do
    end do
     
    ! calculate the radial wavefunctions of atom 2 (if needed)
    if ((is1==is2).and.(ia1==ia2)) then
      do ir = 1, nrmt(is1)
        jr = nrmt(is2)-ir+1
        kr = nrmt(is1)+rpath%nptir+ir-1
        evec(kr) = evecmt(jr)
      end do
    else
      call ylm(-rpath%rvec,input%groundstate%lmaxapw,zylm)
      call wavefmt(input%groundstate%lradstep,input%groundstate%lmaxapw,&
      &    is2,ia2,Gkset%ngk(1,ik),apwalm,evecsv(:,ib,1),lmmaxapw,evecmtlm)
      call zgemv('T',lmmaxapw,nrmt(is2),zone,evecmtlm,lmmaxapw,zylm,1, &
      &    zzero,evecmt,1)
      do ir = 1, nrmt(is2)
        jr = nrmt(is2)-ir+1
        kr = nrmt(is1)+rpath%nptir+ir-1
        evec(kr) = evecmt(jr)
      end do
    end if
      
    deallocate(apwalm)
    deallocate(zylm)
    deallocate(evecmtlm)
    deallocate(evecmt)
    deallocate(evecsv)
      
    return
end subroutine
!EOC      

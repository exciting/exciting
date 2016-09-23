!BOP
!
!!ROUTINE: plotevecmix
!
!!INTERFACE:

subroutine plotevecmix(iq,ik,jk,ib1,ib2,atom1,atom2)

!!DESCRIPTION:
!
!This subroutine calculates the real space representation of the product
!of two eigenvectors using the expansion in mixed basis functions:
!
!\begin{equation}
!\Psi_{n\vec{k}}(\vec{r})\Psi^*_{m\vec{k}-\vec{q}}(\vec{r})=\sum\limits_i
!M^i_{nm}(\vec{k},\vec{q})\tilde{\chi}_i^{\vec{q}}(\vec{r})
!\end{equation}
!In the line joining the two given atoms for ploting. 

!!USES:
    use modmain
    use modgw
    use mod_rpath

!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: iq  ! q-point index
    integer(4), intent(in) :: ik  ! The k-point for which the (L)APW+lo function is ploted
    integer(4), intent(in) :: jk  ! The k-point for which the (L)APW+lo function is ploted
    integer(4), intent(in) :: ib1 ! THe band index of the function
    integer(4), intent(in) :: ib2 ! THe band index of the function
    integer(4), intent(in) :: atom1 ! the atom used as origin
    integer(4), intent(in) :: atom2 ! the atom used as final position
      
!!LOCAL VARIABLES:
    integer(4) :: i, j
    integer(4) :: ia1, ia2, is1, is2, ias1, ias2
    integer(4) :: igq, jgq
    integer(4) :: im, imix, irm
    integer(4) :: ir, jr, kr
    integer(4) :: l, lm, m
    real(8) :: t1, ri(3)
    real(8) :: phs, phsat
    complex(8) :: eph
    complex(8), allocatable :: evecmix(:)
    complex(8), allocatable :: zylm(:)
    character(len=64) :: filename

!EOP
!BOC
    call boxmsg(6,'-','PLOTEVECMIX')
      
    write(*,*) 'Parameters:'
    write(*,*) '1 k-point number (iik): ', ik
    write(*,*) '2 k-point number (jjk): ', jk
    write(*,*) '1 band index (ib1): ', ib1
    write(*,*) '2 band index (ib2): ', ib2
    write(*,*) 'atom 1 (at1): ', atom1
    write(*,*) 'atom 2 (at2): ', atom2
    write(*,*)

    if ((atom1<1).or.(atom1.gt.natmtot)) stop 'atom1 is wrong'
    if ((atom2<1).or.(atom2.gt.natmtot)) stop 'atom2 is wrong'

    ! Set the name of the output file
5   format('evmix-',i4,'-',i4,'-',i4,'-',i4,'-',i4,'-',i4,'.out')
    write(filename,5) ik, jk, ib1, ib2, atom1, atom2
    call str_strip(filename)
    open(unit=72,file=filename,status='unknown')

    is1 = rpath%atom1(1)
    ia1 = rpath%atom1(2)
    ias1 = idxas(ia1,is1)
    
    is2 = rpath%atom2(1)
    ia2 = rpath%atom2(2)
    ias2 = idxas(ia2,is2)
    
    allocate(zylm((maxbigl+1)*(maxbigl+1)))
    allocate(evecmix(rpath%nptot))
    evecmix(:) = zzero
    
    ! MT 1
    phsat = 2.0d0*pi*(kqset%vql(1,iq)*atposl(1,ia1,is1)+ &
                      kqset%vql(2,iq)*atposl(2,ia1,is1)+ &
                      kqset%vql(3,iq)*atposl(3,ia1,is1))
    eph = cmplx(cos(phsat),sin(phsat),8)
    
    call ylm(rpath%rvec,maxbigl,zylm)
    
    im = 0
    do irm = 1, nmix(ias1)
      l = bigl(irm,ias1)  
      do m = -l, l
        im = im+1
        imix = locmixind(ias1,im)
        lm = idxlm(l,m)
        do ir = 1, nrmt(is1)
          t1 = 1.0d0/spr(ir,is1)
          evecmix(ir) = evecmix(ir)+minmmat(imix,ib1,ib2)*zylm(lm)* &
          &             cmplx(umix(ir,irm,ias1)*t1,0.0d0,8)*eph
        end do
      end do
    end do

    ! IR
    do igq = 1, Gqset%ngk(1,iq)
      imix = locmatsiz+igq
      do jgq = 1, Gqset%ngk(1,iq)
        phsat = 2.0d0*pi*(Gqset%vgkl(1,jgq,1,iq)*atposl(1,ia1,is1)+ &
        &                 Gqset%vgkl(2,jgq,1,iq)*atposl(2,ia1,is1)+ &
        &                 Gqset%vgkl(2,jgq,1,iq)*atposl(3,ia1,is1))
        do i = nrmt(is1)+1, nrmt(is1)+rpath%nptir-1
          ri(1:3) = rpath%rvec(1:3)/rpath%rlen*rpath%r(i,1)
          phs = phsat+(Gqset%vgkc(1,jgq,1,iq)*ri(1)+ &
          &            Gqset%vgkc(2,jgq,1,iq)*ri(2)+ &
          &            Gqset%vgkc(3,jgq,1,iq)*ri(3))
          eph = cmplx(cos(phs),sin(phs),8)/sqrt(omega)
          evecmix(i) = evecmix(i)+minmmat(imix,ib1,ib2)*conjg(sgi(jgq,igq))*eph
        end do 
      enddo  
    enddo
    
    ! MT 2
    phsat = 2.0d0*pi*(kqset%vql(1,iq)*atposl(1,ia2,is2)+ &
                      kqset%vql(2,iq)*atposl(2,ia2,is2)+ &
                      kqset%vql(3,iq)*atposl(3,ia2,is2))
    eph = cmplx(cos(phsat),sin(phsat),8)
    
    call ylm(-1.d0*rpath%rvec,maxbigl,zylm)

    im = 0
    do irm = 1, nmix(ias2)
      l = bigl(irm,ias2)  
      do m = -l, l
        im = im+1
        imix = locmixind(ias2,im)
        lm = idxlm(l,m)
        do ir = 1, nrmt(is2)
          jr = nrmt(is2)-ir+1
          kr = nrmt(is1)+rpath%nptir+ir-1
          t1 = 1.0d0/spr(jr,is2)
          evecmix(kr) = evecmix(kr)+minmmat(imix,ib1,ib2)*zylm(lm)* &
          &             cmplx(umix(jr,irm,ias2)*t1,0.0d0,8)*eph
        end do
      end do
    end do    
    
    ! Output
    do i = 1, rpath%nptot
      write(72,'(4g18.10)') rpath%r(i,1), real(evecmix(i)), &
      &                     aimag(evecmix(i)), abs(evecmix(i))   
    end do

    deallocate(evecmix)
    deallocate(zylm)
      
    close(72)
      
    return
end subroutine
!EOC

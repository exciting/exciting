
subroutine calc_mb_functions(iq,npt,wfmb)

    use modmain
    use modgw
    use mod_rpath
    
    implicit none
    integer, intent(in) :: iq
    integer, intent(in) :: npt
    complex(8), intent(out) :: wfmb(matsizmax,npt)
    
    integer :: is1, ia1, ias1, is2, ia2, ias2
    integer :: igq, jgq, ir, jr, kr, l, m, lm, imix, im, irm
    real(8) :: t1, ri(3), phs, phsat
    complex(8) :: eph
    
    complex(8), allocatable :: zylm(:)

    wfmb(:,:) = 0.d0
    
    ! calculate the MB functions
    is1 = rpath%atom1(1)
    ia1 = rpath%atom1(2)
    ias1 = idxas(ia1,is1)
    
    is2 = rpath%atom2(1)
    ia2 = rpath%atom2(2)
    ias2 = idxas(ia2,is2)
    
    ! MT 1
    phsat = 2.0d0*pi*(kqset%vql(1,iq)*atposl(1,ia1,is1)+ &
                      kqset%vql(2,iq)*atposl(2,ia1,is1)+ &
                      kqset%vql(3,iq)*atposl(3,ia1,is1))
    eph = cmplx(cos(phsat),sin(phsat),8)

    allocate(zylm((maxbigl+1)*(maxbigl+1)))
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
          wfmb(imix,ir) = wfmb(imix,ir)+eph*cmplx(umix(ir,irm,ias1)*t1,0.0d0,8)*zylm(lm)
        end do
      end do
    end do

    ! IR
    call diagsgi(iq)
    do igq = 1, Gqset%ngk(1,iq)
      imix = locmatsiz+igq
      do jgq = 1, Gqset%ngk(1,iq)
        phsat = 2.0d0*pi*(Gqset%vgkl(1,jgq,1,iq)*atposl(1,ia1,is1)+ &
        &                 Gqset%vgkl(2,jgq,1,iq)*atposl(2,ia1,is1)+ &
        &                 Gqset%vgkl(3,jgq,1,iq)*atposl(3,ia1,is1))
        do ir = nrmt(is1)+1, nrmt(is1)+rpath%nptir-1
          ri(1:3) = rpath%rvec(1:3)/rpath%rlen*rpath%r(ir,1)
          phs = phsat+(Gqset%vgkc(1,jgq,1,iq)*ri(1)+ &
          &            Gqset%vgkc(2,jgq,1,iq)*ri(2)+ &
          &            Gqset%vgkc(3,jgq,1,iq)*ri(3))
          eph = cmplx(cos(phs),sin(phs),8)/sqrt(omega)
          wfmb(imix,ir) = wfmb(imix,ir)+eph*conjg(sgi(jgq,igq))
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
          wfmb(imix,kr) = wfmb(imix,kr)+eph*cmplx(umix(jr,irm,ias2)*t1,0.0d0,8)*zylm(lm)
        end do
      end do
    end do
    
if (.false.) then
    if (iq==1) then
    do imix = 1, matsiz
      do ir = 1, npt
        write(87,'(i8,f18.6)') ir, dble(wfmb(imix,ir))
        write(88,'(i8,f18.6)') ir, imag(wfmb(imix,ir))
      end do
      write(87,*)
      write(88,*)
    end do
    end if
end if
    
    return
end subroutine


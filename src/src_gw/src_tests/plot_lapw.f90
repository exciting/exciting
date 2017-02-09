!BOP
!
! !ROUTINE: plotlapw
!
! !INTERFACE:
      subroutine plot_lapw()

! !DESCRIPTION:
!
! This subroutine calculates the real space representation of the
! (L)APW+lo basis functions in the line joining the two given atoms 
! for ploting
!
! !USES:

      use modinput
      use modmain
      use modgw
      use mod_rpath

      implicit none

! !LOCAL VARIABLES:
      integer(4) :: igp ! The G vector of the function to be ploted
      integer(4) :: i, j
      integer(4) :: ia, ia1, ia2
      integer(4) :: ir, jr, kr
      integer(4) :: is, is1, is2
      integer(4) :: np
      integer(4) :: ikp
      
      character(len=64) :: filename

      real(8) :: ri(3)
      real(8) :: kgvec(3), kgvecl(3)
      real(8) :: phs, phsat
      
      complex(8) :: pw
      complex(8), allocatable :: yl(:)
      complex(8), allocatable :: apwalm(:,:,:,:)
      complex(8), allocatable :: evecmtlm(:,:)
      complex(8), allocatable :: evecmt(:)
      complex(8), allocatable :: evecfv(:)
      complex(8), allocatable :: evec(:)      

! 
! !REVISION HISTORY:
!
! Created: 9th. July 2004 by RGA
! Last modified: 16th. Sept 2005 by RGA
! Revisited: 02.05.2011 by DIN
!
!EOP
!BOC

      call boxmsg(6,'-','PLOTLAPW')

      ikp = kset%ik2ikp(input%gw%iik)

      write(*,*) 'Parameters:'
      write(*,*) 'k-point number (iik): ', input%gw%iik
      write(*,*) 'irreducible k-point number (ikp): ', ikp
      write(*,*) 'lower bound for k+g-vector number (igmin): ', input%gw%igmin
      write(*,*) 'upper bound for k+g-vector number (igmax): ', input%gw%igmax
      write(*,*) 'atom 1 (at1): ', input%gw%at1
      write(*,*) 'atom 2 (at2): ', input%gw%at2
      write(*,*)

      if ((input%gw%at1 < 1).or.(input%gw%at1>natmtot)) stop 'atom1 is wrong'
      if ((input%gw%at2 < 1).or.(input%gw%at2>natmtot)) stop 'atom2 is wrong'

      ! Generate real-space path
      call init_rpath(rpath,input%gw%at1,input%gw%at2)
      np = rpath%nptot
      
      is1 = rpath%atom1(1)
      ia1 = rpath%atom1(2)
    
      is2 = rpath%atom2(1)
      ia2 = rpath%atom2(2)

      ! Calculate LAPW Basis Functions
      do igp = input%gw%igmin, input%gw%igmax

        write(filename,'("lapw-",i4,"-",i4,"-",i4,"-",i4,".out")') &
        &     ikp, igp, input%gw%at1, input%gw%at2
        call str_strip(filename)
        open(unit=71,file=filename,status='unknown')
        
        allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
        allocate(yl(lmmaxapw))
        allocate(evecmtlm(lmmaxapw,nrmtmax))
        allocate(evecmt(nrmtmax))
        allocate(evecfv(nmatmax))
        allocate(evec(np))

        ! find the matching coefficients
        call match(ngk(1,ikp),gkc(:,1,ikp),tpgkc(:,:,1,ikp), &
        &          sfacgk(:,:,1,ikp),apwalm)

        ! to calculate just a single basis function
        evecfv(:) = 0.0d0
        evecfv(igp) = 1.0d0

        ! calculate the values of the spherical harmonics for atom 1
        call ylm(rpath%rvec,input%groundstate%lmaxapw,yl)

        ! calculate the radial wavefunctions of atom 1
        call wavefmt(1,input%groundstate%lmaxapw, &
        &            is1,ia1,ngk(1,ikp),apwalm,evecfv,lmmaxapw,evecmtlm)

        ! convert from spherical harmonics to spherical coordinates
        call zgemv('T',lmmaxapw,nrcmt(is1),zone,evecmtlm,lmmaxapw,yl,1, &
        &          zzero,evecmt,1)
        do ir = 1, nrcmt(is1)
           evec(ir) = evecmt(ir)
        enddo

        ! Calculate the phase of the plane waves due to the change of origin
        kgvec(:) = vgkc(:,igp,1,ikp)
        kgvecl(:) = vgkl(:,igp,1,ikp)
       
        phsat = 2.0d0*pi*(kgvecl(1)*atposl(1,ia1,is1)+ &
        &                 kgvecl(2)*atposl(2,ia1,is1)+ &
        &                 kgvecl(3)*atposl(3,ia1,is1))
        
        do ir = nrmt(is1)+1, nrmt(is1)+rpath%nptir-1
          ri(1:3) = rpath%rvec(1:3)/rpath%rlen*rpath%r(ir,1)
          phs = phsat+(kgvec(1)*ri(1)+kgvec(2)*ri(2)+kgvec(3)*ri(3))
          pw = cmplx(dcos(phs),dsin(phs),8)/sqrt(omega)
          evec(ir) = pw
        end do 
     
        ! calculate the radial wavefunctions of atom 2 (if needed)
        if ((is1==is2).and.(ia1==ia2)) then
          do ir = 1, nrmt(is2)
            jr = nrmt(is2)-ir+1
            kr = nrmt(is1)+rpath%nptir+ir-1
            evec(kr) = evecmt(jr)
          end do
        else
          call ylm(-1.d0*rpath%rvec,input%groundstate%lmaxapw,yl)
          call wavefmt(1,input%groundstate%lmaxapw,  &
          &            is2,ia2,ngk(1,ikp),apwalm,evecfv,lmmaxapw,evecmtlm)
          call zgemv('t',lmmaxapw,nrcmt(is2),zone,evecmtlm,lmmaxapw,yl,1, &
          &          zzero,evecmt,1)
          do ir = 1, nrcmt(is2)
            jr = nrmt(is2)-ir+1
            kr = nrmt(is1)+rpath%nptir+ir-1
            evec(kr) = evecmt(jr)
          end do
        end if

        ! write to file
        do i = 1, np
          write(71,*) rpath%r(i,1), real(evec(i)), aimag(evec(i))
        end do
        close(71)
        
        deallocate(apwalm)
        deallocate(yl)
        deallocate(evecmtlm)
        deallocate(evecmt)
        deallocate(evecfv)
        deallocate(evec)        

      enddo ! igp    
      
      return
end subroutine
!EOC

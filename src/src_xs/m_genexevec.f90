module m_genexevec
  use modinput, only: input
  use modmpi
  use modbse
  use modscl
  use m_dzmatmult
  use mod_constants
  use m_writecmplxparts

  implicit none

  contains

    subroutine genexevec(i1, i2, nex, cmat, cpmat, auxvec, evals, rvecp, avecp)

      integer(4), intent(in) :: i1, i2, nex
      complex(8), intent(in) :: cmat(:,:), auxvec(:,:)
      complex(8), allocatable, intent(inout) :: cpmat(:,:)
      real(8), intent(in) :: evals(:)
      complex(8), allocatable, intent(out) :: rvecp(:,:), avecp(:,:)

      integer(4) :: nreq, m, i

      ! Reqested number of solutions
      nreq = i2-i1+1
      if(nreq <= 0 .or. nreq > nex .or. i2<1 .or. i1<1) then 
        write(*,*) "Error(genexevec): Index mismatch."
        write(*,*) "i1,i2,nex", i1,i2,nex
        call terminate
      end if

      ! Check input
      m = size(cmat,1)
      if(m/=size(cmat,2) .or. m /= size(cpmat,1) .or. m /= size(auxvec,1)&
        & .or. m /= size(cpmat, 2) .or. nreq > size(auxvec,2)) then
        write(*,*) "Error(genexevec): Size mismatch."
        write(*,*) "shape(cmat)=",shape(cmat)
        write(*,*) "shape(cpmat)=",shape(cpmat)
        write(*,*) "shape(auxvec)=",shape(auxvec)
        call terminate
      end if

      ! Y
      ! Make (A-B)^{-1/2} |E_\lambda|^{1/2} Z_\lambda
      allocate(avecp(m,nreq))
      !  Aux2 = (A-B)^{-1/2} Z -> Y
      call zgemm('N','N', m, nreq, m, zone, cpmat, m,&
        & auxvec(:,i1:i2), m, zzero, avecp, m)
      !  Aux2_{m,\lambda} = Aux2_{m,\lambda} * |E_\lambda|^{1/2}
      do i=1, nreq
        avecp(:,i) = sqrt(evals(i+i1-1)) * avecp(:,i)
      end do
      deallocate(cpmat)

      ! X
      ! Make (A-B)^{1/2} |E_\lambda|^-{1/2} Z_\lambda
      allocate(rvecp(m,nreq))
      !  Aux1 = (A-B)^{1/2} Z -> X
      call zgemm('N','N', m, nreq, m, zone, cmat, m,&
        & auxvec(:,i1:i2), m, zzero, rvecp, m)
      !  Aux1_{m,\lambda} = Aux1_{m,\lambda} * |E_\lambda|^-{1/2}
      do i=1, nreq
        rvecp(:,i) = 1.0d0/sqrt(evals(i+i1-1)) * rvecp(:,i)
      end do

      ! X^+ = 1/2 (aux1+aux2)
      rvecp(1:m,1:nreq) = 0.5d0 *(rvecp+avecp)
      ! Y^+ = 1/2 (aux1-aux2)
      avecp(1:m,1:nreq) = rvecp - avecp

    end subroutine genexevec

    subroutine gendexevec(i1, i2, nex, dcmat, dcpmat, dauxvec, evals,&
       & drvecp, davecp)

      ! I/O
      integer(4), intent(in) :: i1, i2, nex
      type(dzmat), intent(in) :: dcmat, dauxvec
      type(dzmat), intent(inout) :: dcpmat
      real(8), intent(in) :: evals(:)
      type(dzmat), intent(inout) :: drvecp, davecp 

      ! Local
      integer(4) :: nreq, m, i, context, ig

      character(*), parameter :: thisname = "gendexevec"

      ! Reqested number of solutions
      nreq = i2-i1+1
      if(nreq <= 0 .or. nreq > nex .or. i2<1 .or. i1<1) then 
        if(mpiglobal%rank==0) then 
          write(*,*) "Error(gendexevec): Index mismatch."
          write(*,*) "i1,i2,nex", i1,i2,nex
        end if
        call terminate
      end if

      ! Check contexts
      context = dcmat%context
      if(dcpmat%context /= context .or. dauxvec%context /= context) then
        if(mpiglobal%rank == 0) then 
          write(*,*) "Error(gendexevec): Context mismatch"
          write(*,*) "context of dcmat, dcpmat, dauxvec :"
          write(*,*) dcmat%context, dcpmat%context, dauxvec%context
        end if
        call terminate
      end if
      ! Check matrix sized
      m = dcmat%nrows
      if(m/=dcmat%ncols .or. m/=dcpmat%nrows .or. m/=dauxvec%nrows &
        & .or. m/=dcpmat%ncols .or. nreq > dauxvec%ncols) then
        write(*,*) "Error(gendexevec): Size mismatch."
        write(*,*) "shape(dcmat)=", dcmat%nrows, dcmat%ncols
        write(*,*) "shape(dcpmat)=", dcpmat%nrows, dcpmat%ncols
        write(*,*) "shape(dauxvec)=", dauxvec%nrows, dauxvec%ncols
        call terminate
      end if

      ! Y
      call new_dzmat(davecp, m, nreq, bi2d)
      ! Make (A-B)^{-1/2} |E_\lambda|^{1/2} Z_\lambda
      !  Aux2 = (A-B)^{-1/2} Z -> Y
      call dzmatmult(dcpmat, dauxvec, davecp, n=nreq, jb=i1)
      ! C' no longer needed 
      call del_dzmat(dcpmat)
      !  Aux2_{m,\lambda} = Aux2_{m,\lambda} * |E_\lambda|^{1/2} -> Y
      do i=1, davecp%ncols_loc
        if(davecp%isdistributed) then
          ig = davecp%c2g(i)
        else
          ig = i
        end if
        davecp%za(:,i) = sqrt(evals(ig+i1-1)) * davecp%za(:,i)
      end do
      
      ! X
      call new_dzmat(drvecp, m, nreq, bi2d)
      ! Make (A-B)^{1/2} |E_\lambda|^-{1/2} Z_\lambda
      !  Aux1 = (A-B)^{1/2} Z -> X
      call dzmatmult(dcmat, dauxvec, drvecp, n=nreq, jb=i1)
      !  Aux1_{m,\lambda} = Aux1_{m,\lambda} * |E_\lambda|^-{1/2} -> X
      ! Loop over local column index
      do i=1, drvecp%ncols_loc
        ! Get global column index
        if(drvecp%isdistributed) then 
          ig = drvecp%c2g(i)
        else
          ig = i 
        end if
        drvecp%za(:,i) = 1.0d0/sqrt(evals(ig+i1-1)) * drvecp%za(:,i)
      end do

      ! X^+ = 1/2 (aux1+aux2) 
      drvecp%za(:,:) = 0.5d0 *(drvecp%za(:,:)+davecp%za(:,:))
      ! Y^+ = 1/2 (aux1-aux2) = X^+ - aux2
      davecp%za(:,:) = drvecp%za(:,:)-davecp%za(:,:)

    end subroutine gendexevec
    
end module m_genexevec

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
      complex(8), intent(in) :: cmat(:,:), cpmat(:,:), auxvec(:,:)
      real(8), intent(in) :: evals(:)
      complex(8), intent(out), optional :: rvecp(:,:), avecp(:,:)

      complex(8), allocatable :: aux1(:,:), aux2(:,:)
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

      ! Make (A-B)^{1/2} |E_\lambda|^-{1/2} Z_\lambda
      allocate(aux1(m,nreq))
      !  Aux1 = (A-B)^{1/2} Z
      call zgemm('N','N', m, nreq, m, zone, cmat, m,&
        & auxvec(:,i1:i2), m, zzero, aux1, m)
      !  Aux1_{m,\lambda} = Aux1_{m,\lambda} * |E_\lambda|^-{1/2}
      do i=1, nreq
        aux1(:,i) = 1.0d0/sqrt(evals(i+i1-1)) * Aux1(:,i)
      end do

      ! Make (A-B)^{-1/2} |E_\lambda|^{1/2} Z_\lambda
      allocate(aux2(m,nreq))
      !  Aux2 = (A-B)^{-1/2} Z
      call zgemm('N','N', m, nreq, m, zone, cpmat, m,&
        & auxvec(:,i1:i2), m, zzero, aux2, m)
      !  Aux2_{m,\lambda} = Aux2_{m,\lambda} * |E_\lambda|^{1/2}
      do i=1, nreq
        aux2(:,i) = sqrt(evals(i+i1-1)) * aux2(:,i)
      end do

      if(present(rvecp)) then 
        if(size(rvecp,1) < m .or. size(rvecp,2) < nreq) then 
          write(*,*) "Error(genexevec): Size mismatch. shape(rvecp)=", shape(rvecp)
          call terminate
        end if
        ! X^+ = 1/2 (aux1+aux2)
        rvecp(1:m,1:nreq) = 0.5d0 *(aux1+aux2)
      end if
      if(present(avecp)) then 
        if(size(avecp,1) < m .or. size(avecp,2) < nreq) then 
          write(*,*) "Error(genexevec): Size mismatch. shape(avecp)=", shape(avecp)
          call terminate
        end if
        ! Y^+ = 1/2 (aux1-aux2)
        avecp(1:m,1:nreq) = 0.5d0 *(aux1-aux2)
      end if

      deallocate(aux1,aux2)

    end subroutine genexevec

    subroutine gendexevec(i1, i2, nex, dcmat, dcpmat, dauxvec, evals,&
       & drvecp, davecp)

      ! I/O
      integer(4), intent(in) :: i1, i2, nex
      type(dzmat), intent(in) :: dcmat, dcpmat, dauxvec
      real(8), intent(in) :: evals(:)
      type(dzmat), intent(inout), optional :: drvecp, davecp 

      
      ! Local
      type(dzmat) :: daux1, daux2
      integer(4) :: nreq, m, i, context, ig

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
      ! Check output matrices
      if(present(drvecp)) then 
        if(drvecp%context /= context) then 
          if(mpiglobal%rank == 0) then 
            write(*,*) "Error(gendexevec): Context mismatch"
            write(*,*) "context of dcmat, drvecp :"
            write(*,*) dcmat%context, drvecp%context
          end if
          call terminate
        end if
        if(drvecp%nrows /= m .or. drvecp%ncols /= nreq) then
          if(mpiglobal%rank == 0) then 
            write(*,*) "Error(gendexevec): Size mismatch"
            write(*,*) "shape(dcmat)=", dcmat%nrows, dcmat%ncols
            write(*,*) "shape(drvecp)=", drvecp%nrows, drvecp%ncols
          end if
          call terminate
        end if
      end if
      if(present(davecp)) then 
        if(davecp%context /= context) then 
          if(mpiglobal%rank == 0) then 
            write(*,*) "Error(gendexevec): Context mismatch"
            write(*,*) "context of dcmat, davecp :"
            write(*,*) dcmat%context, davecp%context
          end if
          call terminate
        end if
        if(davecp%nrows /= m .or. davecp%ncols /= nreq) then
          if(mpiglobal%rank == 0) then 
            write(*,*) "Error(gendexevec): Size mismatch"
            write(*,*) "shape(dcmat)=", dcmat%nrows, dcmat%ncols
            write(*,*) "shape(davecp)=", davecp%nrows, davecp%ncols
          end if
          call terminate
        end if
      end if
      
      ! Make (A-B)^{1/2} |E_\lambda|^-{1/2} Z_\lambda
      call new_dzmat(daux1, m, nreq, bi2d)
      !  Aux1 = (A-B)^{1/2} Z
      call dzmatmult(dcmat, dauxvec, daux1, n=nreq, jb=i1)
      !  Aux1_{m,\lambda} = Aux1_{m,\lambda} * |E_\lambda|^-{1/2}
      ! Loop over local column index
      do i=1, daux1%ncols_loc
        ! Get global column index
        ig = daux1%c2g(i)
        daux1%za(:,i) = 1.0d0/sqrt(evals(ig+i1-1)) * daux1%za(:,i)
      end do

      ! Make (A-B)^{-1/2} |E_\lambda|^{1/2} Z_\lambda
      call new_dzmat(daux2, m, nreq, bi2d)
      !  Aux2 = (A-B)^{-1/2} Z
      call dzmatmult(dcpmat, dauxvec, daux2, n=nreq, jb=i1)
      !  Aux2_{m,\lambda} = Aux2_{m,\lambda} * |E_\lambda|^{1/2}
      do i=1, daux2%ncols_loc
        ig = daux2%c2g(i)
        daux2%za(:,i) = sqrt(evals(ig+i1-1)) * daux2%za(:,i)
      end do

      if(present(drvecp)) then 
        ! X^+ = 1/2 (aux1+aux2)
        drvecp%za(:,:) = 0.5d0 *(daux1%za(:,:)+daux2%za(:,:))
      end if
      if(present(davecp)) then 
        ! Y^+ = 1/2 (aux1-aux2)
        davecp%za(:,:) = 0.5d0 *(daux1%za(:,:)-daux2%za(:,:))
      end if

      call del_dzmat(daux1)
      call del_dzmat(daux2)

    end subroutine gendexevec
    
end module m_genexevec

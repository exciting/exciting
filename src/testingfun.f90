subroutine testingfun
  use modinput
  use modmpi
  use modscl
  use m_dzmatmult

  implicit none

  type(dzmat) :: amat, bmat, cmat
  integer(4) :: matsize, iloop, jloop, nstart, rb, cb, bl
  real(8) :: t0, t1

  write(*,*) "Hello, testing fun here!"


  !   Make square'ish process grid (context 0)
  call setupblacs(mpiglobal, 'grid', bi2d)

  if(mpiglobal%rank == 0) then 
    call printblacsinfo(bi2d)
  end if
  call barrier(mpiglobal)

  nstart = 2**10
  bl = 8


  do jloop=0, 8

    rb = bl*2**jloop
    cb = rb
    
    if(mpiglobal%rank == 0) then 
        write(*,*) "===="
        write(*,*) "rb x cb = ", rb, cb
        write(*,*) "===="
    end if
    
    do iloop=0, 7

      matsize = nstart*2**iloop

      ! Make matrices
      call new_dzmat(amat, matsize, matsize, bi2d, rblck=rb, cblck=cb, fzero=.false.)
      call new_dzmat(bmat, matsize, matsize, bi2d, rblck=rb, cblck=cb, fzero=.false.)
      call new_dzmat(cmat, matsize, matsize, bi2d, rblck=rb, cblck=cb, fzero=.false.)

      if(mpiglobal%rank == 0) then 
        write(*,*) "matsize = ", matsize
        write(*,*) "gb = ", dble(matsize)**2*8/1024**3
        write(*,*) "3*gb = ", dble(matsize)**2*3*8/1024**3
        write(*, '("  local matrix shape ",i6," x",i6)')&
          & amat%nrows_loc, amat%ncols_loc
        write(*,*) "local gb = ", dble(amat%nrows_loc)**2*8/1024**3
        write(*,*) "local 3*gb = ", dble(amat%nrows_loc)**2*3*8/1024**3
        write(*,*) "----"
      end if
      call barrier(mpiglobal)

      ! Multiply
      call timesec(t0)
      call dzmatmult(amat, bmat, cmat)
      call timesec(t1)

      if(mpiglobal%rank == 0) then 
        write(*, '("Timing (in seconds)	   :", f12.3)') t1 - t0
      end if
      
      ! Clear matrices 
      call del_dzmat(amat)
      call del_dzmat(bmat)
      call del_dzmat(cmat)

    end do

  end do

  call exitblacs(bi2d)

  if(mpiglobal%rank == 0) then 
    write(*,*) "All good."
  end if

  call terminate

end subroutine testingfun

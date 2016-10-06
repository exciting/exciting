module m_writecmplxparts

  implicit none

  contains

    subroutine writecmplxparts(fbasename, remat, immat, ik1, ik2)
      use m_getunit
      character(*), intent(in) :: fbasename
      real(8), intent(in) :: remat(:,:)
      real(8), intent(in), optional :: immat(:,:)
      integer(4), intent(in), optional :: ik1, ik2

      integer(4) :: un, a1, a2, n, m
      character(256) :: fname, tmp1, tmp2, tmp3

      n = size(remat,1)
      m = size(remat,2)

      tmp1 = ''
      if(present(ik1)) then
        write(tmp1, '(I8)') ik1
      end if
      tmp2 = ''
      if(present(ik2)) then
        write(tmp2, '(I8)') ik2
      end if
      tmp3= trim(adjustl(tmp1))//trim(adjustl(tmp2))
      if(present(ik1) .or. present(ik2)) then
        tmp3 = "_"//trim(adjustl(tmp3))
      end if

      fname =''
      write(fname,'("Re_",a,a,".OUT")') trim(adjustl(fbasename)),trim(adjustl(tmp3))

      call getunit(un)
      open(unit=un, file=fname, action='write', status='replace')
      do a1=1,n
        write(un,'(SP,E23.16)', advance='no') remat(a1,1)
        do a2=2,m
          write(un, '(SP,1x,E23.16)', advance='no') remat(a1,a2)
        end do
        write(un,*)
      end do
      close(un)

      if(present(immat)) then
        fname =''
        write(fname,'("Im_",a,a,".OUT")') trim(adjustl(fbasename)),trim(adjustl(tmp3))

        call getunit(un)
        open(unit=un, file=fname, action='write', status='replace')
        do a1=1,n
          write(un,'(SP,E23.16)', advance='no') immat(a1,1)
          do a2=2,m
            write(un, '(SP,1x,E23.16)', advance='no') immat(a1,a2)
          end do
          write(un,*)
        end do
        close(un)

      end if

    end subroutine writecmplxparts

end module m_writecmplxparts

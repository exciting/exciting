module m_writecmplxparts

  implicit none

  contains

    subroutine writecmplxparts(fbasename, remat, immat, ik1, ik2, revec, imvec, veclen, dirname)
      use m_getunit
      character(*), intent(in) :: fbasename
      real(8), intent(in), optional :: remat(:,:), immat(:,:)
      real(8), intent(in), optional :: revec(*), imvec(*)
      integer(4), intent(in), optional :: ik1, ik2, veclen
      character(*), intent(in), optional :: dirname

      integer(4) :: un, a1, a2, n, m
      character(256) :: fname, tmp1, tmp2, tmp3, frmt, frmtnoa, dname, syscommand

      if(present(dirname)) then 
        dname = trim(adjustl(dirname))//'/'
        syscommand = 'test ! -e '//trim(adjustl(dname))//' && mkdir '//trim(adjustl(dname))
        call system(trim(adjustl(syscommand)))
      else
        dname = ''
      end if

      frmt = '(SP,E23.16)'
      frmtnoa = '(SP,1x,E23.16)'

      if(present(remat)) then 
        n = size(remat,1)
        m = size(remat,2)
      end if

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
      fname = trim(adjustl(dname))//trim(adjustl(fname))

      if(present(remat)) then 
        call getunit(un)
        open(unit=un, file=fname, action='write', status='replace')
        do a1=1,n
          write(un, fmt=frmt, advance='no') remat(a1,1)
          do a2=2,m
            write(un, fmt=frmtnoa, advance='no') remat(a1,a2)
          end do
          write(un,*)
        end do
        close(un)
      else if(present(revec) .and. present(veclen)) then
        call getunit(un)
        open(unit=un, file=fname, action='write', status='replace')
        do a1=1,veclen
          write(un, fmt=frmt) revec(a1)
        end do
        close(un)
      end if

      if(present(immat)) then
        call getunit(un)
        fname =''
        write(fname,'("Im_",a,a,".OUT")') trim(adjustl(fbasename)),trim(adjustl(tmp3))
        fname = trim(adjustl(dname))//trim(adjustl(fname))

        open(unit=un, file=fname, action='write', status='replace')
        do a1=1,n
          write(un,fmt=frmt, advance='no') immat(a1,1)
          do a2=2,m
            write(un, fmt=frmtnoa, advance='no') immat(a1,a2)
          end do
          write(un,*)
        end do
        close(un)
      else if(present(imvec) .and. present(veclen)) then
        call getunit(un)
        fname =''
        write(fname,'("Im_",a,a,".OUT")') trim(adjustl(fbasename)),trim(adjustl(tmp3))
        fname = trim(adjustl(dname))//trim(adjustl(fname))
        open(unit=un, file=fname, action='write', status='replace')
        do a1=1,veclen
          write(un,fmt=frmt) imvec(a1)
        end do
        close(un)
      end if

    end subroutine writecmplxparts

end module m_writecmplxparts

Subroutine writebevec
      Use modmain
      Use modmpi
      Use modxs
      Use m_getunit
      Implicit None

      Real (8), Allocatable :: beval (:)
      Complex (8), Allocatable :: bevec (:, :)
      Real (8) :: bevec_ksum, bevec1
      Integer :: ex_min, ex_max, ex_min2, ex_max2, iv, ic, iknr
      integer ::  nrnst1, nrnst3 
      integer :: s1, s2, hamsiz, un, iostat
      integer :: ist, jst
      Logical :: exist
      Character (256) :: lambda

      call init0
      call init1

      if (.not.associated(input%xs%writeexcitons)) &
      &  input%xs%writeexcitons => getstructwriteexcitons (emptynode)

       ex_min2 = input%xs%writeexcitons%MinNumberExcitons
       ex_max2 = input%xs%writeexcitons%MaxNumberExcitons

      !read EXCCOEFF.bin where BSE eigenvectors are stored
       inquire(file='EXCCOEFF.bin', exist=exist)
       if ( (.not. exist) .and. (rank==0) ) then
         write(*,*)
         write(*,'("Error(bse): file EXCCOEFF.bin does not exist!")')
         write(*,*)
         stop
       else
         Call getunit (un)
         open(un,File='EXCCOEFF.bin', Action='READ',Form='UNFORMATTED', IOstat=iostat)
         if ( (iostat.ne.0) .and. (rank==0) ) then
           write(*,*) iostat
           write(*,'("Error(bse): error reading EXCCOEFF.bin")')
           write(*,*)
           stop
         end if
         read(un) ex_min, ex_max, nkptnr, istl3, sta1, sta2, nrnst1, nrnst3, hamsiz 

         Allocate (beval(hamsiz), bevec(hamsiz,ex_min:ex_max))

         do s1 = ex_min, ex_max 
            read(un) beval(s1), bevec(1:hamsiz,s1)
         end do
         Close(un)
       end if


       if (  (ex_min2 .lt. ex_min) .or. &
          &  (ex_min2 .gt. ex_max) .or. &
          &  (ex_max2 .lt. ex_min) .or. &
          &  (ex_max2 .gt. ex_max) .or. &
          &  (ex_min2 .gt.  ex_max2) ) then
          write(*,*)
          write(*,'("Error(bse): wrong range of exciton indices: ", 2I5)') &
                & ex_min2, ex_max2
          write(*,'("Error(bse): range of exciton indices must be within the stored range of EXCCOEFF.bin: ", 2I5)') &
                & ex_min, ex_max 
          write(*,*)
          stop
        end if


         !write bin

!        Write(*,*) 'writebevec'
! 
!        Write (*,'("ex_min, ex_max, nkptnr, nsta1, nrnst1, nrnst3, hamsiz", 7I8 )') ex_min, ex_max, nkptnr, nsta1, nrnst1, nrnst3, hamsiz
!        do s1 = 1, ex_max-ex_min+1
!           write(*,'("s1, beval(s1), bevec(1:5,s1) ", I8, 33F14.8)') s1, beval(s1), bevec(1:16,s1)
!           write(*,*)
!        end do


      !write ASCII output of Abs(BSE eigenvector)^2 
      Do s1 = ex_min2, ex_max2
      Call getunit (un)
      Write (lambda, '("_LAMBDA",i4.4)') s1
      Open (Unit=un, File=trim('BEVEC'//trim(lambda)//'.OUT'), Form='formatted', Action='write')
            Do iknr = 1, nkptnr
               Do iv = 1, nrnst1
                  Do ic = 1, nrnst3
                     s2 = hamidx (iv, ic, iknr, nrnst1, nrnst3)
                     bevec1 = bevec(s2,s1) * conjg(bevec(s2,s1))
                     Write (un, '(I8, 3G18.10, 2I6, G18.10)') iknr, vkl(:, iknr), iv+sta1-1, ic , bevec1
                  End Do
               End Do
            End Do
      Write (un,*)
      Write (un, '("k-point, k-point coordinates, valence band, conduction band, Abs(BSE eigenvector)^2 ")')
      Write (un,*)
      Write (un, '(I6, " : Nr. k-points")') nkptnr
      Write (un, '(I6, " : VBM")') nrnst1+sta1-1
      Write (un, '(I6, " : Nr. valence states")') nrnst1
      Write (un, '(I6, " : Nr. conduction states")') nrnst3
      Write (un,*)
      Close (un)
      End Do


      !write ASCII output of sum of Abs(BSE eigenvector)^2 over all k-points     
      Do s1 = ex_min2, ex_max2
      Call getunit (un)
      Write (lambda, '("_LAMBDA",i4.4)') s1
      Open (Unit=un, File=trim('BEVEC_KSUM'//trim(lambda)//'.OUT'), Form='formatted', Action='write')
         Do iv = 1, nrnst1
             Do ic = 1, nrnst3
                bevec_ksum=0
                Do iknr = 1, nkptnr
                     s2 = hamidx (iv, ic, iknr, nrnst1, nrnst3)
                     bevec_ksum= bevec_ksum + bevec(s2,s1) * conjg(bevec(s2,s1))
                End Do
                Write (un, '(2I8, G18.10)') iv+sta1-1, ic, bevec_ksum
               End Do
            End Do
      Write (un,*)
      Write (un, '("valence band, conduction band, sum of Abs(BSE eigenvector)^2 over all k-points")')
      Write (un,*)
      Write (un, '(I6, " : VBM")') nrnst1+sta1-1
      Write (un, '(I6, " : Nr. valence states")') nrnst1
      Write (un, '(I6, " : Nr. conduction states")') nrnst3
      Write (un,*)
      Close (un)
      End Do

!---------------------------------------------------------------------------
! din: New output file for the bandstructure to be able to post-process it
!---------------------------------------------------------------------------
      do s1 = ex_min2, ex_max2
        call getunit(un)
        write(lambda, '("exciton_evec_",i4.4,".dat")') s1
        open(Unit=un, File=trim(lambda), Form='Formatted', Action='Write')
        ! nkpt total, Nv, iv0, Nc, ic0
        write(un,*) "# ", nkptnr, &
        &                 nrnst1, sta1,  &
        &                 nrnst3, istl3
        do iknr = 1, nkptnr
          do iv = 1, nrnst1
            ist = iv+sta1-1
            do ic = 1, nrnst3
              jst = ic+istl3-1
              s2 = hamidx(iv, ic, iknr, nrnst1, nrnst3)
              write(un, '(I8, 4X, 3F10.6, 2I8, 4X, G18.10)') iknr, vkl(:,iknr), &
              &    ist, jst, dble(bevec(s2,s1)*conjg(bevec(s2,s1)))
            end do
          end do
        end do
        close(un)
      end do ! s1

      Deallocate(beval,bevec)

Contains
!
      Integer Function hamidx (i1, i2, ik, n1, n2)
         Implicit None
         Integer, Intent (In) :: i1, i2, ik, n1, n2
         hamidx = i2 + n2 * (i1-1) + n1 * n2 * (ik-1)
      End Function hamidx

End Subroutine writebevec

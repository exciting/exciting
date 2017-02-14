!
!-----------------------------------------------------------
!
      subroutine flushbuf(fid)
        implicit none
        integer :: fid
        if(fid.eq.6) return 
!If in AIX Fortran2003, flush is a statement 
#ifdef AIX
      flush(fid)
#else
      call flush(fid)
#endif
      end subroutine
!      
!-----------------------------------------------------------
!
      subroutine write_cputime(fid,cput,info)
        integer :: fid
        real(8),intent(in) :: cput
        character(len=*)   :: info
        character(len=30)  :: string
        integer :: i
        
        do i=1,30
           string(i:i)=' '
        end do
        string=info//string
        
        write(fid,10) ADJUSTL(string),cput
    10  format(" CPU Time for    ",A25,F10.2)
        call flushbuf(fid)
      
      end subroutine 
!      
!-----------------------------------------------------------
!
      subroutine init_random_seed()
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))
      
        call system_clock(count=clock)
      
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = seed)
      
        deallocate(seed)
      end subroutine
!      
!-----------------------------------------------------------
!
      subroutine SetName(fn,pref,N)
! Set a file name with the form "<pref>-<N>.<suffix>"
        character(len=*) :: fn,pref
        integer(4) :: N
        integer(4) :: w,fWidth,l
        character  :: ID*10
        external fWidth

        w=fWidth(n)
        
        write(ID,'(I10)') N
        l=LEN_trim(pref)
        if(pref(l:l) .EQ. '-') THEN
           fn=trim(pref)//ID(10-w+1:10)
        else 
           fn=trim(pref)//"-"//ID(10-w+1:10)
        endif
      end subroutine
!      
!-----------------------------------------------------------
!
      function Int2Str(inum)
        implicit none
        integer :: inum
        character(len=10) :: Int2Str
        integer :: w,fWidth
        character :: ID*10
        external fWidth

        Int2Str=' ';
        w=fWidth(inum)
        write(ID,'(I10)') inum
        Int2Str(1:w)=ID(10-w+1:10)
      end function      
!      
!-----------------------------------------------------------
!
      subroutine str2int(str,num,ierr)
        implicit none
        integer:: num,ierr
        character(*)::str

        integer::slen,id,ndig,n
        character::c

        slen=len(str) 
        ndig=0
        num=0
        ierr=0
        do id=slen,1,-1
          c=str(id:id)
          if(c.ne.'') then 
            call char2int(c,n)
            if(n.eq.-1) then 
              ierr=1
              exit 
            else 
              num=num+n*10**ndig
              ndig=ndig+1
            endif 
          endif 
        enddo 
        
        if(c.eq.'-') then 
          num=-num
          ierr=0
        elseif(c.eq.'+') then 
          ierr=0
        endif  

      end subroutine 
!      
!-----------------------------------------------------------
!
      subroutine char2int(c,n)
        implicit none
        integer n
        character c
        character::str(0:9)=(/'0','1','2','3','4','5','6','7','8','9'/)
        integer:: i
         
        n=-1
        do i=0,9
          if(c.eq.str(i)) then 
            n=i
            exit 
          endif 
        enddo 
      end subroutine  
!      
!-----------------------------------------------------------
!
      integer function fWidth(n)
         integer::n
         if(n==0) then
           fWidth=1; return
         endif
         if(n>0) then
            fWidth=int(log10(n*1.0))+1
         else
            fWidth=int(log10(-n*1.0))+2
         endif
      end function 
!      
!-----------------------------------------------------------
!
      subroutine SkipRD(iunit,nlines) 
      implicit none 
      integer*4 iunit,nlines,n
      
      DO N=1,nlines
         READ(iunit,*)
      ENDDO
      end subroutine 
!      
!-----------------------------------------------------------
!
      subroutine DBGMSG(IPRT,msg)
      integer*4 iprt
      character(len=*) msg

      if(iprt .GT.0 ) then
        write(6,*)  trim(msg)
      endif
      end subroutine
!      
!-----------------------------------------------------------
!
      subroutine errmsg(lerr,sname,msg)
      logical,intent(in):: lerr
      character(len=*) sname,msg

      if(lerr) then
        write(6,*) "FATAL ERROR in ",sname
        write(6,*) '----',trim(msg)
        stop 
      endif

      end subroutine
!      
!-----------------------------------------------------------
!
      subroutine errmsg0(ierr,sname,msg)
      integer,intent(in):: ierr
      character(len=*) sname,msg

      if(ierr.ne.0) then
        write(6,*) "FATAL ERROR in ",sname
        write(6,*) '----',trim(msg),"ierr=",ierr
        stop
      endif

      end subroutine
!      
!-----------------------------------------------------------
!
      subroutine wrnmsg(lerr,sname,msg)
      logical,intent(in):: lerr
      character(len=*) sname,msg

      if(lerr) then
        write(6,*) "WARNING in ",sname
        write(6,*) '----',trim(msg)
      endif
      end subroutine
!      
!-----------------------------------------------------------
!
      subroutine LinMSG(fid,c,ttl)
      implicit none
      integer,parameter::MaxLineLength=80
      integer::fid
      character c
      character(len=*) ttl
      character(len=MaxLineLength)::str
      integer::n1,n2,i,length

      write(fid,*) 
      length=len(ttl)
      if(length .GT. MaxLineLength) THEN
         write(fid,'(A80)') ttl
      else

         do i=1,MaxLineLength
           str(i:i)=c
         enddo
         if (length>0) then
           n1=(MaxLineLength-length)/2
           n2=n1+length-1
           str(n1:n2)=ttl(1:length)
           str(n1-1:n1-1)='';
           str(n2+1:n2+1)=''
         end if
         write(fid,'(A80)') str

      endif
      write(fid,*) 
      end subroutine
!      
!-----------------------------------------------------------
!
      subroutine BoxMSG(fid,c,ttl)
      implicit none
      integer,parameter::MaxLineLength=80
      integer:: fid
      character c
      character(len=*) ttl
      character(len=MaxLineLength)::str
      integer::n1,n2,i,length

      write(fid,*) 
      length=len(ttl)
      do i=1,MaxLineLength
        str(i:i)=c
      enddo
      write(fid,'(A80)') str

      if(length .GE. MaxLineLength) THEN
        write(fid,'(A80)') ttl
      else
        str(1:1)=c
        do i=2,MaxLineLength-1
          str(i:i)=''
        enddo
        str(i:i)=c

        n1=(MaxLineLength-length)/2
        n2=n1+length-1
        str(n1:n2)=ttl(1:length)
        str(n1-1:n1-1)='';
        str(n2+1:n2+1)=''
        write(fid,'(A80)') str
      endif
      do i=1,MaxLineLength
        str(i:i)=c
      enddo
      write(fid,'(A80)') str
      write(fid,*) 
      end subroutine
!      
!-----------------------------------------------------------
!
      subroutine str_strip( string )
      ! remove multiple spaces from a string
      implicit none
      integer         :: i, j
      character*(*)   :: string
      character(1024) :: temp
                                             
      do i=1,len(temp)
         temp(i:i)=' '
      end do

      j = 0
      do i = 1, len(string)
         if ( string(i:i) .ne. ' ' ) then
            j = j + 1
            temp(j:j) = string(i:i)             
         end if                              
      end do                                 

      do i=1,len(string)
         string(i:i) = temp(i:i)
      end do   

      end subroutine                                   
!      
!-----------------------------------------------------------
!
      subroutine printbox(fid,c,ttl)
      implicit none
      integer,parameter::MaxLineLength=80
      integer:: fid
      character c
      character(len=*) ttl
      character(len=MaxLineLength)::str
      integer::n1,n2,i,length

      write(fid,*) 
      length=len(ttl)
      do i=1,MaxLineLength
        str(i:i)=c
      enddo
      write(fid,'(A80)') str

      if(length .GE. MaxLineLength) THEN
        write(fid,'(A80)') ttl
      else
        str(1:1) = c
        str(2:2) = ' '
        do i=3,MaxLineLength-1
          str(i:i)=''
        enddo
        str(i:i)=c
        n1=3!(MaxLineLength-length)/2
        n2=n1+length-1
        str(n1:n2)=ttl(1:length)
        str(n1-1:n1-1)='';
        str(n2+1:n2+1)=''
        str(i:i)=c
        write(fid,'(A80)') str
      endif
      do i=1,MaxLineLength
        str(i:i)=c
      enddo
      write(fid,'(A80)') str
!      write(fid,*) 
      end subroutine
!      
!-----------------------------------------------------------
!
      subroutine printline(fid,c)
      implicit none
      integer,parameter::MaxLineLength=80
      integer:: fid
      character c
      character(len=MaxLineLength)::str
      integer :: i

      do i=1,MaxLineLength
         str(i:i)=c
      enddo
      write(fid,'(A80)') str

      end subroutine
!      
!-----------------------------------------------------------
!
      subroutine printtext(fid,c,ttl)
      implicit none
      integer,parameter::MaxLineLength=80
      integer:: fid
      character c
      character(len=*) ttl
      character(len=MaxLineLength)::str
      integer::n1,n2,i,length

      length=len(ttl)
      do i=1,MaxLineLength
        str(i:i)=c
      enddo

      if(length .GE. MaxLineLength) THEN
        write(fid,'(A80)') ttl
      else
        str(1:2)="| "
        do i=3,MaxLineLength-1
          str(i:i)=''
        enddo
        str(i:i)=c
        n1=3!(MaxLineLength-length)/2
        n2=n1+length-1
        str(n1:n2)=ttl(1:length)
        str(n1-1:n1-1)='';
        str(n2+1:n2+1)=''
        str(i:i)=c
        write(fid,'(A80)') str
      endif

      end subroutine

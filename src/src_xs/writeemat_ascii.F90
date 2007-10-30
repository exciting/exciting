
subroutine writeemat_ascii
  use modmain
  use modtddft
  use m_getunit
  use m_getemat
  use m_genfilname
  implicit none
  complex(8), allocatable :: emat(:,:,:)
  complex(8) :: xou,xuo
  character(16) :: f1,f2,f
  character(256) :: filnam
  integer :: un,iq,ik,iv,ic,igq
  real(8) :: vkloff_save(3)

  ! save k-point offset
  vkloff_save = vkloff

  call init0
  call init1
  call init2xs
  call getunit(un)

  ! loop over q-points
  do iq = 1, nqpt
     vkloff(:)=qvkloff(:,iq)
     ! calculate k+q and G+k+q related variables
     call init1td
     if (allocated(xiou)) deallocate(xiou)
     if (allocated(xiuo)) deallocate(xiuo)
     allocate(xiou(nstval,nstcon,ngq(iq)))
     allocate(xiuo(nstcon,nstval,ngq(iq)))
     ! filename for matrix elements file
     call genfilname(basename='EMAT',asc=.true.,iq=iq,filnam=filnam)
     open(un,file=trim(filnam),action='write')
     write(un,'(a)') 'iq,ik,iv,ic,igq,xou,xuo,|xou|^2,|xuo|^2 below'
     ! loop over k-points
     do ik=1,nkpt
        ! read matrix elements of exponential expression
        call genfilname(basename='EMAT',iq=iq,filnam=fnemat)
        call getemat(iq,ik,.true.,trim(fnemat),xiou,xiuo)
        do iv=1,nstval
           f1='v'
           !!!if (ist1.gt.(nstsv-nempty-1)) f1='c'
           do ic=1,nstcon
              f2='c'
              !!!if (ist2.gt.(nstsv-nempty-1)) f2='c'
              f='  '//trim(f1)//'-'//trim(f2)//'  '
              do igq=1,ngq(iq)
                 xou=xiou(iv,ic,igq)
                 xuo=xiuo(ic,iv,igq)
                 write(un,'(5i8,6g18.10)') iq,ik,iv,ic,igq,xou,xuo,&
                      abs(xou)**2,abs(xuo)**2
              end do
           end do
        end do
     end do ! ik
     close(un)
     deallocate(xiou,xiuo)
  end do ! iq

  ! restore offset
  vkloff = vkloff_save
  call genfilname(setfilext=.true.)

end subroutine writeemat_ascii

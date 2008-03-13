
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modifcs
  implicit none
  private

  public :: wien2k_energy_init
  public :: wien2k_energy_fetch
  public :: wien2k_energy_inquire
  public :: wien2k_energy_printinfo

  logical :: isinit_wien2k_energy
  data isinit_wien2k_energy / .false. /

  type t_wien2kfile_energy
     private
     character(256) :: fname
     integer :: nlheader,nkpt,nstmin,nstmax,nvmin,nvmax
     real(8) :: emin,emax
  end type t_wien2kfile_energy

  type(t_wien2kfile_energy) :: wien2kfile_energy

contains
  
  !///////////////////////////////////////////////////////////////////////////
  
  subroutine wien2k_pmat_import()
    implicit none
  end subroutine wien2k_pmat_import

  !///////////////////////////////////////////////////////////////////////////

  subroutine wien2k_energy_init(un,fname,nlheader)
    implicit none
    ! arguments
    integer, intent(in) :: un,nlheader
    character(*), intent(in) :: fname
    ! local variables
    character(3) :: tc3
    character(10) :: tc10
    real(8) :: k(3),en
    integer :: j,ist,nv,ne
    ! initialize
    wien2kfile_energy%fname = trim(fname)
    wien2kfile_energy%nlheader = nlheader
    wien2kfile_energy%nkpt = 0
    wien2kfile_energy%nstmin = 1000000000
    wien2kfile_energy%nstmax = 0
    wien2kfile_energy%nvmin = 1000000000
    wien2kfile_energy%nvmax = 0
    wien2kfile_energy%emin = 1.d9
    wien2kfile_energy%emax = -1.d9
    open(un,file=trim(fname),form='formatted',action='read',status='old')
    ! skip header
    do j=1,nlheader
       read(un,*) tc10
    end do
    do while (.true.)
       ! read k-point (format: '(3e19.12,a10,2i6,f5.1,a3)' )
       read(un,'(3e19.12,a10,2i6)',end=20) k,tc10,nv,ne
write(*,'(3e19.12,a10,2i6,f5.1,a3)') k,tc10,nv,ne
       ! add k-point
       wien2kfile_energy%nkpt = wien2kfile_energy%nkpt + 1
       ! minimum and maximum number of states for k-point
       wien2kfile_energy%nstmin=min(wien2kfile_energy%nstmin,ne)
       wien2kfile_energy%nstmax=max(wien2kfile_energy%nstmax,ne)
       ! minimum and maximum number of vectors for k-point
       wien2kfile_energy%nvmin=min(wien2kfile_energy%nvmin,nv)
       wien2kfile_energy%nvmax=max(wien2kfile_energy%nvmax,nv)
       do j=1,ne
          read(un,*) ist,en
write(*,'(i4,1x,g25.18)') ist,en
          ! minimum and maximum energies for k-point (converted to Hartrees)
          wien2kfile_energy%emin=min(wien2kfile_energy%emin,en/2.d0)
          wien2kfile_energy%emax=max(wien2kfile_energy%emax,en/2.d0)
       end do
    end do
20  continue
    close(un)
    isinit_wien2k_energy=.true.
  end subroutine wien2k_energy_init

  !///////////////////////////////////////////////////////////////////////////

  subroutine wien2k_energy_fetch(un,wvkl,we,nstk,vkl,epslat)
    implicit none
    integer, intent(in) :: un
    real(8), intent(out) :: wvkl(:,:)
    real(8), intent(out) :: we(:,:)
    integer, optional, intent(out) :: nstk(:)
    real(8), optional, intent(in) :: vkl(:,:)
    real(8), optional, intent(in) :: epslat
    ! local variables
    character(10) :: tc10
    integer :: s(1),s1(2),s2(2),iv(3),ik,jk,ist,j,nv,ne
    real(8) :: epsl,vl(3),vlt(3),en
    epsl=1.d-6
    if (present(epslat)) epsl=epslat
    call wien2k_energy_chkinit
    s1=shape(wvkl)
    if (present(vkl)) s2=shape(vkl)
    if ((s1(1).ne.3).or.(s1(2).ne.wien2kfile_energy%nkpt)) then
       write(*,*)
       write(*,'("Error(wien2k_energy_fetch): receiving array for &
            &k-points has wrong shape :",2i9)') s1
       write(*,'(" required shape:",2i9)') 3,wien2kfile_energy%nkpt
       write(*,*)
       stop
    end if
    if (present(vkl)) then
       if ((s2(1).ne.3).or.(s2(2).ne.wien2kfile_energy%nkpt)) then
          write(*,*)
          write(*,'("Error(wien2k_energy_fetch): reference array for &
               &k-points has wrong shape :",2i9)') s2
          write(*,'(" required shape:",2i9)') 3,wien2kfile_energy%nkpt
          write(*,*)
          stop
       end if
    end if
    s1=shape(we)
    if ((s1(1).lt.wien2kfile_energy%nstmax).or.(s1(2).ne. &
         wien2kfile_energy%nkpt)) then
       write(*,*)
       write(*,'("Error(wien2k_energy_fetch): receiving array for &
            &eigenvalues has wrong shape :",2i9)') s1
       write(*,'(" required shape:",2i9)') wien2kfile_energy%nstmax, &
            wien2kfile_energy%nkpt
       write(*,*)
       stop 
    end if
    if (present(nstk)) then
       s=shape(nstk)
       if (s(1).ne.wien2kfile_energy%nkpt) then
          write(*,'("Error(wien2k_energy_fetch): receiving array for &
               &number of eigenvalues has wrong shape :",2i9)') s
          write(*,'(" required shape:",2i9)') wien2kfile_energy%nkpt
          write(*,*)
          stop 
       end if
    end if
    ! initialize output variables
    wvkl(:,:)=0.d0
    we(:,:)=0.d0
    if (present(nstk)) nstk(:)=0
    open(un,file=trim(wien2kfile_energy%fname),form='formatted',action='read', &
         status='old')
    ! skip header
    do j=1,wien2kfile_energy%nlheader
       read(un,*) tc10
    end do
    do ik=1,wien2kfile_energy%nkpt
       read(un,'(3e19.12,a10,2i6)') vl,tc10,nv,ne
       jk=ik
       if (present(vkl)) then
          vlt(:)=vl(:)
          !***call r3frac(epsl,vlt,iv)
          do jk=1,wien2kfile_energy%nkpt
             if (sum(abs(vlt-vkl(:,jk))).lt.epsl) goto 30
          end do
          write(*,'("Error(wien2k_energy_fetch): could not find equivalent &
               &k-point in reference set")')
          write(*,'(" k-point          :",i9)') ik
          write(*,'(" fetched k-point  :",3g18.10)') vl
          write(*,*)
          stop
30        continue
       end if
       wvkl(:,jk)=vl(:)
       if (present(nstk)) then
          nstk(jk)=ne
       end if
       do j=1,ne
          read(un,*) ist,we(j,jk)
          if (ist.ne.j) then
             write(*,'("Error(wien2k_energy_fetch): mismatch in state &
                  &indices")')
             write(*,'(" k-point        :",i9,3g18.10)') ik,vl
             write(*,'(" required index :",3g18.10)') j
             write(*,'(" fetched index  :",3g18.10)') ist
             write(*,*)
             stop
          end if
       end do
    end do
    close(un)
  end subroutine wien2k_energy_fetch

  !///////////////////////////////////////////////////////////////////////////

  subroutine wien2k_energy_inquire(fname,nlheader,nkpt,nstmin,nstmax,nvmin, &
  	nvmax,emin,emax)
    implicit none
    character(*), optional, intent(out) :: fname
    integer, optional, intent(out) :: nlheader,nkpt,nstmin,nstmax,nvmin,nvmax
    real(8), optional, intent(out) :: emin,emax
    call wien2k_energy_chkinit
    if (present(fname)) then
       fname=trim(wien2kfile_energy%fname)
    end if
    if (present(nlheader)) then
       nlheader=wien2kfile_energy%nlheader
    end if
    if (present(nkpt)) then
       nkpt=wien2kfile_energy%nkpt
    end if
    if (present(nstmin)) then
       nstmin=wien2kfile_energy%nstmin
    end if
    if (present(nstmax)) then
       nstmax=wien2kfile_energy%nstmax
    end if
    if (present(nvmin)) then
       nvmin=wien2kfile_energy%nvmin
    end if
    if (present(nvmax)) then
       nvmax=wien2kfile_energy%nvmax
    end if
    if (present(emin)) then
       emin=wien2kfile_energy%emin
    end if
    if (present(emax)) then
       emax=wien2kfile_energy%emax
    end if    
  end subroutine wien2k_energy_inquire

  !///////////////////////////////////////////////////////////////////////////

  subroutine wien2k_energy_printinfo
    implicit none
    character(256) :: fname
    integer :: nlheader,nkpt,nstmin,nstmax,nvmin,nvmax
    real(8) :: emin,emax
    call wien2k_energy_chkinit
    call wien2k_energy_inquire(fname,nlheader,nkpt,nstmin,nstmax,nvmin, &
      nvmax,emin,emax)
    write(*,*)
    write(*,'("Info for WIEN2k energy-file: ",a)') trim(fname)
    write(*,'(" number of lines of header :",i6)') nlheader
    write(*,'(" number of k-points        :",i6)') nkpt
    write(*,'(" minimum number of vectors :",i6)') nvmin
    write(*,'(" maximum number of vectors :",i6)') nvmax
    write(*,'(" minimum number of states  :",i6)') nstmin
    write(*,'(" maximum number of states  :",i6)') nstmax
    write(*,'(" minimum state energy      :",g25.18)') emin
    write(*,'(" maximum state energy      :",g25.18)') emax
    write(*,*)
  end subroutine wien2k_energy_printinfo

  !///////////////////////////////////////////////////////////////////////////

  subroutine wien2k_energy_chkinit
    implicit none
    if (.not.isinit_wien2k_energy) call term_err('energy file not initialized')
  end subroutine wien2k_energy_chkinit

  !///////////////////////////////////////////////////////////////////////////
 
  subroutine wien2k_energy_finit
    implicit none
    ! initialize
    wien2kfile_energy%fname = ''
    wien2kfile_energy%nlheader = 0
    wien2kfile_energy%nkpt = 0
    wien2kfile_energy%nstmin = 0
    wien2kfile_energy%nstmax = 0
    wien2kfile_energy%nvmin = 0
    wien2kfile_energy%nvmax = 0
    wien2kfile_energy%emin = 1.d9
    wien2kfile_energy%emax = -1.d9    
  end subroutine wien2k_energy_finit
 
  !///////////////////////////////////////////////////////////////////////////

  subroutine term_err(str)
    implicit none
    character(*), intent(in) :: str
    write(*,'(a)') trim(str)
    stop 'invoked by: term_err'
  end subroutine term_err

end module modifcs



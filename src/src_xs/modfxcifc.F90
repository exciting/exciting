
! Copyright (C) 2004-2007 S. Sagmeister, J. K. Dewhurst, S. Sharma and 
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modfxcifc
  implicit none
contains

!BOP
! !ROUTINE: fxcifc
! !INTERFACE:
  subroutine fxcifc(fxctype,w,wgrid,iq,ng,alrc,alrcd,blrcd,lrccoef,fxcg)
    use m_fxc_lrc
    use m_fxc_lrcd
!    use m_fxc_lrcmodel
    use m_fxc_alda
    use m_fxc_bse_ma03
! !INPUT/OUTPUT PARAMETERS:
!   fxctype : type of exchange-correlation functional (in,integer)
! !DESCRIPTION:
!   Interface to the exchange-correlation kernel routines. This makes it 
!   relatively
!   simple to add new functionals which do not necessarily depend only on
!   all input parameters. Based upon the routine {\tt modxcifc}.
!
! !REVISION HISTORY:
!   Created October 2007 (Sagmeister)
!EOP
!BOC
    implicit none
    ! mandatory arguments
    integer, intent(in) :: fxctype
    ! optional arguments
    complex(8), optional, intent(in) :: w
    complex(8), optional, intent(in) :: wgrid(:)
    integer, optional, intent(in) :: iq
    integer, optional, intent(in) :: ng
    real(8), optional, intent(in) :: alrc
    real(8), optional, intent(in) :: alrcd
    real(8), optional, intent(in) :: blrcd
    real(8), optional, intent(in) :: lrccoef(:,:)
    complex(8), optional, intent(out) :: fxcg(:,:)
    ! local variables
    ! automatic arrays
    select case(abs(fxctype))
    case(0)
       ! RPA case fxc is zero
       if (present(fxcg).and.(present(ng))) then
          fxcg=(0.d0,0.d0)
       else
          goto 10
       end if
    case(1)
       ! static long-range kernel without local field effects
       ! L. Reining, Phys. Rev. Lett. 88, 06404 (2002)
       if (present(fxcg).and.(present(ng)).and.(present(alrc))) then
          call fxc_lrc(ng,.false.,alrc,fxcg)
       else
          goto 10
       end if
    case(2)
       ! static long-range kernel including local field effects
       ! L. Reining, Phys. Rev. Lett. 88, 06404 (2002)
       if (present(fxcg).and.(present(ng)).and.(present(alrc))) then
          call fxc_lrc(ng,.true.,alrc,fxcg)
       else
          goto 10
       end if
    case(3)
       ! dynamical long-range kernel without local field effects
       ! L. Reining, Phys. Rev. B 72, 125203 (2005)
       if (present(fxcg).and.(present(ng)).and.(present(alrcd)) &
            .and.(present(w)).and.(present(blrcd))) then
          call fxc_lrcd(ng,.false.,alrcd,blrcd,w,fxcg)
       else
          goto 10
       end if
    case(4)
       ! static long-range kernel including local field effects
       ! L. Reining, Phys. Rev. B 72, 125203 (2005)
       if (present(fxcg).and.(present(ng)).and.(present(alrcd)) &
            .and.(present(w)).and.(present(blrcd))) then
          call fxc_lrcd(ng,.true.,alrcd,blrcd,w,fxcg)
       else
          goto 10
       end if
    case(5)
       ! ALDA kernel, [Reference]
       if (present(fxcg).and.(present(ng)).and.present(iq)) then
          call fxc_alda(iq,ng,fxcg)
       else
          goto 10
       end if
    case(6)
       ! LRC frequency dependent model (S. Sagmeister)
       if (present(fxcg).and.(present(ng)).and.(present(w))&
            .and.present(lrccoef)) then
!          call fxc_lrcmodel(ng,.true.,lrccoef,w,fxcg)
       else
          goto 10
       end if
    case(7)
       ! xc-kernel derived from the Bethe-Salpeter equation
       ! no local field effects
       ! A. Marini, Phys. Rev. Lett. 91, 256402 (2003)
       if (present(fxcg).and.(present(ng)).and.(present(wgrid))) then
          call fxc_bse_ma03(ng,.false.,wgrid,fxcg)
       else
          goto 10
       end if
    case(8)
       ! xc-kernel derived from the Bethe-Salpeter equation
       ! inclusion of local field effects
       ! A. Marini, Phys. Rev. Lett. 91, 256402 (2003)
       if (present(fxcg).and.(present(ng)).and.(present(wgrid))) then
          call fxc_bse_ma03(ng,.true.,wgrid,fxcg)
       else
          goto 10
       end if
    case default
       write(*,*)
       write(*,'("Error(fxcifc): fxctype not defined : ",I8)') fxctype
       write(*,*)
       stop
    end select
    return
    ! error treatment
10  continue
    write(*,*)
    write(*,'("Error(fxcifc): missing arguments for exchange-correlation &
         &kernel type ",I5)') fxctype
    write(*,*)
    stop
  end subroutine fxcifc
!EOC

!BOP
! !ROUTINE: getfxcdata
! !INTERFACE:
  subroutine getfxcdata(fxctype,fxcdescr,fxcspin)
! !INPUT/OUTPUT PARAMETERS:
!   fxctype  : type of exchange-correlation functional (in,integer)
!   fxcdescr : description of functional (out,character(256))
!   fxcspin  : spin treatment (out,integer)
!   fxcgrad  : gradient treatment (out,integer)
! !DESCRIPTION:
!   Returns data on the exchange-correlation functional labelled by
!   {\tt fxctype}. The character array {\tt fxctype} contains a short
!   description
!   of the functional including journal references. The variable 
!   {\tt fxcspin} is
!   set to 1 or 0 for spin-polarised or -unpolarised functionals,
!   respectively.
!
! !REVISION HISTORY:
!   Created October 2007 (Sagmeister)
!EOP
!BOC
    implicit none
    integer, intent(in) :: fxctype
    character(256), intent(out) :: fxcdescr
    integer, intent(out) :: fxcspin
    select case(abs(fxctype))
    case(0)
       fxcdescr='xc-kernel set to zero (RPA case)'
       ! spin-polarisation not required
       fxcspin=-1
       return     
    case(1)
       fxcdescr='long-range xc-kernel, no local field effects'
       ! spin-polarisation not required
       fxcspin=-1
       return
    case(2)
       fxcdescr='long-range xc-kernel, including local field effects'
       ! spin-polarisation not required
       fxcspin=0
       return
    case(3)
       fxcdescr='dynamical long-range xc-kernel, no local field effects'
       ! spin-polarisation not required
       fxcspin=-1
       return
    case(4)
       fxcdescr='dynamical long-range xc-kernel, including local field effects'
       ! spin-polarisation not required
       fxcspin=0
       return
    case(5)
       fxcdescr='ALDA kernel, including local field effects'
       ! spin-polarisation not required
       fxcspin=0
       return
    case(7)
       fxcdescr='BSE kernel, A. Marini, Phys. Rev. Lett. 91, 256402 (2003)'
       ! spin-polarisation not required
       fxcspin=0
       return
    case default
       write(*,*)
       write(*,'("Error(getfxcdata): fxctype not defined : ",I8)') fxctype
       write(*,*)
       stop
    end select
  end subroutine getfxcdata
!EOC

end module modfxcifc


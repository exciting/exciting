! Copyright (C) 2007-2010 D. Nabok, P. Puschnig and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine read_input()

  use param
  implicit none
  integer :: i, iostat
    
! Reading of integration parameters
  read(*,*)
  read(*,'(a)') phifile
  phifile=adjustl(phifile)
  i=index(phifile,' ')
  phifile=trim(phifile(1:i))//'/kernel.dat'
  call read_phi()
  
  read(*,*) vdWDF_version
  read(*,*) nrx,nry,nrz
  read(*,*) epsrel
  read(*,*) epsabs
  read(*,*) key1  
  read(*,*) maxeval
! Reading name of the density file
  read(*,*); read(*,*)
  read(*,*) densunits
  read(*,'(a)') xsffile

  xsffile=adjustl(xsffile)
  i=index(xsffile,' ')
  xsffile=xsffile(1:i)
  call read_densities()

! Reading name of the gradient density file (optional)
  read(*,'(a)',IOSTAT=iostat) xsfgradients
  if (iostat==0) then
      xsfgradients=adjustl(xsfgradients)
      i=index(xsfgradients,' ')
      xsfgradients=xsfgradients(1:i)
  else
      xsfgradients=''
  end if
  call read_gradients()
  
  write(*,*)
  write(*,*) '**********************  noloco  **********************' 
  write(* ,'(a,a)')     '   vdW-DF kernel file            : ', trim(phifile)
  write(* ,'(a,a)')     '   vdW-DF type                   : ', vdWDF_version
  write(*,'(a,i2,2x,i2,2x,i2)')&
                        '   supercell size                : ', nrx,nry,nrz
  write(* ,'(a,D12.3)') '   relative integration accuracy : ', epsrel
  write(* ,'(a,D12.3)') '   absolute integration accuracy : ', epsabs
  write(*,*)
  write(* ,'(a,a)')     '   Density file                  : ', trim(xsffile)
  if (LEN_TRIM(xsfgradients)>0) then
    write(* ,'(a,a)')     '   Density gradients file        : ', trim(xsfgradients)
  end if

end subroutine read_input

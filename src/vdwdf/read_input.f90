subroutine read_input(phifile,xsffile)

  use param
  
  implicit none
  character*80    :: phifile
  character*80    :: xsffile
    
  write(*,*)
  write(*,*) '**********************  OUTPUT OF PROGRAM vdWDF.x  **********************' 
! Reading of integration parameters
  read(*,*)
  read(*,*) phifile
  read(*,*) vdWDF_version
  read(*,*) nrx,nry,nrz
  read(*,*) epsrel
  read(*,*) epsabs
  
  write(* ,'(a,a)')     '   vdW-DF kernel file            : ', phifile
  write(* ,'(a,a)')     '   vdW-DF type                   : ', vdWDF_version
  write(*,'(a,i2,2x,i2,2x,i2)')&
                        '   supercell size                : ', nrx,nry,nrz
  write(* ,'(a,D12.3)') '   relative integration accuracy : ', epsrel
  write(* ,'(a,D12.3)') '   absolute integration accuracy : ', epsabs
  write(*,*)
! Reading name of the density file
  read(*,*); read(*,*)
  read(*,*) xsffile
  write(* ,'(a,a)')     '   Density file                  : ', xsffile
! Reading additional integrational parameters
  read(*,*); read(*,*); read(*,*)
  read(*,*) flags
  read(*,*) mineval
  read(*,*) maxeval
  read(*,*) key1
  read(*,*) key2
  read(*,*) key3
  read(*,*) maxpass
  read(*,*) border
  read(*,*) maxchisq
  read(*,*) mindeviation
  
end subroutine read_input


program xcinfo
  use xc_f03_lib_m

  implicit none

  TYPE(xc_f03_func_t) :: xc_func
  TYPE(xc_f03_func_info_t) :: xc_info
  TYPE(xc_f03_func_reference_t) :: xc_ref
  integer :: i, number, func_id
  character(len=120) :: s1, s2

  ! Query libxc for the ID number of PBE correlation
  func_id = xc_f03_functional_get_number("GGA_C_PBE")

  call xc_f03_func_init(xc_func, func_id, XC_UNPOLARIZED)
  xc_info = xc_f03_func_get_info(xc_func)

  select case(xc_f03_func_info_get_kind(xc_info))
  case(XC_EXCHANGE)
    write(*, '(a)') 'Exchange'
  case(XC_CORRELATION)
    write(*, '(a)') 'Correlation'
  case(XC_EXCHANGE_CORRELATION)
    write(*, '(a)') 'Exchange-correlation'
  case(XC_KINETIC)
    write(*, '(a)') 'Kinetic'
  end select

  s1 = xc_f03_func_info_get_name(xc_info)
  select case(xc_f03_func_info_get_family(xc_info))
  case (XC_FAMILY_LDA);       write(s2,'(a)') "LDA"
  case (XC_FAMILY_HYB_LDA);   write(s2,'(a)') "Hybrid LDA"
  case (XC_FAMILY_GGA);       write(s2,'(a)') "GGA"
  case (XC_FAMILY_HYB_GGA);   write(s2,'(a)') "Hybrid GGA"
  case (XC_FAMILY_MGGA);      write(s2,'(a)') "MGGA"
  case (XC_FAMILY_HYB_MGGA);  write(s2,'(a)') "Hybrid MGGA"
  case (XC_FAMILY_LCA);       write(s2,'(a)') "LCA"
  end select
  write(*, '(4a)') trim(s1), ' (', trim(s2), ')'

  do i=0, XC_MAX_REFERENCES-1
     ! Get info on reference number i; can't reuse the loop variable
     ! due to Fortran restrictions for intent(inout)
     number = i
     xc_ref = xc_f03_func_info_get_references(xc_info, number)
     if( number < 0 ) then
        exit
     end if
    write(*, '(a,i1,2a)') '[', i+1, '] ', trim(xc_f03_func_reference_get_ref(xc_ref))
  end do

  call xc_f03_func_end(xc_func)

end program xcinfo

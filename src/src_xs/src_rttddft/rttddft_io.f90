!> Module to expose the IO methods and variables required in RT-TDDFT
module rttddft_io
  use rttddft_file_formats
  use rttddft_file_names
  use rttddft_io_formatted
  use rttddft_io_unformatted
  
  implicit none

  private

  ! procedures
  public :: close_file_etot, &
            close_file_info, &
            close_file_nexc, &
            close_file_timing, &
            close_files_vector_fields, &
            copy_files, &
            delete_jpa_files, &
            delete_pmat_binary_file, &
            delete_pmat_mt_binary_file, &
            delete_wavefunction_binary_file, &
            file_pmat_exists, &
            file_pmat_mt_exists, &
            get_filename_pmat, &
            get_filename_pmat_mt, &
            open_file_etot, &
            open_file_info, &
            open_files_vector_fields, &
            open_file_nexc, &
            open_file_timing, &
            read_phases, &
            read_pmat, &
            read_pmat_mt, &
            read_state_Ehrenfest_MD, &
            read_vector_field, &
            read_wavefunction, &
            write_density_to_file, &
            write_eigenvalues, &
            write_file_info, &
            write_file_info_header,&
            write_file_info_fill_line_with_char, &
            write_nexc, &
            write_occupations, &
            write_phases, &
            write_pmat, &
            write_pmat_mt, &
            write_projection_coefficients, &
            write_state_Ehrenfest_MD, &
            write_timing, &
            write_total_energy, &
            write_vector_field, &
            write_wavefunction

  ! types
  public :: file_handler

  ! kinds and enums
  public :: binary, &
            groundstate, &
            hdf5, &
            restart_format, &
            t, &
            t_minus_dt

end module

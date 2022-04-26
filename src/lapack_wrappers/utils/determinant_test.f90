!> Unit tests for determinant
module determinant_test
  use precision, only: dp
  use constants, only: zzero
  use modmpi, only: mpiinfo
  use unit_test_framework, only: unit_test_type
  use math_utils, only: all_close
  use determinant, only: determinant_LU

  implicit none

  private
  public :: determinant_test_driver

  contains

  !> Run tests for determinant
  subroutine determinant_test_driver(mpiglobal, kill_on_failure)
    !> mpi environment
    type(mpiinfo), intent(in) :: mpiglobal
    !> Kill the program upon failure of an assertion
    logical, intent(in), optional :: kill_on_failure

    !> Test report object
    type(unit_test_type) :: test_report
    !> Number of assertions
    integer, parameter :: n_assertions = 6

    call test_report%init(n_assertions, mpiglobal)

    call test_determinant(test_report)

    if (present(kill_on_failure)) then
      call test_report%report('determinant', kill_on_failure)
    else
      call test_report%report('determinant')
    end if

    call test_report%finalise()
  end subroutine determinant_test_driver

  !> Test determinant
  subroutine test_determinant(test_report)
    !> Test report object
    type(unit_test_type), intent(inout) :: test_report
    !> test input
    real(dp) :: reference_2x2_real, A_real(5, 5), A_2x2_real(2, 2)
    complex(dp) :: reference_2x2_complex, A_complex(5, 5), A_2x2_complex(2, 2)

    A_2x2_real = transpose(reshape([1.00_dp, -3.893_dp, &
                                    0.23_dp,  2.12323_dp], [2, 2]))

    reference_2x2_real = A_2x2_real(1, 1) * A_2x2_real(2, 2) - A_2x2_real(1, 2) * A_2x2_real(2, 1)
    call test_report%assert(all_close(determinant_LU(A_2x2_real), reference_2x2_real), &
               'Test determinant explicitly for a 2 x 2 real matrix.')

    A_real = transpose(reshape([ 48.09009230_dp,  36.23203237_dp, -29.97433087_dp,   7.00464893_dp,  -1.57860058_dp, &
                                 -3.63775329_dp, -10.09285910_dp, -11.26465777_dp, -47.85961607_dp, -37.48520036_dp, &
                                -12.19145971_dp,   0.57774986_dp, -28.75819209_dp,  23.07833972_dp, -27.63199486_dp, &
                                -40.22832003_dp, -27.62653442_dp,  28.07067103_dp,   9.71508023_dp,  12.97875268_dp, &
                                 19.90281295_dp, -20.22667993_dp,   7.61766909_dp,  37.54706795_dp, 32.9927136_dp ], &
                                 [5, 5]))

    call test_report%assert(all_close(determinant_LU(A_real), -3401729.351128587964922_dp, 1e-8_dp), &
               'Test determinant for real a bijective matrix (the deteminant is not zero).')

    A_real = transpose(reshape([ 48.09009230_dp,  36.23203237_dp, 0.0_dp,   7.00464893_dp,  -1.57860058_dp, &
                                 -3.63775329_dp, -10.09285910_dp, 0.0_dp, -47.85961607_dp, -37.48520036_dp, &
                                -12.19145971_dp,   0.57774986_dp, 0.0_dp,  23.07833972_dp, -27.63199486_dp, &
                                -40.22832003_dp, -27.62653442_dp, 0.0_dp,   9.71508023_dp,  12.97875268_dp, &
                                 19.90281295_dp, -20.22667993_dp, 0.0_dp,  37.54706795_dp, 32.9927136_dp ], &
                                 [5, 5]))

    call test_report%assert(all_close(determinant_LU(A_real), 0.0_dp), &
                            'Test determinant for real a non-bijective matrix (the determinant is zero). &
                            Expected: 0.0')

    A_2x2_complex = transpose(reshape(cmplx([  1.00_dp, -3.893_dp, &
                                             -23.032_dp, 5.893_dp], &

                                            [ -2.321_dp, 12.123_dp, &
                                               5.321_dp, -0.123_dp]), [2, 2]))

    reference_2x2_complex = A_2x2_complex(1, 1) * A_2x2_complex(2, 2) - A_2x2_complex(1, 2) * A_2x2_complex(2, 1)
    call test_report%assert(all_close(determinant_LU(A_2x2_real), reference_2x2_real), &
               'Test determinant explicitly for a 2 x 2 complex matrix.')

    A_complex = transpose(reshape(cmplx( &
            [ -44.52110712_dp, -48.17481048_dp,    -19.66887_dp,   2.95379091_dp, -32.48344889_dp, &
              -48.40959198_dp,   36.6255884_dp, -40.95025262_dp,  -9.47597016_dp,  21.96138633_dp, &
              -23.89308839_dp, -46.60346009_dp,  29.51933189_dp,  46.80034723_dp, -26.99063883_dp, &
              -28.34722129_dp, -15.28477719_dp,  -9.98834923_dp,  19.58476894_dp, -45.2623859_dp , &
               25.18235434_dp,  45.25070169_dp,  39.07849915_dp,  15.93610182_dp, -33.5371513_dp ], &

            [  14.86243624_dp,  30.70582514_dp,   44.7784307_dp, -46.66873007_dp, -24.13516175_dp, &
               32.30310610_dp,  33.04884469_dp,  20.86135191_dp, -42.08274516_dp,  48.39385979_dp, &
               11.45829782_dp,  36.55565904_dp,  34.21616482_dp,  29.01076295_dp, -40.37086456_dp, &
               32.39639731_dp, -45.57484939_dp,  27.23970665_dp,  -3.86243919_dp,  22.13109274_dp, &
                 5.8363281_dp,   9.02281896_dp,   -9.2413642_dp, -18.39173734_dp, -27.88252968_dp], &
            kind=dp), [5, 5]))

    call test_report%assert(all_close(determinant_LU(A_complex), &
            cmplx(403357472.347355604171753_dp, 676258729.01441454887390137_dp, kind=dp), 1e-6_dp), &
                           'Test determinant for a complex bijective matrix (the deteminant is not zero).')

    A_complex = transpose(reshape(cmplx( &
            [ -44.52110712_dp, -48.17481048_dp, 0.0_dp,   2.95379091_dp, -32.48344889_dp, &
              -48.40959198_dp,   36.6255884_dp, 0.0_dp,  -9.47597016_dp,  21.96138633_dp, &
              -23.89308839_dp, -46.60346009_dp, 0.0_dp,  46.80034723_dp, -26.99063883_dp, &
              -28.34722129_dp, -15.28477719_dp, 0.0_dp,  19.58476894_dp, -45.2623859_dp , &
               25.18235434_dp,  45.25070169_dp, 0.0_dp,  15.93610182_dp, -33.5371513_dp ], &

            [  14.86243624_dp,  30.70582514_dp, 0.0_dp, -46.66873007_dp, -24.13516175_dp, &
               32.30310610_dp,  33.04884469_dp, 0.0_dp, -42.08274516_dp,  48.39385979_dp, &
               11.45829782_dp,  36.55565904_dp, 0.0_dp,  29.01076295_dp, -40.37086456_dp, &
               32.39639731_dp, -45.57484939_dp, 0.0_dp,  -3.86243919_dp,  22.13109274_dp, &
                 5.8363281_dp,   9.02281896_dp, 0.0_dp, -18.39173734_dp, -27.88252968_dp], &
            kind=dp), [5, 5]))
    call test_report%assert(all_close(determinant_LU(A_complex), zzero), &
                            'Test determinant for a complex non-bijective matrix (the determinant is zero)')
  end subroutine test_determinant

end module determinant_test
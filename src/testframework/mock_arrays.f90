!> Mock arrays for testing.
module mock_arrays
  use precision, only: dp, sp

  implicit none

  private
  public :: fill_array, value_map_real_rank1, &
                        value_map_real_rank2, &
                        value_map_real_rank3, &
                        value_map_complex_rank2, &
                        value_map_complex_rank3, &
                        value_map_complex_rank4


! double arrays

  !> Real vecor with 5 elements
  real(dp), public :: real_vector_5(5) = [ 31.18072641_dp,  -4.72689264_dp,  48.67919055_dp, -41.66440826_dp, &
                                          -40.21374565_dp]
  !> Real vector with 7 elements
  real(dp), public :: real_vector_7(7) = [-26.54377445_dp, -17.22046184_dp,   4.83382119_dp,  -9.47731147_dp, &
                                            0.03664525_dp,   3.04736258_dp,  17.84092835_dp]
  !> Real 5 x 5 matrix
  real(dp), public :: real_matrix_5x5(5, 5) = transpose(reshape([-40.17392854_dp, -29.62435155_dp,  -0.82734985_dp, -10.30943673_dp, &
                                                                  48.46185434_dp, &
                                                                  47.31299703_dp,   6.77807354_dp,  11.67496062_dp, -35.14904297_dp, &
                                                                  14.48758125_dp, &
                                                                 -25.42093406_dp,  23.68830969_dp, -23.00607822_dp,   12.7168522_dp, &
                                                                  22.11271185_dp, &
                                                                  26.66413809_dp, -33.56388754_dp,  42.48611478_dp,   -1.0718386_dp, &
                                                                   -7.1753077_dp , &
                                                                  16.33945082_dp,   39.7244291_dp,  -25.5748391_dp,    39.604497_dp, &
                                                                  33.70929988_dp], [5, 5]))
  !> Real 5 x 7 matrix
  real(dp), public :: real_matrix_5x7(5, 7) = transpose(reshape([   9.06173884_dp, -16.40962765_dp,  15.66983533_dp,  29.14170303_dp, &
                                                                   30.01592387_dp, -33.14855289_dp,   9.77584235_dp, &
                                                                  -16.41480408_dp, -13.45202387_dp, -35.18713482_dp,  -44.0682099_dp , &
                                                                  -34.35152851_dp, -47.16065796_dp,   4.17329693_dp, &
                                                                   47.71787082_dp,  28.31043534_dp,  17.57815792_dp,  -29.0436189_dp, &
                                                                    8.28924983_dp,   2.12933902_dp,    7.8501177_dp , &
                                                                   29.27877003_dp,  -2.88508885_dp, -46.28881163_dp,  18.08679224_dp, &
                                                                   -3.27836804_dp, -27.90913228_dp,  18.47728991_dp, &
                                                                   33.05916149_dp, -35.07317596_dp,  -5.39633473_dp, -39.78075197_dp, &
                                                                   45.95078973_dp, -12.78879693_dp, -15.22917256_dp ], [7, 5]))
  !> Real 7 x 5 matrix
  real(dp), public :: real_matrix_7x5(7, 5) = transpose(reshape([  -2.89185458_dp,  19.07126067_dp,  15.01048731_dp, -38.24203168_dp, &
                                                                   29.64731979_dp,  33.46992354_dp,   9.31958249_dp,  46.69694407_dp, &
                                                                  -39.98383303_dp,  -43.6151068_dp,  -27.3546948_dp,   8.55055018_dp, &
                                                                  -46.92833888_dp,  -5.20754552_dp, -38.24928213_dp, -16.91100668_dp, &
                                                                  -27.22135749_dp,  21.90649908_dp, -16.83657893_dp, -13.16367346_dp, &
                                                                  -34.82759992_dp,  -9.14037716_dp,  45.38223981_dp,  10.98111997_dp, &
                                                                    6.23074262_dp, -20.34205625_dp,  28.38064005_dp,   2.76094858_dp, &
                                                                   -3.15745602_dp,   3.16632168_dp, -16.68867829_dp, -18.46512746_dp, &
                                                                     3.1480205_dp,  18.97778663_dp,  47.67943281_dp,   2.45030698_dp, &
                                                                   18.89052178_dp], [5, 7]))
  !> Real symmetric 5 x 5 matrix
  real(dp), public :: real_symmetric_matrix_5x5(5, 5) = transpose(reshape([  47.3120746_dp,  26.77011639_dp,  28.32428222_dp,  27.30299212_dp, &
                                                                            -34.00535404_dp, &
                                                                             26.77011639_dp,  -25.1535981_dp, -45.97519099_dp,  35.42787447_dp, &
                                                                             28.19186526_dp, &
                                                                             28.32428222_dp, -45.97519099_dp, -46.98635297_dp, -20.04984215_dp, &
                                                                             24.58347956_dp, &
                                                                             27.30299212_dp,  35.42787447_dp, -20.04984215_dp,   0.25961623_dp, &
                                                                              0.69501014_dp, &
                                                                            -34.00535404_dp,  28.19186526_dp,  24.58347956_dp,   0.69501014_dp, &
                                                                            -28.44893406_dp], [5, 5]))
  !> Real 5 x 5 orthogonal matrix
  real(dp), public :: real_orthogonal_matrix_5x5(5, 5) = transpose(reshape([ 0.63806971_dp, -0.63857815_dp, -0.34365517_dp,  0.00522444_dp, -0.25876403_dp, &
                                                                            -0.26043362_dp,  0.3436803_dp , -0.72155634_dp,  0.11296424_dp, -0.52976759_dp, & 
                                                                             0.06634032_dp,  0.09696107_dp, -0.49973653_dp, -0.63693877_dp,  0.575126_dp  , &
                                                                             0.52753549_dp,  0.53139645_dp,  0.28488724_dp, -0.44699968_dp, -0.39793802_dp, &
                                                                             0.49228716_dp,  0.4269873_dp , -0.17424187_dp,  0.61782852_dp,  0.40405803_dp], [5, 5]))
  !> Real 5 x 5 full rank matrix
  real(dp), public :: real_full_rank_matrix_5x5(5, 5) = transpose(reshape([  48.0900923_dp,  36.23203237_dp, -29.97433087_dp,   7.00464893_dp, &
                                                                            -1.57860058_dp, &
                                                                            -3.63775329_dp,  -10.0928591_dp, -11.26465777_dp, -47.85961607_dp, &
                                                                           -37.48520036_dp, &
                                                                           -12.19145971_dp,   0.57774986_dp, -28.75819209_dp,  23.07833972_dp, &
                                                                           -27.63199486_dp, &
                                                                           -40.22832003_dp, -27.62653442_dp,  28.07067103_dp,   9.71508023_dp, &
                                                                            12.97875268_dp, &
                                                                            19.90281295_dp, -20.22667993_dp,   7.61766909_dp,  37.54706795_dp, &
                                                                             32.9927136_dp ], [5, 5]))
  !> Real 5 x 7 matrix with rank 4 
  real(dp), public :: real_rank_4_matrix_5x7(5, 7) = transpose(reshape([  9.06173884_dp, -16.40962765_dp,  15.66983533_dp,  29.14170303_dp, &
                                                                                 0.0_dp,          0.0_dp,          0.0_dp, &
                                                                        -16.41480408_dp, -13.45202387_dp, -35.18713482_dp, -44.0682099_dp , &
                                                                                 0.0_dp,          0.0_dp,          0.0_dp, &
                                                                         47.71787082_dp,  28.31043534_dp,  17.57815792_dp, -29.0436189_dp , &
                                                                                 0.0_dp,          0.0_dp,          0.0_dp, &
                                                                         29.27877003_dp,  -2.88508885_dp, -46.28881163_dp,  18.08679224_dp, &
                                                                                 0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                 0.0_dp,          0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                 0.0_dp,          0.0_dp,          0.0_dp], [7, 5]))
  !> Real 7 x 5 matrix with rank 4 
  real(dp), public :: real_rank_4_matrix_7x5(7, 5) = transpose(reshape([  9.06173884_dp, -16.40962765_dp,  15.66983533_dp,  29.14170303_dp, &
                                                                                 0.0_dp, &
                                                                        -16.41480408_dp, -13.45202387_dp, -35.18713482_dp, -44.0682099_dp , &
                                                                                 0.0_dp, &
                                                                         47.71787082_dp,  28.31043534_dp,  17.57815792_dp, -29.0436189_dp , &
                                                                                 0.0_dp, &
                                                                         29.27877003_dp,  -2.88508885_dp, -46.28881163_dp,  18.08679224_dp, &
                                                                                 0.0_dp, &
                                                                                 0.0_dp,          0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                 0.0_dp, &
                                                                                 0.0_dp,          0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                 0.0_dp, &
                                                                                 0.0_dp,          0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                 0.0_dp], [5, 7]))





! double complex arrays

  !> Complex vector with 5 elements
  complex(dp), public :: complex_vector_5(5) = cmplx([ -2.30866203_dp, -30.44677196_dp,   36.2118037_dp, -30.5659898_dp, &
                                                       47.58889362_dp], &

                                                     [ 12.73245688_dp,  -43.5729583_dp,  27.03688433_dp,  49.7539840_dp, &
                                                      -41.57934484_dp], kind=dp)
  !> Complex vector with 7 elements
  complex(dp), public :: complex_vector_7(7) = cmplx([ 25.62674667_dp,  -14.9142924_dp,  -0.49922657_dp,-20.83173308_dp, &
                                                       21.61559333_dp,  38.52653728_dp,   6.02442546_dp], &

                                                               [ 19.06604381,  10.34856432, -12.13899351,  16.90726442, &
                                                                 16.16283436,   3.91299653,  -9.94204067], kind=dp)
  !> Complex 5 x 5 matrix
  complex(dp), public :: complex_matrix_5x5(5, 5) = transpose(reshape(cmplx([ -40.17392854_dp, -29.62435155_dp,  -0.82734985_dp, -10.30943673_dp, &
                                                                               48.46185434_dp, &
                                                                               47.31299703_dp,   6.77807354_dp,  11.67496062_dp, -35.14904297_dp, &
                                                                               14.48758125_dp, &
                                                                              -25.42093406_dp,  23.68830969_dp, -23.00607822_dp,   12.7168522_dp, &
                                                                               22.11271185_dp, &
                                                                               26.66413809_dp, -33.56388754_dp,  42.48611478_dp,   -1.0718386_dp, &
                                                                                -7.1753077_dp , &
                                                                               16.33945082_dp,   39.7244291_dp,  -25.5748391_dp ,   39.604497_dp, &
                                                                               33.70929988_dp], &
                                                                      
                                                                            [  -3.62805399_dp, 30.17124227_dp,  -3.54995097_dp, -23.27306593_dp, &
                                                                               32.75861251_dp, &
                                                                               27.65944794_dp, 45.37960003_dp, -42.83317035_dp, -41.52218814_dp, &
                                                                                6.15550826_dp, &
                                                                               45.40326741_dp, 42.39384644_dp, -36.45239616_dp, -18.38045087_dp, &
                                                                                9.20916614_dp, &
                                                                                1.59569587_dp, 15.34311375_dp,   12.5354502_dp,  30.52183751_dp, &
                                                                                8.82668351_dp, &
                                                                               46.82891375_dp, 41.46778084_dp,  -7.79986042_dp,  -3.34967487_dp, &
                                                                              -22.63124669_dp], kind=dp), [5, 5]))
  !> Complex 5 x 7 matrix
  complex(dp), public :: complex_matrix_5x7(5, 7) = transpose(reshape(cmplx([ -18.7814739_dp, -35.75636149_dp, -25.89864277_dp,  11.04508789_dp, &
                                                                             -20.79065069_dp,  47.42124332_dp,  13.17763894_dp, &
                                                                              33.75023754_dp,  49.97417692_dp,  -0.48275953_dp,  43.98236553_dp, &
                                                                               1.36495302_dp,   1.01455123_dp,  45.98544819_dp, &
                                                                              21.48351648_dp,   1.63997878_dp,  16.31513844_dp,  29.47589814_dp, &
                                                                               2.67889982_dp,   3.08151508_dp, -21.94527983_dp, &
                                                                              40.49274295_dp, -28.86521605_dp, -49.94130354_dp,  48.91554628_dp, &
                                                                             -22.78898301_dp,   9.67631696_dp, -26.67369113_dp, &
                                                                             -17.20282601_dp, -19.55147567_dp,  18.65011868_dp,  22.76807154_dp, &
                                                                               41.5186408_dp,   -6.0275415_dp,  43.86113634_dp], &

                                                                            [-28.97154895_dp, -30.01164974_dp,  37.37357186_dp,  15.58097948_dp, &
                                                                               6.61000788_dp, -37.94313822_dp,  40.90130461_dp, &
                                                                              -33.8514265_dp,  34.40405017_dp, -44.22889464_dp,  41.54835073_dp, &
                                                                               5.29344049_dp,  36.55668618_dp,   0.28776443_dp, &
                                                                              27.21430366_dp,  34.61978756_dp,  32.66960231_dp, -36.56589004_dp, &
                                                                              18.11818186_dp,  26.94329308_dp,  -0.50794621_dp, &
                                                                             -43.59929192_dp,  34.88338453_dp,   -4.6433614_dp,   9.07841908_dp, &
                                                                             -43.20518396_dp, -18.72251584_dp,  35.96143653_dp, &
                                                                             -16.61908417_dp,  48.22907402_dp,  22.75478376_dp,  19.73577893_dp, &
                                                                              22.80786931_dp, -35.12437444_dp,  -1.31678731_dp], kind=dp), [7, 5]))
  !> Complex 7 x 5 matrix
  complex(dp), public :: complex_matrix_7x5(7, 5) = transpose(reshape(cmplx([-15.76457636_dp, -13.10633076_dp,  24.50742622_dp,  27.78197405_dp, &
                                                                              16.85583405_dp,   40.0532989_dp,  -0.79731881_dp,  10.38030932_dp, &
                                                                              23.18130197_dp,  30.15582799_dp,  29.54712995_dp, -46.36497795_dp, &
                                                                              -2.25686993_dp,   0.85856371_dp, -17.66400024_dp,  -7.45825453_dp, &
                                                                             -41.26169756_dp, -37.50047737_dp,  32.49023951_dp,  36.32448571_dp, &
                                                                             -45.41645606_dp,  12.93965366_dp, -32.02159024_dp,   8.55410596_dp, &
                                                                               1.65255132_dp, -31.28010187_dp,  30.85471657_dp,    7.3187356_dp, &
                                                                             -17.95181164_dp, -38.85225614_dp, -38.25951835_dp, -22.60543684_dp, &
                                                                             -12.69575022_dp,  14.96294671_dp, -33.74281834_dp, -36.73461298_dp, &
                                                                               4.77888767_dp], &
                                                                            
                                                                            [-29.56785948_dp,  -29.3596185_dp, -26.06947443_dp,  49.95145734_dp, &
                                                                             -39.33928642_dp,   4.23622163_dp,  -7.43201071_dp,  38.28895233_dp, &
                                                                               3.59232427_dp, -30.81677821_dp, -27.49969413_dp,  -20.1539012_dp, &
                                                                               -6.9841385_dp,  -5.90651352_dp,   38.0935866_dp,  15.65833432_dp, &
                                                                              -1.76539909_dp, -33.92278221_dp, -21.03924414_dp,   2.04761883_dp, &
                                                                             -48.57899107_dp,  -1.29368174_dp, -36.97271557_dp,   5.60701216_dp, &
                                                                              -4.40070446_dp,  12.45335491_dp,  48.97361474_dp, -41.35112939_dp, &
                                                                             -31.25811089_dp, -11.36623968_dp,     8.885927_dp,  17.71068798_dp, &
                                                                              15.11203334_dp,   1.67623585_dp,  43.41324191_dp, -10.89571735_dp, &
                                                                             -34.18643191_dp], kind=dp), [5, 7]))
  !> Complex hermitian 5 x 5 matrix
  complex(dp), public :: complex_hermitian_matrix_5x5(5, 5) = transpose(reshape(cmplx([-49.72690415_dp,  -14.84768354_dp,  10.05665554_dp,  29.36660616_dp, &
                                                                                        -20.0201864_dp,  -14.84768354_dp, -15.53149339_dp, -11.64317351_dp, &
                                                                                         3.02310325_dp,  -16.29220719_dp,  10.05665554_dp, -11.64317351_dp, &
                                                                                        47.84923206_dp,   -5.72778535_dp,  13.16564869_dp,  29.36660616_dp, &
                                                                                         3.02310325_dp,   -5.72778535_dp,  38.45072588_dp,  21.26837488_dp, &
                                                                                        -20.0201864_dp , -16.29220719_dp,  13.16564869_dp,  21.26837488_dp, &
                                                                                         42.3303856_dp ], &

                                                                                      [         0.0_dp,   5.21297468_dp,  -0.68828716_dp,  -18.61860468_dp, &
                                                                                       -23.20913956_dp,  -5.21297468_dp,          0.0_dp,   31.94011008_dp, &
                                                                                         -2.2229102_dp,  30.48337625_dp,   0.68828716_dp,  -31.94011008_dp, &
                                                                                                0.0_dp, -10.65927965_dp,   8.38358035_dp,   18.61860468_dp, &
                                                                                          2.2229102_dp,  10.65927965_dp,          0.0_dp,   34.28279495_dp, &
                                                                                        23.20913956_dp, -30.48337625_dp,  -8.38358035_dp,  -34.28279495_dp, &
                                                                                                0.0_dp], kind=dp), [5, 5]))
  !> Complex unitary 5 x 5 matrix
  complex(dp), public :: complex_unitary_matrix_5x5(5, 5) = transpose(reshape(cmplx([-0.33365437_dp, -0.21756134_dp,  0.07501474_dp, -0.28585794_dp,  0.10201445_dp, &
                                                                                     -0.4171256_dp ,  0.1862851_dp , -0.09447265_dp, -0.32270405_dp, -0.10475249_dp, &
                                                                                     -0.24272685_dp,  0.15928472_dp, -0.50208627_dp,  0.0213966_dp ,  0.61711886_dp, &
                                                                                     -0.22716718_dp,  0.14406867_dp,  0.25321109_dp, -0.29045102_dp,  0.49170423_dp, &
                                                                                      0.18721289_dp,  0.11869668_dp, -0.36368011_dp, -0.30958696_dp, -0.26682812_dp], &
                                                                                    
                                                                                    [  0.47787822_dp, -0.17924166_dp, -0.18438518_dp, -0.66834485_dp,  0.04914282_dp, &
                                                                                      -0.28302054_dp,  0.58705115_dp,  0.21042146_dp, -0.18848737_dp, -0.40340262_dp, &
                                                                                       0.18424181_dp, -0.0643122_dp ,  0.42917659_dp,  0.14321723_dp,  0.1988545_dp , &
                                                                                      -0.31417612_dp, -0.3084001_dp , -0.46322849_dp,  0.26201446_dp, -0.24564024_dp, &
                                                                                      -0.35775856_dp, -0.61785542_dp,  0.24480553_dp, -0.2492159_dp , -0.14065227_dp], kind=dp), [5, 5]))
  !> Complex full rank 5 x 5 matrix
  complex(dp), public :: complex_full_rank_matrix_5x5(5, 5) = transpose(reshape(cmplx([-44.52110712_dp, -48.17481048_dp,    -19.66887_dp,   2.95379091_dp, &
                                                                                       -32.48344889_dp, &
                                                                                       -48.40959198_dp,   36.6255884_dp, -40.95025262_dp,  -9.47597016_dp, &
                                                                                        21.96138633_dp, &
                                                                                       -23.89308839_dp, -46.60346009_dp,  29.51933189_dp,  46.80034723_dp, &
                                                                                       -26.99063883_dp, &
                                                                                       -28.34722129_dp, -15.28477719_dp,  -9.98834923_dp,  19.58476894_dp, &
                                                                                        -45.2623859_dp , &
                                                                                        25.18235434_dp,  45.25070169_dp,  39.07849915_dp,  15.93610182_dp, &
                                                                                        -33.5371513_dp ], &
                                                                                      
                                                                                      [ 14.86243624_dp,  30.70582514_dp,  44.7784307_dp, -46.66873007_dp, &
                                                                                       -24.13516175_dp, &
                                                                                         32.3031061_dp,  33.04884469_dp, 20.86135191_dp, -42.08274516_dp, &
                                                                                        48.39385979_dp, &
                                                                                        11.45829782_dp,  36.55565904_dp, 34.21616482_dp,  29.01076295_dp, &
                                                                                       -40.37086456_dp, &
                                                                                        32.39639731_dp, -45.57484939_dp, 27.23970665_dp,  -3.86243919_dp, &
                                                                                        22.13109274_dp, &
                                                                                          5.8363281_dp,   9.02281896_dp,  -9.2413642_dp, -18.39173734_dp, &
                                                                                       -27.88252968_dp], kind=dp), [5, 5]))

  !> Real 5 x 7 matrix with rank 4 
  complex(dp), public :: complex_rank_4_matrix_5x7(5, 7) = transpose(reshape(cmplx([ -18.7814739_dp, -35.75636149_dp, -25.89864277_dp,  11.04508789_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                     33.75023754_dp,  49.97417692_dp,  -0.48275953_dp,  43.98236553_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                     21.48351648_dp,   1.63997878_dp,  16.31513844_dp,  29.47589814_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                     40.49274295_dp, -28.86521605_dp, -49.94130354_dp,  48.91554628_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp], &

                                                                                   [-28.97154895_dp, -30.01164974_dp,  37.37357186_dp,  15.58097948_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                     -33.8514265_dp,  34.40405017_dp, -44.22889464_dp,  41.54835073_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                     27.21430366_dp,  34.61978756_dp,  32.66960231_dp, -36.56589004_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                    -43.59929192_dp,  34.88338453_dp,  -4.6433614_dp,   9.07841908_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp], kind=dp), [7, 5]))
                                                                                             
  !> Real 7 x 5 matrix with rank 4 
  complex(dp), public :: complex_rank_4_matrix_7x5(7, 5) = transpose(reshape(cmplx([ -18.7814739_dp, -35.75636149_dp, -25.89864277_dp,  11.04508789_dp, &
                                                                                             0.0_dp, &
                                                                                     33.75023754_dp,  49.97417692_dp,  -0.48275953_dp,  43.98236553_dp, &
                                                                                             0.0_dp, &
                                                                                     21.48351648_dp,   1.63997878_dp,  16.31513844_dp,  29.47589814_dp, &
                                                                                             0.0_dp, &
                                                                                     40.49274295_dp, -28.86521605_dp, -49.94130354_dp,  48.91554628_dp, &
                                                                                             0.0_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                             0.0_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                             0.0_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                             0.0_dp], &

                                                                                   [-28.97154895_dp, -30.01164974_dp,  37.37357186_dp,  15.58097948_dp, &
                                                                                             0.0_dp, &
                                                                                     -33.8514265_dp,  34.40405017_dp, -44.22889464_dp,  41.54835073_dp, &
                                                                                             0.0_dp, &
                                                                                     27.21430366_dp,  34.61978756_dp,  32.66960231_dp, -36.56589004_dp, &
                                                                                             0.0_dp, &
                                                                                    -43.59929192_dp,  34.88338453_dp,  -4.6433614_dp,   9.07841908_dp, &
                                                                                             0.0_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                             0.0_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                             0.0_dp, &
                                                                                             0.0_dp,          0.0_dp,          0.0_dp,          0.0_dp, &
                                                                                             0.0_dp], kind=dp), [5, 7]))

                                                                                
        !> Fills a given array with (arbitrary) numbers. Helpful for I/O tests, 
        !> where the acutal content of the array is not relevant.
        interface fill_array
                procedure ::  fill_real_array_rank0,&
                        fill_real_array_rank1,& 
                        fill_real_array_rank2,&
                        fill_real_array_rank3,&
                        fill_complex_array_rank2,&
                        fill_complex_array_rank3,&
                        fill_complex_array_rank4                
        end interface
        
        
 contains 
                                                                                        
        !> Given its index, this function maps
        !> one element of a real double-precision
        !> array \( \mathbf{a} \) of rank 1 according to the arbitrary expression
        !> \[
        !>    a_{m} = 1.5m   
        !>                                      \].
        pure function value_map_real_rank1(m) result(b)
                integer(sp), intent(in) :: m 
                real(dp) :: b
                b = 1.5_dp * m
        end function

        !> Given its indices, this function maps
        !> one element of a real double-precision
        !> array \( \mathbf{a} \) of rank 2 according to the arbitrary expression
        !> \[
        !>    a_{mn} = m + 2n 
        !>                                      \].
        pure function value_map_real_rank2(m, n) result(b)
                integer(sp), intent(in) :: m, n
                real(dp) :: b
                b = 1_dp*m + 2_dp*n
        end function

        !> Given its indices, this function maps
        !> one element of a real double-precision 
        !> array \( \mathbf{a} \) of rank 3 according to the arbitrary expression
        !> \[
        !>    a_{mnl} = m + 2(n+l) 
        !>                                      \].
        pure function value_map_real_rank3(m, n, l) result(b)
                integer(sp), intent(in) :: m, n, l
                real(dp) :: b
                b = 1._dp*m + 2_dp*(n+l)
        end function
       
        !> Given its indices, this function maps
        !> one element of a complex double-precision
        !> array \( \mathbf{a} \) of rank 2 according to the arbitrary expression
        !> \[
        !>    a_{mn} = m + 2n \cdot i   
        !>                                      \].
        pure function value_map_complex_rank2(m, n) result(b)
                integer(sp), intent(in) :: m, n
                complex(dp) :: b
                b = cmplx(m, 2*n, dp)
        end function

        !> Given its indices, this function maps
        !> one element of a complex double-precision 
        !> array \( \mathbf{a} \) of rank 3 according to the arbitrary expression
        !> \[
        !>    a_{mnl} = m + 2(n+l) \cdot i   
        !>                                      \].
        pure function value_map_complex_rank3(m, n, l) result(b)
                integer(sp), intent(in) :: m, n, l
                complex(dp) :: b
                b = cmplx(m, 2*(n+l), dp)
        end function

        !> Given its indices, this function maps
        !> one element of a complex double-precision 
        !> array \( \mathbf{a} \) of rank 4 according to the arbitrary expression
        !> \[
        !>    a_{mnlk} = m + l + 2nk \cdot i   
        !>                                      \].
        pure function value_map_complex_rank4(m, n, l, k) result(b)
                integer(sp), intent(in) :: m, n, l, k
                complex(dp) :: b
                
                b = cmplx(m + l, 2*n*k, dp)
        end function


        !> Fills a real double-precision scalar \( \mathbf{a} \)
        !> according to a value map.
        subroutine fill_real_array_rank0(array, value_map)
                real(dp), intent(out) :: array

                interface
                pure function value_map(a) result(b)
                        use precision, only: dp
                        integer, intent(in) :: a
                        real(dp) :: b
                end function value_map
                end interface

                array = value_map(1)

        end subroutine fill_real_array_rank0
        
        !> Fills a real double-precision array \( \mathbf{a} \) of rank 1
        !> according to a value map.
        subroutine fill_real_array_rank1(array, value_map)
                real(dp), intent(out) :: array(:)
                integer(sp) :: n

                interface
                        pure function value_map(a) result(b)
                                use precision, only: dp
                                integer, intent(in) :: a
                                real(dp) :: b
                        end function value_map
                end interface

                do n = 1, size(array, dim=1)
                        array(n) = value_map(n)
                end do

        end subroutine fill_real_array_rank1

        !> Fills a complex double-precision array \( \mathbf{a} \) of rank 2
        !> according to a value map.
        subroutine fill_real_array_rank2(array, value_map)
                real(dp), intent(out) :: array(:, :)
                integer(sp) :: n, m

                interface
                        pure function value_map(a, b) result(c)
                                use precision, only: dp
                                integer, intent(in) :: a, b
                                real(dp) :: c
                        end function value_map
                end interface


                do m = 1, size(array, dim=1)
                        do n = 1, size(array, dim=2)
                                array(m, n) = value_map(m, n)
                        end do
                end do

        end subroutine fill_real_array_rank2


        !> Fills a complex double-precision array \( \mathbf{a} \) of rank 3
        !> according to a value map.
        subroutine fill_real_array_rank3(array, value_map)

                real(dp), intent(out) :: array(:, :, :)
                integer(sp) :: m, n, l
                interface
                        pure function value_map(a, b, c) result(d)
                                use precision, only: dp
                                integer, intent(in) :: a, b, c
                                real(dp) :: d
                        end function value_map
                end interface

                do m = 1, size(array, dim=1)
                        do n = 1, size(array, dim=2)
                                do l = 1, size(array, dim=3)
                                        array(m, n, l) = value_map(m, n, l)
                                end do
                        end do
                end do

        end subroutine fill_real_array_rank3
                                                    
        
        !> Fills a complex double-precision array \( \mathbf{a} \) of rank 2
        !> according to a value map.
        subroutine fill_complex_array_rank2(array, value_map)
                complex(dp), intent(out) :: array(:, :)
                integer(sp) :: n, m

                interface
                        pure function value_map(a, b) result(c)
                                use precision, only: dp
                                integer, intent(in) :: a, b
                                complex(dp) :: c
                        end function value_map
                end interface


                do m = 1, size(array, dim=1)
                        do n = 1, size(array, dim=2)
                                array(m, n) = value_map(m, n)
                        end do
                end do

        end subroutine fill_complex_array_rank2



        !> Fills a complex double-precision array \( \mathbf{a} \) of rank 3
        !> according to a value map.
        subroutine fill_complex_array_rank3(array, value_map)

                complex(dp), intent(out) :: array(:, :, :)
                integer(sp) :: m, n, l
                interface
                        pure function value_map(a, b, c) result(d)
                                use precision, only: dp
                                integer, intent(in) :: a, b, c
                                complex(dp) :: d
                        end function value_map
                end interface

                do m = 1, size(array, dim=1)
                        do n = 1, size(array, dim=2)
                                do l = 1, size(array, dim=3)
                                        array(m, n, l) = value_map(m, n, l)
                                end do
                        end do
                end do

        end subroutine fill_complex_array_rank3


        !TODO(Max) Issue 139. Fix Order of Loops in mock_arrays
        !> Fills a complex double-precision array \( \mathbf{a} \) of rank 4
        !> according to a value map.
        subroutine fill_complex_array_rank4(array, value_map)

                complex(dp), intent(out) :: array(:, :, :, :)
                integer(sp) :: m, n, l, k

                !> Allows to 
                interface
                        pure function value_map(a, b, c, d) result(e)
                                use precision, only: dp
                                integer, intent(in) :: a, b, c, d
                                complex(dp) :: e
                        end function value_map
                end interface

                do m = 1, size(array, dim=1)
                        do n = 1, size(array, dim=2)
                                do l = 1, size(array, dim=3)
                                        do k = 1, size(array, dim=4)
                                                array(m, n, l, k) = &
                                                 value_map(m, n, l, k)
                                        end do
                                end do
                        end do
                end do

        end subroutine fill_complex_array_rank4

end module mock_arrays

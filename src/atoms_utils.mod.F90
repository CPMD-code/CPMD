MODULE atoms_utils
  USE adat,                            ONLY: atwt,&
                                             covrad,&
                                             defrag,&
                                             elem,&
                                             mm_core_raggio,&
                                             nelcon
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: atoms

CONTAINS

  ! ==================================================================
  SUBROUTINE atoms
    ! ==--------------------------------------------------------------==
    ! ==  COVALENT RADII                                              ==
    ! ==  ATOMIC MASSES                                               ==
    ! ==  RAGGIO                                                      ==
    ! ==  ELEMENTS                                                    ==
    ! ==  ELECTRON CONFIGURATIONS                                     ==
    ! ==  mm_core_raggio2 in au
    ! ==  covrad2 in angstrom
    ! ==--------------------------------------------------------------==
    INTEGER, DIMENSION(4, 7, 99), PARAMETER :: nelcon2 = RESHAPE((/1,0,0,0,0,0&
      ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0&
      ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0&
      ,0,0,0,0,0,0,2,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0&
      ,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,2,0,0,0,0&
      ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,3,0,0,0,0,0,0,0,0,0,0,0,0&
      ,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0&
      ,0,0,2,0,0,0,2,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6&
      ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,1,0,0,0,0,0&
      ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0&
      ,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0&
      ,0,0,2,6,0,0,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,3&
      ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,4,0,0,0,0,0,0,0,0&
      ,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0&
      ,0,0,2,0,0,0,2,6,0,0,2,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6&
      ,0,0,2,6,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,0,0,2,0&
      ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,1,0,2,0,0,0,0,0,0,0,0,0&
      ,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,2,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0&
      ,0,0,2,6,0,0,2,6,3,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6&
      ,5,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,5,0,2,0,0,0,0,0&
      ,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,6,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0&
      ,0,0,2,0,0,0,2,6,0,0,2,6,7,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6&
      ,0,0,2,6,8,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,1,&
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,0,0,0,0,0,0,0,0&
      ,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
      2,0,0,0,2,6,0,0,2,6,10,0,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0&
      ,2,6,10,0,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,4,0&
      ,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,5,0,0,0,0,0,0,0,0,&
      0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0&
      ,0,0,2,6,0,0,2,6,10,0,2,6,0,0,1,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,&
      6,10,0,2,6,0,0,2,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,1,0,&
      2,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,2,0,2,0,0,0,0,0,0,0&
      ,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,4,0,1,0,0,0,0,0,0,0,0,0,0,0,2,0,0,&
      0,2,6,0,0,2,6,10,0,2,6,5,0,1,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,&
      10,0,2,6,6,0,1,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,7,0,1,&
      0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,8,0,1,0,0,0,0,0,0,0,0&
      ,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0&
      ,2,6,0,0,2,6,10,0,2,6,10,0,1,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,&
      10,0,2,6,10,0,2,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,0,&
      2,1,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,0,2,2,0,0,0,0,0,&
      0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,0,2,3,0,0,0,0,0,0,0,0,0,0,2,0,&
      0,0,2,6,0,0,2,6,10,0,2,6,10,0,2,4,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,&
      6,10,0,2,6,10,0,2,5,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,&
      0,2,6,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,0,2,6,0,0,1,0,&
      0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,0,2,6,0,0,2,0,0,0,0,0,0,0,2,&
      0,0,0,2,6,0,0,2,6,10,0,2,6,10,0,2,6,1,0,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,&
      2,6,10,0,2,6,10,2,2,6,0,0,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,&
      10,3,2,6,0,0,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,4,2,6,0,0,2&
      ,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,5,2,6,0,0,2,0,0,0,0,0,0,0&
      ,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,6,2,6,0,0,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0&
      ,0,2,6,10,0,2,6,10,7,2,6,0,0,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,&
      6,10,7,2,6,1,0,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,9,2,6,0,0&
      ,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,10,2,6,0,0,2,0,0,0,0,0,&
      0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,11,2,6,0,0,2,0,0,0,0,0,0,0,2,0,0,0,2&
      ,6,0,0,2,6,10,0,2,6,10,12,2,6,0,0,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10&
      ,0,2,6,10,13,2,6,0,0,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,&
      2,6,0,0,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,1,0,2,0,0&
      ,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,2,0,2,0,0,0,0,0,0,0,2,&
      0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,3,0,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0&
      ,2,6,10,0,2,6,10,14,2,6,4,0,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6&
      ,10,14,2,6,5,0,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,6,&
      0,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,7,0,2,0,0,0,0,0&
      ,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,9,0,1,0,0,0,0,0,0,0,2,0,0,0,&
      2,6,0,0,2,6,10,0,2,6,10,14,2,6,10,0,1,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,&
      10,0,2,6,10,14,2,6,10,0,2,0,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,&
      14,2,6,10,0,2,1,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,10,0,&
      2,2,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,10,0,2,3,0,0,0,0,&
      0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,10,0,2,4,0,0,0,0,0,0,2,0,0,0,&
      2,6,0,0,2,6,10,0,2,6,10,14,2,6,10,0,2,5,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,&
      10,0,2,6,10,14,2,6,10,0,2,6,0,0,0,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,&
      14,2,6,10,0,2,6,0,0,1,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,10,0,&
      2,6,0,0,2,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,10,0,2,6,1,0,2,0,&
      0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,10,0,2,6,2,0,2,0,0,0,2,0,0,0,&
      2,6,0,0,2,6,10,0,2,6,10,14,2,6,10,2,2,6,1,0,2,0,0,0,2,0,0,0,2,6,0,0,2,6,&
      10,0,2,6,10,14,2,6,10,3,2,6,1,0,2,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,&
      14,2,6,10,4,2,6,1,0,2,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,10,6,&
      2,6,0,0,2,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,10,7,2,6,0,0,2,0,&
      0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10,14,2,6,10,7,2,6,1,0,2,0,0,0,2,0,0,0,&
      2,6,0,0,2,6,10,0,2,6,10,14,2,6,10,9,2,6,0,0,2,0,0,0,2,0,0,0,2,6,0,0,2,6,&
      10,0,2,6,10,14,2,6,10,10,2,6,0,0,2,0,0,0,2,0,0,0,2,6,0,0,2,6,10,0,2,6,10&
      ,14,2,6,10,11,2,6,0,0,2,0,0,0/),(/4,7,99/))
    REAL(real_8), DIMENSION(99), PARAMETER :: atwt2 = (/1.00797_real_8,  &
      4.0026_real_8,    6.939_real_8,   9.0122_real_8,   10.811_real_8,&
      12.01115_real_8,14.0067_real_8,  15.9994_real_8,  18.9984_real_8,   &
      20.183_real_8,22.9898_real_8,  24.312_real_8,  26.9815_real_8,   &
      28.086_real_8,  30.9738_real_8,32.064_real_8,   35.453_real_8,   &
      39.948_real_8,   39.102_real_8,   40.080_real_8,44.956_real_8,   &
      47.900_real_8,   50.942_real_8,   51.996_real_8,   54.938_real_8,&
      55.847_real_8,   58.933_real_8,   58.710_real_8,   63.540_real_8,   &
      65.370_real_8,69.720_real_8,   72.590_real_8,   74.922_real_8,   &
      78.960_real_8,   79.909_real_8,83.800_real_8,   85.470_real_8,   &
      87.620_real_8,   88.905_real_8,   91.220_real_8,92.906_real_8,   &
      95.940_real_8,   98.000_real_8,  101.070_real_8,  102.905_real_8,&
      106.400_real_8,  107.870_real_8,  112.400_real_8,  114.820_real_8,  &
      118.690_real_8,121.750_real_8,  127.600_real_8,  126.904_real_8,  &
      131.300_real_8,  132.905_real_8,137.340_real_8,  138.910_real_8,  &
      140.120_real_8,  140.907_real_8,  144.240_real_8,147.000_real_8,  &
      150.350_real_8,  151.960_real_8,  157.250_real_8,  158.924_real_8,&
      162.500_real_8,  164.930_real_8,  167.260_real_8,  168.934_real_8,  &
      173.040_real_8,174.970_real_8,  178.490_real_8,  180.948_real_8,  &
      183.850_real_8,  186.200_real_8,190.200_real_8,  192.200_real_8,  &
      195.090_real_8,  196.967_real_8,  200.590_real_8,204.370_real_8,  &
      207.190_real_8,  208.980_real_8,  210.000_real_8,  210.000_real_8,&
      222.000_real_8, 250.0_real_8, 250.0_real_8,250.0_real_8,250.0_real_8,&
      250.0_real_8,250.0_real_8,250.0_real_8,250.0_real_8,250.0_real_8,&
      250.0_real_8,250.0_real_8,250.0_real_8,250.0_real_8/)
    REAL(real_8), DIMENSION(99), PARAMETER :: covrad2 = (/0.32_real_8,  &
      0.93_real_8,  1.23_real_8,  0.90_real_8,  0.82_real_8,  0.77_real_8,&
      0.75_real_8,  0.73_real_8,  0.72_real_8,  0.71_real_8,  1.54_real_8,  &
      1.36_real_8,1.18_real_8,  1.11_real_8,  1.06_real_8,  1.02_real_8,  &
      0.99_real_8,  0.98_real_8,2.03_real_8,  1.74_real_8,  1.44_real_8,  &
      1.32_real_8,  1.22_real_8,  1.18_real_8,1.17_real_8,  1.17_real_8,  &
      1.16_real_8,  1.15_real_8,  1.17_real_8,  1.25_real_8,1.26_real_8,  &
      1.22_real_8,  1.20_real_8,  1.16_real_8,  1.14_real_8,  1.12_real_8,&
      2.16_real_8,  1.91_real_8,  1.62_real_8,  1.45_real_8,  1.34_real_8,  &
      1.30_real_8,1.27_real_8,  1.25_real_8,  1.25_real_8,  1.28_real_8,  &
      1.34_real_8,  1.48_real_8,1.44_real_8,  1.41_real_8,  1.40_real_8,  &
      1.36_real_8,  1.33_real_8,  1.31_real_8,2.35_real_8,  1.98_real_8,  &
      1.69_real_8,  1.65_real_8,  1.65_real_8,  1.64_real_8,1.63_real_8,  &
      1.62_real_8,  1.85_real_8,  1.61_real_8,  1.59_real_8,  1.59_real_8,&
      1.58_real_8,  1.57_real_8,  1.56_real_8,  1.56_real_8,  1.56_real_8,  &
      1.44_real_8,1.34_real_8,  1.30_real_8,  1.28_real_8,  1.26_real_8,  &
      1.27_real_8,  1.30_real_8,1.34_real_8,  1.49_real_8,  1.48_real_8,  &
      1.47_real_8,  1.46_real_8,  1.46_real_8,1.45_real_8,  0.00_real_8,  &
      0.00_real_8,  0.00_real_8,  0.00_real_8,  1.65_real_8,0.00_real_8,  &
      1.42_real_8,  0.00_real_8,  0.00_real_8,  0.00_real_8,  0.00_real_8,&
      0.00_real_8,  0.00_real_8,  0.00_real_8/)
    REAL(real_8), DIMENSION(99), PARAMETER :: defrag2 = (/1.2_real_8, &
      1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, &
      1.2_real_8, 1.2_real_8,1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, &
      1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8,1.2_real_8, &
      1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, &
      1.2_real_8, 1.2_real_8,1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, &
      1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8,1.2_real_8, &
      1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, &
      1.2_real_8, 1.2_real_8,1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, &
      1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8,1.2_real_8, &
      1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, &
      1.2_real_8, 1.2_real_8,1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, &
      1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8,1.2_real_8, &
      1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, &
      1.2_real_8, 1.2_real_8,1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, &
      1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8,1.2_real_8, &
      1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, 1.2_real_8, &
      1.2_real_8, 1.2_real_8/)
    REAL(real_8), DIMENSION(99), PARAMETER :: mm_core_raggio2 = (/0.50_real_8,&
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8,  0.50_real_8, 0.50_real_8,  &
      0.50_real_8,  0.50_real_8,  0.50_real_8/)
    INTEGER                                  :: i, j
    CHARACTER(len=2), DIMENSION(99), PARAMETER :: el2 = (/' H','He','Li','Be',&
      ' B',' C',' N',' O',' F','Ne','Na','Mg','Al','Si',' P',' S','Cl','Ar',&
      ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge',&
      'As','Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','TC','Ru','Rh','Pd',&
      'Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd',&
      'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',' W',&
      'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',&
      'Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es'/)

! Extended Huckel parameters taken from
! the table collected by S. Alvarez, Dept. de Quimica Inorganica,
! Univ. de Barcelona, Spain
! AK 2005/03/28: unsused
! DATA (EHTPAR2(I,1),I=1,69) 
! *                       /-13.6_real_8,-23.4_real_8, -5.4_real_8,-10.0_real_8,-15.2_real_8,
! *-21.4_real_8,-26.0_real_8,-32.3_real_8,-40.0_real_8,-43.2_real_8, -5.1_real_8, -9.0_real_8,-12.3_real_8,
! *-17.3_real_8,-18.6_real_8,-20.0_real_8,-26.3_real_8, -0.0_real_8, -4.34_real_8, -7.0_real_8, -8.87_real_8,
! * -8.97_real_8, -8.81_real_8, -8.66_real_8, -9.75_real_8, -9.1_real_8, -9.21_real_8, -9.17_real_8,
! *-11.4_real_8,-12.41_real_8,-14.58_real_8,-16.0_real_8,-16.22_real_8,-20.5_real_8,-22.07_real_8,
! * -0.0_real_8, -4.18_real_8, -6.62_real_8, -0.0_real_8, -8.0_real_8,-10.1_real_8, -8.34_real_8,
! *-10.07_real_8,-10.4_real_8, -8.09_real_8, -7.32_real_8, -0.0_real_8, -0.0_real_8,-12.6_real_8,
! *-16.16_real_8,-18.8_real_8,-20.8_real_8,-18.0_real_8, -0.0_real_8, -3.88_real_8, -0.0_real_8, 
! *-7.67_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -4.86_real_8, -0.0_real_8, 
! *-0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8/
! DATA (EHTPAR2(I,1),I=70,99) 
! *            / -5.35_real_8, -6.05_real_8, -0.0_real_8,-10.1_real_8,-8.26_real_8, 
! * -9.36_real_8, -8.17_real_8,-11.36_real_8, -9.077_real_8,-10.92_real_8,-13.68_real_8,-11.6_real_8,
! *-15.7_real_8,-15.19_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -5.39_real_8, -0.0_real_8, -5.5_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -0.0_real_8/
! DATA (EHTPAR2(I,2),I=1,69) 
! *                       / -0.0_real_8, -0.0_real_8, -3.500, -6.0_real_8, -8.5_real_8,
! *-11.4_real_8,-13.4_real_8,-14.8_real_8,-18.1_real_8,-20.0_real_8, -3.0_real_8, -4.5_real_8, -6.5_real_8,
! * -9.2_real_8,-14.0_real_8,-11.0_real_8,-14.2_real_8, -0.0_real_8, -2.73_real_8, -4.0_real_8, -2.75_real_8,
! * -5.44_real_8, -5.52_real_8, -5.24_real_8, -5.89_real_8, -5.32_real_8, -5.29_real_8, -5.15_real_8,
! * -6.06_real_8, -6.53_real_8, -6.75_real_8, -9.0_real_8,-12.16_real_8,-14.4_real_8,-13.1_real_8,
! * -0.0_real_8, -2.6_real_8, -3.92_real_8, -0.0_real_8, -5.4_real_8, -6.86_real_8, -5.24_real_8,
! * -5.4_real_8, -6.87_real_8, -4.57_real_8, -3.75_real_8, -0.0_real_8, -0.0_real_8, -6.19_real_8,
! * -8.32_real_8,-11.7_real_8,-13.2_real_8,-12.7_real_8,
! * -0.0_real_8, -2.49_real_8, -0.0_real_8, -5.01_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -4.86_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8/
! DATA (EHTPAR2(I,2),I=70,99) 
! *                       / -5.35_real_8, -6.05_real_8, -0.0_real_8, -6.86_real_8,
! * -5.17_real_8, -5.96_real_8, -4.81_real_8, -4.5_real_8, -5.475_real_8, -5.55_real_8, -8.47_real_8,
! * -5.8_real_8, -8.0_real_8, -7.79_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -5.39_real_8, -0.0_real_8, -5.5_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -0.0_real_8, -0.0_real_8/
! DATA (EHTPAR2(I,3),I=1,69) 
! *                       / -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -8.510,
! *-10.81_real_8,-11.0_real_8,-11.22_real_8,-11.67_real_8,-12.6_real_8,-13.18_real_8,-13.49_real_8,
! *-14.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! *-0.0_real_8, -0.0_real_8, -0.0_real_8,-10.2_real_8,-12.1_real_8,-10.5_real_8,-12.82_real_8,-14.9_real_8,
! *-12.5_real_8,-12.02_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! *-0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -8.21_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! *-0.0_real_8, -6.06_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! *-0.0_real_8/
! DATA (EHTPAR2(I,3),I=70,99) 
! *                       / -5.21_real_8, -5.12_real_8, -0.0_real_8,-12.1_real_8,
! *-10.37_real_8,-12.66_real_8,-11.84_real_8,-12.17_real_8,-12.59_real_8,-15.076_real_8,-17.5_real_8, 
! *-0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! *-0.0_real_8, -10.11_real_8, -0.0_real_8, -9.19_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! *-0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8/
! DATA (EHTPAR2(I,4),I=1,69) 
! *                       / -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! *-11.28_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8/
! DATA (EHTPAR2(I,4),I=70,99) 
! *                       /-13.86_real_8,-22.4_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -9.64_real_8, -0.0_real_8,-10.62_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8, -0.0_real_8,
! * -0.0_real_8, -0.0_real_8/
! ==--------------------------------------------------------------==

    DO i=1,99
       elem%el(i)=el2(i)
       covrad(i)=covrad2(i)
       mm_core_raggio(i)=mm_core_raggio2(i)
       atwt(i)=atwt2(i)
       defrag(i)=defrag2(i)
       DO j=1,7
          nelcon(1,j,i)=nelcon2(1,j,i)
          nelcon(2,j,i)=nelcon2(2,j,i)
          nelcon(3,j,i)=nelcon2(3,j,i)
          nelcon(4,j,i)=nelcon2(4,j,i)
       ENDDO
    ENDDO
  END SUBROUTINE atoms
  ! ==================================================================

END MODULE atoms_utils

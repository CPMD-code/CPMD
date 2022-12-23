MODULE fftchk_utils
  USE cp_cuda_types,                   ONLY: cp_cuda_env
  USE error_handling,                  ONLY: stopgm
  USE parac,                           ONLY: paral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fftchk
  !public :: fftchk_fftw
  !public :: fftchk_def
  !public :: fftchk_essl

CONTAINS

  ! ==================================================================
  FUNCTION fftchk(m,n)
    ! ==--------------------------------------------------------------==
    ! N < 0 : take the next smaller one
    ! N = 1 : take the next bigger one
    ! N = 2 : take the next bigger even one
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: m, n, fftchk

    INTEGER                                  :: fftchk_cuda_m, fftchk_cuda_p
    LOGICAL                                  :: found

    fftchk=HUGE(0)
    fftchk_cuda_p = HUGE(0)
    fftchk_cuda_m = HUGE(0)

#if defined(__HAS_FFT_DEFAULT)
    fftchk = fftchk_def(m,n)
#elif defined(__HAS_FFT_ESSL)
    fftchk = fftchk_essl(m,n)
#elif defined(__HAS_FFT_HP)
    fftchk = fftchk_def(m,n)
#elif defined(__HAS_FFT_FFTW3)
    fftchk = fftchk_fftw(m,n)
#else
    CALL stopgm("FFTCHK","TABLE OF ROOTS NOT AVAILABLE",&
         __LINE__,__FILE__)
#endif

    found = .TRUE.
    IF( cp_cuda_env%use_fft ) THEN
       fftchk_cuda_m = FFTCHK_CUDA(fftchk,-1)
       fftchk_cuda_p = FFTCHK_CUDA(fftchk, 1)
       found = fftchk_cuda_m == fftchk .OR. fftchk_cuda_p == fftchk       
    ENDIF

    !vw for the moment we stop, we could have interate till some match size is found
    IF( .NOT. found ) CALL stopgm("FFTCHK","no easy match between GPU and CPU input size",&
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
  END FUNCTION fftchk
  ! ==================================================================
  FUNCTION fftchk_fftw(m,n)
    ! ==--------------------------------------------------------------==
    ! N < 0 : take the next smaller one
    ! N = 1 : take the next bigger one
    ! N = 2 : take the next bigger even one
    ! ==--------------------------------------------------------------==
    ! ==================================================================
    ! ==   The LFT table list of acceptable values for                ==
    ! ==   the transform lengths in the FFT (roots 2, 3, 5, 7)        ==
    ! ==================================================================
    INTEGER                                  :: m, n, fftchk_fftw

    INTEGER, PARAMETER                       :: nmx = 316

    INTEGER                                  :: i, m1
    INTEGER, DIMENSION(nmx) :: lft = (/   2,   3,   4,   5,   6,   7,   8,   9&
      ,  10,  12,14,  15,  16,  18,  20,  21,  24,  25,  27,  28,30,  32,  35,&
      36,  40,  42,  45,  48,  49,  50,54,  56,  60,  63,  64,  70,  72,  75, &
      80,  81,84,  90,  96,  98, 100, 105, 108, 112, 120, 125,126, 128, 135, &
      140, 144, 147, 150, 160, 162, 168,175, 180, 189, 192, 196, 200, 210, 216&
      , 224, 225,240, 243, 245, 250, 252, 256, 270, 280, 288, 294,300, 315, &
      320, 324, 336, 343, 350, 360, 375, 378,384, 392, 400, 405, 420, 432, 441&
      , 448, 450, 480,486, 490, 500, 504, 512, 525, 540, 560, 567, 576,588, &
      600, 625, 630, 640, 648, 672, 675, 686, 700,720, 729, 735, 750, 756, 768&
      , 784, 800, 810, 840,864, 875, 882, 896, 900, 945, 960, 972, 980,1000,&
      1008,1024,1029,1050,1080,1120,1125,1134,1152,1176,1200,1215,1225,1250,&
      1260,1280,1296,1323,1344,1350,1372,1400,1440,1458,1470,1500,1512,1536,&
      1568,1575,1600,1620,1680,1701,1715,1728,1750,1764,1792,1800,1875,1890,&
      1920,1944,1960,2000,2016,2025,2048,2058,2100,2160,2187,2205,2240,2250,&
      2268,2304,2352,2400,2401,2430,2450,2500,2520,2560,2592,2625,2646,2688,&
      2700,2744,2800,2835,2880,2916,2940,3000,3024,3072,3087,3125,3136,3150,&
      3200,3240,3360,3375,3402,3430,3456,3500,3528,3584,3600,3645,3675,3750,&
      3780,3840,3888,3920,3969,4000,4032,4050,4096,4116,4200,4320,4374,4375,&
      4410,4480,4500,4536,4608,4704,4725,4800,4802,4860,4900,5000,5040,5103,&
      5120,5145,5184,5250,5292,5376,5400,5488,5600,5625,5670,5760,5832,5880,&
      6000,6048,6075,6125,6144,6174,6250,6272,6300,6400,6480,6561,6615,6720,&
      6750,6804,6860,6912,7000,7056,7168,7200,7203,7290,7350,7500,7560,7680,&
      7776,7840,7875,7938,8000,8064,8100,8192  /)

    fftchk_fftw=0
    m1 = m
    DO i=1,nmx
       IF (lft(i).GE.m1) THEN
          IF (n.LT.0) THEN
             m1=lft(i-1)
             GOTO 10
          ELSEIF (n.EQ.1) THEN
             m1=lft(i)
             GOTO 10
          ELSEIF (n.EQ.2) THEN
             IF (MOD(lft(i),2).EQ.0) THEN
                m1=lft(i)
                GOTO 10
             ENDIF
          ENDIF
       ENDIF
    ENDDO
    GOTO 999
10  CONTINUE
    fftchk_fftw=m1
    RETURN
    ! ==--------------------------------------------------------------==
999 CONTINUE
    IF (paral%io_parent) THEN
       WRITE(6,*) 'FUNCTION FFTCHK   '
       WRITE(6,*) ' THE MINIMAL MESH SIZE IS LARGER THAN ',lft(nmx),&
            ' M=',M
       WRITE(6,*) ' ADD LARGER MESH VALUES IN FFTCHK'
    ENDIF
    CALL stopgm('FFTCHK ',' ',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END FUNCTION fftchk_fftw
  ! ==================================================================
  FUNCTION fftchk_def(m,n)
    ! ==--------------------------------------------------------------==
    ! N < 0 : take the next smaller one
    ! N = 1 : take the next bigger one
    ! N = 2 : take the next bigger even one
    ! ==--------------------------------------------------------------==
    ! ==================================================================
    ! ==   The following table list of acceptable values for          ==
    ! ==   the transform lengths in the FFT (roots 2, 3 and 5)        ==
    ! ==================================================================
    INTEGER                                  :: m, n, fftchk_def

    INTEGER, PARAMETER                       :: nmx = 100

    INTEGER                                  :: i, m1
    INTEGER, DIMENSION(nmx) :: lft = (/   3,   4,   5,   6,   8,   9,  12,  15&
      ,  16,  18,20,  24,  25,  27,  30,  32,  36,  40,  45,  48,54,  60,  64,&
      72,  75,  80,  81,  90,  96, 100,108, 120, 125, 128, 135, 144, 150, 160,&
      162, 180,192, 200, 216, 225, 240, 243, 256, 270, 288, 300,320, 324, 360,&
      375, 384, 400, 405, 432, 450, 480,486, 500, 512, 540, 576, 600, 625, 640&
      , 648, 675,720, 729, 750, 768, 800, 810, 864, 900, 960, 972,1000,1024,&
      1080,1152,1200,1280,1296,1350,1440,1458,1500,1536,1600,1620,1728,1800,&
      1920,1944,2000,2048 /)

    fftchk_def=0
    m1 = m
    DO i=1,nmx
       IF (lft(i).GE.m1) THEN
          IF (n.LT.0) THEN
             m1=lft(i-1)
             GOTO 10
          ELSEIF (n.EQ.1) THEN
             m1=lft(i)
             GOTO 10
          ELSEIF (n.EQ.2) THEN
             IF (MOD(lft(i),2).EQ.0) THEN
                m1=lft(i)
                GOTO 10
             ENDIF
          ENDIF
       ENDIF
    ENDDO
    GOTO 999
10  CONTINUE
    fftchk_def=m1
    RETURN
    ! ==--------------------------------------------------------------==
999 CONTINUE
    IF (paral%io_parent) THEN
       WRITE(6,*) 'FUNCTION FFTCHK   '
       WRITE(6,*) ' THE MINIMAL MESH SIZE IS LARGER THAN ',lft(nmx),&
            ' M=',M
       WRITE(6,*) ' ADD LARGER MESH VALUES IN FFTCHK'
    ENDIF
    CALL stopgm('FFTCHK ',' ',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END FUNCTION fftchk_def
  ! ==================================================================
  FUNCTION fftchk_essl(m,n)
    ! ==--------------------------------------------------------------==
    ! N < 0 : take the next smaller one
    ! N = 1 : take the next bigger one
    ! N = 2 : take the next bigger even one
    ! ==--------------------------------------------------------------==
    ! ==================================================================
    ! ==   The lft table list of acceptable values for                ==
    ! ==   the transform lengths in the FFT is taken from pag. 758    ==
    ! ==   of the ESSL manual (vol. 3)                                ==
    ! ==================================================================
    INTEGER                                  :: m, n, fftchk_essl

    INTEGER, PARAMETER                       :: nmx = 99

    INTEGER                                  :: i, m1
    INTEGER, DIMENSION(nmx) :: lft = (/   2,   4,   6,   8,  10,  12,  14,  16&
      ,  18,  20,22,  24,  28,  30,  32,  36,  40,  42,  44,  48,56,  60,  64,&
      66,  70,  72,  80,  84,  88,  90,96, 110, 112, 120, 126, 128, 132, 140, &
      144, 154,160, 168, 176, 180, 192, 198, 210, 220, 224, 240,252, 256, 264,&
      280, 288, 308, 320, 330, 336, 352,360, 384, 396, 420, 440, 448, 462, 480&
      , 504, 512,528, 560, 576, 616, 630, 640, 660, 672, 704, 720,768, 770, &
      792, 840, 880, 896, 924, 960, 990,1008,1024,1056,1120,1152,1232,1260,&
      1280,1320,1344/)

    fftchk_essl=0
    m1 = m
    DO i=1,nmx
       IF (lft(i).GE.m1) THEN
          IF (n.LT.0) THEN
             m1=lft(i-1)
             GOTO 10
          ELSEIF (n.EQ.1) THEN
             m1=lft(i)
             GOTO 10
          ELSEIF (n.EQ.2) THEN
             IF (MOD(lft(i),2).EQ.0) THEN
                m1=lft(i)
                GOTO 10
             ENDIF
          ENDIF
       ENDIF
    ENDDO
    GOTO 999
10  CONTINUE
    fftchk_essl=m1
    RETURN
    ! ==--------------------------------------------------------------==
999 CONTINUE
    IF (paral%io_parent) THEN
       WRITE(6,*) 'FUNCTION FFTCHK   '
       WRITE(6,*) ' THE MINIMAL MESH SIZE IS LARGER THAN ',lft(nmx),&
            ' M=',M
       WRITE(6,*) ' ADD LARGER MESH VALUES IN FFTCHK'
    ENDIF
    CALL stopgm('FFTCHK ',' ',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END FUNCTION fftchk_essl
  ! ==================================================================

  !     ==================================================================
  FUNCTION FFTCHK_CUDA(M,N)
    !     ==--------------------------------------------------------------==
    !     N < 0 : take the next smaller one
    !     N = 1 : take the next bigger one
    !     N = 2 : take the next bigger even one
    !     ==--------------------------------------------------------------==
    INTEGER                                  :: M, N, FFTCHK_CUDA

    INTEGER, PARAMETER                       :: nmx = 316

    INTEGER                                  :: I, M1
    INTEGER, DIMENSION(nmx) :: LFT = (/  2,   3,   4,   5,   6,   7,   8,   9,&
      10,  12,14,  15,  16,  18,  20,  21,  24,  25,  27,  28,30,  32,  35,  &
      36,  40,  42,  45,  48,  49,  50,54,  56,  60,  63,  64,  70,  72,  75, &
      80,  81,84,  90,  96,  98, 100, 105, 108, 112, 120, 125,126, 128, 135, &
      140, 144, 147, 150, 160, 162, 168,175, 180, 189, 192, 196, 200, 210, 216&
      , 224, 225,240, 243, 245, 250, 252, 256, 270, 280, 288, 294,300, 315, &
      320, 324, 336, 343, 350, 360, 375, 378,384, 392, 400, 405, 420, 432, 441&
      , 448, 450, 480,486, 490, 500, 504, 512, 525, 540, 560, 567, 576,588, &
      600, 625, 630, 640, 648, 672, 675, 686, 700,720, 729, 735, 750, 756, 768&
      , 784, 800, 810, 840,864, 875, 882, 896, 900, 945, 960, 972, 980,1000,&
      1008,1024,1029,1050,1080,1120,1125,1134,1152,1176,1200,1215,1225,1250,&
      1260,1280,1296,1323,1344,1350,1372,1400,1440,1458,1470,1500,1512,1536,&
      1568,1575,1600,1620,1680,1701,1715,1728,1750,1764,1792,1800,1875,1890,&
      1920,1944,1960,2000,2016,2025,2048,2058,2100,2160,2187,2205,2240,2250,&
      2268,2304,2352,2400,2401,2430,2450,2500,2520,2560,2592,2625,2646,2688,&
      2700,2744,2800,2835,2880,2916,2940,3000,3024,3072,3087,3125,3136,3150,&
      3200,3240,3360,3375,3402,3430,3456,3500,3528,3584,3600,3645,3675,3750,&
      3780,3840,3888,3920,3969,4000,4032,4050,4096,4116,4200,4320,4374,4375,&
      4410,4480,4500,4536,4608,4704,4725,4800,4802,4860,4900,5000,5040,5103,&
      5120,5145,5184,5250,5292,5376,5400,5488,5600,5625,5670,5760,5832,5880,&
      6000,6048,6075,6125,6144,6174,6250,6272,6300,6400,6480,6561,6615,6720,&
      6750,6804,6860,6912,7000,7056,7168,7200,7203,7290,7350,7500,7560,7680,&
      7776,7840,7875,7938,8000,8064,8100,8192  /)

!     Variables
!     ==================================================================
!     ==   The following table list of acceptable values for          ==
!     ==   the transform lengths in the FFT (roots 2, 3, 5, 7)        ==
!     ==   Note that CUFFT can also work with other lenghts           ==
!     ==================================================================
!     ==--------------------------------------------------------------==

    FFTCHK_CUDA=0
    M1 = M
    DO I=1,NMX
       IF(LFT(I).GE.M1) THEN
          IF(N.LT.0) THEN
             M1=LFT(I-1)
             GOTO 10
          ELSEIF(N.EQ.1) THEN
             M1=LFT(I)
             GOTO 10
          ELSEIF(N.EQ.2) THEN
             IF(MOD(LFT(I),2).EQ.0) THEN
                M1=LFT(I)
                GOTO 10
             ENDIF
          ENDIF
       ENDIF
    ENDDO
    GOTO 999
10  CONTINUE
    FFTCHK_CUDA=M1
    RETURN
    ! ==--------------------------------------------------------------==
999 CONTINUE
    IF(paral%IO_PARENT) WRITE(6,*) 'FUNCTION FFTCHK     '
    IF(paral%IO_PARENT) WRITE(6,*) ' THE MINIMAL MESH SIZE IS LARGER THAN ',LFT(NMX), ' M=',M
    IF(paral%IO_PARENT) WRITE(6,*) ' ADD LARGER MESH VALUES IN FFTCHK '
    CALL stopgm("FFTCHK_CUDA"," ",&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END FUNCTION FFTCHK_CUDA


END MODULE fftchk_utils

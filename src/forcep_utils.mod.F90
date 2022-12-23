MODULE forcep_utils
  USE fft_maxfft,                      ONLY: maxfftn
  USE isos,                            ONLY: isos1,&
                                             isos3
  USE linres,                          ONLY: td01,&
                                             tshl
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: fpar

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rhoe_psi_size

CONTAINS

  ! ==================================================================
  SUBROUTINE rhoe_psi_size(il_rhoe,il_psi,il_rhoe_1d, il_rhoe_2d, il_psi_1d, il_psi_2d )
    ! ==--------------------------------------------------------------==
    ! == RHOE is for ELECTRONIC DENSITY ARRAY                         ==
    ! == PSI  is ELECTRONIC POTENTIAL and FFT                         ==
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(OUT), OPTIONAL           :: il_rhoe, il_psi, il_rhoe_1d, &
                                                il_rhoe_2d, il_psi_1d, &
                                                il_psi_2d

    INTEGER                                  :: my_il_psi_1d

    IF(PRESENT(il_rhoe   )) il_rhoe    = fpar%nnr1*clsd%nlsd
    IF(PRESENT(il_rhoe_1d)) il_rhoe_1d = fpar%nnr1

    IF(PRESENT(il_rhoe_2d)) il_rhoe_2d = clsd%nlsd
    IF(td01%ns_tri.GT.0.OR.tshl%isc) il_rhoe_2d = 2

    ! FFT
    my_il_psi_1d = maxfftn
    ! For isolated system hip needs complex(8) :: PSI(2*NNR1)
    IF (isos1%tclust) THEN
       IF (isos3%ps_type.EQ.1) THEN
          my_il_psi_1d = 2 * my_il_psi_1d
       ELSEIF (isos3%ps_type.EQ.2) THEN
       ELSEIF (isos3%ps_type.EQ.3) THEN
       ENDIF
    ENDIF
    IF(PRESENT(il_psi   )) il_psi    = my_il_psi_1d* clsd%nlsd
    IF(PRESENT(il_psi_1d)) il_psi_1d = my_il_psi_1d
    IF(PRESENT(il_psi_2d)) il_psi_2d = clsd%nlsd
    IF(td01%ns_tri.GT.0.OR.tshl%isc) il_psi_2d = 2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rhoe_psi_size
  ! ==================================================================

END MODULE forcep_utils

MODULE switch_functionals_utils
  USE cp_xc_utils,                     ONLY: cp_xc,&
                                             cp_xc_functional,&
                                             cp_xc_kernel
  USE error_handling,                  ONLY: stopgm
  USE func,                            ONLY: func1,&
                                             func2,&
                                             func3
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: lr01,&
                                             lr03,&
                                             lrf1,&
                                             lrf2,&
                                             lrf3
  USE system,                          ONLY: cntl
  USE tbxc,                            ONLY: toldcode

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: switch_functionals

CONTAINS

  ! ==================================================================
  SUBROUTINE switch_functionals()
    !
    ! Moved to own module           28.09.2017 M.P.Bircher @ LCBC/EPFL

    CHARACTER(len=*), PARAMETER              :: procedureN = 'switch_functionals'
    INTEGER                                  :: ihold
    LOGICAL                                  :: lhold
    REAL(real_8)                             :: rhold

    !
    ! For the moment, pointer-based switching. Will be replaced by explicit passing (less ugly).
    !
    IF (cntl%use_xc_driver) THEN
       IF (ASSOCIATED(cp_xc,cp_xc_functional)) THEN
          cp_xc => cp_xc_kernel
       ELSEIF (ASSOCIATED(cp_xc,cp_xc_kernel)) THEN
          cp_xc => cp_xc_functional
       ELSE
          CALL stopgm(procedureN,'Invalid association status for cp_xc pointer',&
                      __LINE__,__FILE__)
       ENDIF
       CALL cp_xc%get( tgc=cntl%tgc, tgcx=cntl%tgcx, tgcc=cntl%tgcc, ttau=cntl%ttau, &
                       thybrid=cntl%thybrid, mhfx=func1%mhfx, phfx=func3%phfx, &
                       msrx=func1%msrx, srxa=func2%srxa, cam_alpha=func2%cam_alpha, cam_beta=func2%cam_beta )
    ELSE
       ihold           = lrf1%td_x
       lrf1%td_x       = func1%mfxcx
       func1%mfxcx     = ihold
       ihold           = lrf1%td_c
       lrf1%td_c       = func1%mfxcc
       func1%mfxcc     = ihold
       ihold           = lrf1%td_gx
       lrf1%td_gx      = func1%mgcx
       func1%mgcx      = ihold
       ihold           = lrf1%td_gc
       lrf1%td_gc      = func1%mgcc
       func1%mgcc      = ihold
       ihold           = lrf1%td_hf
       lrf1%td_hf      = func1%mhfx
       func1%mhfx      = ihold
       ihold           = lrf1%td_mtau
       lrf1%td_mtau    = func1%mtau
       func1%mtau      = ihold   
       lhold           = lrf2%td_tgc
       lrf2%td_tgc     = cntl%tgc
       cntl%tgc        = lhold   
       lhold           = lrf2%td_tgcx
       lrf2%td_tgcx    = cntl%tgcx
       cntl%tgcx       = lhold   
       lhold           = lrf2%td_tgcc
       lrf2%td_tgcc    = cntl%tgcc
       cntl%tgcc       = lhold   
       lhold           = lrf2%td_code
       lrf2%td_code    = toldcode
       toldcode        = lhold   
       lhold           = lrf2%td_hybrid
       lrf2%td_hybrid  = cntl%thybrid
       cntl%thybrid    = lhold   
       rhold           = lrf3%tdpxlda
       lrf3%tdpxlda    = func3%pxlda
       func3%pxlda     = rhold   
       rhold           = lrf3%tdpclda
       lrf3%tdpclda    = func3%pclda
       func3%pclda     = rhold   
       rhold           = lrf3%tdpxgc
       lrf3%tdpxgc     = func3%pxgc
       func3%pxgc      = rhold   
       rhold           = lrf3%tdpcgc
       lrf3%tdpcgc     = func3%pcgc
       func3%pcgc      = rhold   
       rhold           = lrf3%tdphfx
       lrf3%tdphfx     = func3%phfx
       func3%phfx      = rhold   
       lhold           = lrf2%td_ttau
       lrf2%td_ttau    = cntl%ttau
       cntl%ttau       = lhold   
    ENDIF

    ! ==--------------------------------------------------------------==
  END SUBROUTINE switch_functionals
  ! ==================================================================
END MODULE switch_functionals_utils

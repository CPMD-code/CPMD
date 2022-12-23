MODULE getfu_utils
  USE cnst,                            ONLY: factem
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE pimd,                            ONLY: eharv,&
                                             pimd3,&
                                             pmars
  USE system,                          ONLY: cntr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: getfu

CONTAINS

  ! ==================================================================
  SUBROUTINE getfu(fion,fstage,stage,iflag)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:,:), &
                                                fstage(:,:,:,:), &
                                                stage(:,:,:,:)
    INTEGER                                  :: iflag

    INTEGER                                  :: i, ia, ip, is, npx
    REAL(real_8)                             :: fbar, fh, omegap, omp2, rnp

    npx=pimd3%np_total
    rnp=REAL(npx,kind=real_8)
    DO ip=1,npx
       eharv(ip) = 0._real_8
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  PERFORM TRANSFORMATION FR --> FU                            ==
    ! ==  FU_1 = \sum_{s=1}^P FR_s (FR IS ASSUMED TO BE DIVIDED BY P) ==
    ! ==  FU_s = [(s-2)/(s-1)]*FU_{s-1} + FR_s                        ==
    ! ==  IFLAG=1 : CALCULATE FSTAGE (INCLUDING HARMONIC FORCE)       ==
    ! ++            AND HARMONIC ENERGY (IN EHARV(IP))                ==
    ! ==        0   CALCULATE ONLY HARMONIC ENERGY (IN EHARV(IP))     ==
    ! ==            WITHOUT CHANGING FSTAGE                           ==
    ! ==--------------------------------------------------------------==
    omegap = SQRT(rnp)*cntr%tempw/factem
    omp2 = omegap*omegap
    IF (iflag.EQ.1) THEN
       DO i=1,3
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                fbar = 0.0_real_8
                DO ip=1,npx
                   fbar = fbar + fion(i,ia,is,ip)
                ENDDO
                fstage(i,ia,is,1) = fbar
                DO ip=2,npx
                   fstage(i,ia,is,ip) = REAL(ip-2,kind=real_8)*fstage(i,ia,is,ip-1)&
                        /REAL(ip-1,kind=real_8)&
                        + fion(i,ia,is,ip)
                ENDDO
                DO ip=2,npx
                   fh = -pmars(is,ip)*omp2*stage(i,ia,is,ip)
                   fstage(i,ia,is,ip) = fstage(i,ia,is,ip) + fh
                   eharv(ip) = eharv(ip) - 0.5_real_8*fh*stage(i,ia,is,ip)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ELSE
       DO i=1,3
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                DO ip=2,npx
                   fh = -pmars(is,ip)*omp2*stage(i,ia,is,ip)
                   eharv(ip) = eharv(ip) - 0.5_real_8*fh*stage(i,ia,is,ip)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getfu
  ! ==================================================================

END MODULE getfu_utils

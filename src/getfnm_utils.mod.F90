MODULE getfnm_utils
  USE cnst,                            ONLY: factem
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE pimd,                            ONLY: eharv,&
                                             pimd3,&
                                             pmars,&
                                             tnmi
  USE system,                          ONLY: cntr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: getfnm

CONTAINS

  ! ==================================================================
  SUBROUTINE getfnm(fion,fstage,stage,iflag)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:,:), &
                                                fstage(:,:,:,:), &
                                                stage(:,:,:,:)
    INTEGER                                  :: iflag

    INTEGER                                  :: i, ia, ip, is, jp, npx
    REAL(real_8)                             :: fh, omegap, omp2, rnp

    npx=pimd3%np_total
    rnp=REAL(npx,kind=real_8)
    DO ip=1,npx
       eharv(ip) = 0._real_8
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  PERFORM TRANSFORMATION FR --> FU                            ==
    ! ==  FU_i = \sum_j FR_j*TNMI(j,i) (FR ASSUMED TO BE DIVIDED BY P)==
    ! ==  IFLAG=1 : CALCULATE FSTAGE (INCLUDING HARMONIC FORCE)       ==
    ! ==            AND HARMONIC ENERGY (IN EHARV(IP))                ==
    ! ==        0   CALCULATE ONLY HARMONIC ENERGY (IN EHARV(IP))     ==
    ! ==            WITHOUT CHANGING FSTAGE                           ==
    ! ==--------------------------------------------------------------==
    omegap = SQRT(rnp)*cntr%tempw/factem
    omp2 = omegap*omegap
    IF (iflag.EQ.1) THEN
       CALL zeroing(fstage)!,3*maxsys%nax*maxsys%nsx*npx)
       DO i=1,3
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                DO ip=1,npx
                   DO jp=1,npx
                      fstage(i,ia,is,ip) = fstage(i,ia,is,ip) +&
                           fion(i,ia,is,jp)*tnmi(jp,ip)
                   ENDDO
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
  END SUBROUTINE getfnm
  ! ==================================================================

END MODULE getfnm_utils

MODULE pinmtrans_utils
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE pimd,                            ONLY: pimd3,&
                                             tnm,&
                                             tnmi
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pinmtrans

CONTAINS

  ! ==================================================================
  SUBROUTINE pinmtrans(prim,stg,itrans)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: prim(:,:,:,:), stg(:,:,:,:)
    INTEGER                                  :: itrans

    INTEGER                                  :: i, ia, ip, is, jp, npx

    npx=pimd3%np_total
    IF (itrans.EQ.0) THEN
       ! ==--------------------------------------------------------------==
       ! ==  PERFORM TRANSFORMATION R --> U                              ==
       ! ==  U_i = \sum_j TNM(i,j)*R_j                                   ==
       ! ==  WHERE TNM IS THE FORWARD NORMAL MODE TRANSFORMATION MATRIX  ==
       ! ==--------------------------------------------------------------==
       CALL zeroing(stg)!,3*maxsys%nax*maxsys%nsx*npx)
       DO i=1,3
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                DO ip=1,npx
                   DO jp=1,npx
                      stg(i,ia,is,ip) = stg(i,ia,is,ip) +&
                           tnm(ip,jp)*prim(i,ia,is,jp)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ELSE
       ! ==--------------------------------------------------------------==
       ! ==  PERFORM INVERSE TRANSFORMATION U --> R                      ==
       ! ==  R_i = \sum_j TNMI(i,j)*U_j                                  ==
       ! ==  WHERE TNMI IS THE INVERSE TRANSFORMATION MATRIX             ==
       ! ==--------------------------------------------------------------==
       CALL zeroing(prim)!,3*maxsys%nax*maxsys%nsx*npx)
       DO i=1,3
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                DO ip=1,npx
                   DO jp=1,npx
                      prim(i,ia,is,ip) = prim(i,ia,is,ip) +&
                           tnmi(ip,jp)*stg(i,ia,is,jp)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE pinmtrans
  ! ==================================================================

END MODULE pinmtrans_utils

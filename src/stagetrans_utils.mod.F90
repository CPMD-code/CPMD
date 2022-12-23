MODULE stagetrans_utils
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE pimd,                            ONLY: pimd3

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: stagetrans

CONTAINS

  ! ==================================================================
  SUBROUTINE stagetrans(prim,stg,itrans)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: prim(:,:,:,:), stg(:,:,:,:)
    INTEGER                                  :: itrans

    INTEGER                                  :: i, ia, ip, is, npx
    REAL(real_8)                             :: rstar

    npx=pimd3%np_total
    IF (itrans.EQ.0) THEN
       ! ==--------------------------------------------------------------==
       ! ==  PERFORM TRANSFORMATION R --> U                              ==
       ! ==  U_1 = R_1, U_s = R_s - R_s^* ,   s=2,...,P                  ==
       ! ==  R_s^* = [(s-1)R_{s+1} + R_1]/s                              ==
       ! ==--------------------------------------------------------------==
       DO i=1,3
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                stg(i,ia,is,1) = prim(i,ia,is,1)
                stg(i,ia,is,npx) = prim(i,ia,is,npx)-prim(i,ia,is,1)
                DO ip=2,npx-1
                   rstar = (REAL(ip-1,kind=real_8)*prim(i,ia,is,ip+1) +&
                        prim(i,ia,is,1))/REAL(ip,kind=real_8)
                   stg(i,ia,is,ip) = prim(i,ia,is,ip) - rstar
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ELSE
       ! ==--------------------------------------------------------------==
       ! ==  PERFORM RECURSIVE INVERSE TRANSFORMATION U --> R            ==
       ! ==  R_1 = U_1, R_s = U_1 + \sum_{t=s}^P [(s-1)/(t-1)]*U_t       ==
       ! ==--------------------------------------------------------------==
       DO i=1,3
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                prim(i,ia,is,1) = stg(i,ia,is,1)
                prim(i,ia,is,npx) = stg(i,ia,is,npx)+stg(i,ia,is,1)
                DO ip=npx-1,2,-1
                   prim(i,ia,is,ip) = stg(i,ia,is,ip) +&
                        REAL(ip-1,kind=real_8)*prim(i,ia,is,ip+1)/REAL(ip,kind=real_8) +&
                        prim(i,ia,is,1)/REAL(ip,kind=real_8)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE stagetrans
  ! ==================================================================

END MODULE stagetrans_utils

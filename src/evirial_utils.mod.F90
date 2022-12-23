MODULE evirial_utils
  USE cnst,                            ONLY: factem
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nose,                            ONLY: glib
  USE pbc_utils,                       ONLY: pbc
  USE pimd,                            ONLY: fionks,&
                                             pimd3
  USE system,                          ONLY: cntr,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: evirial

CONTAINS

  ! ==================================================================
  SUBROUTINE evirial(taup,eviri)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==                REAL IONIC KINETIC ENERGY                     ==
    ! ==                WITH THE VIRIAL ESTIMATOR                     ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: taup(:,:,:,:), eviri

    INTEGER                                  :: ia, ip, is, k
    REAL(real_8)                             :: fglib, rmean(3), rnp, xlm, &
                                                ylm, zlm, xlm_, ylm_, zlm_

    rnp=REAL(pimd3%np_total,kind=real_8)
    eviri=0.0_real_8
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO k=1,3
             rmean(k)=0.0_real_8
             DO ip=1,pimd3%np_total
                rmean(k)=rmean(k)+taup(k,ia,is,ip)
             ENDDO
             rmean(k)=rmean(k)/rnp
          ENDDO
          DO ip=1,pimd3%np_total
             xlm_=taup(1,ia,is,ip)-rmean(1)
             ylm_=taup(2,ia,is,ip)-rmean(2)
             zlm_=taup(3,ia,is,ip)-rmean(3)
             CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
             eviri=eviri-xlm*fionks(1,ia,is,ip)&
                  -ylm*fionks(2,ia,is,ip)&
                  -zlm*fionks(3,ia,is,ip)
          ENDDO
       ENDDO
    ENDDO
    ! note that FIONKS is already divided by NP
    ! GLIB ACCOUNT ONLY FOR DYNAMICAL DEGREES OF FREEDOM IF CONTRAINTS!
    fglib=(glib+3.0_real_8*ions1%nat*(rnp-1))/rnp
    eviri=fglib*cntr%tempw/(2.0_real_8*factem) + eviri/2.0_real_8
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE evirial
  ! ==================================================================

END MODULE evirial_utils

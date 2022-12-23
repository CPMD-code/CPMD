MODULE fharm_utils
  USE cnst,                            ONLY: factem
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE pbc_utils,                       ONLY: pbc
  USE pimd,                            ONLY: ehar,&
                                             ipm1,&
                                             ipp1,&
                                             pimd3,&
                                             pmar
  USE system,                          ONLY: cntr,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fharm

CONTAINS

  ! ==================================================================
  SUBROUTINE fharm(tau0,fion,inp,tforce)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==               THE HARMONIC SPRING ENERGY                     ==
    ! ==               THE HARMONIC SPRING FORCES                     ==
    ! ==                      OF THE IONS                             ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:,:), fion(:,:,:)
    INTEGER                                  :: inp
    LOGICAL                                  :: tforce

    INTEGER                                  :: ia, is
    REAL(real_8)                             :: face, facf, rnp, xlm, xlmm1, &
                                                xlmp1, ylm, ylmm1, ylmp1, &
                                                zlm, zlmm1, zlmp1, xlm_, ylm_, zlm_

    rnp=REAL(pimd3%np_total,kind=real_8)
    ehar=0.0_real_8
    DO is=1,ions1%nsp
       facf = pmar(is)*rnp*(cntr%tempw/factem)**2
       face = facf/2.0_real_8
       DO ia=1,ions0%na(is)
          xlm_=tau0(1,ia,is,inp)-tau0(1,ia,is,ipp1(inp))
          ylm_=tau0(2,ia,is,inp)-tau0(2,ia,is,ipp1(inp))
          zlm_=tau0(3,ia,is,inp)-tau0(3,ia,is,ipp1(inp))
          CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
          xlmp1=xlm
          ylmp1=ylm
          zlmp1=zlm
          IF (tforce) THEN
             xlm_=tau0(1,ia,is,inp)-tau0(1,ia,is,ipm1(inp))
             ylm_=tau0(2,ia,is,inp)-tau0(2,ia,is,ipm1(inp))
             zlm_=tau0(3,ia,is,inp)-tau0(3,ia,is,ipm1(inp))
             CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
             xlmm1=xlm
             ylmm1=ylm
             zlmm1=zlm
             fion(1,ia,is) = fion(1,ia,is) - facf*(xlmm1+xlmp1)
             fion(2,ia,is) = fion(2,ia,is) - facf*(ylmm1+ylmp1)
             fion(3,ia,is) = fion(3,ia,is) - facf*(zlmm1+zlmp1)
          ENDIF
          ehar = ehar+face*( xlmp1**2+ylmp1**2+zlmp1**2 )
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fharm
  ! ==================================================================

END MODULE fharm_utils

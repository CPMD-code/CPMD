MODULE csize_utils
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_max,&
                                             mp_sum
  USE parac,                           ONLY: parai
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             nkpt,&
                                             spar
  USE utils,                           ONLY: zgive

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: csize

CONTAINS

  ! ==================================================================
  SUBROUTINE csize(c2,nstate,gemax,cnorm)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,nstate)
    REAL(real_8)                             :: gemax, cnorm

    CHARACTER(*), PARAMETER                  :: procedureN = 'csize'

    INTEGER                                  :: i, iiabs, nocc, nod, nou
    REAL(real_8)                             :: gd, gu
    REAL(real_8), EXTERNAL                   :: ddot

#if defined(__SR11KIBM)
    REAL(real_8), EXTERNAL                   :: dzamax
#else
    INTEGER, EXTERNAL                        :: izamax
#endif
    ! ==--------------------------------------------------------------==
    nocc=0
    nou=0
    IF (nkpt%ngwk.GT.0) THEN
       IF (cntl%tlsd) THEN
          nou=0
          DO i=1,spin_mod%nsup
             IF (crge%f(i,1).GT.1.e-5_real_8) nou=nou+1
          ENDDO
          IF (nou.GT.0) THEN
#if defined(__SR11KIBM)
             gu=dzamax(nou*nkpt%ngwk,c2,1)
#else
             iiabs=izamax(nou*nkpt%ngwk,c2,1)
             IF (iiabs<1.OR.iiabs>nou*nkpt%ngwk) CALL stopgm(procedureN,&
                  'out of bound',& 
                  __LINE__,__FILE__)
             gu=ABS(zgive(c2,iiabs))
#endif
          ELSE
             gu=0._real_8
          ENDIF
          nod=0
          DO i=spin_mod%nsup+1,nstate
             IF (crge%f(i,1).GT.1.e-5_real_8) nod=nod+1
          ENDDO
          IF (nod.GT.0) THEN
#if defined(__SR11KIBM)
             gd=dzamax(nod*nkpt%ngwk,c2(1,spin_mod%nsup+1),1)
#else
             iiabs=izamax(nod*nkpt%ngwk,c2(1,spin_mod%nsup+1),1)
             IF (iiabs<1.OR.iiabs>nod*nkpt%ngwk) CALL stopgm(procedureN,&
                  'out of bound',& 
                  __LINE__,__FILE__)
             gd=ABS(zgive(c2(1,spin_mod%nsup+1),iiabs))
#endif
          ELSE
             gd=0._real_8
          ENDIF
          gemax=MAX(gu,gd)
          nocc=nod+nou
          cnorm=0.0_real_8
          IF (tkpts%tkpnt) THEN
             DO i=1,nou
                cnorm=cnorm+ddot(2*nkpt%ngwk,c2(1,i),1,c2(1,i),1)
             ENDDO
             DO i=spin_mod%nsup+1,spin_mod%nsup+nod
                cnorm=cnorm+ddot(2*nkpt%ngwk,c2(1,i),1,c2(1,i),1)
             ENDDO
          ELSE
             DO i=1,nou
                cnorm=cnorm+dotp(ncpw%ngw,c2(:,i),c2(:,i))
             ENDDO
             DO i=spin_mod%nsup+1,spin_mod%nsup+nod
                cnorm=cnorm+dotp(ncpw%ngw,c2(:,i),c2(:,i))
             ENDDO
          ENDIF
       ELSE
          nocc=0
          DO i=1,nstate
             IF (crge%f(i,1).GT.1.e-5_real_8) nocc=nocc+1
          ENDDO
#if defined(__SR11KIBM)
          gemax=dzamax(nocc*nkpt%ngwk,c2,1)
#else 
          iiabs=izamax(nocc*nkpt%ngwk,c2,1)
          IF (iiabs<1.OR.iiabs>nocc*nkpt%ngwk) CALL stopgm(procedureN,&
               'out of bound',& 
               __LINE__,__FILE__)
          gemax=ABS(zgive(c2(1,1),iiabs))
#endif
          cnorm=0.0_real_8
          IF (tkpts%tkpnt) THEN
             DO i=1,nocc
                cnorm=cnorm+ddot(2*nkpt%ngwk,c2(1,i),1,c2(1,i),1)
             ENDDO
          ELSE
             DO i=1,nocc
                cnorm=cnorm+dotp(ncpw%ngw,c2(:,i),c2(:,i))
             ENDDO
          ENDIF
       ENDIF
    ELSE
       cnorm=0._real_8
       gemax=0._real_8
    ENDIF
    CALL mp_sum(cnorm,parai%allgrp)
    CALL mp_max(gemax,parai%allgrp)
    cnorm=SQRT(cnorm/REAL(nocc*spar%ngwks,kind=real_8))
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE csize
  ! ==================================================================

END MODULE csize_utils

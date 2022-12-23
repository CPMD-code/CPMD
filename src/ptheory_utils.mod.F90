MODULE ptheory_utils
  USE error_handling,                  ONLY: stopgm
  USE fint,                            ONLY: fint1
  USE hpsi_utils,                      ONLY: give_scr_hpsi,&
                                             hpsi
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE rgs_utils,                       ONLY: rgs_c
  USE system,                          ONLY: fpar,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ptheory
  PUBLIC :: give_scr_ptheory

CONTAINS

  ! ==================================================================
  SUBROUTINE ptheory(c0,cs,sc0,we,adav,vpot,psi,&
       nstate,ikind,nhpsi)
    ! ==--------------------------------------------------------------==
    ! == CHANGES THE STATES ACCORDING TO FIRST-ORDER PERTUBATION      ==
    ! == THEORY                                                       ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: sc0(nkpt%ngwk,*)
    REAL(real_8)                             :: vpot(fpar%nnr1)
    COMPLEX(real_8)                          :: psi(fpar%nnr1)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: adav(nstate,nstate)
    REAL(real_8)                             :: we(nstate)
    COMPLEX(real_8)                          :: cs(nkpt%ngwk,nstate), &
                                                c0(nkpt%ngwk,nstate)
    INTEGER                                  :: ikind, nhpsi

    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) ,&
                                                zzero = (0._real_8,0._real_8)
    INTEGER, PARAMETER                       :: ispin = 1 

    INTEGER                                  :: i, isub, j
    REAL(real_8)                             :: factor

    CALL tiset('   PTHEORY',isub)
    ! ==--------------------------------------------------------------==
    IF (tkpts%tkpnt) THEN
       CALL hpsi(c0,cs,sc0,vpot,psi,nstate,ikind,ispin)
       CALL ovlap_h(nstate,adav,c0,cs)
       DO i=1,nstate
          adav(i,i)=CMPLX(0._real_8,0._real_8,kind=real_8)
          DO j=i+1,nstate
             IF (fint1%ttrot) THEN
                factor=-(EXP(-fint1%betap*we(j))-EXP(-fint1%betap*we(i)))/&
                     (we(j)-we(i))**2
             ELSE
                factor=-1._real_8/(we(j)-we(i))
             ENDIF
             adav(j,i)=adav(j,i)*factor
             adav(i,j)=CONJG(adav(j,i))
          ENDDO
       ENDDO
       CALL zgemm('N','N',nkpt%ngwk,nstate,nstate,zone,c0,nkpt%ngwk,&
            adav,nstate,zzero,cs,nkpt%ngwk)
       DO i=1,nstate
          CALL daxpy(nkpt%ngwk*2,cs(1,i),1,zone,c0(1,i),1)
       ENDDO
       CALL rgs_c(c0,nstate,adav)
    ELSE
       CALL stopgm('PTHEORY',' TKPNT .FALSE.',& 
            __LINE__,__FILE__)
    ENDIF
    nhpsi=nhpsi+1
    CALL tihalt('   PTHEORY',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ptheory
  ! ==================================================================
  SUBROUTINE give_scr_ptheory(lptheory,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lptheory
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    CALL give_scr_hpsi(lptheory,tag,nstate)
    lptheory=MAX(lptheory,nstate)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_ptheory
  ! ==================================================================

END MODULE ptheory_utils

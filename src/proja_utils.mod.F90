MODULE proja_utils
  USE csize_utils,                     ONLY: csize
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE rotate_utils,                    ONLY: rotate
  USE spin,                            ONLY: spin_mod
  USE spsi_utils,                      ONLY: spsi
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: proja

CONTAINS

  ! ==================================================================
  SUBROUTINE proja(c0,c2,sc0,nstate,mproj)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: sc0(ncpw%ngw,nstate), &
                                                c2(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    INTEGER                                  :: mproj

    CHARACTER(*), PARAMETER                  :: procedureN = 'proja'

    INTEGER                                  :: i, ierr, isub, length
    REAL(real_8), ALLOCATABLE                :: aux(:)

    CALL tiset('     PROJA',isub)
    IF (mproj.EQ.0) THEN
       IF (pslo_com%tivan) THEN
          length=MAX(nstate,ncpw%ngw*2)
       ELSE
          length=nstate
       ENDIF
    ELSEIF (mproj.EQ.2) THEN
       IF (pslo_com%tivan) THEN
          length=MAX(nstate*nstate,ncpw%ngw*2)
       ELSE
          length=nstate*nstate
       ENDIF
    ENDIF
    ALLOCATE(aux(length),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (mproj.EQ.0) THEN
       ! do nothing
    ELSEIF (mproj.EQ.1) THEN
       ! ==------------------------------------------------------------==
       ! ==  Projection C2(I) <= C2(I) - <C0(I)|C2(I)> SC0(I)          ==
       ! ==------------------------------------------------------------==
       IF (pslo_com%tivan) THEN
          CALL dcopy(2*ncpw%ngw*nstate,c0(1,1),1,sc0(1,1),1)
          CALL spsi(nstate,sc0)
          DO i=1,nstate
             aux(i)=dotp(ncpw%ngw,c0(:,i),c2(:,i))
          ENDDO
          CALL mp_sum(aux,nstate,parai%allgrp)
          DO i=1,nstate
             CALL daxpy(2*ncpw%ngw,-aux(i),sc0(1,i),1,c2(1,i),1)
          ENDDO
       ELSE
          DO i=1,nstate
             aux(i)=dotp(ncpw%ngw,c0(:,i),c2(:,i))
          ENDDO
          CALL mp_sum(aux,nstate,parai%allgrp)
          DO i=1,nstate
             CALL daxpy(2*ncpw%ngw,-aux(i),c0(1,i),1,c2(1,i),1)
          ENDDO
       ENDIF
    ELSEIF (mproj.EQ.2) THEN
       ! ==------------------------------------------------------------==
       ! ==  Projection A(I) <= A(I) - SUM{J} <B(J)|A(I)> SB(J)        ==
       ! ==------------------------------------------------------------==
       IF (pslo_com%tivan) THEN
          CALL dcopy(2*ncpw%ngw*nstate,c0(1,1),1,sc0(1,1),1)
          CALL spsi(nstate,sc0)
          CALL ovlap(nstate,aux,c0,c2)
          CALL mp_sum(aux,nstate*nstate,parai%allgrp)
          CALL rotate(-1.0_real_8,sc0,1.0_real_8,c2,aux,nstate,2*ncpw%ngw,&
               cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
       ELSE
          CALL ovlap(nstate,aux,c0,c2)
          CALL mp_sum(aux,nstate*nstate,parai%allgrp)
          CALL rotate(-1.0_real_8,c0,1.0_real_8,c2,aux,nstate,2*ncpw%ngw,&
               cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
       ENDIF
    ELSE
       CALL stopgm('PROJA','MPROJ OPTION',& 
            __LINE__,__FILE__)
    ENDIF
    CALL csize(c2,nstate,gemax,cnorm)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('     PROJA',isub)
    RETURN
  END SUBROUTINE proja
  ! ==================================================================

END MODULE proja_utils

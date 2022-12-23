MODULE gettrans_utils
  USE cnst,                            ONLY: ry
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: nua,&
                                             nub,&
                                             td01,&
                                             td03,&
                                             urot
  USE mp_interface,                    ONLY: mp_sum
  USE orbrot_utils,                    ONLY: get_u
  USE parac,                           ONLY: parai,&
                                             paral
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gettrans
  PUBLIC :: printtrans

CONTAINS

  ! ==================================================================
  SUBROUTINE gettrans(transition,c1,cvir,nstate,nvir,nap,&
       nbp,nlr)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: transition(:,:,:)
    COMPLEX(real_8)                          :: cvir(:,:)
    INTEGER                                  :: nstate, nvir, nap, nbp, nlr
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate,nlr)

    CHARACTER(*), PARAMETER                  :: procedureN = 'gettrans'

    INTEGER                                  :: i1, i2, idamax, ierr, is, iv, &
                                                iz, j, nz
    LOGICAL                                  :: debug, ostda
    REAL(real_8)                             :: val, vlast
    REAL(real_8), ALLOCATABLE                :: peig(:), pmat(:,:), &
                                                trans(:,:), trot(:,:), &
                                                umat(:,:)

! ==================================================================

    ostda=(td03%tda.AND.td01%msubta.GT.0.AND..NOT.td03%tdacanon)
    debug=.FALSE.
    IF (ostda) THEN
       nz=nua
    ELSE
       nz=nstate
    ENDIF
    ALLOCATE(trans(nvir, nz),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zeroing(trans)!,nvir*nz)
    DO is=1,nlr
       IF (cntl%tlsd) THEN
          DO iv=1,nap
             DO j=1,spin_mod%nsup
                trans(iv,j)=SQRT(2._real_8)*dotp(ncpw%ngw,c1(:,j,is),cvir(:,iv))
             ENDDO
          ENDDO
          DO iv=nap+1,nap+nbp
             DO j=spin_mod%nsup+1,nstate
                trans(iv,j)=SQRT(2._real_8)*dotp(ncpw%ngw,c1(:,j,is),cvir(:,iv))
             ENDDO
          ENDDO
       ELSE
          DO iv=1,nvir
             DO j=1,nstate
                trans(iv,j)=dotp(ncpw%ngw,c1(:,j,is),cvir(:,iv))
             ENDDO
          ENDDO
       ENDIF
       CALL mp_sum(trans,nvir*nstate,parai%allgrp)
       IF (ostda) THEN
          ALLOCATE(trot(nvir, nz),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(umat(nua, nub),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(pmat(nua, nub),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          ALLOCATE(peig(nub),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)

          CALL get_u(urot(1,1,is),umat,pmat,peig,nua,nub)
          CALL dcopy(nvir*nz,trans,1,trot,1)
          CALL zeroing(trans)!,nvir*nstate)
          DO iv=1,nvir
             DO iz=1,nz
                DO j=1,nstate
                   trans(iv,j)=trans(iv,j)+umat(j,iz)*trot(iv,iz)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       vlast=0._real_8
       DO iz=1,4
          j=idamax(nvir*nz,trans,1)
          i1=MOD(j-1,nvir)+1
          i2=(j-1)/nvir+1
          val=ABS(trans(i1,i2))
          IF (val.LT.0.2_real_8*vlast ) GOTO 100
          IF (val.LT.0.1_real_8) GOTO 100
          transition(1,iz,is)=i2
          transition(2,iz,is)=i1
          transition(3,iz,is)=NINT(10000._real_8*val)
          trans(i1,i2)=0._real_8
          vlast=val
       ENDDO
100    CONTINUE
    ENDDO
    ! ==--------------------------------------------------------------==
    DEALLOCATE(trans,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF (ostda) THEN
       DEALLOCATE(trot,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(umat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(pmat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(peig,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gettrans
  ! ==================================================================
  SUBROUTINE printtrans(transition,eigv,nstate,nvir,npa,npb,nlr)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: transition(3,4,20)
    REAL(real_8)                             :: eigv(*)
    INTEGER                                  :: nstate, nvir, npa, npb, nlr

    CHARACTER(len=5), DIMENSION(2)           :: state = (/"ALPHA","BETA "/)
    CHARACTER(len=7), DIMENSION(2) :: stype = (/"HOMO - ","TOP -  "/)
    INTEGER                                  :: i, io, is, istat, it, iv, &
                                                ixst, nz
    LOGICAL                                  :: ostda

    ostda=(td03%tda.AND.td01%msubta.GT.0.AND..NOT.td03%tdacanon)
    IF (paral%io_parent) THEN
       IF (td03%treorder.OR.td03%tdlocal) THEN
          ixst=2
       ELSE
          ixst=1
       ENDIF
       IF (ostda) THEN
          nz=nua
       ELSE
          nz=nstate
       ENDIF
       DO is=1,nlr
          it=0
          DO i=1,4
             IF (transition(1,i,is).EQ.0) GOTO 100
             it=it+1
          ENDDO
100       CONTINUE
          WRITE(6,'(//,A,I3,T42,A,F10.3,A)') " STATE=",is,"EIGENVALUE=",&
               2._real_8*ry*eigv(is)," eV"
          IF (cntl%tlsd) THEN
             DO i=1,it
                io=transition(1,i,is)
                IF (io.GT.spin_mod%nsup) THEN
                   istat=2
                   io=nstate-io
                ELSE
                   istat=1
                   io=spin_mod%nsup-io
                ENDIF
                iv=transition(2,i,is)
                IF (iv.GT.npa) iv=iv-npa
                WRITE(6,'(T3,A,2X,A,T24,F7.3,T41,A,I3,A,A,I3)')&
                     " TRANSITION ",state(istat),(REAL(transition(3,i,&
                     is),kind=real_8)/10000._real_8) **2,stype(ixst),io," --> ","LUMO + ",&
                     iv-1
             ENDDO
          ELSE
             DO i=1,it
                WRITE(6,'(T10,A,T24,F7.3,T41,A,I3,A,A,I3)')&
                     " TRANSITION " ,&
                     (REAL(transition(3,i,is),kind=real_8)/10000._real_8)**2,stype(ixst),&
                     nz-transition(1,i,is)," --> ","LUMO + ",&
                     transition(2,i,is)-1
             ENDDO
          ENDIF
       ENDDO
       WRITE(6,'(/,A,T53,F10.4,A)') " CHKSUM(EIGENVALUES) = ",&
            2.0_real_8*ry*SUM(eigv(1:nlr)),' eV'
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE printtrans
  ! ==================================================================

END MODULE gettrans_utils

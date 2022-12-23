MODULE exterp_utils
  USE efld,                            ONLY: extf
  USE error_handling,                  ONLY: stopgm
  USE extpotmod,                       ONLY: extpot
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_recv,&
                                             mp_send,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             parap,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: exterp

CONTAINS

  ! ==================================================================
  SUBROUTINE exterp
    ! ==--------------------------------------------------------------==
    ! Variables

    CHARACTER(*), PARAMETER                  :: procedureN = 'exterp'

    INTEGER                                  :: ierr, isub, nn

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF (cntl%texpot) THEN
       nn=fpar%kr1*fpar%kr2s*fpar%kr3s
       ! nn=nr1*nr2s*nr3s
       ALLOCATE(extpot(nn),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL read_pot(extpot)

       IF (cntl%texadd) THEN
          extf => extpot
       ENDIF

       IF (paral%parent) CALL prmem('    EXTERP')
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE exterp
  ! ==================================================================
  SUBROUTINE read_pot(extpot)
    ! ==--------------------------------------------------------------==
    REAL(real_8) :: extpot(fpar%kr1,fpar%kr2s,fpar%kr3s)

    CHARACTER(*), PARAMETER                  :: procedureN = 'read_pot'

    INTEGER                                  :: i, ierr, ip, ipp, ipx, ir, &
                                                ir2, irr, kk, msgid
    REAL(real_8), ALLOCATABLE                :: pscr(:,:)

    kk=fpar%kr2s*spar%nr3s
    ! if(parent) open(99,file='extpot.grid')
    IF (paral%io_parent) THEN
       OPEN(99,file='extpot.unfo.grid',form='unformatted')
       WRITE (6,'(A)')&
            ' READING EXTERNAL POTENTIAL FROM extpot.unfo.grid'
       IF (cntl%texadd) WRITE (6,'(A)')&
            ' EXTERNAL POTENTIAL ADDED TO FORCES ACTING ON IONS'
       IF (.NOT.isos1%tclust) WRITE (6,'(A)')&
            ' WARNING! APPLYING FIELD IN REAL-SPACE FOR A PERIODIC SYSTEM!'
    ENDIF
    CALL zeroing(extpot)!,kr1*kr2s*kr3s)
    ALLOCATE(pscr(fpar%kr2s,kk/fpar%kr2s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO ir=1,spar%nr1s
       IF (paral%io_parent) THEN
          CALL zeroing(pscr)!,kk)
          DO ir2=1,spar%nr2s     ! read y
             READ(99)(pscr(ir2,i),i=1,spar%nr3s)! read z
          ENDDO
          ! READ(99) (PSCR(I),I=1,KK)
       ENDIF
       DO ipp=1,parai%nproc
          IF (ir.GE.parap%nrxpl(ipp-1,1).AND.ir.LE.parap%nrxpl(ipp-1,2)) THEN
             ipx=ipp
             GOTO 10
          ENDIF
       ENDDO
       CALL stopgm('RD30POT','ERROR IN NRXPL',& 
            __LINE__,__FILE__)
10     CONTINUE
       ip=parap%pgroup(ipx)
       CALL mp_sync(parai%allgrp)
       IF (ip.NE.parai%source) THEN
          IF (paral%parent) THEN
             msgid=1
             !msglen=kk * 8
             CALL mp_send(pscr,kk,ip,msgid,parai%allgrp)
          ELSEIF (parai%me.EQ.ip) THEN
             msgid=1
             !msglen=kk * 8
             CALL mp_recv(pscr,kk,parap%pgroup(1),msgid,parai%allgrp)
          ENDIF
          CALL mp_sync(parai%allgrp)
       ENDIF
       IF (ip.EQ.parai%me) THEN
          irr=ir-parap%nrxpl(parai%mepos,1)+1
          CALL dcopy(kk,pscr,1,extpot(irr,1,1),fpar%kr1)
       ENDIF
    ENDDO
    DEALLOCATE(pscr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (paral%io_parent) CLOSE(99)
    RETURN
  END SUBROUTINE read_pot
  ! ==================================================================

END MODULE exterp_utils

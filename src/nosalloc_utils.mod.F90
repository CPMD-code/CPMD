MODULE nosalloc_utils
  USE bsym,                            ONLY: bsfac
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_app,&
                                             fo_old,&
                                             fo_verb
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_revert
  USE mm_input,                        ONLY: lqmmm
  USE mp_interface,                    ONLY: mp_bcast
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE nose,                            ONLY: &
       cafescl, cafesini, cafesinr, cafesnse, cafestmp, chainer, dofsl, eta, &
       etadot, etap, etap1, etap1dot, etapdot, etapm, etapmdot, etc, etcdot, &
       gkt, lctmemb, lcttab, loct, mapdof, ncafesgrp, nchx, ndfnt, nlctmbm, &
       nosl, qnosp, qnospc, qnospm, tcafes, tempwr, wnosep2r, wnosepr, gkt1
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: grandparent,&
                                             pimd1,&
                                             pimd3,&
                                             supergroup,&
                                             supersource,&
                                             np_high,&
                                             np_low
  USE prmem_utils,                     ONLY: prmem
  USE system,                          ONLY: cntl,&
                                             maxsys
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nosalloc

CONTAINS

  ! ==================================================================
  SUBROUTINE nosalloc
    ! ==--------------------------------------------------------------==
    ! ==  ALLOCATE NOSE ARRAYS                                        ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'nosalloc'

    INTEGER                                  :: i, ia, ierr, is, isub, k, nn, &
                                                nr, nrepl, ip
    LOGICAL                                  :: ferror, ionode, oldstatus

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    ! initialize oldstatus to keep overly eager compilers to optimize it away.
    oldstatus=.FALSE.
    IF (lqmmm%qmmm)CALL mm_dim(mm_go_mm,oldstatus)

    IF (cntl%tpath.OR.tmw) THEN
       ionode=grandparent
       IF (cntl%tpath)nrepl=pimd3%np_total
       IF (tmw)nrepl=mwi%nwalk
       CALL mp_bcast(nosl%tultra,supersource,supergroup)
       CALL mp_bcast(nosl%tmnose,supersource,supergroup)
    ELSE
       ionode=paral%parent
       nrepl=1
       CALL mp_bcast(nosl%tultra,parai%source,parai%allgrp)
       CALL mp_bcast(nosl%tmnose,parai%source,parai%allgrp)
    ENDIF
    ! keep track of the size that we can bcast data while restarting
    nosl%nose_nrepl = nrepl
    IF (cntl%tnosee) THEN
       ! NN: Two thermostats for HS & BS wavefunction
       ALLOCATE(eta(nchx,bsfac*nrepl),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(etadot(nchx,bsfac*nrepl),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tnosep) THEN
       IF (nosl%tultra) THEN
          ALLOCATE(gkt(2,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(etap(maxsys%nsx,nchx,nrepl),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(etapdot(maxsys%nsx,nchx,nrepl),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(qnosp(maxsys%nsx,nchx,nrepl),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ! (GM) MAXIMUM NUMBER OF LOCAL THERMOSTATS == maxsys%nax*maxsys%nsx
       ELSEIF (loct%tloct) THEN
          ALLOCATE(gkt(2,maxsys%nax*maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(etap(maxsys%nax*maxsys%nsx,nchx,nrepl),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(etapdot(maxsys%nax*maxsys%nsx,nchx,nrepl),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(qnosp(maxsys%nax*maxsys%nsx,nchx,nrepl),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(lcttab(maxsys%nsx,maxsys%nax),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(lctmemb(2,maxsys%nax*maxsys%nsx,loct%nloct),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(nlctmbm(loct%nloct),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(dofsl(loct%nloct),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSEIF (nosl%tmnose) THEN
          ALLOCATE(gkt(2,3*maxsys%nax*maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(etapm(3*maxsys%nax*maxsys%nsx,nchx,nrepl),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(etapmdot(3*maxsys%nax*maxsys%nsx,nchx,nrepl),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(qnospm(3*maxsys%nax*maxsys%nsx,nchx,nrepl),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          nr=MIN(nrepl,2)
          nn=(3*maxsys%nax*maxsys%nsx*nchx*nr)
          ALLOCATE(ndfnt(3*maxsys%nax*maxsys%nsx,nn/(3*maxsys%nax*maxsys%nsx)),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(mapdof(3,3*maxsys%nax*maxsys%nsx,nn/(3*maxsys%nax*maxsys%nsx)),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          IF (tcafes) THEN
             ALLOCATE(tempwr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(wnosepr(3,maxsys%nax,maxsys%nsx),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(wnosep2r(3,maxsys%nax,maxsys%nsx),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(chainer(3,maxsys%nax,maxsys%nsx),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(cafescl(maxsys%nax,maxsys%nsx),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(cafestmp(ncafesgrp),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(cafesnse(ncafesgrp),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ! -----------------------------------------------------------------    
             IF (paral%parent) THEN
                IF (paral%io_parent)&
                     CALL fileopen(177,'CAFES_IN',fo_old,ferror)
                ! CAFES_IN not present use values from input
                IF (ferror) GOTO 113
                IF (paral%io_parent)&
                     WRITE(6,'(/,1X,A,/,1X,A,I4)')&
                     'USING CAFES GROUP DEFINITIONS FROM CAFES_IN FILE:',&
                     'NUMBER OF CAFES GROUPS:',ncafesgrp
                DO is=1,ions1%nsp
                   DO ia=1,ions0%na(is)
                      IF (paral%io_parent)&
                           READ(177,err=112,END=112,fmt=*) tempwr(1,ia,is),&
                           wnosepr(1,ia,is), cafescl(ia,is)
                      tempwr(2,ia,is)=tempwr(1,ia,is)
                      tempwr(3,ia,is)=tempwr(1,ia,is)
                      wnosepr(2,ia,is)=wnosepr(1,ia,is)
                      wnosepr(3,ia,is)=wnosepr(1,ia,is)
                      IF (cafescl(ia,is).GT.ncafesgrp) THEN
                         CALL stopgm('NOSALLOC','INCONSISTENT CAFES GROUPS',& 
                              __LINE__,__FILE__)
                      ENDIF
                   ENDDO
                ENDDO
                IF (paral%io_parent)&
                     CALL fileclose(177)
                GOTO 114
112             CALL stopgm('NOSALLOC','PROBLEMS READING CAFES_IN',& 
                     __LINE__,__FILE__)
113             CONTINUE
                k=0
                IF (paral%io_parent)&
                     WRITE(6,'(/,1X,A,/,1X,A,I4)')&
                     'USING CAFES PARAMETERS FROM THE INPUT FILE',&
                     'NUMBER OF CAFES GROUPS:',ncafesgrp
                IF (paral%io_parent)&
                     WRITE(6,122) 'GROUP', '1ST', 'LAST', 'TEMP', 'FREQ'
                DO i=1,ncafesgrp
                   IF (paral%io_parent)&
                        WRITE(6,121) i,cafesini(1,i),cafesini(2,i),&
                        cafesinr(1,i),cafesinr(2,i)
                ENDDO
                DO is=1,ions1%nsp
                   DO ia=1,ions0%na(is)
                      k=k+1
                      DO i=1,ncafesgrp
                         IF ((cafesini(1,i).LE.k).AND. (cafesini(2,i).GE.k))&
                              THEN
                            cafescl(ia,is)=i
                            tempwr(1,ia,is)=cafesinr(1,i)
                            wnosepr(1,ia,is)=cafesinr(2,i)
                         ENDIF
                      ENDDO
                      tempwr(2,ia,is)=tempwr(1,ia,is)
                      tempwr(3,ia,is)=tempwr(1,ia,is)
                      wnosepr(2,ia,is)=wnosepr(1,ia,is)
                      wnosepr(3,ia,is)=wnosepr(1,ia,is)
                   ENDDO
                ENDDO
114             CONTINUE
             ENDIF
             CALL mp_bcast(cafescl,SIZE(cafescl),parai%source,parai%allgrp)
             CALL mp_bcast(tempwr,SIZE(tempwr),parai%source,parai%allgrp)
             CALL mp_bcast(wnosepr,SIZE(wnosepr),parai%source,parai%allgrp)
          ENDIF
          ! Addtionally allocate arrays for a single NH chain if CMD
          IF (cntl%tpath.AND.cntl%tpimd.AND.pimd1%tcentro) THEN
             DO ip=np_low,np_high
                IF (ip.NE.1) CYCLE
                ALLOCATE(gkt1(2,maxsys%nsx),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ALLOCATE(etap1(nchx,nrepl),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ALLOCATE(etap1dot(nchx,nrepl),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ALLOCATE(qnospc(nchx,nrepl),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
             ENDDO
          ENDIF
       ELSE
          ALLOCATE(gkt(2,maxsys%nsx),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(etap1(nchx,nrepl),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(etap1dot(nchx,nrepl),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(qnospc(nchx,nrepl),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    IF (cntl%tnosec) THEN
       ALLOCATE(etc(nchx,nrepl),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(etcdot(nchx,nrepl),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF (cntl%tnosep.OR.cntl%tnosee.OR.cntl%tnosec) THEN
       IF (ionode) CALL prmem('  NOSALLOC')
    ENDIF

    IF (tcafes.AND.paral%parent) THEN
       IF (paral%io_parent)&
            CALL fileopen(171,'CAFES',fo_app+fo_verb,ferror)
    ENDIF

    IF (lqmmm%qmmm)CALL mm_dim(mm_revert,oldstatus)

120 FORMAT(i12,i12,f12.2,f12.2,f12.4 )
121 FORMAT(i12,i12,i12,f12.2,f12.2)
122 FORMAT(a12,a12,a12,a12,a12)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE nosalloc
  ! ==================================================================

END MODULE nosalloc_utils

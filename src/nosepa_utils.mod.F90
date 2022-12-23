MODULE nosepa_utils
  USE cnst,                            ONLY: factem,&
                                             pi
  USE cotr,                            ONLY: cotc0,&
                                             lskcor,&
                                             ntcnst
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE mm_dimmod,                       ONLY: cpat,&
                                             cpsp,&
                                             gratom,&
                                             mmdim,&
                                             nam,&
                                             nat_cpmd,&
                                             solsolv,&
                                             solvvv
  USE mm_input,                        ONLY: lqmmm
  USE nose,                            ONLY: &
       dofsl, dofst, dofsu, dtsuz, gckt, gkt, lctmemb, lctrng, lcttab, loct, &
       loctpin, mapdof, msuzx, ncalls, ncdof, nchain, nchc, nche, nchx, &
       ndfnt, nedof, nit, nlctmbm, nosl, ntherm, qnoscc, qnosee, qnosp, &
       qnospc, qnospm, tcafes, tempwr, wnosec, wnosee, wnosep, wnosep2r, &
       wnosepr, gkt1, tnosepc
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: grandparent,&
                                             pimd1,&
                                             pimd2,&
                                             pimd3,&
                                             flnm
  USE prcp,                            ONLY: prcpl
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             maxsp,&
                                             maxsys,&
                                             spar
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nosepa
  !public :: noscnst
  !public :: masscnst
  !public :: ilcttab
  !public :: ilctmemb
  !public :: loccnst
  !public :: locgrcns

CONTAINS

  ! ==================================================================
  SUBROUTINE nosepa(np1,np2)
    ! ==--------------------------------------------------------------==
    ! == WARNING: USE DELT_ELEC FOR TIMESTEP                          ==
    ! ==          SO WE ASSUME THAT DELT_ELEC == DELT_IONS FOR cntl%tnosep ==
    ! ==          (CHECK IN CONTROL_TEST)                             ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: np1, np2

    REAL(real_8), PARAMETER                  :: wntau = 7.26e-7_real_8

    INTEGER                                  :: i, ia, iat, imatch, intt, ip, &
                                                ipp, is, j, k, l, m, n1, nc, &
                                                no, ntrue
    REAL(real_8) :: ddof, drnp, grscnst, psuz(msuzx,5), pyosh, w_yosh6_0, &
      w_yosh6_1, w_yosh6_2, w_yosh6_3, w_yosh7_0, w_yosh7_1, w_yosh7_2, &
      w_yosh7_3, w_yosh7_4, w_yosh8_0, w_yosh8_1, w_yosh8_2, w_yosh8_3, &
      w_yosh8_4, w_yosh8_5, w_yosh8_6, w_yosh8_7, wdum, wnosec2, wnosee2, &
      wnosep2, xnt
    LOGICAL                                  :: tmnose_save

! ==--------------------------------------------------------------==
! ==  SET THE PARAMETERS OF THE NOSE-HOOVER THERMOSTATS           ==
! ==--------------------------------------------------------------==

    wnosee=cntr%wnose0*wntau
    wnosep=cntr%wnosp0*wntau
    wnosec=cntr%wnosc0*wntau
    ! constraint contribution for each DOF of the gromos solvent molecules
    grscnst=0.0_real_8
    IF (lqmmm%qmmm.AND.(solvvv%nram_gr.GT.0)) THEN
       grscnst=REAL(solvvv%ncons_gr,kind=real_8)/REAL(solvvv%nram_gr,kind=real_8)/3.0_real_8
    ENDIF
    ! set parameters for all CAFES chains.
    IF (loct%tloct) THEN
       CALL ilcttab
       CALL ilctmemb
       CALL zeroing(dofsl)!,loct%nloct)
    ENDIF
    IF (tcafes) THEN
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO k=1,3
                wnosepr(k,ia,is)=wnosepr(k,ia,is)*wntau
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    nit=cnti%nit0
    nchain=cnti%nchp
    nche=cnti%nchs
    nchc=cnti%nchb
    ncalls=cnti%ncalls0
    ! ..ELECTRON NOSE PARAMETERS
    IF (cntl%tnosee) THEN
       n1=crge%n*(2*spar%ngws-1)
       nc=crge%n*crge%n
       IF (cntl%nonort) nc=crge%n
       ntrue=(n1-nc)
       ! if NEDOF<0 use true number of DOFs
       IF ( cntr%nedof0.LE.0._real_8 ) THEN
          nedof=NINT(ntrue*ABS(cntr%nedof0))
       ELSE
          nedof=NINT(crge%n*cntr%nedof0)
       ENDIF
       wnosee2 = wnosee*wnosee
       qnosee(1) = 2.0_real_8*cntr%ekinw/wnosee2
       DO i=2,nche
          qnosee(i) = 2._real_8*cntr%ekinw/wnosee2/REAL(nedof,kind=real_8)
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,I10)')&
            ' NOSEPA| TOTAL # OF ELECTRONIC DEGREES OF FREEDOM: ',ntrue
       IF (paral%io_parent)&
            WRITE(6,'(A,T56,I10)')&
            ' NOSEPA| USED # OF ELECTRONIC DEGREES OF FREEDOM : ',nedof
    ENDIF
    ! ..DEGREES OF FREEDOM
    dofst=0.0_real_8
    CALL zeroing(dofsu)!,maxsp)
    iat=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          DO k=1,3
             IF (lskcor(k,iat).NE.0) THEN
                dofst=dofst+1._real_8
                dofsu(is)=dofsu(is)+1._real_8
                ! .. subtract lost DOFs due to gromos solvent constraints 
                ! ES this we should do at a later stage, otherwise
                ! ES the if block after 'DDOF=abs(3._real_8*real(NAT,kind=real_8)-DOFST)'
                ! ES will not work
                IF (lqmmm%qmmm) THEN
                   IF (iat.GT.solsolv%nrpt) THEN
                      ! ES                  DOFST=DOFST-GRSCNST
                      dofsu(is)=dofsu(is)-grscnst
                   ENDIF
                ENDIF
                IF (loct%tloct) THEN
                   IF (lcttab(is,ia).GT.0)  THEN
                      dofsl(lcttab(is,ia))=dofsl(lcttab(is,ia))+1._real_8
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    IF (nosl%tultra) CALL noscnst
    IF (loct%tloct)  CALL loccnst
    ddof=ABS(3._real_8*REAL(ions1%nat,kind=real_8)-dofst)
    IF (.NOT.(cntl%tpath.AND.cntl%tpimd)) THEN
       IF (ddof.LT.1.e-2_real_8) THEN
          dofst=dofst-3._real_8
          IF (isos1%tisos.AND..NOT.lqmmm%qmmm) THEN
             IF (ions1%nat.EQ.1) THEN
                dofst=dofst
             ELSEIF (ions1%nat.EQ.2) THEN
                dofst=dofst-2._real_8
             ELSE
                dofst=dofst-3._real_8
             ENDIF
          ENDIF
       ENDIF
    ELSE
       ! ..Momenta are not conserved in path integral runs
       ! No DOFST-3.0
       drnp=SQRT(REAL(pimd3%np_total,kind=real_8))
       ! DM1
       ! ..Classical test
       IF (pimd1%testpi) THEN
          IF (ddof.LT.1.e-2_real_8) THEN
             dofst=dofst-3._real_8
             ! DM1        IF(TISOS) DOFST=DOFST-3._real_8
             IF (isos1%tisos) THEN
                IF (ions1%nat.EQ.1) THEN
                   dofst=dofst
                ELSEIF (ions1%nat.EQ.2) THEN
                   dofst=dofst-2._real_8
                ELSE
                   dofst=dofst-3._real_8
                ENDIF
             ENDIF
             ! DM2
          ENDIF
       ENDIF
       ! DM2
    ENDIF

    ! consider (CPMD) constraints.
    dofst=dofst-REAL(cotc0%mcnstr,kind=real_8)

    IF (lqmmm%qmmm.AND.(nosl%tultra.OR.(nosl%tmnose.AND..NOT.tcafes))&
         .AND.(solsolv%nsolv*solvvv%ncons_gr.NE.0)) THEN
       CALL stopgm('NOSEPA',&
            'NO SOLVENT CONSTRAINTS WITH NOSE ULTRA/MASSIVE',& 
            __LINE__,__FILE__)
    ELSE
       ! consider (classical) solvent constraints.
       ! these are already handled in Gromos and thus not in MCNSTR.
       IF (lqmmm%qmmm.AND.(solvvv%nram_gr*solvvv%ncons_gr.GT.0)) THEN
          dofst=dofst-REAL(solsolv%nsolv/solvvv%nram_gr*solvvv%ncons_gr,kind=real_8)
          IF (loct%tloct) CALL locgrcns
       ENDIF
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,'(A,T56,I10)')&
         ' NOSEPA| USED # OF NUCLEAR DEGREES OF FREEDOM :    ',&
         INT(dofst)
    ! 
    ! ..ION NOSE PARAMETERS
    DO ip=np1,np2
       IF (cntl%tnosep) THEN
          IF (cntl%tpath.AND.cntl%tpimd) THEN
             wnosep2 = (drnp*cntr%tempw/factem)**2
             IF ((pimd1%tpinm.OR.pimd1%tstage).AND.ip.EQ.1) wnosep2 = wnosep*wnosep
             IF (pimd1%tringp.AND.pimd1%tpinm.AND.ip.GT.1) wnosep2 = wnosep2*flnm(ip)
             IF (pimd1%tcentro.AND.(pimd1%tpinm.OR.pimd1%tstage)&
                 .AND.ip.GT.1) wnosep2 = wnosep2*pimd2%facstage
             ipp=MIN(ip,2)
             ! DM1
             ! ..Classical test
             IF (pimd1%testpi) THEN
                wnosep2 = wnosep*wnosep
                ipp=1
             ENDIF
          ELSE
             wnosep2 = wnosep*wnosep
             ipp=1
          ENDIF
          ! switch off massive NHCs of centroids (ip=1) for CMD
          tmnose_save=nosl%tmnose
          IF (cntl%tpath.AND.cntl%tpimd.AND.&
             pimd1%tcentro.AND.ip.EQ.1) nosl%tmnose=.FALSE.
          IF (nosl%tultra) THEN
             ! ..EACH ATOM TYPE HAS A SEPERATE CHAIN
             IF (ipp.GT.1) THEN
                DO is=1,ions1%nsp
                   gkt(ipp,is) = 3._real_8*ions0%na(is)*cntr%tempw/factem
                   qnosp(is,1,ip) = gkt(ipp,is)/wnosep2
                   DO l=2,nchain
                      qnosp(is,l,ip) = gkt(ipp,is)/wnosep2/(3._real_8*ions0%na(is))
                   ENDDO
                ENDDO
             ELSE
                DO is=1,ions1%nsp
                   gkt(ipp,is) = dofsu(is)*cntr%tempw/factem
                   qnosp(is,1,ip) = gkt(ipp,is)/wnosep2
                   DO l=2,nchain
                      qnosp(is,l,ip) = gkt(ipp,is)/wnosep2/dofsu(is)
                   ENDDO
                ENDDO
             ENDIF
          ELSEIF (loct%tloct) THEN
             IF (ipp.GT.1) THEN
                DO is=1,loct%nloct
                   wnosep=loctpin(2,is)*wntau
                   wnosep2=wnosep*wnosep
                   gkt(ipp,is) = 3._real_8*nlctmbm(is)*loctpin(1,is)/factem
                   qnosp(is,1,ip) = gkt(ipp,is)/wnosep2
                   DO l=2,nchain
                      qnosp(is,l,ip) = gkt(ipp,is)/wnosep2/&
                           (3._real_8*nlctmbm(is))
                   ENDDO
                ENDDO
             ELSE
                DO is=1,loct%nloct
                   IF (dofsl(is).LE.1.0e-6) THEN
                      CALL stopgm('NOSEPARA','LOCAL DOF NOT GT THAN 0',& 
                           __LINE__,__FILE__)
                   ENDIF
                   wnosep=loctpin(2,is)*wntau
                   wnosep2=wnosep*wnosep
                   gkt(ipp,is) = dofsl(is)*loctpin(1,is)/factem
                   qnosp(is,1,ip) = gkt(ipp,is)/wnosep2
                   DO l=2,nchain
                      qnosp(is,l,ip) = gkt(ipp,is)/wnosep2/dofsl(is)
                   ENDDO
                ENDDO
             ENDIF
          ELSEIF (nosl%tmnose) THEN
             IF (tcafes) THEN
                ntherm(ipp)=3*ions1%nat
                intt=0
                DO is=1,ions1%nsp
                   DO ia=1,ions0%na(is)
                      DO k=1,3
                         wnosep2r(k,ia,is)=wnosepr(k,ia,is)**2
                         intt=intt+1
                         ndfnt(intt,ipp)=1
                         mapdof(1,intt,ipp)=k
                         mapdof(2,intt,ipp)=ia
                         mapdof(3,intt,ipp)=is
                         gkt(ipp,intt) = tempwr(k,ia,is)/factem
                         qnospm(intt,1,ip) = gkt(ipp,intt)/wnosep2r(k,ia,is)
                         DO l=2,nchain
                            qnospm(intt,l,ip) = gkt(ipp,intt)/wnosep2r(k,ia,is)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                CALL masscnst(ip,ipp,wnosep2)
             ENDIF
          ELSE
             IF (ipp.GT.1) THEN
                gkt(ipp,1) = 3._real_8*ions1%nat*cntr%tempw/factem
                qnospc(1,ip) = gkt(ipp,1)/wnosep2
                DO l=2,nchain
                   qnospc(l,ip) = gkt(ipp,1)/wnosep2/(3._real_8*ions1%nat)
                ENDDO
             ELSE
                gkt(ipp,1) = dofst*cntr%tempw/factem
                qnospc(1,ip) = gkt(ipp,1)/wnosep2
                DO l=2,nchain
                   qnospc(l,ip) = gkt(ipp,1)/wnosep2/dofst
                ENDDO
             ENDIF
          ENDIF
          IF (ALLOCATED(gkt1)) CALL dcopy(2*maxsys%nsx,gkt(1,1),1,gkt1(1,1),1)
          nosl%tmnose=tmnose_save
       ENDIF
    ENDDO
    IF (cntl%tnosec) THEN
       ncdof = 9
       IF (prcpl%tisot) ncdof = 1
       IF (prcpl%tzflex) ncdof = 1
       gckt = REAL(ncdof,kind=real_8)*cntr%tempc/factem
       wnosec2 = wnosec*wnosec
       qnoscc(1) = REAL(ncdof,kind=real_8)*cntr%tempc/(factem*wnosec2)
       DO l=2,nchc
          qnoscc(l) = cntr%tempc/(factem*wnosec2)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cntl%tpath.AND.cntl%tpimd.AND.cntl%tnosep) THEN
       IF (grandparent) THEN
          IF (paral%io_parent)&
               WRITE(6,*)
          IF (paral%io_parent)&
               WRITE(6,'(A,A)') ' *********************************',&
               '*******************************'
          IF (paral%io_parent)&
               WRITE(6,'(A,A)') ' THERMOSTATS: CHANGE IONIC NOSE',&
               ' FREQUENCIES IN PATH INTEGRAL CASE'
          IF (pimd1%tpinm.OR.pimd1%tstage) THEN
             IF (pimd1%tcentro) THEN
                IF ((.NOT.nosl%tmnose).AND.paral%io_parent)&
                     WRITE(6,'(/,A,A,/)')&
                     ' WARNING FROM NOSEPA: YOU SHOULD USE MASSIVE',&
                     ' THERMOSTATTING FOR CENTROID PIMD ! !'

                IF (paral%io_parent)&
                     WRITE(6,'(A,A)') ' FOR STAGING OR NORMAL MODE PROPAGATOR',&
                     ' WITH CENTROID DYNAMICS'
                IF (.NOT.tnosepc) THEN
                   IF (paral%io_parent)&
                     WRITE(6,'(A,A)')&
                     '    NO THERMOSTAT ON FIRST BEAD  IP=1 '
                ELSE 
                   wdum=wnosep/(wntau*2._real_8*pi)
                   IF (paral%io_parent)&
                     WRITE(6,'(A,F9.2,A)')&
                     '    CHARACTERISTIC FREQUENCY FOR IP=1 :',wdum,' CM**-1 (SINGLE)'
                ENDIF
                wdum=(drnp*cntr%tempw/factem*SQRT(pimd2%facstage))/(wntau*2._real_8*pi)
                IF (paral%io_parent)&
                     WRITE(6,'(A,F9.2,A)')&
                     '    CHARACTERISTIC FREQUENCY FOR IP>1 :',wdum,' CM**-1 '
             ELSE IF (pimd1%tringp) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(A,A)') ' FOR NORMAL MODE PROPAGATOR',&
                     ' WITH RING-POLYMER DYNAMICS'
                IF (.NOT.tnosepc) THEN
                   IF (paral%io_parent)&
                     WRITE(6,'(A,A)')&
                     '    NO THERMOSTAT ON FIRST BEAD  IP=1 '
                ELSE
                   wdum=wnosep/(wntau*2._real_8*pi)
                   IF (paral%io_parent)&
                     WRITE(6,'(A,F9.2,A)')&
                     '    CHARACTERISTIC FREQUENCY FOR IP=1 :',wdum,' CM**-1 '
                ENDIF
                DO ip=2,pimd3%np_total
                   wdum=(drnp*cntr%tempw/factem*SQRT(flnm(ip)))/(wntau*2._real_8*pi)
                   IF (paral%io_parent) &
                   WRITE(6,'(A,I3,A,F9.2,A)') &
                   '    CHARACTERISTIC FREQUENCY FOR IP=',ip ,':',wdum,' CM**-1 '
                ENDDO
             ELSE
                IF (paral%io_parent)&
                     WRITE(6,'(A)') ' FOR STAGING OR NORMAL MODE PROPAGATOR'
                wdum=wnosep/(wntau*2._real_8*pi)
                IF (paral%io_parent)&
                     WRITE(6,'(A,F12.2,A)')&
                     '    CHARACTERISTIC FREQUENCY FOR IP=1 :',wdum,' CM**-1 '
                wdum=(drnp*cntr%tempw/factem)/(wntau*2._real_8*pi)
                IF (paral%io_parent)&
                     WRITE(6,'(A,F12.2,A)')&
                     '    CHARACTERISTIC FREQUENCY FOR IP>1 :',wdum,' CM**-1 '
             ENDIF
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(A)') ' FOR PRIMITIVE PROPAGATOR'
             wdum=(drnp*cntr%tempw/factem)/(wntau*2._real_8*pi)
             IF (paral%io_parent)&
                  WRITE(6,'(A,F12.2,A)')&
                  '    CHARACTERISTIC FREQUENCY FOR IP=1 :',wdum,' CM**-1 '
             wdum=(drnp*cntr%tempw/factem)/(wntau*2._real_8*pi)
             IF (paral%io_parent)&
                  WRITE(6,'(A,F12.2,A)')&
                  '    CHARACTERISTIC FREQUENCY FOR IP>1 :',wdum,' CM**-1 '
             IF (pimd1%tringp.AND..NOT.tnosepc) THEN
                IF (paral%io_parent)&
                    WRITE(6,'(A)')&
                    '    WARNING FROM NOSEPA: CENTROIDOFF IS MEANINGFUL ONLY FOR NORMAL MODE PROPAGATOR'
             ENDIF
          ENDIF
          IF (paral%io_parent)&
               WRITE(6,'(A,A)') ' *********************************',&
               '*******************************'
          IF (paral%io_parent)&
               WRITE(6,*)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (nchain.GT.nchx) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' NOSE-HOOVER CHAIN TOO LONG (IONS) ',nchain,nchx
       CALL stopgm('NOSEPA',' ',& 
            __LINE__,__FILE__)
    ENDIF
    IF (nche.GT.nchx) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' NOSE-HOOVER CHAIN TOO LONG (ELECTRONS) ',nche,nchx
       CALL stopgm('NOSEPA',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ..THE FOLLOWING IS RELATED TO SUZUKI-NOSE INTEGRATION
    DO no = 2,msuzx
       xnt = 1._real_8/(2._real_8*REAL(no,kind=real_8) - 1._real_8)
       DO i = 1,5
          psuz(no,i) = 1._real_8/(4._real_8 - 4._real_8**xnt)
          IF (MOD(i,3).EQ.0)&
               psuz(no,i) = 1._real_8 - 4._real_8/(4._real_8 - 4._real_8**xnt)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ..MATCH FLAG FOR NUMBER OF TIME STEPS
    imatch=-1
    ! ==--------------------------------------------------------------==
    IF (ncalls.EQ.3) THEN
       pyosh = 1.0_real_8/(2.0_real_8 - 2.0_real_8**(1._real_8/3._real_8))
       dtsuz(1) = pyosh*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(2) = -(2._real_8**(1._real_8/3._real_8))*pyosh*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(3) = pyosh*cntr%delt_elec/REAL(nit,kind=real_8)
       ncalls = 3
       imatch=1
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (ncalls.EQ.5) THEN
       k=0
       DO i=1,5
          k=k+1
          dtsuz(k) = psuz(2,i)*cntr%delt_elec/REAL(nit,kind=real_8)
       ENDDO
       imatch=1
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (ncalls.EQ.7) THEN
       w_yosh6_1 =  0.784513610477560_real_8
       w_yosh6_2 =  0.235573213359357_real_8
       w_yosh6_3 = -0.117767998417887e1_real_8
       w_yosh6_0 = 1.0_real_8 - 2.0_real_8*(w_yosh6_1+w_yosh6_2+w_yosh6_3)
       dtsuz(1) = w_yosh6_3*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(2) = w_yosh6_2*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(3) = w_yosh6_1*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(4) = w_yosh6_0*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(5) = w_yosh6_1*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(6) = w_yosh6_2*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(7) = w_yosh6_3*cntr%delt_elec/REAL(nit,kind=real_8)
       imatch=1
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (ncalls.EQ.9) THEN
       w_yosh7_1 = 0.192_real_8
       w_yosh7_2 =  0.554910818409783619692725006662999_real_8
       w_yosh7_3 =  0.124659619941888644216504240951585_real_8
       w_yosh7_4 =  -0.843182063596933505315033808282941_real_8
       w_yosh7_0 = 1.0_real_8 - 2.0_real_8*(w_yosh7_1+w_yosh7_2+&
            w_yosh7_3+w_yosh7_4)
       dtsuz(1) = w_yosh7_4*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(2) = w_yosh7_3*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(3) = w_yosh7_2*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(4) = w_yosh7_1*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(5) = w_yosh7_0*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(6) = w_yosh7_1*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(7) = w_yosh7_2*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(8) = w_yosh7_3*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(9) = w_yosh7_4*cntr%delt_elec/REAL(nit,kind=real_8)
       imatch=1
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (ncalls.EQ.15) THEN
       w_yosh8_1 =  0.102799849391985_real_8
       w_yosh8_2 = -0.196061023297549e1_real_8
       w_yosh8_3 =  0.193813913762276e1_real_8
       w_yosh8_4 = -0.158240635368243_real_8
       w_yosh8_5 = -0.144485223686048e1_real_8
       w_yosh8_6 =  0.253693336566229_real_8
       w_yosh8_7 =  0.914844246229740_real_8
       w_yosh8_0 = 1.0_real_8 - 2.0_real_8*(w_yosh8_1+w_yosh8_2+w_yosh8_3&
            +w_yosh8_4+w_yosh8_5+w_yosh8_6&
            +w_yosh8_7)
       dtsuz(1) = w_yosh8_7*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(2) = w_yosh8_6*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(3) = w_yosh8_5*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(4) = w_yosh8_4*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(5) = w_yosh8_3*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(6) = w_yosh8_2*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(7) = w_yosh8_1*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(8) = w_yosh8_0*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(9) = w_yosh8_1*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(10) = w_yosh8_2*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(11) = w_yosh8_3*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(12) = w_yosh8_4*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(13) = w_yosh8_5*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(14) = w_yosh8_6*cntr%delt_elec/REAL(nit,kind=real_8)
       dtsuz(15) = w_yosh8_7*cntr%delt_elec/REAL(nit,kind=real_8)
       imatch=1
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (ncalls.EQ.25) THEN
       k=0
       DO j=1,5
          DO i=1,5
             k=k+1
             dtsuz(k) = psuz(2,i)*psuz(3,j)*cntr%delt_elec/REAL(nit,kind=real_8)
          ENDDO
       ENDDO
       imatch=1
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (ncalls.EQ.125) THEN
       k=0
       DO l=1,5
          DO j=1,5
             DO i=1,5
                k=k+1
                dtsuz(k) = psuz(2,i)*psuz(3,j)*psuz(4,l)*&
                     cntr%delt_elec/REAL(nit,kind=real_8)
             ENDDO
          ENDDO
       ENDDO
       imatch=1
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (ncalls.EQ.625) THEN
       k=0
       DO m=1,5
          DO l=1,5
             DO j=1,5
                DO i=1,5
                   k=k+1
                   dtsuz(k) = psuz(2,i)*psuz(3,j)*psuz(4,l)*psuz(5,m)*&
                        cntr%delt_elec/REAL(nit,kind=real_8)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       imatch=1
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (imatch.LT.0) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ncalls,' SUZUKI/YOSHIDA STEPS NOT PROGRAMMED.'
       IF (paral%io_parent)&
            WRITE(6,*) 'SORRY, DUDE'
       IF (paral%io_parent)&
            WRITE(6,*) 'POSSIBLE CHOICES FOR NCALLS ARE '
       IF (paral%io_parent)&
            WRITE(6,*) '35 7 915 25125 625'
       CALL stopgm('NOSEPA',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nosepa
  ! ==================================================================
  SUBROUTINE noscnst
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'noscnst'

    INTEGER                                  :: i, ia, iat, ierr, is, is1, &
                                                is2, is3, is4
    INTEGER, ALLOCATABLE                     :: n1(:)
    REAL(real_8)                             :: ddc(maxsp)

! ==--------------------------------------------------------------==

    IF (cotc0%mcnstr.EQ.0) RETURN
    ALLOCATE(n1(3*maxsys%nax*ions1%nsp+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    iat=0
    n1(1)=0
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          iat=iat+1
          n1(iat+1)=is
       ENDDO
       ddc(is)=0._real_8
    ENDDO
    DO i=1,cotc0%mcnstr
       is1=n1(ntcnst(2,i)+1)
       ddc(is1)=ddc(is1)+1._real_8
       is2=n1(ntcnst(3,i)+1)
       is3=n1(ntcnst(4,i)+1)
       is4=n1(ntcnst(5,i)+1)
       IF (is1.NE.is2.AND.is2.NE.0) GOTO 100
       IF (is1.NE.is3.AND.is3.NE.0) GOTO 100
       IF (is1.NE.is4.AND.is4.NE.0) GOTO 100
    ENDDO
    DO is=1,ions1%nsp
       dofsu(is)=dofsu(is)-ddc(is)
    ENDDO
    DEALLOCATE(n1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
100 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' CONSTRAINTS ONLY ON EQUAL ATOM TYPES ALLOWED '
    IF (paral%io_parent)&
         WRITE(6,*) ' WITH THIS TYPE OF NOSE THERMOSTATS  '
    CALL stopgm('NOSCNST','NOSE + CONSTRAINTS',& 
         __LINE__,__FILE__)
  END SUBROUTINE noscnst
  ! ==================================================================
  SUBROUTINE masscnst(ip,ipp,wnosep2)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ip, ipp
    REAL(real_8)                             :: wnosep2

    CHARACTER(*), PARAMETER                  :: procedureN = 'masscnst'

    INTEGER                                  :: i, ia, iat, id, ierr, intt, &
                                                ipn, is, k, l, nt
    INTEGER, ALLOCATABLE                     :: nxx(:)

    IF (cotc0%mcnstr.EQ.0.OR.ipp.GT.1) THEN
       ntherm(ipp)=3*ions1%nat
       intt=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             DO k=1,3
                intt=intt+1
                ndfnt(intt,ipp)=1
                mapdof(1,intt,ipp)=k
                mapdof(2,intt,ipp)=ia
                mapdof(3,intt,ipp)=is
                gkt(ipp,intt) = cntr%tempw/factem
                qnospm(intt,1,ip) = gkt(ipp,intt)/wnosep2
                DO l=2,nchain
                   qnospm(intt,l,ip) = gkt(ipp,intt)/wnosep2
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ELSE
       ALLOCATE(nxx(3*maxsys%nax*maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(nxx)!,3*maxsys%nax*maxsys%nsx)
       ! ..mark all atoms involved in constraints
       DO i=1,cotc0%mcnstr
          IF (ntcnst(2,i).NE.0) nxx(ntcnst(2,i))=1
          IF (ntcnst(3,i).NE.0) nxx(ntcnst(3,i))=1
          IF (ntcnst(4,i).NE.0) nxx(ntcnst(4,i))=1
          IF (ntcnst(5,i).NE.0) nxx(ntcnst(5,i))=1
       ENDDO

       nt=0
       DO ia=1,ions1%nat
          nt=nt+nxx(ia)
       ENDDO
       ! ..assume that atoms involved in constraints have no fixed coordinates
       ndfnt(1,ipp)=3*nt
       gkt(ipp,1)=(3*nt-cotc0%mcnstr)*cntr%tempw/factem
       qnospm(1,1,ip) = gkt(ipp,1)/wnosep2
       DO l=2,nchain
          qnospm(1,l,ip) = gkt(ipp,1)/wnosep2/(3*nt-cotc0%mcnstr)
       ENDDO
       ! ..other degrees of freedom
       ntherm(ipp)=1
       id=0
       iat=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             iat=iat+1
             IF (nxx(iat).EQ.0) THEN
                DO k=1,3
                   ntherm(ipp)=ntherm(ipp)+1
                   ipn=ntherm(ipp)
                   ndfnt(ipn,ipp)=1
                   intt=ipn+ndfnt(1,ipp)-1
                   mapdof(1,intt,ipp)=k
                   mapdof(2,intt,ipp)=ia
                   mapdof(3,intt,ipp)=is
                   gkt(ipp,ipn) = cntr%tempw/factem
                   qnospm(ipn,1,ip) = gkt(ipp,ipn)/wnosep2
                   DO l=2,nchain
                      qnospm(ipn,l,ip) = gkt(ipp,ipn)/wnosep2
                   ENDDO
                ENDDO
             ELSE
                DO k=1,3
                   id=id+1
                   mapdof(1,id,ipp)=k
                   mapdof(2,id,ipp)=ia
                   mapdof(3,id,ipp)=is
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       DEALLOCATE(nxx,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE masscnst
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE ilcttab
    ! ==--------------------------------------------------------------==
    ! == BUILD MAP CPMD NO -> LOCAL THERMOSTAT NO.                    ==
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: ia, ias, ica, idxlct, ile, &
                                                isp, istart, istop

! ==--------------------------------------------------------------==

    CALL zeroing(lcttab)!,maxsys%nsx*maxsys%nax)
    DO ile=1,loct%nlocrng
       idxlct=lctrng(1,ile)
       istart=lctrng(2,ile)
       istop= lctrng(3,ile)
       IF (istart.GT.istop) CALL stopgm('LOCAL THERMOSTATS',&
            'RANGE START GREATER THAN RANGE STOP',& 
            __LINE__,__FILE__)
       DO ia=istart,istop
          ica=NAT_cpmd(ia)
          isp=cpsp(ia)
          ias=cpat(ia)
          IF (lcttab(isp,ias).NE.0) THEN
             IF (lcttab(isp,ias).EQ.idxlct) THEN
                IF (paral%io_parent)&
                     WRITE(6,'(1X,A,I6,1X,A)')&
                     'LOCAL TEMPERATURE WARNING: ATOM',ia,&
                     'BELONGS TO MORE THAN ONE THERMOSTAT RANGE'
             ELSE
                CALL stopgm('LOCAL TEMPERATURE',&
                     'FOUND ATOM BELONGING TO TWO THERMOSTATS',& 
                     __LINE__,__FILE__)
             ENDIF
          ELSE
             lcttab(isp,ias)=idxlct
          ENDIF

       ENDDO

    ENDDO
    ! DO ISP=1,NSP
    ! DO IAS=1,NA(ISP)
    ! WRITE(6,*) 'LCTTAB(',ISP,',',IAS,') GR_no',
    ! *           gratom((ISP-1)*maxsys%nax+IAS),' LOCTST', LCTTAB(ISP,IAS)
    ! ENDDO
    ! ENDDO

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ilcttab
  ! ==================================================================
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE ilctmemb
    ! ==--------------------------------------------------------------==
    ! == LIST ATOMS BELONGING TO THERMOSTAT NO.                       ==
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: ia, iag, idxlct, is, k

! ==--------------------------------------------------------------==

    CALL zeroing(lctmemb)!,2*maxsys%nsx*maxsys%nax*loct%nloct)
    CALL zeroing(nlctmbm)!,loct%nloct)

    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          idxlct=lcttab(is,ia)
          IF (idxlct.GT.0) THEN
             nlctmbm(idxlct)=nlctmbm(idxlct)+1
             lctmemb(1,nlctmbm(idxlct),idxlct)=is
             lctmemb(2,nlctmbm(idxlct),idxlct)=ia
          ENDIF
       ENDDO
    ENDDO
    ! DEBUG 
    DO idxlct=1,loct%nloct
       ! WRITE(6,*)PARENT,'ATOMS IN THERMOSTATNO.',IDXLCT,NLCTMBM(IDXLCT)
       DO k=1,nlctmbm(idxlct)
          is=lctmemb(1,k,idxlct)
          ia=lctmemb(2,k,idxlct)
          iag=gratom((is-1)*maxsys%nax+ia)
          ! WRITE(6,'(I5)') IAG
       ENDDO
       ! PRINT *
    ENDDO


    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ilctmemb
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE loccnst
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'loccnst'

    INTEGER                                  :: i, ia, iat, ierr, il1, il2, &
                                                il3, il4, is
    INTEGER, ALLOCATABLE                     :: n1(:)
    REAL(real_8)                             :: ddc(loct%nloct)

! ==--------------------------------------------------------------==

    IF (cotc0%mcnstr.EQ.0) RETURN
    ALLOCATE(n1(3*maxsys%nax*mmdim%nspm+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ddc)!,loct%nloct)
    iat=0
    n1(1)=0
    DO is=1,mmdim%nspm
       DO ia=1,NAm(is)
          iat=iat+1
          n1(iat+1)=lcttab(is,ia)
       ENDDO
    ENDDO
    DO i=1,cotc0%mcnstr
       il1=n1(ntcnst(2,i)+1)
       ddc(il1)=ddc(il1)+1._real_8
       il2=n1(ntcnst(3,i)+1)
       il3=n1(ntcnst(4,i)+1)
       il4=n1(ntcnst(5,i)+1)
       IF (il1.NE.il2.AND.il2.NE.0) GOTO 200
       IF (il1.NE.il3.AND.il3.NE.0) GOTO 200
       IF (il1.NE.il4.AND.il4.NE.0) GOTO 200
    ENDDO
    DO i=1,loct%nloct
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,F8.1,1X,A,I4)') 'FOUND ',ddc(i)&
            ,'CPMD CONSTRAINTS FOR THERMOSTAT',i
       dofsl(i)=dofsl(i)-ddc(i)
    ENDDO
    DEALLOCATE(n1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
200 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' CONSTRAINTS ONLY ON EQUAL LOCAL THERMOSTATS ALLOWED '
    IF (paral%io_parent)&
         WRITE(6,*) ' WITH THIS TYPE OF NOSE THERMOSTATS  '
    CALL stopgm('NOSCNST','NOSE + CONSTRAINTS',& 
         __LINE__,__FILE__)
  END SUBROUTINE loccnst
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE locgrcns
    ! ==--------------------------------------------------------------==
    ! CORRECT LOCAL DOF FOR SOLVENT CONSTRAINTS

    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'locgrcns'

    INTEGER                                  :: i, ia, iat, idxlct, ierr, is, &
                                                isl, lastidx
    INTEGER, ALLOCATABLE                     :: n1(:)
    REAL(real_8)                             :: ddc(loct%nloct)

! ==--------------------------------------------------------------==

    ALLOCATE(n1(maxsys%nax*mmdim%nspm+1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ddc)!,loct%nloct)

    iat=0
    n1(1)=0
    DO is=1,mmdim%nspm
       DO ia=1,NAm(is)
          iat=iat+1
          n1(iat)=lcttab(is,ia)
          ! WRITE(6,*)is,ia,iat,n1(iat)
       ENDDO
    ENDDO

    DO isl=0,solsolv%nsolv/solvvv%nram_gr-1
       lastidx=0
       DO i=1,solvvv%nram_gr
          iat=NAT_cpmd(solsolv%nrpt+isl*solvvv%nram_gr+i)
          idxlct=n1(iat)
          ! WRITE(6,*)isl,i,NRPT+ISL*NRAM_gr+I,iat,idxlct
          IF (i.GT.1.AND.lastidx.NE.idxlct) THEN
             ! ATOMS OF THE SOLVENT DO NOT BELONG TO THE SAME THERMOSTAT
             GOTO 201
          ENDIF
          lastidx=idxlct
       ENDDO
       ddc(idxlct)= ddc(idxlct)+REAL(solvvv%ncons_gr,kind=real_8)
    ENDDO

    DO i=1,loct%nloct
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,F8.1,1X,A,I4)') 'FOUND ',ddc(i),&
            'SOLVENT CONSTRAINTS FOR THERMOSTAT',i
       dofsl(i)=dofsl(i)-ddc(i)
    ENDDO
    DEALLOCATE(n1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
201 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*)solsolv%nrpt,isl,i,iat,solsolv%nrpt+isl*solvvv%nram_gr+i,idxlct,lastidx
    IF (paral%io_parent)&
         WRITE(6,*) ' SOLVENT MOLECULES ONLY IN  EQUAL LOCAL THERMOSTATS '
    IF (paral%io_parent)&
         WRITE(6,*) ' ALLOWED. CHECK YOUR INPUT'
    CALL stopgm('NOSCNST','NOSE + SOLVENT CONSTRAINTS',& 
         __LINE__,__FILE__)
  END SUBROUTINE locgrcns
  ! ==================================================================

END MODULE nosepa_utils

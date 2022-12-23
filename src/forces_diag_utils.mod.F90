MODULE forces_diag_utils
  USE andp,                            ONLY: rin0,&
                                             rmix,&
                                             rout0
  USE andr,                            ONLY: andr2
  USE calc_alm_utils,                  ONLY: calc_alm
  USE cdft_utils,                      ONLY: cdft_w,&
                                             v_pred,&
                                             vupdate
  USE cdftmod,                         ONLY: cdftci,&
                                             cdftcom,&
                                             cdftlog
  USE cnstfc_utils,                    ONLY: restfc
  USE cotr,                            ONLY: cotr007,&
                                             gsrate,&
                                             resval,&
                                             resval_dest
  USE ehpsi_utils,                     ONLY: set_b2l
  USE elct,                            ONLY: crge
  USE elct2,                           ONLY: tfixo
  USE ener,                            ONLY: ener_com
  USE fint,                            ONLY: fint1
  USE forcedr_driver,                  ONLY: forcedr
  USE k_updwf_utils,                   ONLY: give_scr_kupdwf,&
                                             k_updwf
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_walltime
  USE mm_extrap,                       ONLY: numcold
  USE mp_interface,                    ONLY: mp_bcast
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE rhoofr_c_utils,                  ONLY: rhoofr_c
  USE rhoofr_utils,                    ONLY: rhoofr
  USE rinitwf_driver,                  ONLY: rinitwf
  USE rnlsm_utils,                     ONLY: rnlsm
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE rrandd_utils,                    ONLY: rrandd
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: store1
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpbl,&
                                             nkpt
  USE testex_utils,                    ONLY: testex
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt_ions
  USE updrho_utils,                    ONLY: give_scr_updrho,&
                                             updrho
  USE updwf_utils,                     ONLY: give_scr_updwf,&
                                             updwf
  USE wrener_utils,                    ONLY: wrprint_wfopt
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: forces_diag
  PUBLIC :: updiag
  PUBLIC :: dens_for_diag
  PUBLIC :: give_scr_forces_diag

CONTAINS

  ! ==================================================================
  SUBROUTINE forces_diag(nstate,c0,c2,cr,sc0,cscr,vpp,eigv,&
       rhoe,psi,&
       tau0,velp,taui,fion,ifcalc,&
       irec,tfor,tinfo)
    ! ==--------------------------------------------------------------==
    ! == ENTER WITH THE IONIC POSITIONS IN TAU0 AND AN INPUT GUESS    ==
    ! == FOR THE DENSITY IN RIN0, THIS ROUTINE RETURNS THE CONVERGED  ==
    ! == DENSITY, WAVEFUNCTIONS AND IONIC FORCES.                     ==
    ! ==--------------------------------------------------------------==
    ! == TAU0: ATOMIC POSITION                                        ==
    ! == VELP and TAUI are necessary for RESTART FILE                 ==
    ! == FION: IONIC FORCES                                           ==
    ! == IFCALC: total number of iterations                           ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(:,:,:), &
                                                c2(nkpt%ngwk,nstate), cr(*), &
                                                sc0(nkpt%ngwk,nstate), cscr(*)
    REAL(real_8)                             :: vpp(ncpw%ngw), eigv(*), &
                                                rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: tau0(:,:,:), velp(:,:,:), &
                                                taui(:,:,:), fion(:,:,:)
    INTEGER                                  :: ifcalc, irec(:)
    LOGICAL                                  :: tfor, tinfo

    CHARACTER(*), PARAMETER                  :: procedureN = 'forces_diag'

    INTEGER                                  :: cdfti, i, infr, isub, nhpsi
    LOGICAL                                  :: calste_save, convtest, diisc, &
                                                pcgc, pcgminc, tsdec
    REAL(real_8)                             :: detot, dummy(2), etot0, tcpu, &
                                                thl(2), time1, time2, vener

!    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate,1), &

    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ropt_mod%convwf=.FALSE.
    ! Skip stress tensor calculation during wavefunction optimization
    calste_save=ropt_mod%calste
    ropt_mod%calste=.FALSE.
    etot0=0.0_real_8
    IF (cntl%tdiag.AND.cntl%tlanc) THEN
       CALL set_b2l()
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cntl%tdiag) THEN
       IF (.NOT.tinfo.AND.paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(T9,A)')&
               'INFR       DRHOMAX        ETOT       DETOT   HPSI    TCPU'
       ENDIF
    ELSE
       IF (.NOT.tinfo.AND.paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(T9,A)')&
               'INFR         GEMAX            ETOT         DETOT     TCPU'
       ELSEIF (paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,'(2A)') ' NFI      GEMAX       CNORM',&
               '           ETOT        DETOT      TCPU'
       ENDIF
       ropt_mod%sdiis=.TRUE.
       ropt_mod%spcg=.TRUE.
    ENDIF

    ! ==================================================================
    ! ==                        MAIN  LOOP  (cntl%cdft)                    ==
    ! ==================================================================
    IF (cntl%cdft) THEN
       IF (cdftlog%recw)CALL cdft_w(rhoe,tau0,dummy)
       IF (cdftlog%tpred)CALL v_pred()
       DO cdfti=1,cdftci%cdft_end
          DO infr=1,cnti%nomore_iter
             time1=m_walltime()
             ifcalc=ifcalc+1
             ! Diagonalization
             IF (cntl%tdiag) THEN
                CALL updrho(c0,c2,cr,sc0,cscr,vpp,tau0,fion,eigv,&
                     rhoe,psi,&
                     nstate,.FALSE.,tinfo,.FALSE.,infr,thl,nhpsi)
             ELSE
                ! Switching to cntl%pcg MINIMIZE
                IF (cntl%tpcgfic.AND.cdfti.EQ.1) THEN
                   diisc=cntl%diis
                   pcgc=cntl%pcg
                   tsdec=cntl%tsde
                   pcgminc=cntl%pcgmin
                   cntl%diis=.FALSE.
                   cntl%pcg=.TRUE.
                   cntl%tsde=.FALSE.
                   cntl%pcgmin=.TRUE.
                ENDIF
                IF (tkpts%tkpnt) THEN
                   CALL k_updwf(c0,c2,sc0,tau0,fion,cr,cscr,vpp,eigv,&
                        rhoe,psi,nstate,.FALSE.,infr)
                ELSE
                   CALL updwf(c0(:,:,1),c2,sc0,tau0,fion,cr,cscr,vpp,eigv,&
                        rhoe,psi,nstate,.FALSE.,.TRUE.)
                ENDIF
                ! Switching back to the selected optimizer
                IF (cntl%tpcgfic.AND.cdfti.EQ.1) THEN
                   cntl%diis=diisc
                   cntl%pcg=pcgc
                   cntl%tsde=tsdec
                   cntl%pcgmin=pcgminc
                ENDIF
             ENDIF
             ! Printout the evolution of the iterative optimization
             IF (paral%parent) THEN
                detot=ener_com%etot-etot0
                vener=cdftcom%cdft_v(1)*(cdftcom%cdft_nc+cdftcom%vgrad(1))
                IF (infr.EQ.1) detot=0.0_real_8
                time2=m_walltime()
                tcpu=(time2-time1)*0.001_real_8
                IF (tinfo) THEN
                   i=-1
                   IF (cntl%tdiag) i=infr
                   CALL wrprint_wfopt(eigv,crge%f,ener_com%amu,nstate,ener_com%etot+vener,etot0+vener,&
                        tcpu,gemax,cnorm,thl,.FALSE.,ifcalc,i)
                ELSEIF (cntl%tdiag) THEN
                   IF (paral%io_parent)&
                        WRITE(6,&
                        '(T9,I4,T14,1PE12.3,T27,0PF12.6,T39,1PE12.3,'//&
                        'T53,0PF5.2,T59,0PF7.2) ')&
                        infr,gemax,ener_com%etot+vener,detot,nhpsi/REAL(nstate,kind=real_8),tcpu
                ELSE
                   IF (paral%io_parent)&
                        WRITE(6,&
                        '(T9,I4,T15,1PE12.3,T31,0PF12.6,T45,1PE12.3,T59,0PF7.2) ')&
                        infr,gemax,ener_com%etot+vener,detot,tcpu
                ENDIF
                etot0=ener_com%etot
             ENDIF
             IF (MOD(ifcalc,store1%istore).EQ.0)&
                  CALL zhwwf(2,irec,c0,c2,nstate,eigv,tau0,velp,taui,iteropt%nfi)
             ! Check to break the loop.
             IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
             IF (soft_com%exsoft) GOTO 200   ! Soft Exit (no force calculation)
             IF (ropt_mod%convwf) GOTO 101   ! Convergence of wavefunctions
          ENDDO
101       CONTINUE
          ropt_mod%convwf=.FALSE.
          CALL vupdate(c0,psi,rhoe,convtest)
          ropt_mod%sdiis=.TRUE.
          ropt_mod%spcg=.TRUE.
          ! IF(PARENT)CALL WRENER
          IF (convtest) THEN
             GOTO 100
          ENDIF
          IF (cdftlog%reswf) THEN
             IF (paral%io_parent)&
                  WRITE(6,*) "Re-initialising Wavefunction"
             CALL rinitwf(c0,c2,sc0,crge%n,tau0,taui,rhoe,psi)
          ENDIF
       ENDDO
    ELSE
       ! ==================================================================
       ! ==                        MAIN  LOOP (NORMAL)                   ==
       ! ==================================================================
       DO infr=1,cnti%nomore_iter
          time1=m_walltime()
          ifcalc=ifcalc+1
          ! Diagonalization
          IF (cntl%tdiag) THEN
             CALL updrho(c0,c2,cr,sc0,cscr,vpp,tau0,fion,eigv,&
                  rhoe,psi,&
                  nstate,.FALSE.,tinfo,.FALSE.,infr,thl,nhpsi)
          ELSE
             CALL updwf(c0(:,:,1),c2,sc0,tau0,fion,cr,cscr,vpp,eigv,&
                  rhoe,psi,nstate,.FALSE.,.TRUE.)
          ENDIF
          ! Printout the evolution of the iterative optimization
          IF (paral%parent) THEN
             detot=ener_com%etot-etot0
             IF (infr.EQ.1) detot=0.0_real_8
             time2=m_walltime()
             tcpu=(time2-time1)*0.001_real_8
             IF (tinfo) THEN
                i=-1
                IF (cntl%tdiag) i=infr
                CALL wrprint_wfopt(eigv,crge%f,ener_com%amu,nstate,ener_com%etot,etot0,&
                     tcpu,gemax,cnorm,thl,.FALSE.,ifcalc,i)
             ELSEIF (cntl%tdiag) THEN
                IF (paral%io_parent)&
                     WRITE(6,&
                     '(T9,I4,T14,1PE12.3,T27,0PF12.6,T39,1PE12.3,'//&
                     'T53,0PF5.2,T59,0PF7.2) ')&
                     infr,gemax,ener_com%etot,detot,nhpsi/REAL(nstate,kind=real_8),tcpu
             ELSE
                IF (paral%io_parent)&
                     WRITE(6,&
                     '(T9,I4,T15,1PE12.3,T31,0PF12.6,T45,1PE12.3,T59,0PF7.2) ')&
                     infr,gemax,ener_com%etot,detot,tcpu
             ENDIF
             etot0=ener_com%etot
          ENDIF
          IF (MOD(ifcalc,store1%isctore).EQ.0)&
               CALL zhwwf(2,irec,c0,c2,nstate,eigv,tau0,velp,taui,iteropt%nfi)
          ! Check to break the loop.
          IF (.NOT.soft_com%exsoft) CALL testex(soft_com%exsoft)
          IF (soft_com%exsoft) GOTO 200   ! Soft Exit (no force calculation)
          IF (ropt_mod%convwf) GOTO 100
          ! ==--------------------------------------------------------------==
          ! special case for ASPC extrapolator with CSTEPS flag. once we have 
          ! a full wavefunction history, we don't converge out but use the
          ! specified number of corrector steps. 
          ! for the first wavefunction during init we have TINFO=.TRUE. .
          ! that one is not taken from the history and may no be converged,
          ! we have to converge that one in full.
          IF (cntl%taspc.AND.(cnti%naspc.GT.0)) THEN
             IF (.NOT.tinfo) THEN
                IF ((infr.GE.cnti%naspc).AND.(numcold.GE.cnti%mextra)) GOTO 100
             ENDIF
          ENDIF
          ! ==--------------------------------------------------------------==
       ENDDO
    ENDIF
    ! ==================================================================
    ! ==                      END OF MAIN  LOOP                       ==
    ! ==================================================================
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' WARNING! DENSITY NOT CONVERGED !!!'
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' WARNING! DENSITY NOT CONVERGED !!!'
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' WARNING! DENSITY NOT CONVERGED !!!'
100 CONTINUE
    ! Compute stress tensor if requested
    ropt_mod%calste=calste_save
    ! ==            DIAGONALIZATION FOLLOWING OPTIMIZATION            ==
    IF (cntl%tdiagopt) THEN
       cntl%tdiag=.TRUE.
       tfixo=.TRUE.
       CALL dens_for_diag(c0,rhoe,psi(:,1),nstate)
       CALL updrho(c0,c2,cr,sc0,cscr,vpp,tau0,fion,eigv,&
            rhoe,psi,&
            nstate,tfor,tinfo,ropt_mod%calste,infr,thl,nhpsi)
       ropt_mod%sdiis=.TRUE.
       cntl%tdiag=.FALSE.
    ENDIF
    ! ==    IF WE NEED FORCES AND DON'T HAVE THEM YET, GET THEM NOW   ==
    IF (tfor) THEN
       IF (cntl%tdiag) THEN
          ! Diagonalization
          CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
          CALL updrho(c0,c2,cr,sc0,cscr,vpp,tau0,fion,eigv,&
               rhoe,psi,&
               nstate,.TRUE.,tinfo,ropt_mod%calste,infr,thl,nhpsi)
          ! Forces from geometrical restraints
          IF (paral%parent) THEN
             DO i=1,cotr007%mrestr
                resval(i)=resval(i)+gsrate(i)*dt_ions
                IF (resval_dest(i).NE.-999._real_8) THEN
                   IF (gsrate(i).GT.0._real_8.AND.resval(i).GT.resval_dest(i))&
                        resval(i)=resval_dest(i) ! increase
                   IF (gsrate(i).LT.0._real_8.AND.resval(i).LT.resval_dest(i))&
                        resval(i)=resval_dest(i) ! decrease
                ENDIF
             ENDDO
             CALL restfc(tau0,fion)
          ENDIF
       ELSE
          CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
          IF (tkpts%tkpnt) THEN
             infr=infr+1
             ropt_mod%convwf=.FALSE.
             CALL k_updwf(c0,c2,sc0,tau0,fion,cr,cscr,vpp,eigv,&
                  rhoe,psi,nstate,.TRUE.,infr)
             ropt_mod%convwf=.TRUE.
             ! Forces from geometrical restraints
             IF (paral%parent) THEN
                DO i=1,cotr007%mrestr
                   resval(i)=resval(i)+gsrate(i)*dt_ions
                   IF (resval_dest(i).NE.-999._real_8) THEN
                      IF (gsrate(i).GT.0._real_8.AND.resval(i).GT.resval_dest(i))&
                           resval(i)=resval_dest(i) ! increase
                      IF (gsrate(i).LT.0._real_8.AND.resval(i).LT.resval_dest(i))&
                           resval(i)=resval_dest(i) ! decrease
                   ENDIF
                ENDDO
                CALL restfc(tau0,fion)
             ENDIF
          ELSE
             CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,tau0,fion,eigv,&
                  nstate,1,.TRUE.,tfor)
          ENDIF
          ! AK: NOTE: Restraints are already treated in FORCEDR()
       ENDIF
       !CALL mp_bcast(fion,3*maxsys%nax*maxsys%nsx,parai%source,parai%allgrp)
       CALL mp_bcast(fion,3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
    ENDIF
200 CONTINUE
    ! Reset calste correctly even for soft exit
    ropt_mod%calste=calste_save
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE forces_diag
  ! ==================================================================
  SUBROUTINE updiag(c0,c2,cr,sc0,cscr,vpp,tau0,fion,eigv,&
       rhoe,psi,nstate,ifcalc,tfor,tinfo,&
       time1,etot0)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: cr(*), cscr(*)
    REAL(real_8)                             :: vpp(ncpw%ngw), tau0(:,:,:), &
                                                fion(:,:,:), eigv(*), &
                                                rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: sc0(nkpt%ngwk,nstate), &
                                                c2(nkpt%ngwk,nstate), &
                                                c0(nkpt%ngwk,nstate)
    INTEGER                                  :: ifcalc
    LOGICAL                                  :: tfor, tinfo
    REAL(real_8)                             :: time1, etot0

    INTEGER                                  :: nhpsi
    LOGICAL                                  :: tdodiag
    REAL(real_8)                             :: detot, tcpu, thl(2), time2

!(fpar%nnr1,clsd%nlsd)
! ==--------------------------------------------------------------==
! Determine if we do diagonalization this time

    tdodiag = .FALSE.
    IF (cnti%nperid.GT.0 .AND. cnti%nrperi.EQ.0) THEN
       tdodiag=(MOD(ifcalc,cnti%nperid).EQ.0)
    ELSE IF (cnti%nperid.EQ.0 .AND. cnti%nrperi.GT.0) THEN
       tdodiag=(iteropt%ndisrs.GE.cnti%nrperi)
       IF (tdodiag) iteropt%ndisrs=0
    ELSE IF (cnti%nperid.GT.0 .AND. cnti%nrperi.GT.0) THEN
       tdodiag=(iteropt%ndisrs.GE.cnti%nrperi .AND.&
            (iteropt%ndistp.GE.cnti%nperid.OR.iteropt%ndisrs.GT.cnti%nrperi))
       IF (tdodiag) THEN
          iteropt%ndisrs=0
          iteropt%ndistp=0
       ENDIF
    ELSEIF (ropt_mod%convwf .AND. cnti%nperid.EQ.0 .AND. cnti%nrperi.EQ.0) THEN
       tdodiag = .TRUE.
    ENDIF
    ! Yes, we do
    IF (tdodiag) THEN
       ! ampre=2.0e-4_real_8
       ! if (parent) write (6,*) 'randomizing WF, amplitude: ', ampre
       ! call rrane(c0,c2,nstate)
       ! sdiis=.true.
       ! endif
       ! if(.false.) then
       time2=m_walltime()
       tcpu=(time2-time1)*0.001_real_8
       detot=ener_com%etot-etot0
       IF (.NOT.ropt_mod%convwf) THEN
          IF (tinfo) THEN
             CALL wrprint_wfopt(eigv,crge%f,ener_com%amu,nstate,ener_com%etot,etot0,&
                  tcpu,gemax,cnorm,thl,.FALSE.,ifcalc,ifcalc)
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,&
                  '(T9,I4,T15,1PE12.3,T31,0PF12.6,T45,1PE12.3,T59,0PF7.2) ')&
                  ifcalc,gemax,ener_com%etot,detot,tcpu
          ENDIF
       ENDIF
       time1=m_walltime()
       etot0=ener_com%etot
       cntl%tdiag=.TRUE.
       tfixo=.TRUE.
       ! Get density as required by UPDRHO
       CALL dens_for_diag(c0,rhoe,psi(:,1),nstate)
       ifcalc=ifcalc+1
       ! Full diagonalization
       CALL updrho(c0,c2,cr,sc0,cscr,vpp,tau0,fion,eigv,&
            rhoe,psi,&
            nstate,tfor,tinfo,ropt_mod%calste,ifcalc,thl,nhpsi)
       ! Calculate energies the usual way and init cntl%diis if required
       ropt_mod%sdiis=.TRUE.
       iteropt%ndisrs=0
       cntl%tdiag=.FALSE.
       ropt_mod%convwf=.FALSE.
       IF (tinfo.AND.paral%io_parent)&
            WRITE(6,'(2A)') ' NFI      GEMAX',&
            '       CNORM           ETOT        DETOT      TCPU'
       IF (.NOT.tinfo.AND.paral%io_parent)&
            WRITE(6,'(T9,A)') 'INFR',&
            '         GEMAX            ETOT         DETOT     TCPU'
       IF (tfor) CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
       ! Re-init energy arrays, cntl%diis, etc.
       CALL updwf(c0,c2,sc0,tau0,fion,cscr,cr,vpp,eigv,&
            rhoe,psi,nstate,tfor,.TRUE.)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE updiag
  ! ==================================================================
  SUBROUTINE dens_for_diag(c0,rhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8) :: c0(nkpt%ngwk,crge%n,nkpt%nkpnt)
    REAL(real_8)                             :: rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate

    INTEGER                                  :: ikind, ikpt, kbeg, kend, &
                                                kinc, nkpoint

    IF (pslo_com%tivan) THEN
       CALL inq_swap(kbeg,kend,kinc)
       DO ikpt=kbeg,kend,kinc
          nkpoint=nkpbl(ikpt)
          IF (tkpts%tkblock) CALL rkpt_swap(c0,nstate,ikpt,&
               'HGKP HGKM MASKGW EIGKR TWNL C0')
          DO ikind=1,nkpoint
             CALL rnlsm(c0(:,:,ikind),nstate,&
                  ikpt,ikind,.FALSE.)
          ENDDO
       ENDDO
    ENDIF
    IF (tkpts%tkpnt) THEN
       CALL rhoofr_c(c0,rhoe,psi,nstate)
    ELSE
       CALL rhoofr(c0(:,:,1),rhoe,psi,nstate)
    ENDIF
    CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe,1,rin0,1)
    ! DENSITY IS NOW THERE, DO MIXING STUFF
    IF (andr2%trand) CALL rrandd(rin0,andr2%amprd)
    CALL zeroing(rout0)!,nnr1*clsd%nlsd)
    CALL zeroing(rmix)!,nnr1*clsd%nlsd)
    ! NL-Projector Overlap Matrix
    IF (fint1%ttrot) CALL calc_alm
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dens_for_diag
  ! ==================================================================
  SUBROUTINE give_scr_forces_diag(lforces_diag,tag,nstate,tfor)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lforces_diag
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: tfor

    INTEGER                                  :: ltscr

! ==--------------------------------------------------------------==

    IF (cntl%tdiag) THEN
       CALL give_scr_updrho(lforces_diag,tag,nstate,&
            tfor,tfor.AND.cntl%tpres)
    ELSE
       IF (tkpts%tkpnt) THEN
          CALL give_scr_kupdwf(lforces_diag,tag,nstate,tfor)
       ELSE
          CALL give_scr_updwf(lforces_diag,tag,nstate,tfor)
       ENDIF
    ENDIF
    IF (tfor) THEN
       ltscr=3*maxsys%nax*maxsys%nsx
    ELSE
       ltscr=0
    ENDIF
    lforces_diag=MAX(lforces_diag,ltscr)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_forces_diag
  ! ==================================================================


END MODULE forces_diag_utils

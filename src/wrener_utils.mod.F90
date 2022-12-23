MODULE wrener_utils
  USE adat,                            ONLY: elem
  USE cdftmod,                         ONLY: cdftcom,&
                                             cdftlog
  USE clas,                            ONLY: tclas
  USE cnst,                            ONLY: au_kb,&
                                             pi,&
                                             ry
  USE cnstpr_utils,                    ONLY: cnstpr
  USE conv,                            ONLY: nac,&
                                             nbc
  USE coor,                            ONLY: velp
  USE ener,                            ONLY: chrg,&
                                             ener_c,&
                                             ener_com,&
                                             ener_d,&
                                             tenergy_ok
  USE fint,                            ONLY: fint1
  USE hubbardu,                        ONLY: hubbu  
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: rk,&
                                             wk
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_flush
  USE metr,                            ONLY: metr_com
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_qm,&
                                             mm_revert
  USE mm_input,                        ONLY: eqm_r,&
                                             lqmmm
  USE parac,                           ONLY: paral
  USE pslo,                            ONLY: pslo_com
  USE response_pmod,                   ONLY: ener1,&
                                             response1
  USE rmas,                            ONLY: rmass
  USE ropt,                            ONLY: ropt_mod
  USE spin,                            ONLY: amudo,&
                                             amuup,&
                                             lspin2,&
                                             spin_mod
  USE store_types,                     ONLY: &
       cprint, iprint_coor, iprint_eband, iprint_ebogo, iprint_ecas, &
       iprint_eeig, iprint_eext, iprint_egc, iprint_ehee, iprint_ehep, &
       iprint_ehii, iprint_ehsic, iprint_eht, iprint_eigen, iprint_ekin, &
       iprint_elec1, iprint_elec2, iprint_enl, iprint_entropy, iprint_epen, &
       iprint_epseu, iprint_eself, iprint_esr, iprint_etddft, iprint_etot1, &
       iprint_etot2, iprint_evdw, iprint_exc, iprint_force, iprint_info, &
       iprint_vxc, iprint_ehub
  USE strs,                            ONLY: paiu
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             nkpt,&
                                             parm
  USE temps,                           ONLY: tempcm,&
                                             tempqm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vdwcmod,                         ONLY: vdwr
  USE wrgeo_utils,                     ONLY: wrgeo,&
                                             wrgeof,&
                                             wrgeox

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wrprint_wfopt
  PUBLIC :: wrprint_geopt
  PUBLIC :: wrprint_md
  PUBLIC :: wrprint_pr
  PUBLIC :: wrener_init
  PUBLIC :: wrener
  PUBLIC :: wreigen
  PUBLIC :: wrstress
  PUBLIC :: wrcell

CONTAINS

  ! ==================================================================
  SUBROUTINE wrprint_wfopt(eigv,f,amu,nstate,&
       etot,etot0,tcpu,gemax,cnorm,thl,&
       convergence,nfi,infi)
    ! ==--------------------------------------------------------------==
    ! == PRINT FOR WAVEFUNCTIONS OPTIMIZATION                         ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: amu
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate,nkpt%nkpts), &
                                                eigv(nstate,nkpt%nkpts), &
                                                etot, etot0, tcpu, gemax, &
                                                cnorm, thl(2)
    LOGICAL                                  :: convergence
    INTEGER                                  :: nfi, infi

    LOGICAL                                  :: engpri, tend
    REAL(real_8)                             :: detot, tthl

    IF (.NOT.paral%io_parent) RETURN
    detot=etot-etot0
    IF (infi.EQ.1) detot=0.0_real_8
    engpri=MOD(infi-1,cprint%iprint_step).EQ.0.AND.infi.GE.0
    tend=infi.EQ.1.OR.infi.EQ.cnti%nomore.OR.convergence
    IF (cntl%tdiag) THEN
       IF (ok_print(cprint%iprint(iprint_eigen),engpri,tend))&
            CALL wreigen(eigv,f,amu,nstate)
       IF (ok_print(cprint%iprint(iprint_info),engpri,tend)) CALL wrener
       tthl=thl(1)
       IF (cntl%tlsd) tthl=0.5_real_8*(thl(1)+thl(2))
       IF (paral%io_parent)&
            WRITE(6,'(" ==",60("-"),"==")')
       IF (paral%io_parent)&
            WRITE(6,'(A12,I7,T25,A6,F14.6,T48,A5,F10.2,A3) ')&
            ' ==     NFI=',nfi,  ' ETOT=',etot, 'TCPU=',tcpu,' =='
       IF (paral%io_parent)&
            WRITE(6,'(A12,1PE10.3,T25,A6,1PE14.3,T48,A5,1PE10.3,A3) ')&
            ' == DRHOMAX=',gemax,'DETOT=',detot,' THL=',tthl,' =='
       IF (paral%io_parent)&
            WRITE(6,'(" ==",60("-"),"==")')
    ELSE
       IF ((infi.EQ.0).AND.paral%io_parent)&
            WRITE(6,'(2A)') ' NFI      GEMAX       CNORM',&
            '           ETOT        DETOT      TCPU'
       IF (ABS(cnorm).LT.1.e-20_real_8) cnorm=0._real_8
       IF ( (nfi.EQ.1.AND.infi.GE.0) ) THEN
          IF (cntl%tfrho_upw)  THEN
             IF (paral%io_parent)&
                  WRITE(6,'(2A)') ' NFI      GEMAX       CNORM',&
                  '        H-TRACE       DTRACE      TCPU'
          ELSE
             CALL wrener
             IF (paral%io_parent)&
                  WRITE(6,'(2A)') ' NFI      GEMAX       CNORM',&
                  '           ETOT        DETOT      TCPU'
          ENDIF
       ENDIF
       IF (paral%io_parent) THEN
          WRITE(6,'(I4,1PE11.3,1PE12.3,0PF15.6,1PE13.3,0PF10.2) ')&
               nfi,gemax,cnorm,etot,detot,tcpU
       ENDIF
       IF (.NOT.(nfi.EQ.1.AND.infi.GE.0) .AND.&
            ok_print(cprint%iprint(iprint_info),engpri,.FALSE.)) THEN
          IF (cntl%tfrho_upw)  THEN
             IF (paral%io_parent)&
                  WRITE(6,'(2A)') ' NFI      GEMAX       CNORM',&
                  '        H-TRACE       DTRACE      TCPU'
          ELSE
             CALL wrener
             IF (paral%io_parent)&
                  WRITE(6,'(2A)') ' NFI      GEMAX       CNORM',&
                  '           ETOT        DETOT      TCPU'
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    etot0=etot
    CALL m_flush(6)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wrprint_wfopt
  ! ==================================================================
  SUBROUTINE wrprint_geopt(eigv,f,amu,nstate,tau0,fion,etot,etot0,&
       tcpu,gemax,gnmax,dxmax,gnorm,cnmax,dstre,&
       dsmax,convergence,ifcalc,infi0)
    ! ==--------------------------------------------------------------==
    ! == PRINT FOR GEOMETRY OPTIMIZATION                              ==
    ! ==--------------------------------------------------------------==
    ! == EIGV(NSTATE,NKPTS) Eigenvalues                               ==
    ! == F(NSTATE,NKPTS)    Occupation numbers                        ==
    ! == AMU                Chemical potential                        ==
    ! == NSTATE             Number of states                          ==
    ! == TAU0(3,maxsys%nax,maxsys%nsx)    Atomic positions                          ==
    ! == FION(3,maxsys%nax,maxsys%nsx)    Atomic forces                             ==
    ! == ETOT               Total energy                              ==
    ! == ETOT0              Reference total energy                    ==
    ! == TCPU               CPU time                                  ==
    ! == GEMAX              Max. change in wavefunction components    ==
    ! ==                    (or for cntl%tdiag, max. change for density)   ==
    ! == GNMAX              Max. force components                     ==
    ! == DXMAX              Max. atomic position change               ==
    ! == GNORM              Norm of Forces                            ==
    ! == CNMAX              Max. Constraint components                ==
    ! == DSTRE              Sum Abs(HTFOR)                            ==
    ! == DSMAX              Max(HTFOR)                                ==
    ! == CONVERGENCE        .TRUE. if convergence                     ==
    ! == IFCALC             Total number of iterations                ==
    ! == INFI0              Number of steps for geometry optimisation ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: amu
    INTEGER                                  :: nstate
    REAL(real_8) :: f(nstate,nkpt%nkpts), eigv(nstate,nkpt%nkpts), &
      tau0(:,:,:), fion(:,:,:), etot, etot0, tcpu, gemax, gnmax, dxmax, &
      gnorm, cnmax, dstre, dsmax
    LOGICAL                                  :: convergence
    INTEGER                                  :: ifcalc, infi0

    CHARACTER(*), PARAMETER                  :: procedureN = 'wrprint_geopt'

    INTEGER                                  :: isub
    LOGICAL                                  :: engpri0, tend
    REAL(real_8)                             :: detot

    IF (.NOT.paral%parent) RETURN
    IF (.NOT.cntl%tdiag.AND.infi0.EQ.0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(2A)') ' NFI      GEMAX       CNORM',&
            '           ETOT        DETOT      TCPU'
       RETURN
    ENDIF

    CALL tiset(procedureN,isub)

    detot=etot-etot0
    IF (infi0.EQ.1) detot=0.0_real_8
    engpri0=MOD(infi0,cprint%iprint_step).EQ.0
    tend=infi0.EQ.cnti%nomore.OR.convergence
    IF (cntl%tdiag) THEN
       IF (ok_print(cprint%iprint(iprint_eigen),engpri0,tend))&
            CALL wreigen(eigv,f,amu,nstate)
    ENDIF
    IF (ok_print(cprint%iprint(iprint_info),engpri0,tend)) CALL wrener
    IF (cprint%iprint(iprint_force).GT.0) THEN
       CALL wrgeof(tau0,fion)
       CALL wrgeox(tau0)
    ELSEIF (cprint%iprint(iprint_coor).NE.-1) THEN
       CALL cnstpr
       CALL wrgeo(tau0)
       CALL wrgeox(tau0)
    ENDIF
    IF (ropt_mod%calste) CALL wrstress
    IF (cntl%tprcp)  CALL wrcell
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("*"))')
    IF (paral%io_parent)&
         WRITE(6,'(A,I6,A,I7,A)') ' *** TOTAL STEP NR.',ifcalc,&
         '           GEOMETRY STEP NR.',infi0,'  ***'
    IF (dxmax.EQ.0._real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,1PE14.6,T40,A7,0PF14.6,A5)')&
            ' *** GNMAX=',gnmax,'  ETOT=',etot,'  ***'
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(A,1PE14.6,A,1PE8.2,T36,A11,0PF14.6,A5)')&
            ' *** GNMAX=',gnmax,' [',dxmax,']     ETOT=',etot,'  ***'
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,'(A,1PE14.6,T40,A7,1PE14.3,A5)')&
         ' *** GNORM=',gnorm,' DETOT=',detot,'  ***'
    IF (paral%io_parent)&
         WRITE(6,'(A,1PE14.6,T40,A7,0PF14.2,A5)')&
         ' *** CNSTR=',cnmax,'  TCPU=',tcpu,'  ***'
    IF ((cntl%tprcp).AND.paral%io_parent)&
         WRITE(6,'(A,1PE14.6,T40,A7,1PE14.3,A5)')&
         ' *** DSTRE=',dstre,'[DSMAX=',dsmax,'] ***'
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("*"))')
    ! ==--------------------------------------------------------------==
    etot0=etot
    CALL m_flush(6)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE wrprint_geopt
  ! ==================================================================
  SUBROUTINE wrprint_md(eigv,f,amu,nstate,tau0,fion,&
       ekinc,tempp,etot,econs,eham,disa,&
       tcpu,convergence,nfi,infi)
    ! ==--------------------------------------------------------------==
    ! == PRINT FOR MOLECULAR DYNAMICS                                 ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: amu
    INTEGER                                  :: nstate
    REAL(real_8) :: f(nstate,nkpt%nkpts), eigv(nstate,nkpt%nkpts), &
      tau0(:,:,:), fion(:,:,:), ekinc, tempp, etot, econs, eham, disa, tcpu
    LOGICAL                                  :: convergence
    INTEGER                                  :: nfi, infi

    LOGICAL                                  :: engpri, tend

    IF (.NOT.paral%io_parent) RETURN
    IF (cntl%tdiag) THEN
       tend=(infi.EQ.0).OR.(infi.EQ.cnti%nomore).OR.convergence
       engpri=MOD(nfi-1,cprint%iprint_step).EQ.0
       IF (infi.NE.0) THEN
          WRITE(6,'(A,A)')&
               '   NFI   TEMPP        EKS       ECLASSIC',&
               '         DIS           TCPU'
          ! mb-> DISA can become VERY big in low density/high temperature
          ! mb-> systems. Format set to E14.6
          WRITE(6,'(T2,I10,F8.1,2F14.6,E14.6,F9.2)')&
               nfi,tempp,etot,econs,disa,tcpU
          IF (.NOT.cntl%cdft)THEN
             WRITE(3,'(T2,I10,F8.1,2F14.6,E14.6,F9.2)')&
                  nfi,tempp,etot,econs,disa,tcpU
          ELSE
             WRITE(3,'(T2,I10,F8.1,2F14.6,E14.6,F9.2,F12.6,F12.6)')&
                  nfi,tempp,etot,econs,disa,tcpu,cdftcom%cdft_v,cdftcom%vgrad

          ENDIF
       ENDIF
       IF (ok_print(cprint%iprint(iprint_info),engpri,tend)) THEN
          WRITE(6,'(/,1X,64("*"))')
          WRITE(6,'(" *",T20,A,I4,T65,A)')&
               ' RESULTS OF STEP',infi,'*'
          WRITE(6,'(1X,64("*"))')
       ENDIF
       IF (ok_print(cprint%iprint(iprint_eigen),engpri,tend))&
            CALL wreigen(eigv,f,amu,nstate)
       IF (ok_print(cprint%iprint(iprint_info),engpri,tend)) CALL wrener
       IF (ok_print(cprint%iprint(iprint_force),engpri,tend)) THEN
          CALL wrgeof(tau0,fion)
       ELSEIF (ok_print(cprint%iprint(iprint_coor),engpri,tend)) THEN
          CALL cnstpr
          CALL wrgeo(tau0)
       ENDIF
       IF (ok_print(cprint%iprint(iprint_info),engpri,tend)) THEN
          WRITE(6,'(1X,64("*"),/)')
       ENDIF
    ELSE
       IF (infi.EQ.cnti%nomore) GOTO 100
       tend=infi.EQ.1.OR.infi.EQ.cnti%nomore.OR.convergence
       engpri=cprint%tprint.AND.MOD(nfi-1,cprint%iprint_step).EQ.0
       IF (ok_print(cprint%iprint(iprint_force),engpri,tend)) THEN
          CALL wrgeof(tau0,fion)
          IF (cntl%bsymm) WRITE(6,*)
       ELSEIF (ok_print(cprint%iprint(iprint_coor),engpri,tend)) THEN
          CALL cnstpr
          CALL wrgeo(tau0)
       ENDIF
       IF (ok_print(cprint%iprint(iprint_info),engpri,tend)) THEN
          CALL cnstpr
          IF (.NOT.cntl%bsymm) CALL wrener
          IF (.NOT.(infi.EQ.cnti%nomore.OR.convergence)) THEN
             IF (cntl%tmdbo) THEN
                WRITE(6,'(T6,A,A,11X,A,6X,A,A,A)')&
                     '  NFI','   TEMPP','EKS','ECLASSIC',&
                     '         DIS','    TCPU'
             ELSE IF (lqmmm%qmmm) THEN
                WRITE(6,'(T6,A,A,A,11X,A,6X,A,10X,A,A,A,A)')&
                     '  NFI','    EKINC','   TEMPP','EKS','ECLASSIC','EHAM'&
                     ,'           EQM','         DIS','    TCPU'
             ELSE
                WRITE(6,'(T6,A,A,A,11X,A,6X,A,10X,A,A,A)')&
                     '  NFI','    EKINC','   TEMPP','EKS','ECLASSIC','EHAM',&
                     '         DIS','    TCPU'
             ENDIF
          ENDIF
       ELSE IF (cntl%tmdfile) THEN
          IF (cntl%tmdbo) THEN
             WRITE(6,'(T6,A,A,11X,A,6X,A,A,A)')&
                  '  NFI','   TEMPP','EKS','ECLASSIC',&
                  '         DIS','    TCPU'
          ELSE
             WRITE(6,'(T6,A,A,A,11X,A,6X,A,10X,A,A,A)')&
                  '  NFI','    EKINC','   TEMPP','EKS','ECLASSIC','EHAM',&
                  '         DIS','    TCPU'
          ENDIF
       ENDIF
       ! mb-> DISA can become VERY big in low density/high temperature
       ! mb-> systems. Format set to  E10.3 and E14.6
100    CONTINUE
       IF (cntl%tmdbo) THEN
          WRITE(6,'(I10,F8.1,F14.5,F14.5,E12.3,F8.2)')&
               nfi,tempp,etot,econs,disa,tcpU
       ELSE IF (lqmmm%qmmm) THEN
          WRITE(6,'(I10,F9.5,F8.1,4F14.5,E12.3,F8.2)')&
               nfi,ekinc,tempp,etot,econs,eham,eqm_r%eqm,disa,tcpU
       ELSE
          WRITE(6,'(I10,F9.5,F8.1,F14.5,F14.5,F14.5,E12.3,F8.2)')&
               nfi,ekinc,tempp,etot,econs,eham,disa,tcpU
       ENDIF
       IF(lqmmm%qmmm) THEN 
          WRITE(3,'(I10,F12.8,F9.3,F18.10,F18.10,F18.10,F18.10,E15.6,F9.2)')&
               nfi,ekinc,tempp,etot,econs,eham,eqm_r%eqm,disa,tcpU
       ELSE 
          IF (.NOT.cntl%cdft)THEN
             WRITE(3,'(I10,F12.8,F9.3,F18.10,F18.10,F18.10,E15.6,F9.2)')&
                  nfi,ekinc,tempp,etot,econs,eham,disa,tcpU
          ELSE
             WRITE(3,'(I10,F12.8,F9.3,F18.10,F18.10,F18.10,E15.6,F9.2,F12.6,F12.6)')&
                  nfi,ekinc,tempp,etot,econs,eham,disa,tcpu,cdftcom%cdft_v(1),&
                  cdftcom%vgrad(1)
          ENDIF
       ENDIF
       IF (tclas)WRITE(15,'(I10,3(F9.3,3X))') nfi,tempcm,tempqm,tempp
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL m_flush(3)
    CALL m_flush(6)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wrprint_md
  ! ==================================================================
  SUBROUTINE wrprint_pr(ekinc,ekinh,tempp,etot,econs,&
       eham,disa,tau0,fion,tcpu,nfi,infi)
    ! ==--------------------------------------------------------------==
    ! == PRINT FOR MOLECULAR DYNAMICS                                 ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: ekinc, ekinh, tempp, etot, &
                                                econs, eham, disa, &
                                                tau0(:,:,:), fion(:,:,:), tcpu
    INTEGER                                  :: nfi, infi

    LOGICAL                                  :: engpri, tend

    IF (.NOT.paral%parent) RETURN
    tend=infi.EQ.1.OR.infi.EQ.cnti%nomore
    engpri=cprint%tprint.AND.MOD(nfi-1,cprint%iprint_step).EQ.0
    IF (ok_print(cprint%iprint(iprint_force),engpri,tend)) THEN
       CALL wrgeof(tau0,fion)
    ELSEIF (ok_print(cprint%iprint(iprint_coor),engpri,tend)) THEN
       CALL cnstpr
       CALL wrgeo(tau0)
    ENDIF
    IF (ok_print(cprint%iprint(iprint_info),engpri,tend)) THEN
       CALL wrener
       IF ((infi.NE.cnti%nomore).AND.paral%io_parent)&
            WRITE(6,'(A,A)')&
            '       NFI    EKINC    EKINH   TEMPP           EKS      ',&
            'ECLASSIC          EHAM         DIS     TCPU'
    ENDIF
    ! mb-> DISA can become VERY big in low density/high temperature
    ! mb-> systems. Format set to E10.3 and E14.6
    IF (paral%io_parent)&
         WRITE(6,'(I10,F9.5,F9.5,F8.1,3F14.5,E12.3,F9.2)')&
         nfi,ekinc,ekinh,tempp,etot,econs,eham,disa,tcpu
    IF (paral%io_parent)&
         WRITE(3,'(I10,F12.8,F12.8,F10.3,3F18.10,E15.6,F9.2)')&
         nfi,ekinc,ekinh,tempp,etot,econs,eham,disa,tcpu
    ! ==--------------------------------------------------------------==
    CALL m_flush(6)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wrprint_pr
  ! ==================================================================
  SUBROUTINE wrener_init
    ! ==--------------------------------------------------------------==
    ! == INITIALISATION OF IPRINT ARRAY FOR OUTPUT                    ==
    ! ==--------------------------------------------------------------==
    ! Variables
    ! ==--------------------------------------------------------------==
    ! General setup
    ! If unspecified -> display
    IF (cprint%iprint(iprint_info ).NE.-1) cprint%iprint(iprint_info ) =  1
    IF (cprint%iprint(iprint_eht  ).NE.-1) cprint%iprint(iprint_eht  ) =  1
    IF (cprint%iprint(iprint_eself).NE.-1) cprint%iprint(iprint_eself) =  1
    IF (cprint%iprint(iprint_esr  ).NE.-1) cprint%iprint(iprint_esr  ) =  1
    ! If unspecified -> do not display
    IF (cprint%iprint(iprint_ehep ).EQ. 0) cprint%iprint(iprint_ehep ) = -1
    IF (cprint%iprint(iprint_ehee ).EQ. 0) cprint%iprint(iprint_ehee ) = -1
    IF (cprint%iprint(iprint_ehii ).EQ. 0) cprint%iprint(iprint_ehii ) = -1
    IF (cprint%iprint(iprint_elec1).EQ. 0) cprint%iprint(iprint_elec1) = -1
    IF (cprint%iprint(iprint_elec2).EQ. 0) cprint%iprint(iprint_elec2) = -1
    IF (cprint%iprint(iprint_etot1).EQ. 0) cprint%iprint(iprint_etot1) = -1
    IF (cprint%iprint(iprint_ehsic).EQ. 0) cprint%iprint(iprint_ehsic) = -1 ! cmb_ssic
    IF (cprint%iprint(iprint_ehub ).EQ. 0) cprint%iprint(iprint_ehub ) = -1 
    ! Special setup
    IF (cntl%tdiag) THEN
       ! For diagonalisation
       ! If unspecified -> display
       IF (cntl%tfint) THEN
          IF (cprint%iprint(iprint_eeig ).EQ.0) cprint%iprint(iprint_eeig ) =  1
          IF (cprint%iprint(iprint_eband).EQ.0) cprint%iprint(iprint_eband) = -1
       ELSE
          IF (cprint%iprint(iprint_eeig ).EQ.0) cprint%iprint(iprint_eeig ) = -1
          IF (cprint%iprint(iprint_eband).EQ.0) cprint%iprint(iprint_eband) =  1
       ENDIF
       IF (cprint%iprint(iprint_exc  ).EQ.0) cprint%iprint(iprint_exc  ) =  1
       IF (cprint%iprint(iprint_vxc  ).EQ.0) cprint%iprint(iprint_vxc  ) =  1
       ! If unspecified -> do not display
       IF (cprint%iprint(iprint_ekin ).EQ.0) cprint%iprint(iprint_ekin ) = -1
       IF (cprint%iprint(iprint_epseu).EQ.0) cprint%iprint(iprint_epseu) = -1
       IF (cprint%iprint(iprint_enl  ).EQ.0) cprint%iprint(iprint_enl  ) = -1
       IF (cprint%iprint(iprint_entropy).EQ.0) cprint%iprint(iprint_entropy) = -1
       IF (cprint%iprint(iprint_etot2).EQ. 0) cprint%iprint(iprint_etot2) = -1
    ELSE
       ! For Car-Parrinello scheme
       ! If unspecified -> display
       IF (cprint%iprint(iprint_ekin ).EQ.0) cprint%iprint(iprint_ekin ) =  1
       IF (cprint%iprint(iprint_epseu).EQ.0) cprint%iprint(iprint_epseu) =  1
       IF (cprint%iprint(iprint_enl  ).EQ.0) cprint%iprint(iprint_enl  ) =  1
       IF (cprint%iprint(iprint_exc  ).EQ.0) cprint%iprint(iprint_exc  ) =  1
       ! If unspecified -> do not display
       IF (cprint%iprint(iprint_vxc  ).EQ.0) cprint%iprint(iprint_vxc  ) = -1
       ! Never displayed
       cprint%iprint(iprint_eband) = -1
       cprint%iprint(iprint_eeig)  = -1
       cprint%iprint(iprint_entropy) = -1
       cprint%iprint(iprint_ebogo) = -1
       cprint%iprint(iprint_etot2) = -1
    ENDIF
    ! Gradient correction
    IF (cntl%tgc) THEN
       IF (cprint%iprint(iprint_egc  ).EQ.0) cprint%iprint(iprint_egc  ) =  1
    ELSE
       IF (cprint%iprint(iprint_egc  ).EQ.0) cprint%iprint(iprint_egc  ) = -1
    ENDIF
    ! Bogoliubov correction
    IF (fint1%tbogo) THEN
       IF (cprint%iprint(iprint_ebogo).EQ.0) cprint%iprint(iprint_ebogo) =  1
    ELSE
       IF (cprint%iprint(iprint_ebogo).EQ.0) cprint%iprint(iprint_ebogo) = -1
    ENDIF
    ! CAS 22 Correction
    IF (lspin2%tcas22) THEN
       IF (cprint%iprint(iprint_ecas).EQ.0) cprint%iprint(iprint_ecas) =  1
    ELSE
       IF (cprint%iprint(iprint_ecas).EQ.0) cprint%iprint(iprint_ecas) = -1
    ENDIF
    ! PENALTY FUNCTION 
    IF (lspin2%tpenal) THEN
       IF (cprint%iprint(iprint_epen).EQ.0) cprint%iprint(iprint_epen) =  1
    ELSE
       IF (cprint%iprint(iprint_epen).EQ.0) cprint%iprint(iprint_epen) = -1
    ENDIF
    ! cntl%tddft
    IF (cntl%tddft) THEN
       IF (cprint%iprint(iprint_etddft).EQ.0) cprint%iprint(iprint_etddft) =  1
    ELSE
       IF (cprint%iprint(iprint_etddft).EQ.0) cprint%iprint(iprint_etddft) = -1
    ENDIF
    ! Energy of external field
    IF (cntl%texadd) THEN
       IF (cprint%iprint(iprint_eext).EQ.0) cprint%iprint(iprint_eext) = 1
    ELSE
       IF (cprint%iprint(iprint_eext).EQ.0) cprint%iprint(iprint_eext) = -1
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wrener_init
  ! ==================================================================
  SUBROUTINE wrener
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(len=1)                         :: defe
    CHARACTER(len=12)                        :: defelectro, defelectro1, &
                                                defelectro2
    CHARACTER(len=14)                        :: defenergy, defenergy1, &
                                                defenergy2
    CHARACTER(len=16)                        :: defenergy3
    INTEGER                                  :: is
    LOGICAL                                  :: statusdummy
    REAL(real_8)                             :: cdft_ener, elec1, elec2, &
                                                etot1, etot2

! ==--------------------------------------------------------------==
! Initialisation of ouput

    CALL wrener_init
    ! ==--------------------------------------------------------------==
    IF (.NOT.tenergy_ok) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,8(" WARNING"))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,A)') " VXC METHOD: TOTAL ENERGY NOT CORRECT! "
       IF (paral%io_parent)&
            WRITE(6,'(8(" WARNING"))')
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,'(/,A,/,A,T49,F17.10,/,A,T49,F17.10)')&
         ' TOTAL INTEGRATED ELECTRONIC DENSITY',&
         '    IN G-SPACE =',chrg%csumg,&
         '    IN R-SPACE =',chrg%csumr
    IF (cntl%tlsd) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,T53,F13.6)')&
            ' TOTAL INTEGRATED SPIN DENSITY',chrg%csums
       IF (paral%io_parent)&
            WRITE(6,'(A,T53,F13.6)')&
            ' TOTAL INTEGRATED ABSOLUTE VALUE OF SPIN DENSITY',chrg%csumsabs
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL mm_dim(mm_go_qm,statusdummy)
    IF (pslo_com%tivan) THEN
       IF (paral%io_parent)&
            WRITE(6,2)
       DO is=1,ions1%nsp
          IF (paral%io_parent)&
               WRITE(6,3) elem%el(ions0%iatyp(is)),ions0%na(is),chrg%vdbchg(is)
       ENDDO
    ENDIF
2   FORMAT(//' VANDERBILT AUGMENTATION CHARGES (MEAN VALUE PER ATOM)',&
         /,' ATOM TYPE    NR. OF ATOMS        CHARGE ')
3   FORMAT(5x,a3,8x,i4,8x,f12.3)
    CALL mm_dim(mm_revert,statusdummy)
    ! ==--------------------------------------------------------------==
    IF (cprint%iprint(iprint_info).EQ.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,*)
       RETURN
    ENDIF
    ! ==--------------------------------------------------------------==
    ! Definition of ELEC1 (electrostatic energy if not cntl%tdiag).
    IF (cntl%tsic) THEN                           ! cmb_ssic
       defelectro1='(A-S+R+SIC)'
       elec1 =  ener_com%ehep + ener_com%esr - ener_com%eself + ener_com%ehsic
    ELSE                                    ! cmb_ssic
       defelectro1='(E1=A-S+R)'
       elec1 =  ener_com%ehep + ener_com%esr - ener_com%eself
    ENDIF
    ! Definition of ELEC2 (electrostatic energy if cntl%tdiag).
    defelectro2='(E2=I-H-S+R)'
    ! Ion-Ion interaction - e-e (double counting) [cntl%tdiag]
    elec2 =  ener_com%ehii-ener_com%ehee + ener_com%esr - ener_com%eself
    ! Definition of TOTAL ENERGY for UPDWF
    IF (cntl%tfint) THEN
       etot1 = ener_com%ekin + elec1 + ener_com%epseu + ener_com%enl + ener_com%exc + &
               ener_com%entropy + hubbu%ehub
       defenergy1='(K+E1+L+N+X+T)'
       IF(cntl%thubb) defenergy1='(K+E1+L+N+X+U+T)'
    ELSE
       etot1 = ener_com%ekin + elec1 + ener_com%epseu + ener_com%enl + ener_com%exc &
             + hubbu%ehub
       defenergy1='(K+E1+L+N+X)'
       IF(cntl%thubb) defenergy1='(K+E1+L+N+X+U)'
    ENDIF
    ! Definition of TOTAL ENERGY for cntl%tdiag
    etot2 = ener_com%eeig+(ener_com%exc-ener_com%vxc)+ener_com%eht+ener_com%ebogo
    IF (cntl%tfint) THEN
       defe='F'
    ELSE
       defe='B'
    ENDIF
    IF (fint1%tbogo) THEN
       IF (paral%io_parent)&
            WRITE(defenergy2,'("(",A1,"+E2+X-V+O)")') defe
    ELSE
       IF (paral%io_parent)&
            WRITE(defenergy2,'("(",A1,"+E2+X-V)")') defe
    ENDIF
    IF (lqmmm%qmmm) THEN
      IF(cntl%thubb)THEN 
         defenergy3='(K+E1+L+N+X+U+Q+M)'
      ELSE 
         defenergy3='(K+E1+L+N+X+Q+M)'
      ENDIF
    ENDIF
    IF (cntl%cdft)THEN
       cdft_ener=cdftcom%cdft_v(1)*(cdftcom%cdft_nc+cdftcom%vgrad(1))
       IF (cdftlog%tcall)cdft_ener=cdft_ener+cdftcom%cdft_v(2)*(cdftcom%cdft_ns+cdftcom%vgrad(2))
    ENDIF
    IF (cntl%tdiag) THEN
       defelectro=defelectro2
       defenergy=defenergy2
    ELSE
       IF (lspin2%tlse .AND. (lspin2%tcas22.OR.lspin2%tpenal)) THEN
          defenergy='(K+E+L+N+X+C)'
       ELSE
          defenergy=defenergy1
       ENDIF
       defelectro=defelectro1
    ENDIF
    IF (response1%response_running) THEN
       IF (paral%io_parent)&
            WRITE(6,8) ener_com%etot,ener1%eh0,ener1%eloc1,ener1%enl1,ener1%eht1,ener1%exc1,ener1%elag1
       RETURN
8      FORMAT(/ '   TOTAL SECOND ORDER  ENERGY = ',t41,f20.8,' A.U.'/&
            'ZERO ORDER HAMILTONIAN ENERGY = ',t41,f20.8,' A.U.'/&
            '     BARE PERTURBATION ENERGY = ',t41,f20.8,' A.U.'/&
            ' N-L BARE PERTURBATION ENERGY = ',t41,f20.8,' A.U.'/&
            ' HARTREE  PERTURBATION ENERGY = ',t41,f20.8,' A.U.'/&
            '   X-C    PERTURBATION ENERGY = ',t41,f20.8,' A.U.'/&
            'ZERO ORDER LAGR-MULT.  ENERGY = ',t41,f20.8,' A.U.'/)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (cprint%iprint(iprint_info).NE.-1) THEN
       IF (lqmmm%qmmm) THEN
          IF (paral%io_parent)&
               WRITE(6,'(/,1X,A,T25,A,T41,F20.8,A)')&
               defenergy3,'TOTAL ENERGY =',ener_com%etot,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,T22,A,T41,F20.8,A)')&
               defenergy,'TOTAL QM ENERGY =',eqm_r%eqm,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,T19,A,T41,F20.8,A)')&
               '(Q)','TOTAL QM/MM ENERGY =',eqm_r%eqmmm,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,T22,A,T41,F20.8,A)')&
               '(M)','TOTAL MM ENERGY =',eqm_r%emm,' A.U.'
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,T22,A,T41,F20.8,A)')&
               '   ','DIFFERENCE      =',&
               ener_com%etot-(eqm_r%eqm+eqm_r%eqmmm+eqm_r%EMM),' A.U.'
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(/,1X,A,T25,A,T41,F20.8,A)')&
               defenergy,'TOTAL ENERGY =',ener_com%etot,' A.U.'
       ENDIF
    ENDIF
    IF (cprint%iprint(iprint_etot1).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T24,A,T41,F20.8,A)')&
            defenergy1,'ENERGY SUM(1) =',etot1,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_etot2).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T24,A,T41,F20.8,A)')&
            defenergy2,'ENERGY SUM(2) =',etot2,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_ekin).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T23,A,T41,F20.8,A)')&
            '(K)','KINETIC ENERGY =',ener_com%ekin,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_eeig).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T15,A,T41,F20.8,A)')&
            '(F)','ELECTRONIC FREE ENERGY =',ener_com%eeig,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_eband).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T26,A,T41,F20.8,A)')&
            '(B)','BAND ENERGY =',ener_com%eband,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_entropy).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T26,A,T41,F20.8,A)')&
            '(T)','-KT*ENTROPY =',ener_com%entropy,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_eht).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T17,A,T41,F20.8,A)')&
            defelectro,'ELECTROSTATIC ENERGY =',ener_com%eht,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_elec1).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T17,A,T41,F20.8,A)')&
            defelectro1,'ELECTROSTATIC ENERGY =',elec1,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_elec2).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T17,A,T41,F20.8,A)')&
            defelectro2,'ELECTROSTATIC ENERGY =',elec2,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_ehep).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T11,A,T41,F20.8,A)')&
            '(A)','(E+I)-(E+I) HARTREE ENERGY =',ener_com%ehep,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_ehii).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T11,A,T41,F20.8,A)')&
            '(I)','(PSEUDO CHARGE I-I) ENERGY =',ener_com%ehii,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_ehee).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T17,A,T41,F20.8,A)')&
            '(H)','(E-E) HARTREE ENERGY =',ener_com%ehee,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_eself).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T32,A,T41,F20.8,A)')&
            '(S)','ESELF =',ener_com%eself,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_esr).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T34,A,T41,F20.8,A)')&
            '(R)','ESR =',ener_com%esr,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_epseu).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T9,A,T41,F20.8,A)')&
            '(L)','LOCAL PSEUDOPOTENTIAL ENERGY =',ener_com%epseu,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_enl).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T11,A,T41,F20.8,A)')&
            '(N)','N-L PSEUDOPOTENTIAL ENERGY =',ener_com%enl,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_exc).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T10,A,T41,F20.8,A)')&
            '(X)','EXCHANGE-CORRELATION ENERGY =',ener_com%exc,' A.U.'
    ENDIF
    IF (cntl%cdft) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T10,A,T41,F20.8,A)')&
            '(W)','CHARGE-CONSTRAINT ENERGY =',cdft_ener,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_vxc).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T10,A,T41,F20.8,A)')&
            '(V)','EXCHANGE-CORRELATION POTEN. =',ener_com%vxc,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_egc).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T11,A,T41,F20.8,A)')&
            ' ','GRADIENT CORRECTION ENERGY =',ener_com%egc,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_ebogo).NE.-1.AND.ener_com%ebogo.NE.0._real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T9,A,T41,F20.8,A)')&
            '(O)','BOGOLIUBOV CORRECTION ENERGY =',ener_com%ebogo,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_evdw).NE.-1.AND.vdwr%evdw.NE.0._real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T9,A,T41,F20.8,A)')&
            ' ','                  VDW ENERGY =',vdwr%evdw,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_etddft).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T12,A,T41,F20.6,A)')&
            '(E)','EXCITATION ENERGY(TDDFT) =',ener_com%etddft*2._real_8*ry,' e.V.'
    ENDIF
    IF (cprint%iprint(iprint_eext).NE.-1) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T12,A,T41,F20.8,A)')&
            '(E)','ENERGY OF EXTERNAL FIELD =',ener_com%eext,' A.U.'
    ENDIF
    IF (cprint%iprint(iprint_ehsic).NE.-1) THEN   ! cmb_ssic
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T13,A,T41,F20.8,A)')&
            '(H)','(E-E) HARTREE SIC ENERGY =',ener_com%ehsic,' A.U.'
    ENDIF                                 ! cmb_ssic
    IF(cprint%iprint(iprint_ehub).NE.-1) THEN   ! cmb_ssic
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T12,A,T41,F20.8,A)')&
            '(U)','HUBBARD CORRECTION ENERGY =',hubbu%ehub,' A.U.'
    ENDIF                                 ! cmb_ssic
    IF (.NOT.tenergy_ok) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,8(" WARNING"))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,A)') " VXC METHOD: TOTAL ENERGY NOT CORRECT! "
       IF (paral%io_parent)&
            WRITE(6,'(8(" WARNING"))')
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,*)
    ! ------------------------------------------------
    ! write QMMM energy components if necessary
    ! ------------------------------------------------
#if defined (__QMECHCOUPL)
    CALL mm_cpmd_wrenergy
#endif
    ! ==--------------------------------------------------------------==
    IF (lspin2%tlse .AND. lspin2%tcas22) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A)') " TWO STATE (CAS22) HAMILTONIAN "
       IF (paral%io_parent)&
            WRITE(6,'(T10,A,T41,F20.8,A)')&
            "    H(A,A)                = ",(ener_c%etot_a)," A.U."
       IF (paral%io_parent)&
            WRITE(6,'(T10,A,T41,F20.8,A)')&
            "    H(B,B)                = ",(ener_d%etot_b)," A.U."
       IF (paral%io_parent)&
            WRITE(6,'(T10,A,T41,F20.8,A)')&
            "    H(A,B)                = ",(ener_c%etot_ab)," A.U."
       IF (paral%io_parent)&
            WRITE(6,'(T10,A,T41,F20.8," eV  ")')&
            "    H(B,B) - H(A,A)       = ",(ener_d%etot_b-ener_c%etot_a)*ry*2._real_8
       IF (paral%io_parent)&
            WRITE(6,'(T10,A,T41,F20.8," eV  ")')&
            "    Double Excitation     = ",(ener_c%etot_2-ener_c%etot_a)*ry*2._real_8
       IF (paral%io_parent)&
            WRITE(6,'(T10,A,T41,F20.8," eV  ")')&
            "    S-T Splitting         = ",(ener_d%etot_b-ener_d%etot_t)*ry*2._real_8
       IF (paral%io_parent)&
            WRITE(6,'(T10,A,T41,F20.8," Deg.")')&
            "    ROTATION ANGLE        = ",180._real_8*ener_d%casang/pi
       IF (paral%io_parent)&
            WRITE(6,*)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wrener
  ! ==================================================================
  SUBROUTINE wreigen(eigv,f,amu,nstate)
    ! ==--------------------------------------------------------------==
    ! mb-> Higher precision (needed !) on eigenvalues printout
    REAL(real_8)                             :: amu
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate,nkpt%nkpts), &
                                                eigv(nstate,nkpt%nkpts)

    CHARACTER(*), PARAMETER                  :: procedureN = 'wreigen'

    CHARACTER(len=2)                         :: ac(2)
    INTEGER                                  :: i, ic(0:1), ii, ik, isub, j, &
                                                nleft

    CALL tiset(procedureN,isub)
    ac(1)="  "
    ac(2)="NC"
    ! ==--------------------------------------------------------------==
    IF (cntl%tlsd) THEN
       IF (nac.LT.0) nac=spin_mod%nsup
       IF (nbc.LT.0) nbc=spin_mod%nsdown
    ELSE
       IF (nac.LT.0) nac=nstate
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,'(/,A)') ' EIGENVALUES(EV) AND OCCUPATION:'
    DO ik=1,nkpt%nkpts
       IF (tkpts%tkpnt) THEN
          IF (paral%io_parent)&
               WRITE(6,'(1X,A,I8,4F12.6)')&
               'K POINT:',ik,(rk(ii,ik),ii=1,3),wk(ik)
       ENDIF
       IF (cntl%tlsd) THEN
          IF (paral%io_parent)&
               WRITE(6,'(A)') ' ALPHA STATES:'
          DO i=1,spin_mod%nsup,2
             nleft=MIN(spin_mod%nsup-i,1)
             ic(0)=1
             ic(1)=1
             IF (nac.LT.i) ic(0)=2
             IF (nac.LT.i+1) ic(1)=2
             IF (paral%io_parent)&
                  WRITE(6,'(3X,2(I5,F15.7,A2,2X,F13.8,6X))')&
                  (i+j,eigv(i+j,ik)*(2*ry),ac(ic(j)),&
                  f(i+j,ik),j=0,nleft)
          ENDDO
          IF (paral%io_parent)&
               WRITE(6,'(A)') ' BETA STATES:'
          DO i=spin_mod%nsup+1,nstate,2
             nleft=MIN(nstate-i,1)
             ic(0)=1
             ic(1)=1
             IF (spin_mod%nsup+nbc.LT.i) ic(0)=2
             IF (spin_mod%nsup+nbc.LT.i+1) ic(1)=2
             IF (paral%io_parent)&
                  WRITE(6,'(3X,2(I5,F15.7,A2,2X,F13.8,6X))')&
                  (i+j-spin_mod%nsup,eigv(i+j,ik)*(2*ry),ac(ic(j)),&
                  f(i+j,ik),j=0,nleft)
          ENDDO
       ELSE
          DO i=1,nstate,2
             nleft=MIN(nstate-i,1)
             ic(0)=1
             ic(1)=1
             IF (nac.LT.i) ic(0)=2
             IF (nac.LT.i+1) ic(1)=2
             IF (paral%io_parent)&
                  WRITE(6,'(3X,2(I5,F15.7,A2,2X,F13.8,6X))')&
                  (i+j,eigv(i+j,ik)*(2*ry),ac(ic(j)),f(i+j,ik),j=0,nleft)
          ENDDO
       ENDIF
    ENDDO
    IF (cntl%tlsd) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T43,F20.10,A3)')&
            'ALPHA CHEMICAL POTENTIAL =',amuup*(2*ry),' EV'
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T43,F20.10,A3)')&
            ' BETA CHEMICAL POTENTIAL =',amudo*(2*ry),' EV'
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(1X,A,T43,F20.10,A3)')&
            'CHEMICAL POTENTIAL =',amu*(2*ry),' EV'
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
  END SUBROUTINE wreigen
  ! ==================================================================
  SUBROUTINE wrstress
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: ia, is, j
    REAL(real_8)                             :: fact, out(3,3)

! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE(6,'(/,A)') ' TOTAL STRESS TENSOR (kBar):'
    CALL dcopy(9,paiu,1,out,1)
    ! mb-> loop unrolled for Hitachi SR8K
    DO is=1,ions1%nsp
       fact=rmass%pma(is)
       DO ia=1,ions0%na(is)
          out(1,1)=out(1,1)+fact*velp(1,ia,is)*velp(1,ia,is)
          out(2,1)=out(2,1)+fact*velp(2,ia,is)*velp(1,ia,is)
          out(3,1)=out(3,1)+fact*velp(3,ia,is)*velp(1,ia,is)
          out(1,2)=out(1,2)+fact*velp(1,ia,is)*velp(2,ia,is)
          out(2,2)=out(2,2)+fact*velp(2,ia,is)*velp(2,ia,is)
          out(3,2)=out(3,2)+fact*velp(3,ia,is)*velp(2,ia,is)
          out(1,3)=out(1,3)+fact*velp(1,ia,is)*velp(3,ia,is)
          out(2,3)=out(2,3)+fact*velp(2,ia,is)*velp(3,ia,is)
          out(3,3)=out(3,3)+fact*velp(3,ia,is)*velp(3,ia,is)
       ENDDO
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,'(5X,3(1PE20.8))') ((out(1,j)/parm%omega)*au_kb,j=1,3)
    IF (paral%io_parent)&
         WRITE(6,'(5X,3(1PE20.8))') ((out(2,j)/parm%omega)*au_kb,j=1,3)
    IF (paral%io_parent)&
         WRITE(6,'(5X,3(1PE20.8))') ((out(3,j)/parm%omega)*au_kb,j=1,3)
    IF (paral%io_parent)&
         WRITE(6,*)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wrstress
  ! ==================================================================
  SUBROUTINE wrcell
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER                                  :: i, j

! ==--------------------------------------------------------------==

    IF (paral%io_parent)&
         WRITE(6,'(/,A,F12.4,A)') ' CELL PARAMETER (VOLUME=',parm%omega,'):'
    DO i=1,3
       IF (paral%io_parent)&
            WRITE(6,'(T21,3(F15.8))') (metr_com%ht(i,j),j=1,3)
    ENDDO
    IF (paral%io_parent)&
         WRITE(6,*)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE wrcell
  ! ==================================================================
  FUNCTION ok_print(iprint,engpri,tend)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iprint
    LOGICAL                                  :: engpri, tend, ok_print

! ==--------------------------------------------------------------==

    ok_print=(engpri.AND.(iprint.EQ.1)).OR.&
         (tend.AND.(iprint.NE.-1))
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION ok_print
  ! ==================================================================

END MODULE wrener_utils

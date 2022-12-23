MODULE mm_forces_prop_utils
  USE ehpsi_utils,                     ONLY: set_b2l
  USE ehrenfest_utils,                 ONLY: ehrenfest
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE k_updwf_utils,                   ONLY: k_updwf
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: td01
  USE lr_tddft_utils,                  ONLY: lr_tddft
  USE machine,                         ONLY: m_cputime,&
                                             m_flush
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_go_qm
  USE mm_input,                        ONLY: lqmmm
  USE mp_interface,                    ONLY: mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE ropt,                            ONLY: infi,&
                                             ropt_mod
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE updrho_utils,                    ONLY: updrho
  USE updwf_utils,                     ONLY: updwf
  USE utils,                           ONLY: zclean_k
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_forces_prop

CONTAINS

  ! ==================================================================
  SUBROUTINE mm_forces_prop(nstate,c0,c1,c2,cr,sc0,cscr,vpp,eigv,&
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
    ! EHR

    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(:,:,:), c1(*), c2(:,:), &
                                                cr(*), sc0(nkpt%ngwk,nstate), &
                                                cscr(*)
    REAL(real_8)                             :: vpp(ncpw%ngw), eigv(*), &
                                                rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: tau0(:,:,:), velp(*), &
                                                taui(*), fion(:,:,:)
    INTEGER                                  :: ifcalc, irec(:)
    LOGICAL                                  :: tfor, tinfo

    CHARACTER(*), PARAMETER                  :: procedureN = 'mm_forces_prop'

    COMPLEX(real_8), ALLOCATABLE             :: gde(:)
    INTEGER                                  :: i, ierr, ik, infr, isub, nhpsi
    LOGICAL                                  :: statusdummy, update_pot
    REAL(real_8)                             :: detot, etot0, tcpu, thl(2), &
                                                time1

! Variables
! EHR  Ehrenfest

    CALL tiset(procedureN,isub)
    ik=1
    ! ==--------------------------------------------------------------==
    ! 
    IF (paral%qmnode)THEN
       time1 = m_cputime()
       ! CALL GIVE_SCR_FORCES_DIAG(LFORCES_DIAG,TAG,NSTATE,TFOR)
       ! CALL TEST_SCR('FORCES_DIAG',TAG,LSCR,LFORCES_DIAG)
       ! ==--------------------------------------------------------------==
       ener_com%ecnstr=0.0_real_8
       ropt_mod%convwf=.FALSE.
       etot0=0.0_real_8
       IF (cntl%tdiag.AND.cntl%tlanc) THEN
          CALL set_b2l()
       ENDIF
       ! EHR
       IF (cntl%tmdeh) ALLOCATE(gde(((nkpt%ngwk*nstate+8)*cnti%mdiis)/4),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
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
          ENDIF
          ropt_mod%sdiis=.TRUE.
          ropt_mod%spcg=.TRUE.
       ENDIF
    ENDIF
    ! 
    ! EHR
    update_pot=.TRUE.
    ! ==================================================================
    ! ==                     EHRENFEST DYNAMICS                       ==
    ! ==================================================================
    ! 
    CALL mm_dim(mm_go_qm,statusdummy)
    IF (cntl%tmdeh) THEN
       IF (infi.GT.0) THEN
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(/,1X,A)') 'occupation for prop'
             IF (paral%io_parent)&
                  WRITE(6,'(13F5.1)') (crge%f(i,ik),i=1,nstate)
          ENDIF
          ! call m_flush(6)
          ! write(6,*) 'me' , me
          ! STOP
          CALL ehrenfest(c0,c2,rhoe,psi(:,1),sc0,eigv)
       ENDIF
    ENDIF
    CALL mm_dim(mm_go_mm,statusdummy)
    IF (cntl%tmdeh.AND.paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,1X,A)') 'occupation for forces computation'
       IF (paral%io_parent)&
            WRITE(6,'(13F5.1)') (crge%f(i,ik),i=1,nstate)
       CALL m_flush(6)
    ENDIF
    ! 
    IF (geq0) CALL zclean_k(c0,nstate,ncpw%ngw)
    ! 
    IF (tfor) THEN
       IF (cntl%tdiag) THEN
          ! Diagonalization
          CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
          CALL stopgm("mm_forces_diag","UPDRHO not compatible",& 
               __LINE__,__FILE__)
          CALL updrho(c0,c2,cr,sc0,cscr,vpp,tau0,fion,eigv,&
               rhoe,psi,&
               nstate,.TRUE.,tinfo,ropt_mod%calste,infr,thl,nhpsi)
       ELSE
          IF (.NOT.tkpts%tkpnt) THEN
             CALL updwf(c0(:,:,1),c2,sc0,tau0,fion,cr,cscr,vpp,eigv,&
                  rhoe,psi,nstate,.TRUE.,.TRUE.)
          ELSEIF (tkpts%tkpnt.AND.(.NOT.cntl%tmdeh)) THEN
             ! CALL K_UPDWF(C0,C2,SC0,TAU0,FION,PME,GDE,VPP,EIGV,
             CALL k_updwf(c0,c2,sc0,tau0,fion,cr,gde,vpp,eigv,&
                  rhoe,psi,nstate,.FALSE.,1)
          ELSE
             CALL updwf(c0(:,:,1),c2,sc0,tau0,fion,cr,cscr,vpp,eigv,&
                  rhoe,psi,nstate,.TRUE.,.TRUE.)
          ENDIF
       ENDIF
    ENDIF
    IF (cntl%tddft) THEN
       CALL lr_tddft(c0(:,:,1),c1,c2,sc0,rhoe,psi,tau0,fion,eigv,&
            crge%n,tfor,td01%ioutput)
       ! Now we are ready to calculate the QM/MM forces
       ! The full cntl%tddft QM forces are already in FION
    ENDIF
200 CONTINUE
    IF (cntl%tmdeh) DEALLOCATE(gde,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (lqmmm%qmmm) CALL mp_sync(parai%qmmmgrp)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE mm_forces_prop
  ! ==================================================================

END MODULE mm_forces_prop_utils

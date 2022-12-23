MODULE mm_forces_diag_utils
  USE ehpsi_utils,                     ONLY: set_b2l
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: td01
  USE lr_tddft_utils,                  ONLY: lr_tddft
  USE machine,                         ONLY: m_walltime
  USE mm_extrap,                       ONLY: numcold
  USE mm_input,                        ONLY: clc,&
                                             lqmmm
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE parac,                           ONLY: parai,&
                                             paral
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE store_types,                     ONLY: store1
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             ncpw,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE updrho_utils,                    ONLY: give_scr_updrho,&
                                             updrho
  USE updwf_utils,                     ONLY: give_scr_updwf,&
                                             updwf
  USE wrener_utils,                    ONLY: wrprint_wfopt
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_forces_diag
  PUBLIC :: give_scr_mm_forces_diag

CONTAINS

  ! ==================================================================
  SUBROUTINE mm_forces_diag(nstate,c0,c1,c2,cr,sc0,cscr,vpp,eigv,&
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
    COMPLEX(real_8)                          :: c0(:,:), c1(*), &
                                                c2(nkpt%ngwk,nstate), cr(*), &
                                                sc0(nkpt%ngwk,nstate), cscr(*)
    REAL(real_8)                             :: vpp(ncpw%ngw), eigv(*), &
                                                rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)
    REAL(real_8)                             :: tau0(:,:,:), velp(:,:,:), &
                                                taui(:,:,:), fion(:,:,:)
    INTEGER                                  :: ifcalc, irec(:)
    LOGICAL                                  :: tfor, tinfo

    CHARACTER(*), PARAMETER                  :: procedureN = 'mm_forces_diag'

    INTEGER                                  :: infr, isub, nhpsi
    LOGICAL                                  :: update_pot
    REAL(real_8)                             :: detot, etot0, tcpu, thl(2), &
                                                time1, time2, tthl

! Variables
! ==--------------------------------------------------------------==
! 

    CALL tiset(procedureN,isub)
    IF (paral%qmnode)THEN
       time1 =m_walltime()
       ! ==--------------------------------------------------------------==
       ener_com%ecnstr=0.0_real_8
       ropt_mod%convwf=.FALSE.
       etot0=0.0_real_8
       IF (cntl%tdiag.AND.cntl%tlanc) THEN
          CALL set_b2l()

          IF (clc%classical)THEN
             ener_com%ekin=0._real_8
             ener_com%etot=0._real_8
             CALL zeroing(c2)!,ngw*nstate)
             CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
             GOTO 100
          ENDIF
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
          ENDIF
          ropt_mod%sdiis=.TRUE.
          ropt_mod%spcg=.TRUE.
       ENDIF
    ENDIF
    ! 
    update_pot=.TRUE.
    ! ==================================================================
    ! ==                        MAIN  LOOP                            ==
    ! ==================================================================
    DO infr=1,cnti%nomore_iter
       time1=m_walltime()
       ! calculate the electrostatic potential in the first loop only
       IF (lqmmm%qmmm .AND.infr.GT.1 ) update_pot = .FALSE.
       ifcalc=ifcalc+1
       ! Diagonalization
       IF (cntl%tdiag) THEN
          CALL stopgm("mm_forces_diag","UPDRHO not compatible",& 
               __LINE__,__FILE__)
          CALL updrho(c0,c2,cr,sc0,cscr,vpp,tau0,fion,eigv,&
               rhoe,psi,&
               nstate,.FALSE.,tinfo,.FALSE.,infr,thl,nhpsi)
       ELSE
          CALL updwf(c0,c2,sc0,tau0,fion,cr,cscr,vpp,eigv,&
               rhoe,psi,nstate,.FALSE.,update_pot)
       ENDIF
       ! Printout the evolution of the iterative optimization
       IF (paral%parent) THEN
          detot=ener_com%etot-etot0
          IF (infr.EQ.1) detot=0.0_real_8
          time2=m_walltime()
          tcpu=(time2-time1)*0.001_real_8
          tthl=thl(1)
          IF (cntl%tlsd) tthl=0.5_real_8*(thl(1)+thl(2))
          IF (tinfo) THEN
             CALL wrprint_wfopt(eigv,crge%f,ener_com%amu,nstate,ener_com%etot,etot0,&
                  tcpu,gemax,cnorm,thl,.FALSE.,ifcalc,infr)
          ELSEIF (cntl%tdiag) THEN
             IF (paral%io_parent)&
                  WRITE(6,&
                  '(T9,I4,T15,1PE12.3,T27,0PF12.6,T39,1PE12.3,'//&
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
       IF (paral%qmnode) THEN
          IF (MOD(ifcalc,store1%isctore).EQ.0)&
               CALL zhwwf(2,irec,c0,c2,nstate,eigv,tau0,velp,taui,iteropt%nfi)
       ENDIF
       ! Check to break the loop.
       IF (lqmmm%qmmm)CALL mp_bcast(ropt_mod%convwf,parai%qmmmsource,parai%qmmmgrp)
       IF (ropt_mod%convwf) GOTO 100   ! Convergence of wavefunctions
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
    ! ==================================================================
    ! ==                      END OF MAIN  LOOP                       ==
    ! ==================================================================
300 CONTINUE
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' WARNING! DENSITY NOT CONVERGED !!!'
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' WARNING! DENSITY NOT CONVERGED !!!'
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' WARNING! DENSITY NOT CONVERGED !!!'
100 CONTINUE
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
          CALL updwf(c0,c2,sc0,tau0,fion,cr,cscr,vpp,eigv,&
               rhoe,psi,nstate,.TRUE.,.TRUE.)
       ENDIF
    ENDIF
    IF (cntl%tddft) THEN
       CALL lr_tddft(c0,c1,c2,sc0,rhoe,psi,tau0,fion,eigv,&
            crge%n,tfor,td01%ioutput)
       ! Now we are ready to calculate the QM/MM forces
       ! The full cntl%tddft QM forces are already in FION
    ENDIF
200 CONTINUE
    IF (lqmmm%qmmm) CALL mp_sync(parai%qmmmgrp)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE mm_forces_diag
  ! ==================================================================
  SUBROUTINE give_scr_mm_forces_diag(lforces_diag,tag,nstate,tfor)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lforces_diag
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: tfor

! ==--------------------------------------------------------------==

    IF (cntl%tdiag) THEN
       CALL give_scr_updrho(lforces_diag,tag,nstate,&
            tfor,tfor.AND.cntl%tpres)
    ELSE
       CALL give_scr_updwf(lforces_diag,tag,nstate,tfor)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_mm_forces_diag
  ! ==================================================================

END MODULE mm_forces_diag_utils

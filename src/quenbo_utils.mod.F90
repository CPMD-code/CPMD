MODULE quenbo_utils
  USE coor,                            ONLY: tau0
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_walltime
  USE mm_input,                        ONLY: lqmmm
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxdis,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE updwf_utils,                     ONLY: give_scr_updwf,&
                                             updwf
  USE wrener_utils,                    ONLY: wrprint_wfopt

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: quenbo
  PUBLIC :: give_scr_quenbo

CONTAINS

  ! ==================================================================
  SUBROUTINE quenbo(c0,c2,sc0,taur,rhoe,psi)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(:,:), c2(ncpw%ngw,crge%n), &
                                                sc0(ncpw%ngw,crge%n)
    REAL(real_8)                             :: taur(:,:,:), rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'quenbo'

    COMPLEX(real_8), ALLOCATABLE             :: gde(:), pme(:)
    INTEGER                                  :: ierr, isub, nmin, nnfi
    LOGICAL                                  :: update_pot
    REAL(real_8)                             :: etoto, tcpu, thl(2), time1, &
                                                time2
    REAL(real_8), ALLOCATABLE                :: eigv(:), vpp(:)

    CALL tiset(procedureN,isub)
    IF (paral%qmnode)THEN
       ! ==--------------------------------------------------------------==
       IF (paral%io_parent)&
            WRITE(6,'(/,A,A,/)')&
            ' >>>>>>>> QUENCH SYSTEM TO THE ',&
            'BORN-OPPENHEIMER SURFACE <<<<<<<< '
       ALLOCATE(eigv(crge%n),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (cntl%tsde) THEN
          ALLOCATE(pme(crge%n),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(gde(crge%n),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE IF (cntl%diis) THEN
          ALLOCATE(pme((ncpw%ngw*crge%n+8)*cnti%mdiis),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(gde(((ncpw%ngw*crge%n+8)*cnti%mdiis)/4),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE IF (cntl%pcg) THEN
          ALLOCATE(pme(ncpw%ngw*crge%n),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(gde(crge%n),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          cntl%diis=.TRUE.
          cntl%prec=.TRUE.
          cnti%mdiis=MIN(cnti%mdiis,maxdis)
          ALLOCATE(pme((ncpw%ngw*crge%n+8)*cnti%mdiis),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(gde(((ncpw%ngw*crge%n+8)*cnti%mdiis)/4),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(vpp(ncpw%ngw),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF (paral%parent) CALL prmem('    QUENBO')
       ropt_mod%convwf=.FALSE.
       ropt_mod%sdiis=.TRUE.
       ropt_mod%spcg=.TRUE.
       ropt_mod%modens=.FALSE.
       ropt_mod%engpri=.FALSE.
       IF (.NOT.cntl%tprcp) ropt_mod%calste=.FALSE.
       etoto=0.0_real_8
    ELSE
       nmin = 10
       ALLOCATE(eigv(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(pme(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gde(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vpp(nmin),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==      THE BASIC LOOP FOR WAVEFUNCTION OPTIMIZATION            ==
    ! ==--------------------------------------------------------------==
    iteropt%iinfi=0
    update_pot=.TRUE.
100 CONTINUE
    time1=m_walltime()
    ! UPDATE THE WAVEFUNCTIONS
    ! qmmm: this is for skipping the full electrostatic coupling when IINFI> 0
    ! (the potential at fixed ions does not change)
    IF (lqmmm%qmmm .AND. iteropt%iinfi.GT.0 ) update_pot=.FALSE.
    CALL updwf(c0,c2,sc0,tau0,taur,pme,gde,vpp,eigv,&
         rhoe,psi,crge%n,.FALSE.,update_pot)
    ! PRINTOUT THE EVOLUTION OF THE ITERATIVE OPTIMIZATION
    IF (paral%parent) THEN
       time2=m_walltime()
       tcpu=(time2-time1)*0.001_real_8
       nnfi=iteropt%iinfi+1
       CALL wrprint_wfopt(eigv,crge%f,ener_com%amu,crge%n,ener_com%etot,etoto,tcpu,&
            gemax,cnorm,thl,ropt_mod%convwf,nnfi,iteropt%iinfi)
    ENDIF
    iteropt%iinfi=iteropt%iinfi+1
    IF (lqmmm%qmmm)CALL mp_bcast(ropt_mod%convwf,parai%qmmmsource,parai%qmmmgrp)
    IF (.NOT.ropt_mod%convwf.AND.iteropt%iinfi<cnti%nomore_iter) GOTO 100
    ! ==--------------------------------------------------------------==
    ! ==     END OF MAIN LOOP                                         ==
    ! ==--------------------------------------------------------------==
    DEALLOCATE(pme,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(gde,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(eigv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! .. apsi 061101: Added the 'IF' (suggested by Daniel Sebastiani)
    IF (lqmmm%qmmm) CALL mp_sync(parai%qmmmgrp)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
  END SUBROUTINE quenbo

  ! ==================================================================
  SUBROUTINE give_scr_quenbo(lquenbo,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lquenbo
    CHARACTER(len=30)                        :: tag

    CALL give_scr_updwf(lquenbo,tag,crge%n,.FALSE.)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_quenbo
  ! ==================================================================

END MODULE quenbo_utils

MODULE bswfo_utils
  USE bsym,                            ONLY: autocm,&
                                             bsclcs,&
                                             cnstwgt,&
                                             restbs,&
                                             rtsasb
  USE bsympnt,                         ONLY: fnbs
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE finalp_utils,                    ONLY: finalp
  USE geofile_utils,                   ONLY: geofile
  USE gsize_utils,                     ONLY: gsize
  USE kinds,                           ONLY: real_8
  USE lsforce_utils,                   ONLY: lsforce
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm
  USE norm,                            ONLY: gnmax,&
                                             gnorm
  USE parac,                           ONLY: paral
  USE ropt,                            ONLY: bsnfi,&
                                             iteropt,&
                                             ropt_mod
  USE rwfopt_utils,                    ONLY: rwfopt
  USE setbsstate_utils,                ONLY: setbsstate
  USE setirec_utils,                   ONLY: write_irec
  USE soft,                            ONLY: soft_com
  USE store_types,                     ONLY: restart1,&
                                             store1
  USE system,                          ONLY: cnti,&
                                             maxsys
  USE wrgeo_utils,                     ONLY: wrgeof
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: bs_wfo

CONTAINS

  ! ==================================================================
  SUBROUTINE bs_wfo(c0,c2,sc0,pme,gde,vpp,eigv)
    ! CB: Manages wf optimization with BROKEN SYMMETRY option      
    ! ==--------------------------------------------------------------==
    ! QMMM

    COMPLEX(real_8), INTENT(inout)           :: c0(:,:,:), c2(:,:), sc0(:,:,:)
    COMPLEX(real_8)                          :: pme(:), gde(:)
    REAL(real_8)                             :: vpp(:), eigv(crge%n,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'bs_wfo'

    INTEGER                                  :: bsnmr, ierr, irec(100)
    LOGICAL                                  :: bsconv, bscycle, hsconv, &
                                                oldstatus
    REAL(real_8)                             :: couplj, etotbs, etoths, &
                                                etotls, jwn

! Variables
! ==--------------------------------------------------------------==
! Format for for BS messages

111 FORMAT (/,t2,a,/)
    ! If QMMM, change to MM dimensions
    CALL mm_dim(mm_go_mm,oldstatus)
    ! 
    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fnbs(3*maxsys%nax*maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(taup)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    CALL zeroing(fnbs)!,3*maxsys%nax*maxsys%nsx)
    ! ==--------------------------------------------------------------==
    CALL write_irec(irec)
    bsnmr=cnti%nomore
    ! CB: RWFOPT returns for every writing of RESTART
    cnti%nomore=MIN(store1%istore-1,bsnmr-1)
    bsnfi=0
    ! ==--------------------------------------------------------------==
    restbs=.TRUE.
    IF (paral%io_parent)&
         WRITE(6,111)&
         ' BS_WFO| CALCULATING BROKEN SYMMETRY STATE'
    CALL rwfopt(c0(:,:,:),c2,sc0,pme,gde,vpp,eigv)
    bsnfi=bsnfi+iteropt%nfi
    etotbs=ener_com%etot
    ! 
    IF (paral%parent)CALL dcopy(3*maxsys%nax*maxsys%nsx,fion,1,fnbs,1)
    bsconv=ropt_mod%convwf
    IF (.NOT.bsconv.AND.paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,111)'BS_WFO| UNCONVERGED BROKEN SYMMETRY STATE:'
       CALL finalp(tau0,fion,tau0,eigv)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (restbs) restart1%rwf=.FALSE.
    ! CB: Second wave function (HS)
    bsclcs=2
    ! CB: Setup occupation for HS state
    CALL setbsstate
    IF (.NOT.soft_com%exsoft) THEN
       cnti%nomore=MIN(store1%istore,bsnmr)-bsnfi
    ELSE
       cnti%nomore=1
       IF (paral%io_parent)&
            WRITE(6,111)'BS_WFO| EXIT AFTER ONE MORE CYCLE'
    ENDIF
    IF (paral%io_parent)&
         WRITE(6,111)'BS_WFO| SWITCHING TO HIGH SPIN STATE'
    CALL rwfopt(c0(:,:,2:2),c2,sc0,pme,gde,vpp,eigv(1,2))
    cnti%nomore=MIN(store1%istore,bsnmr)
    bsnfi=bsnfi+iteropt%nfi
    etoths=ener_com%etot
    hsconv=ropt_mod%convwf
    IF (.NOT.hsconv.AND.paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,111)'BS_WFO| UNCONVERGED HIGH SPIN STATE'
       CALL finalp(tau0,fion,tau0,eigv)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (soft_com%exsoft) THEN
       bsconv=.TRUE.
       hsconv=.TRUE.
    ENDIF
    CALL zhwwf(2,irec,c0,c2,crge%n,eigv,tau0,velp,taup,bsnfi)
    restart1%restart=.TRUE.
    restart1%rlate=.TRUE.
    restart1%rwf=.TRUE.
    ! ==--------------------------------------------------------------==
    ! CB: First wave function (HS)
    bsclcs=1
    ! CB: Setup occupation for BS state
    CALL setbsstate
    ropt_mod%convwf=bsconv
    bscycle=.NOT.ropt_mod%convwf.AND.(bsnfi.LT.bsnmr)
    DO WHILE(.NOT.ropt_mod%convwf.AND.(bsnfi.LT.bsnmr))
       IF (paral%io_parent)&
            WRITE(6,111)'BS_WFO| REINITIALIZING BROKEN SYMMETRY STATE'

       CALL rwfopt(c0(:,:,:),c2,sc0,pme,gde,vpp,eigv)
       bsnfi=bsnfi+iteropt%nfi
       CALL zhwwf(2,irec,c0,c2,crge%n,eigv,tau0,velp,taup,bsnfi)
    ENDDO
    bsconv=ropt_mod%convwf
    IF (bscycle) THEN
       etotbs=ener_com%etot
       CALL dcopy(3*maxsys%nax*maxsys%nsx,fion,1,fnbs,1)
       IF (.NOT.bsconv.AND.paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,111)'BS_WFO| UNCONVERGED BROKEN SYMMETRY STATE:'
          CALL finalp(tau0,fion,tau0,eigv)
       ENDIF
    ENDIF
    IF (soft_com%exsoft) hsconv=.TRUE.
    ! ==--------------------------------------------------------------==
    ! CB: Second wave function (HS)
    bsclcs=2
    ! CB: Setup occupation for HS state
    CALL setbsstate
    ropt_mod%convwf=hsconv
    bscycle=.NOT.ropt_mod%convwf.AND.(bsnfi.LT.bsnmr)
    DO WHILE(.NOT.ropt_mod%convwf.AND.(bsnfi.LT.bsnmr))
       IF (paral%io_parent)&
            WRITE(6,111)'BS_WFO| REINITIALZING HIGH SPIN STATE'
       CALL rwfopt(c0(:,:,2:2),c2,sc0,pme,gde,vpp,eigv(1,2))
       bsnfi=bsnfi+iteropt%nfi
       CALL zhwwf(2,irec,c0,c2,crge%n,eigv,tau0,velp,taup,bsnfi)
    ENDDO
    hsconv=ropt_mod%convwf
    IF (bscycle) THEN
       etoths=ener_com%etot
       IF (.NOT.hsconv.AND.paral%parent) THEN
          IF (paral%io_parent)&
               WRITE(6,111)'BS_WFO| UNCONVERGED HIGH SPIN STATE:'
          CALL finalp(tau0,fion,tau0,eigv)
       ENDIF
    ENDIF
    IF (soft_com%exsoft) bsconv=.TRUE.
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (.NOT.bsconv.AND.paral%io_parent)&
            WRITE(6,111)'BS_WFO| BROKEN SYMMETRY STATE NOT CONVERGED! !'
       IF (.NOT.hsconv.AND.paral%io_parent)&
            WRITE(6,111)'BS_WFO| HIGH SPIN STATE NOT CONVERGED! !'
       IF (bsconv.AND.hsconv) THEN
          etotls = (1.0_real_8 + cnstwgt) * etotbs - cnstwgt * etoths
          couplj = (etoths - etotbs) * rtsasb
          jwn=couplj*autocm
          CALL lsforce(fnbs,fion)
          CALL gsize(fion,gnmax,gnorm)
          IF (paral%io_parent)&
               WRITE(6,'(/,T2,A)') 'BROKEN SYMMETRY CALCULATIONS CONVERGED:'
          IF (paral%io_parent)&
               WRITE(6,'(T2,A47,F17.8,A7)')&
               'THE PROJECTED ENERGY OF THE LOW SPIN STATE IS:',&
               etotls, 'a.u.'
          IF (paral%io_parent)&
               WRITE(6,'(T2,A47,F17.8,A7)')&
               'THE COUPLING CONSTANT J IS                   :',&
               couplj, 'a.u.'
          IF (paral%io_parent)&
               WRITE(6,'(T2,A47,F17.3,A7)') '~', jwn, 'cm^-1.'
          IF (paral%io_parent)&
               WRITE(6,'(/,T2,A)') 'BS_WFO| PROJECTED ATOMIC FORCES:'
          CALL geofile(tau0,fion,'WRITE')
          CALL wrgeof(tau0,fion)
          IF (paral%io_parent)&
               WRITE(6,'(A)') ' NUCLEAR GRADIENT: '
          IF (paral%io_parent)&
               WRITE(6,'(2(A,1PE15.5))') '    MAX. COMPONENT =',&
               gnmax,'         NORM =',gnorm
       ENDIF
    ENDIF
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fnbs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE bs_wfo

END MODULE bswfo_utils

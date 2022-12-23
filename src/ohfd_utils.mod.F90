MODULE ohfd_utils
  USE andp,                            ONLY: rin0,&
                                             rmix,&
                                             rout0
  USE canon_utils,                     ONLY: canon,&
                                             give_scr_canon
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE copot_utils,                     ONLY: copot,&
                                             give_scr_copot
  USE dynit_utils,                     ONLY: dynit
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE forcedr_driver,                  ONLY: forcedr
  USE forcep_utils,                    ONLY: rhoe_psi_size
  USE forces_diag_utils,               ONLY: forces_diag,&
                                             give_scr_forces_diag
  USE initrun_driver,                  ONLY: initrun
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_flush
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE nlcc,                            ONLY: corel
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE poin,                            ONLY: rhoo
  USE pslo,                            ONLY: pslo_com
  USE rhoofr_utils,                    ONLY: give_scr_rhoofr,&
                                             rhoofr
  USE rinitwf_utils,                   ONLY: give_scr_rinitwf
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: iteropt,&
                                             ropt_mod
  USE setirec_utils,                   ONLY: read_irec,&
                                             write_irec
  USE soft,                            ONLY: soft_com
  USE spin,                            ONLY: clsd
  USE store_types,                     ONLY: irec_rho,&
                                             restart1
  USE system,                          ONLY: cntr,&
                                             fpar,&
                                             maxsys,&
                                             ncpw
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ohfd
  PUBLIC :: give_scr_ohfd

CONTAINS

  ! ==================================================================
  SUBROUTINE ohfd(c0,cm,cn,c2,sc0,vpp)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(ncpw%ngw,crge%n,1), cm(*), &
                                                cn(*), c2(ncpw%ngw,crge%n), &
                                                sc0(ncpw%ngw,*)
    REAL(real_8)                             :: vpp(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ohfd'

    CHARACTER(len=100)                       :: filen
    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8), ALLOCATABLE             :: cref(:), psi(:,:)
    INTEGER                                  :: i, ierr, ifcalc, ii, iii, &
                                                il_psi_1d, il_psi_2d, &
                                                il_rhoe_1d, il_rhoe_2d, &
                                                irec(100), j, knfi, lscr, nnx
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: ekin1, ekin2, ekincp, ekinh1, &
                                                ekinh2, temp1, temp2
    REAL(real_8), ALLOCATABLE                :: eigm(:), eigp(:), eigv(:), &
                                                hmat(:,:), rhoe(:,:), &
                                                rhom(:,:), rhop(:,:), scr(:)

! Variables
! ==================================================================

    ALLOCATE(taup(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(fion(3,maxsys%nax,maxsys%nsx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(taup)!,SIZE(taup))
    CALL zeroing(fion)!,SIZE(fion))
    ! Memory for densities
    nnx=fpar%nnr1*clsd%nlsd
    ALLOCATE(rin0(fpar%nnr1,clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rout0(fpar%nnr1,clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rmix(fpar%nnr1,clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    rhoo => rin0
    CALL zeroing(rin0)!,nnx)
    CALL zeroing(rout0)!,nnx)
    CALL zeroing(rmix)!,nnx)
    ALLOCATE(rhop(fpar%nnr1,clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(rhom(fpar%nnr1,clsd%nlsd),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(cref(ncpw%ngw*crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.FALSE.
    ALLOCATE(eigv(crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eigp(crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(eigm(crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(hmat(crge%n,crge%n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! SCR ALLOCATION AND PARTITION (SCRATCH ARRAY).
    CALL rhoe_psi_size(il_rhoe_1d=il_rhoe_1d, il_rhoe_2d=il_rhoe_2d,&
         il_psi_1d=il_psi_1d,il_psi_2d=il_psi_2d)
    ALLOCATE(rhoe(il_rhoe_1d,il_rhoe_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi(il_psi_1d,il_psi_2d),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL give_scr_ohfd(lscr,tag)
    ALLOCATE(scr(lscr),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! == INITIALISATION                                               ==
    ! ==--------------------------------------------------------------==
    ! TIME STEP FUNCTIONS
    CALL dynit(ekincp,ekin1,ekin2,temp1,temp2,ekinh1,ekinh2)
    ropt_mod%modens=.FALSE.
    ropt_mod%engpri=.FALSE.
    ropt_mod%calste=.FALSE.
    ! Set IREC to restart file.
    CALL read_irec(irec)
    ! INITIALIZATION OF WAVEFUNCTION AND COORDINATES
    CALL initrun(irec,c0,c2,sc0,rhoe,psi,eigv)
    CALL write_irec(irec)
    ! ==--------------------------------------------------------------==
    CALL phfac(tau0)
    IF (corel%tinlc) CALL copot(rhoe,psi,ropt_mod%calste)
    IF (pslo_com%tivan) THEN
       CALL rnlsm(c0(:,:,1),crge%n,1,1,.FALSE.)
    ENDIF
    IF (irec(irec_rho).EQ.0) THEN
       CALL rhoofr(c0(:,:,1),rhoe,psi(:,1),crge%n)
       CALL dcopy(nnx,rhoe,1,rin0,1)
    ENDIF
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,"==",T25,A,T64,"==")')&
            '   REFERENCE POINT'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
    ENDIF
    ifcalc=0
    CALL forces_diag(crge%n,c0,c2,cm,sc0,cn,vpp,eigv,&
         rhoe,psi,&
         tau0,velp,taup,fion,ifcalc,&
         irec,.FALSE.,.TRUE.)
    CALL dcopy(2*ncpw%ngw*crge%n,c0,1,cref,1)
    CALL zhwwf(2,irec,c0,c2,crge%n,eigv,tau0,tau0,tau0,iteropt%nfi)
    ! ..transform to canonical orbitals and get eigenvalues
    CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,tau0,fion,eigv,&
         crge%n,1,.FALSE.,.FALSE.)
    CALL canon(c0,c2,crge%f,crge%n,eigv)
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
       IF (paral%io_parent)&
            WRITE(6,'(1X,"==",T20,A,T64,"==")')&
            'END OF REFERENCE CALCULATION'
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="),/)')
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == END INITIALIZATION                                           ==
    ! ==--------------------------------------------------------------==
    ! Open the FINDIF file
    filen='FINDIF'
    IF (paral%io_parent)&
         CALL fileopen(24,filen,fo_def,ferror)
    IF (restart1%roh.AND.paral%io_parent) REWIND(24)
    ! ==================================================================
    ! ==   LOOP OVER ORBITAL OCCUPATIONS                              ==
    ! ==================================================================
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/," ",22("*"),A,22("*"),/)') ' FINITE DIFFERENCES '
    ENDIF
    DO knfi=1,crge%n
       CALL mp_sync(parai%allgrp)
       IF (paral%io_parent)&
            WRITE(6,&
            '(" **** STATE=",I7,T24,A,F12.6,T54,A,F7.5)')&
            knfi,'EIGVAL(a.u.)=',eigv(knfi),'DISP=',cntr%fdiff
       IF (restart1%roh.AND.paral%parent) THEN
          IF (paral%io_parent)&
               READ(24,fmt=*,END=300)
          IF (paral%io_parent)&
               READ(24,*) (hmat(i,knfi),i=1,crge%n)
          GOTO 301
300       CONTINUE
          restart1%roh=.FALSE.
       ENDIF
301    CONTINUE
       CALL mp_bcast(restart1%roh,parai%source,parai%allgrp)
       IF (.NOT.restart1%roh) THEN
          crge%f(knfi,1)=crge%f(knfi,1)+cntr%fdiff
          ifcalc=0
          CALL dcopy(nnx,rhoe,1,rhop,1)
          CALL dcopy(2*ncpw%ngw*crge%n,cref,1,c0,1)
          CALL forces_diag(crge%n,c0,c2,cm,sc0,cn,vpp,eigp,&
               rhop,psi,&
               tau0,velp,taup,fion,ifcalc,&
               irec,.FALSE.,.FALSE.)
          CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,tau0,fion,eigp,&
               crge%n,1,.FALSE.,.FALSE.)
          CALL canon(c0,c2,crge%f,crge%n,eigp)
          IF (soft_com%exsoft) GOTO 100
          IF (paral%io_parent)&
               WRITE(6,&
               '(" **** STATE=",I7,T24,A,F12.6,T54,A,F7.5)')&
               knfi,'EIGVAL(a.u.)=',eigv(knfi),'DISP=',-cntr%fdiff
          crge%f(knfi,1)=crge%f(knfi,1)-2._real_8*cntr%fdiff
          ifcalc=0
          CALL dcopy(nnx,rhoe,1,rhom,1)
          CALL dcopy(2*ncpw%ngw*crge%n,cref,1,c0,1)
          CALL forces_diag(crge%n,c0,c2,cm,sc0,cn,vpp,eigm,&
               rhom,psi,&
               tau0,velp,taup,fion,ifcalc,&
               irec,.FALSE.,.FALSE.)
          CALL forcedr(c0(:,:,1),c2,sc0,rhoe,psi,tau0,fion,eigm,&
               crge%n,1,.FALSE.,.FALSE.)
          CALL canon(c0,c2,crge%f,crge%n,eigm)
          crge%f(knfi,1)=crge%f(knfi,1)+cntr%fdiff
          DO i=1,crge%n
             hmat(i,knfi)=0.5_real_8*(eigp(i)-eigm(i))/(2._real_8*cntr%fdiff)
          ENDDO
          IF (paral%parent) THEN
             IF (paral%io_parent)&
                  WRITE(24,'(I10,1PE20.10)') knfi,cntr%fdiff
             IF (paral%io_parent)&
                  WRITE(24,*) (hmat(i,knfi),i=1,crge%n)
             CALL m_flush(24)
          ENDIF
       ENDIF
       IF (soft_com%exsoft) GOTO 100
       ! ==================================================================
       ! ==     END OF MAIN LOOP                                         ==
       ! ==================================================================
    ENDDO
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/,A,A)')&
            ' *********************** ORBITAL HARDNESS',&
            ' ***********************'
       IF (paral%io_parent)&
            WRITE(6,*)
       DO i=1,crge%n,8
          ii=MIN(8,crge%n-i+1)
          IF (paral%io_parent)&
               WRITE(6,'(6X,8I8)') (iii,iii=i,i+ii-1)
          DO j=1,crge%n
             IF (paral%io_parent)&
                  WRITE(6,'(I6,8F8.3)') j,(hmat(j,iii),iii=i,i+ii-1)
          ENDDO
       ENDDO
       IF (paral%io_parent)&
            WRITE(6,'(A,A)')&
            ' ****************************************',&
            '*** *********************'
       IF (paral%io_parent)&
            WRITE(6,*)
       IF (paral%io_parent)&
            CALL fileopen(20,'HARDNESS',fo_def,ferror)
       DO i=1,crge%n
          IF (paral%io_parent)&
               WRITE(20,*) (hmat(i,j),j=1,crge%n)
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(20)
    ENDIF
    ! ==--------------------------------------------------------------==
100 CONTINUE
    IF ((paral%parent).AND.paral%io_parent)&
         CALL fileclose(24)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(taup,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fion,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rin0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rout0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rmix,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cref,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhop,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhom,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(rhoe,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psi,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(scr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(hmat,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ohfd
  ! ==================================================================
  SUBROUTINE give_scr_ohfd(lohfd,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lohfd
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lcanon, lcopot, lforces_diag, &
                                                lrhoofr, lrinitwf, lrnlsm, &
                                                nstate

    nstate=crge%n
    lrnlsm=0
    lcopot=0
    CALL give_scr_rinitwf(lrinitwf,tag,nstate)
    IF (pslo_com%tivan) CALL give_scr_rnlsm(lrnlsm,tag,nstate,.FALSE.)
    CALL give_scr_rhoofr(lrhoofr,tag)
    CALL give_scr_forces_diag(lforces_diag,tag,nstate,.TRUE.)
    CALL give_scr_canon(lcanon,tag,nstate)
    IF (corel%tinlc) CALL give_scr_copot(lcopot,tag)
    lohfd=MAX(lrinitwf,lrnlsm,lrhoofr,lforces_diag,lcopot,lcanon)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_ohfd
  ! ==================================================================

END MODULE ohfd_utils

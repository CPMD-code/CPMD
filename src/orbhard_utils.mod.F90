MODULE orbhard_utils
  USE atwf,                            ONLY: atwp
  USE ddip,                            ONLY: lenbk
  USE elct,                            ONLY: crge
  USE elct2,                           ONLY: tfixo
  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnlalloc,&
                                             fnldealloc
  USE inscan_utils,                    ONLY: inscan
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE linres,                          ONLY: lrhd,&
                                             lrhi
  USE ohfd_utils,                      ONLY: ohfd
  USE ohlr_utils,                      ONLY: ohlr
  USE parac,                           ONLY: parai,&
                                             paral
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             ncpw,&
                                             nkpt
  USE utils,                           ONLY: nxxfun
  USE wann,                            ONLY: wannl
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: orbhard

CONTAINS

  ! ==================================================================
  SUBROUTINE orbhard
    ! ==--------------------------------------------------------------==
    ! ==  ORBITAL HARDNESS MATRIX                                     ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'orbhard'

    COMPLEX(real_8), ALLOCATABLE             :: c0(:), c2(:), cm(:), cn(:), &
                                                sc0(:)
    INTEGER                                  :: ierr, nc0, nc2, ncm, ncn, &
                                                nsc0, nstate
    REAL(real_8), ALLOCATABLE                :: vpp(:)

! ==--------------------------------------------------------------==

    IF ( cntl%tqmmm ) CALL stopgm("ORBHARD","QMMM NOT IMPLEMENTED",& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (tkpts%tkpnt) CALL stopgm('ORBHARD','K-Points not implemented',& 
         __LINE__,__FILE__)
    nstate=crge%n
    ! ..now the orbitals are stored in RESTART.1
    IF (cntl%tsdan) THEN
       CALL ohin
       ! ..linear response code
       nc0=2*ncpw%ngw*nstate+8
       nc2=2*ncpw%ngw*nstate+8
       ncm=nc2
       ncn=1
       nsc0=1
       IF (lrhd%local_orb) THEN
          wannl%twann=.TRUE.
          lenbk=nxxfun(nstate)
          nsc0=MAX(2*lenbk*parai%nproc,2*nkpt%ngwk*nstate)
       ENDIF
       ALLOCATE(c0(nc0),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cm(ncm),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cn(ncn),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c2(nc2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(sc0(nsc0),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vpp(1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! ..FNL and DFNL
       CALL fnlalloc(nstate,.FALSE.,.FALSE.)
       ! 
       CALL ohlr(c0,cm,c2,sc0)
    ELSE
       ! ..finite difference code
       IF (cntl%tdiag) THEN
          tfixo=.TRUE.
          nc0=2*ncpw%ngw*nstate+8
          nc2=2*nkpt%ngwk*MAX(atwp%numaormax,nstate)+8
          IF (cntl%tdavi) THEN
             ncm=2*nkpt%ngwk*cnti%ndavv
             ncn=2*nkpt%ngwk*cnti%ndavv
          ELSE
             IF (cntl%tfrsblk) THEN
                ncm=2*nkpt%ngwk*MAX(nstate,(cnti%nkry_max+1)*cnti%nkry_block)
             ELSE
                ncm=2*nkpt%ngwk*MAX(nstate,cnti%nkry_max*cnti%nkry_block)
             ENDIF
             ncn=1
          ENDIF
          nsc0=1
       ELSE
          nc0=2*ncpw%ngw*nstate+8
          nc2=2*ncpw%ngw*nstate+8
          nsc0=1
          IF (cntl%tsde) THEN
             ncm=1
             ncn=1
          ELSE IF (cntl%diis) THEN
             ncm=(nkpt%ngwk*nstate+8)*cnti%mdiis
             ncn=((nkpt%ngwk*nstate+8)*cnti%mdiis)/4
          ELSE IF (cntl%pcg) THEN
             ncm=2*nkpt%ngwk*nstate
             ncn=1
          ENDIF
       ENDIF
       ALLOCATE(c0(nc0),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cm(ncm),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cn(ncn),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(c2(nc2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(sc0(nsc0),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(vpp(nkpt%ngwk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ! ..FNL and DFNL
       CALL fnlalloc(nstate,.FALSE.,.FALSE.)
       ! 
       CALL ohfd(c0,cm,cn,c2,sc0,vpp)
    ENDIF
    ! 
    DEALLOCATE(c0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(cn,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(c2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(sc0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL fnldealloc(.FALSE.,.FALSE.)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE orbhard
  ! ==================================================================
  SUBROUTINE ohin
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE SECTION &HARNESS &END                ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==     &HARDNESS                                                ==
    ! ==        Options                                               ==
    ! ==     &END                                                     ==
    ! ==--------------------------------------------------------------==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==  LOCALIZE                                                    ==
    ! ==  DIAGONAL [OFF]                                              ==
    ! ==  ORBITALS                                                    ==
    ! ==    numo                                                      ==
    ! ==  REFATOM                                                     ==
    ! ==    refat                                                     ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=80)                        :: line
    INTEGER                                  :: ierr, iunit

! ==--------------------------------------------------------------==
! ==  DEFAULTS                                                    ==
! ==--------------------------------------------------------------==

    IF (.NOT.paral%io_parent) GOTO 9999
    iunit = 5
    lrhd%local_orb=.FALSE.
    lrhd%diagonal=.TRUE.
    lrhi%numo=-1
    lrhi%refat=-1
    ! ==--------------------------------------------------------------==
    ierr=inscan(iunit,'&HARDNESS')
    IF (ierr.NE.0) GOTO 100
    ! ==--------------------------------------------------------------==
10  CONTINUE
    IF (paral%io_parent)&
         READ(iunit,err=99,END=99,fmt='(A80)') line
    IF (INDEX(line,'&END').NE.0) GOTO 100
    IF (INDEX(line,'LOCALIZE').NE.0) THEN
       lrhd%local_orb=.TRUE.
       GOTO 10
    ELSE IF (INDEX(line,'DIAGONAL').NE.0) THEN
       lrhd%diagonal=.TRUE.
       IF (INDEX(line,'OFF').NE.0) lrhd%diagonal=.FALSE.
       GOTO 10
    ELSE IF (INDEX(line,'ORBITAL').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit ,err=99,END=99,fmt=*) lrhi%numo
       GOTO 10
    ELSE IF (INDEX(line,'REFATOM').NE.0) THEN
       IF (paral%io_parent)&
            READ(iunit ,err=99,END=99,fmt=*) lrhi%refat
       GOTO 10
    ELSE
       ! ..Dummy line
       GOTO 10
    ENDIF
    ! ==--------------------------------------------------------------==
99  CONTINUE
    IF (paral%io_parent)&
         WRITE(6,*) ' OHIN  : ERROR IN READING INPUT FILE'
    CALL stopgm('OHIN',' ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
100 CONTINUE
    ! ==--------------------------------------------------------------==
    ! ==  Test of options                                             ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! ==  WRITE INFO TO OUTPUT                                        ==
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         WRITE(6,*)
    IF (paral%io_parent)&
         WRITE(6,'(A)') ' ORBITAL HARDNESS INPUT OPTIONS '
    IF (lrhd%diagonal) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' CALCULATE DIAGONAL ORBITAL HARDNESS'
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' CALCULATE GENERALIZE ORBITAL HARDNESS'
    ENDIF
    IF ((lrhd%local_orb).AND.paral%io_parent)&
         WRITE(6,'(A)') ' USE WANNIER FUNCTION BASIS '
    IF (lrhi%numo.LT.0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A)') ' INCLUDE ALL ORBITALS '
    ELSE
       IF (paral%io_parent)&
            WRITE(6,'(A,T60,I10)') ' NUMBER OF ORBITALS INCLUDED ',lrhi%numo
    ENDIF
    IF ((lrhi%refat.GT.0).AND.paral%io_parent)&
         WRITE(6,'(A,T60,I10)')&
         ' REFERENZ ATOM ',lrhi%refat
    IF (paral%io_parent)&
         WRITE(6,*)
    ! ==--------------------------------------------------------------==
9999 CONTINUE
    CALL mp_bcast_byte(lrhd, size_in_bytes_of(lrhd),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(lrhi, size_in_bytes_of(lrhi),parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ohin
  ! ==================================================================
END MODULE orbhard_utils

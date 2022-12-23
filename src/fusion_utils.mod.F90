MODULE fusion_utils
  USE coor,                            ONLY: tau0,&
                                             velp
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE filnmod,                         ONLY: filbod,&
                                             filn
  USE kinds,                           ONLY: real_8
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_revert
  USE parac,                           ONLY: paral
  USE readsr_utils,                    ONLY: xstring
  USE rv30_utils,                      ONLY: zhrwf
  USE setirec_utils,                   ONLY: read_irec,&
                                             write_irec
  USE shop,                            ONLY: s0_filn,&
                                             s1_filn,&
                                             sh02
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: maxsys,&
                                             ncpw
  USE wv30_utils,                      ONLY: zhwwf
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fusion
  PUBLIC :: separate

CONTAINS

  ! ==================================================================
  SUBROUTINE fusion(c0,cm)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,2*crge%n,1), &
                                                cm(ncpw%ngw,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'fusion'

    CHARACTER(len=20)                        :: filn_back
    INTEGER                                  :: ierr, irec(100), nfi, nstate
    LOGICAL                                  :: oldstatus
    REAL(real_8)                             :: eigv(1)
    REAL(real_8), ALLOCATABLE                :: taui(:,:,:)

! Variables
! McB
! McB
! ==================================================================

    IF ( paral%qmnode ) THEN

       IF ( paral%io_parent ) THEN
          WRITE(6,'(/1X,''MERGING RESTART FILES FOR'', '' SURFACE HOPPING'')')
          IF ( .NOT.(restart1%rco.AND.restart1%rwf) ) THEN
             WRITE(6,'(/1X,'' RESTART OPTION NOT ACTIVATED'')')
             CALL stopgm('FUSION','        ',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDIF

#if defined (__GROMOS)
       CALL mm_dim(mm_go_mm,oldstatus)
#endif
       ALLOCATE(taui(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(taui)!,3*maxsys%nax*maxsys%nsx)

       filn_back = filn
       CALL read_irec(irec)

       ! ...    read S0 state
       filn = s0_filn
       CALL zhrwf(1,irec,c0,cm,sh02%nst_s0,eigv,tau0,velp,taui,nfi)
       IF ( paral%io_parent ) THEN
          WRITE(6,'(/1X,''GROUND STATE: '',I5,'' -'',I5)') 1,sh02%nst_s0
       ENDIF

       ! ...    read S1 state
       filn = s1_filn
       CALL zhrwf(1,irec,c0(:,sh02%nst_s0+1:,:),cm(1,sh02%nst_s0+1),&
            sh02%nst_s1,EIGV,TAU0,VELP,TAUI,NFI)
       IF ( paral%io_parent ) THEN
          WRITE(6,'(/1X,''EXCITED STATE: '',I5,'' -'',I5)')&
               sh02%nst_s0+1,sh02%nst_s0+sh02%nst_s1
       ENDIF

       ! ...    write combined wfn
       CALL write_irec(irec)
       filn = filn_back
       nstate = sh02%nst_s0 + sh02%nst_s1
       nfi = 1
       ! McB... cf. elct.inc
       crge%n = nstate
       CALL zhwwf(2,irec,c0,c0,nstate,eigv,tau0,velp,taui,nfi)

#if defined (__GROMOS)
       CALL mm_dim(mm_revert,oldstatus)
#endif
    ENDIF                     ! qmnode
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fusion
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE separate(c0,cm)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,crge%n,1), &
                                                cm(ncpw%ngw,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'separate'

    CHARACTER(len=20)                        :: filn_back
    INTEGER                                  :: ia, ie, ierr, irec(100), nfi, &
                                                nstate
    LOGICAL                                  :: oldstatus
    REAL(real_8)                             :: eigv(1)
    REAL(real_8), ALLOCATABLE                :: taui(:,:,:)

! Variables
! ==================================================================

    IF ( paral%qmnode ) THEN

#if defined (__GROMOS)
       CALL mm_dim(mm_go_mm,oldstatus)
#endif
       ALLOCATE(taui(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(taui)!,3*maxsys%nax*maxsys%nsx)

       filn_back = filbod
       CALL read_irec(irec)

       ! ...    read fused restart state
       nstate = sh02%nst_s0 + sh02%nst_s1
       crge%n = nstate
       nfi = 1
       CALL zhrwf(1,irec,c0,cm,nstate,eigv,tau0,velp,taui,nfi)
       CALL write_irec(irec)

       ! ...    write S0 
       CALL xstring(s0_filn,ia,ie)
       filbod = s0_filn(ia:ie)//'.'
       IF ( paral%io_parent ) THEN
          WRITE(6,'(/1X,''GROUND STATE: '',I5,'' -'',I5)') 1,sh02%nst_s0
       ENDIF
       ! McB... cf. elct.inc
       crge%n = sh02%nst_s0
       CALL zhwwf(2,irec,c0,cm,sh02%nst_s0,eigv,tau0,velp,taui,nfi)

       ! ...    write S1 
       CALL xstring(s1_filn,ia,ie)
       filbod = s1_filn(ia:ie)//'.'
       IF ( paral%io_parent ) THEN
          WRITE(6,'(/1X,''EXCITED STATE: '',I5,'' -'',I5)')&
               sh02%nst_s0+1,sh02%nst_s0+sh02%nst_s1
       ENDIF
       ! McB... cf. elct.inc
       crge%n = sh02%nst_s1
       CALL zhwwf(2,irec,c0(:,sh02%nst_s0+1:,:),cm(1,sh02%nst_s0+1),&
            sh02%nst_s1,EIGV,TAU0,VELP,TAUI,NFI)

       filbod = filn_back

#if defined (__GROMOS)
       CALL mm_dim(mm_revert,oldstatus)
#endif
    ENDIF                     ! qmnode

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE separate
  ! ==================================================================

END MODULE fusion_utils

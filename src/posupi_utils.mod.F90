MODULE posupi_utils
  USE cnst,                            ONLY: fbohr
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE jacobi_utils,                    ONLY: jacobi
  USE kinds,                           ONLY: real_8
  USE metr,                            ONLY: eps,&
                                             metr_com,&
                                             veps
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE readsr_utils,                    ONLY: readsi,&
                                             readsr
  USE store_types,                     ONLY: rout1
  USE system,                          ONLY: cnti,&
                                             iatpt,&
                                             maxsys
  USE tpar,                            ONLY: dt_ions
  USE zeroing_utils,                   ONLY: zeroing
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: posupi
  PUBLIC :: posupih_iso
  PUBLIC :: posupihshock
  PUBLIC :: posupih
  PUBLIC :: posupif
  PUBLIC :: posuphf
  PUBLIC :: sinhxbx

CONTAINS

  ! ==================================================================
  SUBROUTINE posupi(tau0,taup,velp)
    ! ==--------------------------------------------------------------==
    ! ==  UPDATE OF THE POSITIONS FOR VELOCITY VERLET                 ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), taup(:,:,:), &
                                                velp(:,:,:)

    INTEGER                                  :: ia, iat, is

#if defined(__VECTOR)
    !$omp parallel do private(IA,IS,IAT)
#else
    !$omp parallel do private(IA,IS,IAT) schedule(static)
#endif
    DO iat=1,ions1%nat
       ia=iatpt(1,iat)
       is=iatpt(2,iat)
       taup(1,ia,is)=tau0(1,ia,is)+dt_ions*velp(1,ia,is)
       taup(2,ia,is)=tau0(2,ia,is)+dt_ions*velp(2,ia,is)
       taup(3,ia,is)=tau0(3,ia,is)+dt_ions*velp(3,ia,is)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE posupi
  ! ==================================================================
  SUBROUTINE posupih_iso(tau0,taup,velp)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), taup(:,:,:), &
                                                velp(:,:,:)

    INTEGER                                  :: i, ia, iat, is, j
    REAL(real_8)                             :: arg, s1, s2, se

! Variables
! ==--------------------------------------------------------------==
! ==  UPDATE OF THE POSITIONS CONSTANT PRESSURE                   ==
! ==--------------------------------------------------------------==

    s1 = EXP(dt_ions*veps)
    arg = 0.5_real_8*dt_ions*veps
    se = EXP(arg)
    s2 = se*sinhxbx(arg)
#if defined (__VECTOR)
    !$omp parallel do private(IA,IS,IAT)
#else
    !$omp parallel do private(IA,IS,IAT) schedule(static)
#endif
    DO iat=1,ions1%nat
       ia=iatpt(1,iat)
       is=iatpt(2,iat)
       taup(1,ia,is)=tau0(1,ia,is)*s1+dt_ions*velp(1,ia,is)*s2
       taup(2,ia,is)=tau0(2,ia,is)*s1+dt_ions*velp(2,ia,is)*s2
       taup(3,ia,is)=tau0(3,ia,is)*s1+dt_ions*velp(3,ia,is)*s2
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  UPDATE OF THE BAROSTAT                                      ==
    ! ==--------------------------------------------------------------==
    eps = eps + dt_ions*veps
    ! ==--------------------------------------------------------------==
    ! ==  COMPUTE NEW CELL MATRIX                                     ==
    ! ==--------------------------------------------------------------==
    DO j=1,3
       DO i=1,3
          metr_com%ht(i,j) = metr_com%ht(i,j)*s1
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE posupih_iso
  ! ==================================================================
  SUBROUTINE posupihshock(tau0,taup,velp)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), taup(:,:,:), &
                                                velp(:,:,:)

    INTEGER                                  :: ia, is
    REAL(real_8)                             :: arg, s1, s2, se

! Variables
! ==--------------------------------------------------------------==
! ==  UPDATE OF THE POSITIONS CONSTANT PRESSURE                   ==
! ==--------------------------------------------------------------==

    s1 = EXP(dt_ions*veps)
    arg = 0.5_real_8*dt_ions*veps
    se = EXP(arg)
    s2 = se*sinhxbx(arg)
    !$omp   parallel do private(IA,IS) schedule(dynamic) 
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          taup(1,ia,is)=tau0(1,ia,is)*s1+dt_ions*velp(1,ia,is)*s2
          taup(2,ia,is)=tau0(2,ia,is)+dt_ions*velp(2,ia,is)
          taup(3,ia,is)=tau0(3,ia,is)+dt_ions*velp(3,ia,is)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  COMPUTE NEW CELL MATRIX                                     ==
    ! ==--------------------------------------------------------------==
    metr_com%ht(1,1) = metr_com%ht(1,1)*s1
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE posupihshock
  ! ==================================================================
  SUBROUTINE posupih(tau0,taup,velp)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:), taup(:,:,:), &
                                                velp(:,:,:)

    INTEGER                                  :: i, ia, ierr, is, j
    REAL(real_8) :: arg1, arg2, arg3, dth, e1, e2, e3, hdscr(3,3), &
      htmp1(3,3), hvd(3), hvev(3,3), rtmp1(3), s1, s2, s3, vtmp1(3)

! Variables
! ==--------------------------------------------------------------==
! == DIAGONALIZE HTVEL MATRIX                                     ==
! ==--------------------------------------------------------------==

    DO j=1,3
       DO i=1,3
          hdscr(i,j) = metr_com%htvel(i,j)
       ENDDO
    ENDDO
    CALL jacobi(3,3,hdscr,hvd,hvev,ierr)
    ! ==--------------------------------------------------------------==
    ! ==  UPDATE OF THE POSITIONS CONSTANT PRESSURE                   ==
    ! ==--------------------------------------------------------------==
    dth = 0.5_real_8*dt_ions
    arg1 = hvd(1)*dth
    arg2 = hvd(2)*dth
    arg3 = hvd(3)*dth
    e1 = EXP(hvd(1)*dt_ions)
    e2 = EXP(hvd(2)*dt_ions)
    e3 = EXP(hvd(3)*dt_ions)
    s1 = EXP(arg1)*sinhxbx(arg1)
    s2 = EXP(arg2)*sinhxbx(arg2)
    s3 = EXP(arg3)*sinhxbx(arg3)
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          rtmp1(1) = hvev(1,1)*tau0(1,ia,is) +&
               hvev(2,1)*tau0(2,ia,is) +&
               hvev(3,1)*tau0(3,ia,is)
          rtmp1(2) = hvev(1,2)*tau0(1,ia,is) +&
               hvev(2,2)*tau0(2,ia,is) +&
               hvev(3,2)*tau0(3,ia,is)
          rtmp1(3) = hvev(1,3)*tau0(1,ia,is) +&
               hvev(2,3)*tau0(2,ia,is) +&
               hvev(3,3)*tau0(3,ia,is)
          vtmp1(1) = hvev(1,1)*velp(1,ia,is) +&
               hvev(2,1)*velp(2,ia,is) +&
               hvev(3,1)*velp(3,ia,is)
          vtmp1(2) = hvev(1,2)*velp(1,ia,is) +&
               hvev(2,2)*velp(2,ia,is) +&
               hvev(3,2)*velp(3,ia,is)
          vtmp1(3) = hvev(1,3)*velp(1,ia,is) +&
               hvev(2,3)*velp(2,ia,is) +&
               hvev(3,3)*velp(3,ia,is)
          ! 
          rtmp1(1) = rtmp1(1)*e1 + dt_ions*vtmp1(1)*s1
          rtmp1(2) = rtmp1(2)*e2 + dt_ions*vtmp1(2)*s2
          rtmp1(3) = rtmp1(3)*e3 + dt_ions*vtmp1(3)*s3
          ! 
          taup(1,ia,is) = hvev(1,1)*rtmp1(1) +&
               hvev(1,2)*rtmp1(2) +&
               hvev(1,3)*rtmp1(3)
          taup(2,ia,is) = hvev(2,1)*rtmp1(1) +&
               hvev(2,2)*rtmp1(2) +&
               hvev(2,3)*rtmp1(3)
          taup(3,ia,is) = hvev(3,1)*rtmp1(1) +&
               hvev(3,2)*rtmp1(2) +&
               hvev(3,3)*rtmp1(3)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! ==  UPDATE OF THE CELL MATRIX                                   ==
    ! ==--------------------------------------------------------------==
    DO j=1,3
       htmp1(1,j) = 0._real_8
       htmp1(2,j) = 0._real_8
       htmp1(3,j) = 0._real_8
    ENDDO
    DO j=1,3
       DO i=1,3
          htmp1(i,j) = htmp1(i,j)&
               + hvev(1,i)*metr_com%ht(j,1)&
               + hvev(2,i)*metr_com%ht(j,2)&
               + hvev(3,i)*metr_com%ht(j,3)
       ENDDO
    ENDDO
    DO i=1,3
       htmp1(1,i) = htmp1(1,i)*e1
       htmp1(2,i) = htmp1(2,i)*e2
       htmp1(3,i) = htmp1(3,i)*e3
    ENDDO
    DO j=1,3
       metr_com%ht(1,j) = 0._real_8
       metr_com%ht(2,j) = 0._real_8
       metr_com%ht(3,j) = 0._real_8
    ENDDO
    DO j=1,3
       DO i=1,3
          metr_com%ht(i,j) = metr_com%ht(i,j)&
               + hvev(j,1)*htmp1(1,i)&
               + hvev(j,2)*htmp1(2,i)&
               + hvev(j,3)*htmp1(3,i)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE posupih
  ! ==================================================================
  FUNCTION sinhxbx(x)
    ! ==--------------------------------------------------------------==
    ! ==  COMPUTES SINH(X)/X BY A 10TH ORDER MACLAURIN SERIES USING   ==
    ! ==  HORNER'S METHOD (HORNY'S METHOD?)                           ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x, sinhxbx

    REAL(real_8)                             :: a10, a2, a4, a6, a8, x2

    a2  = 1._real_8/6._real_8
    a4 = 1._real_8/120._real_8
    a6  = 1._real_8/5040._real_8
    a8  = 1._real_8/362880._real_8
    a10 = 1._real_8/39916800._real_8
    x2 = x*x
    sinhxbx = ((((a10*x2 + a8)*x2 + a6)*x2 + a4)*x2 + a2)*x2 + 1._real_8
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION sinhxbx
  ! ==--------------------------------------------------------------==
  ! ==================================================================
  SUBROUTINE posupif(iunit,tau0,taup,velp,erread)
    ! ==--------------------------------------------------------------==
    ! ==  UPDATE OF THE POSITIONS FOR VELOCITY VERLET                 ==
    ! ==  form a previous TRAJECTORY file named TRAJSAVED             ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iunit
    REAL(real_8)                             :: tau0(:,:,:), taup(:,:,:), &
                                                velp(:,:,:)
    LOGICAL                                  :: erread

    CHARACTER(len=2)                         :: cato
    CHARACTER(len=200)                       :: line
    INTEGER                                  :: i, ia, iin, in, iout, is, nf
    INTEGER, SAVE                            :: ifirst = 1

    erread=.FALSE.
    IF (paral%io_parent) THEN
       IF (rout1%xtin) THEN
          IF (ifirst.EQ.1) THEN
             !            CALL fileopen(64,'TRAJSAVED.xyz',fo_old,erread)
             !            IF (erread) CALL stopgm('POSUPIF',' TRAJSAVED.xyz NOT FOUND',& 
             !                 __LINE__,__FILE__)
             ifirst=0
             DO in=1,cnti%nskip+cnti%nsample-1
                READ(iunit,'(A200)',err=99,END=99) line
                READ(iunit,'(A200)',err=99,END=99) line
                DO is=1,ions1%nsp
                   DO ia=1,ions0%na(is)
                      READ(iunit,'(A200)',err=99,END=99) line
                      IF (INDEX(line,'DATA').NE.0)&
                           READ(iunit,'(A80)',err=99,END=99) line
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO in=1,cnti%nsample-1
                READ(iunit,'(A200)',err=99,END=99) line
                READ(iunit,'(A200)',err=99,END=99) line
                DO is=1,ions1%nsp
                   DO ia=1,ions0%na(is)
                      READ(iunit,'(A200)',err=99,END=99) line
                      IF (INDEX(line,'DATA').NE.0)&
                           READ(iunit,'(A200)',err=99,END=99) line
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          READ(iunit,'(A200)',err=99,END=99) line
          READ(iunit,'(A200)',err=99,END=99) line
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                READ(iunit,*,err=99,END=99)&
                     cato,taup(1,ia,is),taup(2,ia,is),taup(3,ia,is)
             ENDDO
          ENDDO
          CALL dscal(3*ions1%nsp*maxsys%nax,fbohr,taup,1)
          CALL zeroing(velp)!,3*ions1%nsp*maxsys%nax)
       ELSE
          IF (ifirst.EQ.1) THEN
             !            CALL fileopen(64,'TRAJSAVED',fo_old,erread)
             !            IF (erread) CALL stopgm('POSUPIF',' TRAJSAVED NOT FOUND',& 
             !                 __LINE__,__FILE__)
             ifirst=0
             DO in=1,cnti%nskip+cnti%nsample-1
                DO is=1,ions1%nsp
                   DO ia=1,ions0%na(is)
                      READ(iunit,'(A200)',err=99,END=99) line
                      IF (INDEX(line,'DATA').NE.0)&
                           READ(iunit,'(A200)',err=99,END=99) line
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO in=1,cnti%nsample-1
                DO is=1,ions1%nsp
                   DO ia=1,ions0%na(is)
                      READ(iunit,'(A200)',err=99,END=99) line
                      IF (INDEX(line,'DATA').NE.0)&
                           READ(iunit,'(A80)',err=99,END=99) line
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                READ(iunit,'(A200)',err=99,END=99) line
                IF (INDEX(line,'DATA').NE.0)&
                     READ(iunit,'(A200)',err=99,END=99) line
                iin=1
                CALL readsi(line,iin,iout,nf,erread)
                IF (erread) GOTO 99
                DO i=1,3
                   iin=iout+1
                   CALL readsr(line,iin,iout,taup(i,ia,is),erread)
                ENDDO
                IF (erread) GOTO 99
                DO i=1,3
                   iin=iout+1
                   CALL readsr(line,iin,iout,velp(i,ia,is),erread)
                ENDDO
                IF (erread) GOTO 99
             ENDDO
          ENDDO
       ENDIF
    ENDIF ! if (io_parent)

    GOTO 100

99  CONTINUE
    erread=.TRUE.
100 CONTINUE

    ! bcast data that each cp_groups get the data
    CALL mp_bcast(tau0,3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
    CALL mp_bcast(taup,3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
    CALL mp_bcast(velp,3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
    CALL mp_bcast(erread,parai%io_source,parai%cp_grp)

    RETURN
  END SUBROUTINE posupif
  ! ==================================================================
  SUBROUTINE posuphf(iunit,ht,htvel,erread)
    ! ==--------------------------------------------------------------==
    ! ==  UPDATE OF THE CELL VECTORS FROM A PREVIOUS TRAJECTORY NAMED ==
    ! ==  CELLSAVED READ FROM UNIT OF IUNIT                           ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iunit
    REAL(real_8)                             :: ht(3,3), htvel(3,3)
    LOGICAL                                  :: erread

    CHARACTER(len=200)                       :: line
    INTEGER                                  :: i, iin, iout, j, n
    INTEGER, SAVE                            :: ifirst = 1

! Variables
! ==--------------------------------------------------------------==

    erread=.FALSE.
    IF (paral%io_parent) THEN
       IF (ifirst==1) THEN
          ifirst=0
          DO n=1,cnti%nskip+cnti%nsample-1
             READ(iunit,'(A200)',err=99,END=99) line
             DO i=1,3
                READ(iunit,'(A200)',err=99,END=99) line
                IF (INDEX(line,'DATA').NE.0)&
                     READ(iunit,'(A200)',err=99,END=99) line
             ENDDO
          ENDDO
       ELSE
          DO n=1,cnti%nsample-1
             READ(iunit,'(A200)',err=99,END=99) line
             DO i=1,3
                READ(iunit,'(A200)',err=99,END=99) line
                IF (INDEX(line,'DATA').NE.0)&
                     READ(iunit,'(A200)',err=99,END=99) line
             ENDDO
          ENDDO
       ENDIF
       READ(iunit,'(A200)',err=99,END=99) line
       DO i=1,3
          READ(iunit,'(A200)',err=99,END=99) line
          IF(INDEX(line,'DATA').NE.0)&
               READ(iunit,'(A200)',err=99,END=99) line
          iin=1
          DO j=1,3
             CALL readsr(line,iin,iout,ht(i,j),erread)
             iin=iout+1
          ENDDO
          IF (erread) GOTO 99
          DO j=1,3
             CALL readsr(line,iin,iout,htvel(i,j),erread)
             iin=iout+1
          ENDDO
          IF (erread) GOTO 99
       ENDDO
    ENDIF
    GOTO 100
99  CONTINUE
    erread=.TRUE.
100 CONTINUE
    ! bcast data that each cp_groups get the data
    CALL mp_bcast(ht,SIZE(ht),parai%io_source,parai%cp_grp)
    CALL mp_bcast(htvel,SIZE(htvel),parai%io_source,parai%cp_grp)
    CALL mp_bcast(erread,parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE posuphf
  ! ==================================================================

END MODULE posupi_utils

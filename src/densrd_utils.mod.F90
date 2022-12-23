MODULE densrd_utils
  USE error_handling,                  ONLY: stopgm
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_old,&
                                             fo_ufo
  USE gvec,                            ONLY: gvec_com
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_4,&
                                             real_8
  USE mp_interface,                    ONLY: mp_recv,&
                                             mp_send,&
                                             mp_sync
  USE parac,                           ONLY: parai,&
                                             paral
  USE response_pmod,                   ONLY: response1
  USE system,                          ONLY: mapgp,&
                                             ncpw,&
                                             parap,&
                                             spar

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: densrd

CONTAINS

  ! ==================================================================
  SUBROUTINE densrd(rg,filin)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: rg(ncpw%nhg)
    CHARACTER(len=*)                         :: filin

    INTEGER                                  :: i, ig, len, nhgsr, nr1sr, &
                                                nr2sr, nr3sr
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: gcutr, gcutwr

! Variables

#if defined(__NEC)
    COMPLEX(real_8), ALLOCATABLE :: rgsx(:)
    COMPLEX(real_8), ALLOCATABLE :: rg2(:)

#else
    COMPLEX(real_4), ALLOCATABLE :: rgsx(:)
    COMPLEX(real_4), ALLOCATABLE :: rg2(:)
#endif

    INTEGER, ALLOCATABLE :: mapw(:)

    REAL(real_8) :: x1
    INTEGER :: n1,nhgx,ipp,ip,msgid,msglen,nbytes
    INTEGER :: ierr
    CHARACTER(*),PARAMETER :: procedureN='densrd'
    ! ==--------------------------------------------------------------==
    ! ==  READS THE DENSITY FROM FILE                                 ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       len=2*spar%nhgs/2
    ELSE
       len=2*ncpw%nhg/2
    ENDIF
    ALLOCATE(rgsx(len),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (paral%parent) THEN
       ferror=.FALSE.
       IF (paral%io_parent)CALL fileopen(12,filin,fo_old+fo_ufo,ferror)
       IF (ferror) GOTO 30
       IF (paral%io_parent) THEN
          READ(12,END=20,err=30)
          READ(12,END=20,err=30)
          READ(12,END=20,err=30) nr1sr,nr2sr,nr3sr
          READ(12,END=20,err=30) gcutr,nhgsr,gcutwr
       ENDIF
       IF (nr1sr.NE.spar%nr1s) CALL stopgm('DENSRD','NUMBER OF PLANE WAVES',& 
            __LINE__,__FILE__)
       IF (nr2sr.NE.spar%nr2s) CALL stopgm('DENSRD','NUMBER OF PLANE WAVES',& 
            __LINE__,__FILE__)
       IF (nr3sr.NE.spar%nr3s) CALL stopgm('DENSRD','NUMBER OF PLANE WAVES',& 
            __LINE__,__FILE__)
       IF (nhgsr.NE.spar%nhgs) CALL stopgm('DENSRD','NUMBER OF PLANE WAVES',& 
            __LINE__,__FILE__)
       IF ((ABS(gcutr-gvec_com%gcut)+ABS(gcutwr-gvec_com%gcutw)).GT.1.e-12_real_8)&
            CALL stopgm(' DENSRD','PLANE WAVES CUTOFFS',& 
            __LINE__,__FILE__)
       IF (paral%io_parent) THEN
          READ(12,END=20,err=30)
          READ(12,END=20,err=30)
       ENDIF
       ! For Anatoles reference density: the number (!) of atoms in the
       ! reference molecule may be different:
       IF (paral%io_parent.AND..NOT.response1%tdummyatom) THEN
          DO i=1,ions1%nat
             READ(12,END=20,err=30)
          ENDDO
       ENDIF
       IF (paral%io_parent)&
            READ(12,END=20,err=30) (rgsx(ig),ig=1,spar%nhgs)
       IF (paral%io_parent)&
            CALL fileclose(12)
    ENDIF

    IF (paral%parent) THEN
       x1=REAL(spar%nhgs,kind=real_8)/REAL(parai%nproc,kind=real_8)*1.4_real_8
       n1=MIN(spar%nhgs,NINT(x1))*2/2
       ALLOCATE(rg2(n1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       nhgx=NINT(REAL(ncpw%nhg,kind=real_8)*1.4_real_8)
       nhgx=MIN(spar%nhgs,nhgx)
       nhgx=nhgx + 1
       ALLOCATE(mapw(nhgx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DO ipp=1,parai%nproc
       ip=parap%pgroup(ipp)
       msgid=ip
       IF (paral%parent) THEN
          IF (ip.EQ.parai%me) THEN
             DO ig=1,ncpw%nhg
                rg(ig)=rgsx(mapgp(ig))
             ENDDO
          ELSE
             msgid=1
             !msglen = parap%sparm(1,ipp-1) * 8/irat
             CALL mp_recv(mapw,parap%sparm(1,ipp-1),ip,msgid,parai%allgrp)
             DO ig=1,parap%sparm(1,ipp-1)
                rg2(ig)=rgsx(mapw(ig))
             ENDDO
             msgid=2
             !msglen = parap%sparm(1,ipp-1) * 8
             CALL mp_send(rg2,parap%sparm(1,ipp-1),ip,msgid,parai%allgrp)
          ENDIF
       ELSE
          IF (ip.EQ.parai%me) THEN
             msgid=1
             !msglen = nhg * 8/irat
             CALL mp_send(mapgp,ncpw%nhg,parap%pgroup(1),msgid,parai%allgrp)
             msgid=2
             !msglen = nhg * 8
             CALL mp_recv(rgsx,ncpw%nhg,parap%pgroup(1),msgid,parai%allgrp)
             DO ig=1,ncpw%nhg
                rg(ig)=rgsx(ig)
             ENDDO
          ENDIF
       ENDIF
       CALL mp_sync(parai%allgrp)
    ENDDO


    IF (paral%parent) THEN
       IF (paral%io_parent) WRITE(6,*) ' WAVEFUNCTION READ FROM FILE ',filin
       DEALLOCATE(rg2,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(mapw,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    DEALLOCATE(rgsx,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
20  CONTINUE
    CALL stopgm('DENSRD','END OF FILE',& 
         __LINE__,__FILE__)
30  CONTINUE
    CALL stopgm('DENSRD','READ ERROR',& 
         __LINE__,__FILE__)
  END SUBROUTINE densrd
  ! ==================================================================

END MODULE densrd_utils

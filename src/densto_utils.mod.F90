MODULE densto_utils
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen

  IMPLICIT NONE

  PRIVATE

!!$public :: densto
!!$
!!$contains


END MODULE densto_utils


#ifdef  BIG_STO
! ==================================================================
SUBROUTINE densto(rg,tau0,filin)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_send, mp_recv
  USE parac, ONLY : paral,parai
  USE fileopen_utils, ONLY : fileopen
  USE fileopen_utils, ONLY : fileclose
  USE system, ONLY: parap
  IMPLICIT NONE
  COMPLEX(real_8)                            :: rg(nhg)
  REAL(real_8)                               :: tau0(3,maxsys%nax,maxsys%nsx)
  CHARACTER(len=*)                           :: filin

  INTEGER                                    :: ia, ig, is, len
  LOGICAL                                    :: ferror

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
  INTEGER :: n1,nhgx,ipp,ip,msgid,msglen,nbytes,ipw,ierr,isub
  CHARACTER(*),PARAMETER::procedureN='densto'

  CALL tiset(procedureN,isub)

  ! ==--------------------------------------------------------------==
  ! ==  STORE THE DENSITY ON FILE                                   ==
  ! ==--------------------------------------------------------------==
  IF (paral%parent) THEN
     len=2*nhgs/2
  ELSE
     len=2*nhg/2
  ENDIF
  IF (.NOT.paral%parent) THEN
     ALLOCATE(rgsx(len),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
  ENDIF
  IF (paral%parent) THEN
     x1=REAL(nhgs,kind=real_8)/REAL(parai%nproc,kind=real_8)*1.4_real_8
     n1=MIN(nhgs,NINT(x1))
     ALLOCATE(rg2(n1),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)
     nhgx=NINT(REAL(nhg,kind=real_8)*1.4_real_8)
     nhgx=MIN(nhgs,nhgx)
     nhgx=nhgx + 1
     ALLOCATE(mapw(nhgx),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
          __LINE__,__FILE__)

     IF (paral%io_parent) THEN
        CALL fileopen(12,filin,fo_def+fo_ufo,ferror)
        CALL fileopen(15,"MAP",fo_def+fo_ufo,ferror)
        IF (ferror) GOTO 30

        WRITE(12,err=30) parm%ibrav
        WRITE(12,err=30) cell_com%celldm
        WRITE(12,err=30) nr1s,nr2s,nr3s
        WRITE(12,err=30) gvec_com%gcut,nhgs,gvec_com%gcutw
        WRITE(12,err=30) ions1%nsp
        WRITE(12,err=30) (ions0%na(is),is=1,ions1%nsp)
     ENDIF
     ! For Anatoles reference density: the number of atoms in the
     ! reference molecule may be different
     IF (.NOT.response1%tdummyatom_ref) THEN
        DO is=1,ions1%nsp
           DO ia=1,ions0%na(is)
              IF (paral%io_parent)&
                   WRITE(12,err=30) tau0(1,ia,is)-clsaabox%mm_c_trans(1),&
                   tau0(2,ia,is)-clsaabox%mm_c_trans(2),&
                   tau0(3,ia,is)-clsaabox%mm_c_trans(3),ions0%iatyp(is)
           ENDDO
        ENDDO
     ENDIF

  ENDIF
  DO ipp=1,parai%nproc
     ip=parap%pgroup(ipp)
     msgid = ip
     IF (paral%parent) THEN
        IF (ip.EQ.parai%me) THEN
           DO ig=1,nhg
              rg2(ig)=rg(ig)
              mapw(ig)=mapgp(ig)
           ENDDO
        ELSE
           msgid=1
           !msglen = parap%sparm(1,ipp-1) * 8
           CALL mp_recv(rg2,parap%sparm(1,ipp-1),ip,msgid,parai%allgrp)
           msgid=2
           !msglen = parap%sparm(1,ipp-1) * 8/irat
           CALL mp_recv(mapw,parap%sparm(1,ipp-1),ip,msgid,parai%allgrp)
        ENDIF
        IF (paral%io_parent) THEN
           DO ipw=1,parap%sparm(1,ipp-1)
              WRITE(12,err=30) rg2(ipw)
              WRITE(15,err=30) mapw(ipw)
           ENDDO
        ENDIF
     ELSE
        IF (ip.EQ.parai%me) THEN
           DO ig=1,nhg
              rgsx(ig)=rg(ig)
           ENDDO
           msgid=1
           !msglen = nhg * 8
           CALL mp_send(rgsx,nhg,parap%pgroup(1),msgid,parai%allgrp)
           msgid=2
           !msglen = nhg * 8/irat
           CALL mp_send(mapgp,nhg,parap%pgroup(1),msgid,parai%allgrp)
        ENDIF
     ENDIF
  ENDDO
  IF (paral%io_parent) THEN
     WRITE(12,err=30) (clsaabox%mm_c_trans(ig),ig=1,3)! origin of density block
     CALL fileclose(12)
     CALL fileclose(15)
     WRITE(6,'(A,A)') ' DENSITY WRITTEN TO FILE ',filin
     DEALLOCATE(rg2,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
     DEALLOCATE(mapw,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
  ENDIF
  IF (.NOT.paral%parent) THEN
     DEALLOCATE(rgsx,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
          __LINE__,__FILE__)
  ENDIF
  CALL tihalt(procedureN,isub)
  RETURN
  ! ==--------------------------------------------------------------==
30 CONTINUE
  CALL stopgm('DENSTO','WRITE ERROR',& 
       __LINE__,__FILE__)
END SUBROUTINE densto
! ==================================================================
#else
! ==================================================================
SUBROUTINE densto(rg,tau0,filin)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_send, mp_recv
  USE system , ONLY:mapgp,maxsys,ncpw,parm,spar,parap
  USE parac, ONLY : paral,parai
  USE cell , ONLY:cell_com
  USE gvec , ONLY:gvec_com
  USE ions , ONLY:ions0,ions1
  USE fileopenmod , ONLY:fo_def,fo_ufo
  USE fileopen_utils, ONLY : fileopen
  USE fileopen_utils, ONLY : fileclose
  USE response_pmod , ONLY:response1
  USE mm_dimmod , ONLY:clsaabox
  IMPLICIT NONE
  COMPLEX(real_8)                            :: rg(ncpw%nhg)
  REAL(real_8)                               :: tau0(3,maxsys%nax,maxsys%nsx)
  CHARACTER(len=*)                           :: filin

  CHARACTER(*), PARAMETER                    :: procedureN = 'densto'

  INTEGER                                    :: ia, ierr, ig, ip, ipp, ipw, &
                                                is, isub, len, msgid, n1, nhgx
  INTEGER, ALLOCATABLE                       :: mapw(:)
  LOGICAL                                    :: ferror
  REAL(real_8)                               :: x1

#if defined(__NEC)
  COMPLEX(real_8), ALLOCATABLE :: rgsx(:)
  COMPLEX(real_8), ALLOCATABLE :: rg2(:)
#else
  COMPLEX(real_4), ALLOCATABLE :: rgsx(:)
  COMPLEX(real_4), ALLOCATABLE :: rg2(:)
#endif
  CALL tiset(procedureN,isub)
  ! ==--------------------------------------------------------------==
  ! ==  STORE THE DENSITY ON FILE                                   ==
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
     x1=REAL(spar%nhgs,kind=real_8)/REAL(parai%nproc,kind=real_8)*1.4_real_8
     n1=MIN(spar%nhgs,NINT(x1))
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
     msgid = ip
     IF (paral%parent) THEN
        IF (ip.EQ.parai%me) THEN
           DO ig=1,ncpw%nhg
              rg2(ig)=rg(ig)
              mapw(ig)=mapgp(ig)
           ENDDO
        ELSE
           msgid=1
           !msglen = parap%sparm(1,ipp-1) * 8
           CALL mp_recv(rg2,parap%sparm(1,ipp-1),ip,msgid,parai%allgrp)
           msgid=2
           !msglen = parap%sparm(1,ipp-1) * 8/irat
           CALL mp_recv(mapw,parap%sparm(1,ipp-1),ip,msgid,parai%allgrp)
        ENDIF
        DO ipw=1,parap%sparm(1,ipp-1)
           rgsx(mapw(ipw))=rg2(ipw)
        ENDDO
     ELSE
        IF (ip.EQ.parai%me) THEN
           DO ig=1,ncpw%nhg
              rgsx(ig)=rg(ig)
           ENDDO
           msgid=1
           !msglen = nhg * 8
           CALL mp_send(rgsx,ncpw%nhg,parap%pgroup(1),msgid,parai%allgrp)
           msgid=2
           !msglen = nhg * 8/irat
           CALL mp_send(mapgp,ncpw%nhg,parap%pgroup(1),msgid,parai%allgrp)
        ENDIF
     ENDIF
  ENDDO
  IF (paral%parent.AND.paral%io_parent) THEN
     CALL fileopen(12,filin,fo_def+fo_ufo,ferror)
     IF (ferror) GOTO 30
     WRITE(12,err=30) parm%ibrav
     WRITE(12,err=30) cell_com%celldm
     WRITE(12,err=30) spar%nr1s,spar%nr2s,spar%nr3s
     WRITE(12,err=30) gvec_com%gcut,spar%nhgs,gvec_com%gcutw
     WRITE(12,err=30) ions1%nsp
     WRITE(12,err=30) (ions0%na(is),is=1,ions1%nsp)

     ! For Anatoles reference density: the number of atoms in the
     ! reference molecule may be different
     IF (.NOT.response1%tdummyatom_ref) THEN
        DO is=1,ions1%nsp
           DO ia=1,ions0%na(is)
              WRITE(12,err=30) tau0(1,ia,is)-clsaabox%mm_c_trans(1),&
                   tau0(2,ia,is)-clsaabox%mm_c_trans(2),&
                   tau0(3,ia,is)-clsaabox%mm_c_trans(3),ions0%iatyp(is)
           ENDDO
        ENDDO
     ENDIF

     WRITE(12,err=30) (rgsx(ig),ig=1,spar%nhgs)
     WRITE(12,err=30) (clsaabox%mm_c_trans(ig),ig=1,3)! origin of density block
     CALL fileclose(12)
     WRITE(6,'(A,A)') ' DENSITY WRITTEN TO FILE ',filin
  ENDIF
  DEALLOCATE(rgsx,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
       __LINE__,__FILE__)
  CALL tihalt(procedureN,isub)
  RETURN
  ! ==--------------------------------------------------------------==
30 CONTINUE
  CALL stopgm('DENSTO','WRITE ERROR',& 
       __LINE__,__FILE__)
END SUBROUTINE densto
! ==================================================================
#endif

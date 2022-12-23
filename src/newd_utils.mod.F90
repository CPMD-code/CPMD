MODULE newd_utils
  USE cppt,                            ONLY: gk,&
                                             inyh
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: nlps_com
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE qvan2_utils,                     ONLY: qvan2
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE system,                          ONLY: cntl,&
                                             iatpe,&
                                             ipept,&
                                             maxsys,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: newd
  !public :: newd_worker
  !public :: cftemp
  PUBLIC :: give_scr_newd

CONTAINS

  ! ==================================================================
  ! forwarding subroutine, to be able to use cray pointers with OpenMP
  ! the real stuff is in NEWD_WORKER.
  ! ==================================================================
  SUBROUTINE newd(fnla,deeq,f,vpot,fion,nstate,tfor)
    REAL(real_8)                             :: fnla(:,:,:), &
                                                deeq(ions1%nat,maxsys%nhxs,*)
    COMPLEX(real_8)                          :: vpot(ncpw%nhg)
    REAL(real_8)                             :: fion(:,:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate)
    LOGICAL                                  :: tfor

    CHARACTER(*), PARAMETER                  :: procedureN = 'newd'

    COMPLEX(real_8), ALLOCATABLE             :: qg(:,:), vtmp(:,:)
    INTEGER                                  :: ierr, il_qg, il_vtmp, il_ylm, &
                                                isub, nhh
    REAL(real_8), ALLOCATABLE                :: ylm(:,:)

! cray pointers
! ==--------------------------------------------------------------==

    CALL tiset('      NEWD',isub)
    ! SCR partition
    IF (cntl%bigmem) THEN
       nhh = (maxsys%nhxs*(maxsys%nhxs+1))/2
       il_vtmp=ncpw%nhg*6
       il_qg  =nhh*ncpw%nhg*2
       il_ylm =maxsys%nax*(4+nhh)
    ELSE
       il_vtmp=ncpw%nhg*2
       il_qg  =ncpw%nhg*2
       il_ylm =maxsys%nax
    ENDIF

    ! allocate as 1D, because dims are not important in this subroutine
    ALLOCATE(vtmp(il_vtmp,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(qg(il_qg,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ylm(il_ylm,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL newd_worker(fnla,deeq,f,vpot,fion,nstate,tfor,vtmp,qg,ylm)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(vtmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(qg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(ylm,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('      NEWD',isub)
    RETURN
  END SUBROUTINE newd
  ! ==================================================================
  SUBROUTINE newd_worker(fnla,deeq,f,vpot,fion,nstate,tfor,vtmp,qg,&
       ylm)
    ! ==--------------------------------------------------------------==
    ! ==         THIS ROUTINE CALCULATES THE ARRAY DEEQ               ==
    ! ==                                                              ==
    ! ==  DEEQ_I,IJ = OMEGA( V(G=0) Q_I,IJ(G=0) +                     ==
    ! ==       2 SUM_G> RE[ V*(G) Q_I,IJ(G) E^-IG.R_I ] )             ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fnla(:,:,:), deeq(ions1%nat,&
                                                maxsys%nhxs,*), f(*)
    COMPLEX(real_8)                          :: vpot(*)
    REAL(real_8)                             :: fion(:,:,:)
    INTEGER                                  :: nstate
    LOGICAL                                  :: tfor
    COMPLEX(real_8)                          :: vtmp(ncpw%nhg,*), &
                                                qg(ncpw%nhg,*)
    REAL(real_8)                             :: ylm(maxsys%nax,*)

    INTEGER                                  :: i, ia, ig, ii, ijv, is, isa, &
                                                isa0, iv, jv, k, nhh
    REAL(real_8)                             :: ftmp, otr, rhovan, tom
    REAL(real_8), EXTERNAL                   :: dotp

    IF (cntl%bigmem) THEN
       isa0=0
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is)) THEN
             ijv=0
             DO iv=1,nlps_com%ngh(is)
                DO jv=iv,nlps_com%ngh(is)
                   ijv=ijv+1
                   CALL qvan2(iv,jv,is,qg(1,ijv))
#if defined(__VECTOR)
                   !$omp parallel do private(IG)
#else
                   !$omp parallel do private(IG) schedule(static)
#endif
                   DO ig=1,ncpw%nhg
                      qg(ig,ijv)=CONJG(qg(ig,ijv))*vpot(ig)
                   ENDDO
                ENDDO
             ENDDO
             nhh=(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2
             tom=2.0_real_8*parm%omega
             CALL dgemm('T','N',ions0%na(is),nhh,2*ncpw%nhg,tom,eigrb(1,isa0+1),2*&
                  ncpw%nhg,qg(1,1),2*ncpw%nhg,0.0_real_8,ylm(1,5),maxsys%nax)
             ijv=0
             DO iv=1,nlps_com%ngh(is)
                DO jv=iv,nlps_com%ngh(is)
                   ijv=ijv+1
                   ! TOM=2.0_real_8*OMEGA
                   ! CALL DGEMV('T',2*NHG,NA(IS),TOM,EIGRB(1,ISA0+1),2*NHG,
                   ! *                     QG(1,IJV),1,0.0_real_8,DEEQ(ISA0+1,IV,JV),1)
                   !$omp parallel do private(IA,ISA) shared(DEEQ)
                   DO ia=1,ions0%na(is)
                      isa=isa0+ia
                      deeq(isa,iv,jv)=ylm(ia,4+ijv)
                      IF (geq0) deeq(isa,iv,jv)=deeq(isa,iv,jv)-parm%omega*REAL(&
                           qg(1,ijv)*eigrb(1,isa))
                      deeq(isa,jv,iv)=deeq(isa,iv,jv)
                   ENDDO
                ENDDO
             ENDDO
             IF (tfor) THEN
                ijv=0
                DO iv=1,nlps_com%ngh(is)
                   DO jv=iv,nlps_com%ngh(is)
                      ijv=ijv+1
                      CALL zeroing(ylm(:,1))!,maxsys%nax)
                      IF (cntl%tfdist) THEN
                         !$omp parallel do private(IA,ISA,II,I) shared(YLM)
                         DO ia=1,ions0%na(is)
                            isa=isa0+ia
                            IF (iatpe(isa).EQ.parai%mepos) THEN
                               ii=isa-ipept(1,parai%mepos)+1
                               DO i=1,nstate
                                  ylm(ia,1)=ylm(ia,1)+f(i)*fnla(ii,iv,i)*&
                                       fnla(ii,jv,i)
                               ENDDO
                            ENDIF
                         ENDDO
                         CALL mp_sum(ylm,ions0%na(is),parai%allgrp)
                      ELSE
                         !$omp parallel do private(IA,ISA,II,I) shared(YLM)
                         DO ia=1,ions0%na(is)
                            isa=isa0+ia
                            DO i=1,nstate
                               ylm(ia,1)=ylm(ia,1)+f(i)*fnla(isa,iv,i)*&
                                    fnla(isa,jv,i)
                            ENDDO
                         ENDDO
                      ENDIF
#if defined(__VECTOR)
                      !$omp parallel do private(IG) shared(QG)
#else
                      !$omp parallel do private(IG) shared(QG) schedule(static)
#endif
                      DO ig=1,ncpw%nhg
                         vtmp(ig,1)=qg(ig,ijv)*CMPLX(0.0_real_8,-gk(1,ig),kind=real_8)
                         vtmp(ig,2)=qg(ig,ijv)*CMPLX(0.0_real_8,-gk(2,ig),kind=real_8)
                         vtmp(ig,3)=qg(ig,ijv)*CMPLX(0.0_real_8,-gk(3,ig),kind=real_8)
                      ENDDO
                      IF (ions0%na(is).GT.2) THEN
                         isa=isa0+1
                         CALL dgemm('T','N',ions0%na(is),3,2*ncpw%nhg,2._real_8,eigrb(1,isa),&
                              2*ncpw%nhg,vtmp,2*ncpw%nhg,0.0_real_8,ylm(1,2),maxsys%nax)
                         IF (geq0)CALL dger(ions0%na(is),3,-1.0_real_8,eigrb(1,isa),2*&
                              ncpw%nhg,vtmp,2*ncpw%nhg,ylm(1,2),maxsys%nax)
                      ELSE
                         DO ia=1,ions0%na(is)
                            isa=isa0+ia
                            ylm(ia,2)=dotp(ncpw%nhg,vtmp(1,1),eigrb(1,isa))
                            ylm(ia,3)=dotp(ncpw%nhg,vtmp(1,2),eigrb(1,isa))
                            ylm(ia,4)=dotp(ncpw%nhg,vtmp(1,3),eigrb(1,isa))
                         ENDDO
                      ENDIF
                      ftmp=1.0_real_8
                      IF (iv.NE.jv) ftmp=2.0_real_8
                      DO ia=1,ions0%na(is)
                         otr=ftmp*parm%omega*parm%tpiba*ylm(ia,1)
                         fion(1,ia,is)=fion(1,ia,is)+otr*ylm(ia,2)
                         fion(2,ia,is)=fion(2,ia,is)+otr*ylm(ia,3)
                         fion(3,ia,is)=fion(3,ia,is)+otr*ylm(ia,4)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
          isa0=isa0+ions0%na(is)
       ENDDO
    ELSE
       ! ==--------------------------------------------------------------==
       isa0=0
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is)) THEN
             ijv=0
             DO iv=1,nlps_com%ngh(is)
                DO jv=iv,nlps_com%ngh(is)
                   ijv=ijv+1
                   CALL qvan2(iv,jv,is,qg)
#if defined(__VECTOR)
                   !$omp parallel do private(IG) shared(QG)
#else
                   !$omp parallel do private(IG) shared(QG) schedule(static)
#endif
                   DO ig=1,ncpw%nhg
                      qg(ig,1)=CONJG(qg(ig,1))*vpot(ig)
                   ENDDO
                   DO ia=1,ions0%na(is)
                      isa=isa0+ia
#if defined(__VECTOR)
                      !$omp parallel do private(IG) shared(VTMP)
#else
                      !$omp parallel do private(IG) shared(VTMP) schedule(static)
#endif
                      DO ig=1,ncpw%nhg
                         vtmp(ig,1)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                              ei3(isa,inyh(3,ig))
                      ENDDO
                      deeq(isa,iv,jv)=parm%omega*dotp(ncpw%nhg,qg(1,1),vtmp(1,1))
                      deeq(isa,jv,iv)=deeq(isa,iv,jv)
                      IF (tfor) THEN
                         rhovan=0.0_real_8
                         IF (cntl%tfdist) THEN
                            IF (iatpe(isa).EQ.parai%mepos) THEN
                               ii=isa-ipept(1,parai%mepos)+1
                               !$omp parallel do private(I) reduction(+:RHOVAN)
                               DO i=1,nstate
                                  rhovan=rhovan+f(i)*fnla(ii,iv,i)*fnla(ii,jv,i)
                               ENDDO
                            ENDIF
                            CALL mp_sum(rhovan,parai%allgrp)
                         ELSE
                            !$omp parallel do private(I) reduction(+:RHOVAN)
                            DO i=1,nstate
                               rhovan=rhovan+f(i)*fnla(isa,iv,i)*fnla(isa,jv,i)
                            ENDDO
                         ENDIF
                         otr=parm%omega*parm%tpiba*rhovan
                         DO k=1,3
                            CALL cftemp(ncpw%nhg,qg(1,1),vtmp(1,1),gk,k,ftmp)
                            IF (iv.NE.jv) ftmp=2.0_real_8*ftmp
                            fion(k,ia,is)=fion(k,ia,is)+ftmp*otr
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          isa0=isa0+ions0%na(is)
       ENDDO
    ENDIF
    CALL mp_sum(deeq,ions1%nat*maxsys%nhxs*maxsys%nhxs,parai%allgrp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE newd_worker
  ! ==================================================================
  SUBROUTINE give_scr_newd(lnewd,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lnewd
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: nhh

    IF (cntl%bigmem) THEN
       nhh = (maxsys%nhxs*(maxsys%nhxs+1))/2
       lnewd = 2*ncpw%nhg*nhh + 6*ncpw%nhg + maxsys%nax*(4+nhh) + 100
    ELSE
       lnewd = 2*ncpw%nhg + 6*ncpw%nhg + maxsys%nax + 100
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_newd
  ! ==================================================================


END MODULE newd_utils

! ==================================================================
SUBROUTINE cftemp(nhg,qg,vtmp,gk,k,ftmp)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: nhg
  REAL(real_8)                               :: qg(2,nhg), vtmp(2,nhg), &
                                                gk(3,nhg)
  INTEGER                                    :: k
  REAL(real_8)                               :: ftmp

  INTEGER                                    :: ig

! Variables
! ==--------------------------------------------------------------==
! G=0 TERM IS ZERO !

  ftmp=0.0_real_8
#if defined(__VECTOR)
  !$omp parallel do private(IG) reduction(+:FTMP)
#else
  !$omp parallel do private(IG) reduction(+:FTMP) schedule(static)
#endif
  DO ig=1,nhg
     ftmp=ftmp+gk(k,ig)*(vtmp(1,ig)*qg(2,ig)-vtmp(2,ig)*qg(1,ig))
  ENDDO
  ftmp=2.0_real_8*ftmp
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE cftemp

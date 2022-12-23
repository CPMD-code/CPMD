#include "cpmd_global.h"

MODULE opeigr_utils
  USE ddip,                            ONLY: lenbk,&
                                             ngg1,&
                                             ngg2,&
                                             ngg3,&
                                             ngwmax
  USE error_handling,                  ONLY: stopgm
  USE gvec,                            ONLY: gvec_com
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  USE numpw_utils,                     ONLY: gkpwf
  USE parac,                           ONLY: parai
  USE sphe,                            ONLY: gcutka,&
                                             gcutwmax,&
                                             gcutwmin,&
                                             tsphere
  USE spin,                            ONLY: lspin2,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             parap,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: determ
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: opeigr
  PUBLIC :: opapp
  PUBLIC :: give_scr_opeigr
  !public :: getngg

CONTAINS

  ! ==================================================================
  SUBROUTINE opeigr(c0,c2,sc0,nstate,mapful,mapcol,ddmat,&
       nop1,nop2,dd)
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    ! Arguments
    INTEGER :: nstate,mapful(*),mapcol(*),nop1,nop2
    COMPLEX(real_8) :: c0(*),c2(*),sc0(*),ddmat(nstate*nstate),dd
    ! Variables
    COMPLEX(real_8) :: dds
    COMPLEX(real_8), ALLOCATABLE :: aux(:,:)
    INTEGER :: isub,ierr
    CHARACTER(*),PARAMETER::procedureN='OPEIGR'
    ! ==--------------------------------------------------------------==
    CALL tiset('    OPEIGR',isub)
    CALL zeroing(ddmat)!,nstate*nstate)
    IF (lenbk.LE.0) CALL stopgm('OPEIGR','PARALLEL DIPOLE DYNAMICS '&
         // 'NOT CORRECTLY PROGRAMMED FOR THIS JOB TYPE.',& 
         __LINE__,__FILE__)
    CALL opeigr_para(c0,c2,sc0,nstate,mapful,mapcol,ddmat,nop1,&
         nop2)
    IF (nop2.EQ.0) THEN
       ALLOCATE(aux(nstate,nstate),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL dcopy(2*nstate*nstate,ddmat,1,aux,1)
       CALL determ(aux,nstate,nstate,dd)
       IF (lspin2%tlse) THEN
          CALL dcopy(2*nstate*nstate,ddmat,1,aux,1)
          CALL determ(aux,nstate,nstate-2,dds)
          dd=dd*dds
       ENDIF
       DEALLOCATE(aux,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ELSE
       dd=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
    ENDIF
    CALL tihalt('    OPEIGR',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  CONTAINS
    ! ==================================================================
    SUBROUTINE opeigr_seri(c0,c2,nstate,mapful,ddmat,nop1,nop2)
      ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*), c2(ncpw%ngw,*)
    INTEGER                                  :: nstate, mapful(2,*)
    COMPLEX(real_8)                          :: ddmat(nstate,*)
    INTEGER                                  :: nop1, nop2

    CHARACTER(*), PARAMETER                  :: procedureN = 'opeigr_seri'

    COMPLEX(real_8), ALLOCATABLE             :: aux(:)
    INTEGER                                  :: i, ig

      ALLOCATE(aux(ngg1*ngg2*ngg3),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
           __LINE__,__FILE__)
      DO i=1,nstate
         CALL zeroing(aux)!,ngg1*ngg2*ngg3)
         !$omp parallel do private(IG)
         DO ig=1,spar%ngws
            aux(mapful(1,ig))=c0(ig,i)
            aux(mapful(2,ig))=CONJG(c0(ig,i))
         ENDDO
         CALL opapp(aux,ngg1,ngg2,ngg3,nop1)
         CALL opapp(aux,ngg1,ngg2,ngg3,nop2)
         !$omp parallel do private(IG)
         DO ig=1,spar%ngws
            c2(ig,i)=aux(mapful(1,ig))
         ENDDO
      ENDDO
      IF (cntl%tlsd) THEN
         CALL zgemm('C','N',spin_mod%nsup,spin_mod%nsup,ncpw%ngw,&
              CMPLX(1._real_8,0._real_8,kind=real_8),c0,ncpw%ngw,c2,&
              ncpw%ngw,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat,nstate)
         CALL zgemm('C','N',spin_mod%nsdown,spin_mod%nsdown,ncpw%ngw,&
              CMPLX(1._real_8,0._real_8,kind=real_8),c0(1,&
              spin_mod%nsup+1),ncpw%ngw,c2(1,spin_mod%nsup+1),ncpw%ngw,&
              CMPLX(0._real_8,0._real_8,kind=real_8),ddmat(spin_mod%nsup+1,&
              spin_mod%nsup+1),nstate)
      ELSE
         CALL zgemm('C','N',nstate,nstate,ncpw%ngw,&
              CMPLX(1._real_8,0._real_8,kind=real_8),c0,ncpw%ngw,&
              c2,ncpw%ngw,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat,nstate)
      ENDIF
      ! 
      DO i=1,nstate
         CALL zeroing(aux)!,ngg1*ngg2*ngg3)
         !ocl novrec
         !$omp parallel do private(IG)
         DO ig=1,spar%ngws
            aux(mapful(1,ig))=c0(ig,i)
            aux(mapful(2,ig))=CONJG(c0(ig,i))
         ENDDO
         CALL opapp(aux,ngg1,ngg2,ngg3,nop1)
         CALL opapp(aux,ngg1,ngg2,ngg3,nop2)
         !$omp parallel do private(IG)
         DO ig=1,spar%ngws
            c2(ig,i)=aux(mapful(2,ig))
         ENDDO
         c2(1,i)=CMPLX(0._real_8,0._real_8,kind=real_8)
      ENDDO
      IF (cntl%tlsd) THEN
         CALL zgemm('T','N',spin_mod%nsup,spin_mod%nsup,ncpw%ngw, &
              CMPLX(1._real_8,0._real_8,kind=real_8),c0,ncpw%ngw,c2,&
              ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),ddmat,nstate)
         CALL zgemm('T','N',spin_mod%nsdown,spin_mod%nsdown,ncpw%ngw,&
              CMPLX(1._real_8,0._real_8,kind=real_8),c0(1,&
              spin_mod%nsup+1),ncpw%ngw,c2(1,spin_mod%nsup+1),ncpw%ngw,&
              CMPLX(1._real_8,0._real_8,kind=real_8),ddmat(spin_mod%nsup+1,&
              spin_mod%nsup+1),nstate)
      ELSE
         CALL zgemm('T','N',nstate,nstate,ncpw%ngw,&
              CMPLX(1._real_8,0._real_8,kind=real_8),c0,ncpw%ngw,&
              c2,ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),ddmat,nstate)
      ENDIF
      DEALLOCATE(aux,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
           __LINE__,__FILE__)
      ! ==--------------------------------------------------------------==
      RETURN
    END SUBROUTINE opeigr_seri
    ! ==================================================================
    SUBROUTINE opeigr_para(c0,c2,sc0,nstate,mapful,mapcol,ddmat,&
         nop1,nop2)
      ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)
    COMPLEX(real_8), TARGET                  :: c2(ncpw%ngw,*), sc0(lenbk,*)
    INTEGER                                  :: nstate, mapful(2,*), &
                                                mapcol(ngwmax,*)
    COMPLEX(real_8)                          :: ddmat(nstate,*)
    INTEGER                                  :: nop1, nop2

    COMPLEX(real_8), ALLOCATABLE             :: aux(:), c2f(:,:), c2s(:,:), &
                                                ctr(:,:)
    INTEGER                                  :: i, ig, igii, ii, ip, n1, &
                                                nggp, nn, nnst

      nnst=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
      ALLOCATE(c2s(lenbk,parai%nproc),c2f(spar%ngws,nnst),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)
      ALLOCATE(ctr(spar%ngws,2*nnst),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
           __LINE__,__FILE__)
      ! 
      DO ip=1,parai%nproc
         nn=parap%nst12(ip-1,2)-parap%nst12(ip-1,1)+1
         n1=parap%nst12(ip-1,1)
         CALL dcopy(2*ncpw%ngw*nn,c0(1,n1),1,c2s(1,ip),1)
      ENDDO
      CALL my_trans(c2s,sc0,16*lenbk,1)
      DO ip=1,parai%nproc
         nggp=parap%sparm(3,ip-1)
#ifdef __VECTOR 
         !$omp parallel do private(IGII,II,IG)
         DO igii=1,nnst*nggp
            ii=1+(igii-1)/nggp
            ig=1+MOD(igii-1,nggp)
            c2f(mapcol(ig,ip),ii)=sc0(igii,ip)
         ENDDO
#else 
         !$omp parallel do private(II,IG,IGII) __COLLAPSE2
         DO ii=1,nnst
            DO ig=1,nggp
               igii=(ii-1)*nggp+ig
               c2f(mapcol(ig,ip),ii)=sc0(igii,ip)
            ENDDO
         ENDDO
#endif 
      ENDDO
      ALLOCATE(aux(ngg1*ngg2*ngg3),STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
           __LINE__,__FILE__)
      DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
         ii=i-parap%nst12(parai%mepos,1)+1
         CALL zeroing(aux)!,ngg1*ngg2*ngg3)
#ifdef __NEC 
         !CDIR NODEP
#endif
         !ocl novrec
         !$omp parallel do private(IG)
         DO ig=1,spar%ngws
            aux(mapful(1,ig))=c2f(ig,ii)
            aux(mapful(2,ig))=CONJG(c2f(ig,ii))
         ENDDO
         CALL opapp(aux,ngg1,ngg2,ngg3,nop1)
         CALL opapp(aux,ngg1,ngg2,ngg3,nop2)
         !$omp parallel do private(IG)
         DO ig=1,spar%ngws
            ctr(ig,ii)=aux(mapful(1,ig))
            ctr(ig,nnst+ii)=aux(mapful(2,ig))
         ENDDO
         ctr(1,nnst+ii)=CMPLX(0._real_8,0._real_8,kind=real_8)
      ENDDO
      DEALLOCATE(aux,c2f,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
           __LINE__,__FILE__)
      DO ip=1,parai%nproc
         nggp=parap%sparm(3,ip-1)
#ifdef __VECTOR 
         !$omp parallel do private(IGII,II,IG)
         DO igii=1,nnst*nggp
            ii=1+(igii-1)/nggp
            ig=1+MOD(igii-1,nggp)
            c2s(igii,ip)=ctr(mapcol(ig,ip),ii)
         ENDDO
#else 
         !$omp parallel do private(II,IG,IGII) __COLLAPSE2
         DO ii=1,nnst
            DO ig=1,nggp
               igii=(ii-1)*nggp+ig
               c2s(igii,ip)=ctr(mapcol(ig,ip),ii)
            ENDDO
         ENDDO
#endif 
      ENDDO
      CALL my_trans(c2s,sc0,16*lenbk,1)
      DO ip=1,parai%nproc
         nn=parap%nst12(ip-1,2)-parap%nst12(ip-1,1)+1
         n1=parap%nst12(ip-1,1)
         CALL dcopy(2*ncpw%ngw*nn,sc0(1,ip),1,c2(1,n1),1)
      ENDDO
      IF (cntl%tlsd) THEN
         CALL zgemm('C','N',spin_mod%nsup,spin_mod%nsup,ncpw%ngw,&
              CMPLX(1._real_8,0._real_8,kind=real_8),c0,ncpw%ngw,c2,&
              ncpw%ngw,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat,nstate)
         CALL zgemm('C','N',spin_mod%nsdown,spin_mod%nsdown,ncpw%ngw,&
              CMPLX(1._real_8,0._real_8,kind=real_8),c0(1,&
              spin_mod%nsup+1),ncpw%ngw,c2(1,spin_mod%nsup+1),ncpw%ngw,&
              CMPLX(0._real_8,0._real_8,kind=real_8),ddmat(spin_mod%nsup+1,&
              spin_mod%nsup+1),nstate)
      ELSE
         CALL zgemm('C','N',nstate,nstate,ncpw%ngw,&
              CMPLX(1._real_8,0._real_8,kind=real_8),c0,ncpw%ngw,&
              c2,ncpw%ngw,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat,nstate)
      ENDIF
      ! 
      DO ip=1,parai%nproc
         nggp=parap%sparm(3,ip-1)
#ifdef __VECTOR 
         !$omp parallel do private(IGII,II,IG)
         DO igii=1,nnst*nggp
            ii=1+(igii-1)/nggp
            ig=1+MOD(igii-1,nggp)
            c2s(igii,ip)=ctr(mapcol(ig,ip),nnst+ii)
         ENDDO
#else 
         !$omp parallel do private(II,IG,IGII) __COLLAPSE2
         DO ii=1,nnst
            DO ig=1,nggp
               igii=(ii-1)*nggp+ig
               c2s(igii,ip)=ctr(mapcol(ig,ip),nnst+ii)
            ENDDO
         ENDDO
#endif 
      ENDDO
      CALL my_trans(c2s,sc0,16*lenbk,1)
      DEALLOCATE(c2s,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)

      DO ip=1,parai%nproc
         nn=parap%nst12(ip-1,2)-parap%nst12(ip-1,1)+1
         n1=parap%nst12(ip-1,1)
         CALL dcopy(2*ncpw%ngw*nn,sc0(1,ip),1,c2(1,n1),1)
      ENDDO
      IF (cntl%tlsd) THEN
         CALL zgemm('T','N',spin_mod%nsup,spin_mod%nsup,ncpw%ngw,&
              CMPLX(1._real_8,0._real_8,kind=real_8),c0,ncpw%ngw,c2,&
              ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),ddmat,nstate)
         CALL zgemm('T','N',spin_mod%nsdown,spin_mod%nsdown,ncpw%ngw,&
              CMPLX(1._real_8,0._real_8,kind=real_8),c0(1,&
              spin_mod%nsup+1),ncpw%ngw,c2(1,spin_mod%nsup+1),ncpw%ngw,&
              CMPLX(1._real_8,0._real_8,kind=real_8),ddmat(spin_mod%nsup+1,&
              spin_mod%nsup+1),nstate)
      ELSE
         CALL zgemm('T','N',nstate,nstate,ncpw%ngw,&
              CMPLX(1._real_8,0._real_8,kind=real_8),c0,ncpw%ngw,&
              c2,ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),ddmat,nstate)
      ENDIF
      CALL mp_sum(ddmat,nstate*nstate,parai%allgrp)
      DEALLOCATE(ctr,STAT=ierr)
      IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
      ! ==--------------------------------------------------------------==
      RETURN
    END SUBROUTINE opeigr_para
  END SUBROUTINE opeigr
  ! ==================================================================
  SUBROUTINE opapp(aux,ngg1,ngg2,ngg3,nop)
    ! ==--------------------------------------------------------------==
    ! == Performs some operation of translation in function of NOP    ==
    ! == NOP == 1 -> AUX(x-1,:,:) = AUX(x,:,:)                        ==
    ! == NOP == 2 -> AUX(:,y-1,:) = AUX(:,y,:)                        ==
    ! == NOP == 3 -> AUX(:,:,z-1) = AUX(:,:,z)                        ==
    ! == NOP == -1 or -2 or -3 are the inverse operations             ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ngg1, ngg2, ngg3
    COMPLEX(real_8)                          :: aux(ngg1,ngg2,ngg3)
    INTEGER                                  :: nop

    INTEGER                                  :: ix, iy, iz

    IF (nop.EQ.1) THEN
       DO ix=2,ngg1
          CALL zcopy(ngg2*ngg3,aux(ix,1,1),ngg1,aux(ix-1,1,1),ngg1)
       ENDDO
    ELSEIF (nop.EQ.2) THEN
       DO iz=1,ngg3
          DO iy=2,ngg2
             CALL zcopy(ngg1,aux(1,iy,iz),1,aux(1,iy-1,iz),1)
          ENDDO
       ENDDO
    ELSEIF (nop.EQ.3) THEN
       DO iz=2,ngg3
          CALL zcopy(ngg1*ngg2,aux(1,1,iz),1,aux(1,1,iz-1),1)
       ENDDO
    ELSEIF (nop.EQ.(-1)) THEN
       DO ix=ngg1,2,-1
          CALL zcopy(ngg2*ngg3,aux(ix-1,1,1),ngg1,aux(ix,1,1),ngg1)
       ENDDO
    ELSEIF (nop.EQ.(-2)) THEN
       DO iz=1,ngg3
          DO iy=ngg2,2,-1
             CALL zcopy(ngg1,aux(1,iy-1,iz),1,aux(1,iy,iz),1)
          ENDDO
       ENDDO
    ELSEIF (nop.EQ.(-3)) THEN
       DO iz=ngg3,2,-1
          CALL zcopy(ngg1*ngg2,aux(1,1,iz-1),1,aux(1,1,iz),1)
       ENDDO
    ELSEIF (nop.EQ.0) THEN
    ELSE
       CALL stopgm('OPAPP','WRONG NOP',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE opapp
  ! ==================================================================
  SUBROUTINE give_scr_opeigr(lopeigr,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lopeigr
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: l1, l2

    l1=2*nstate*nstate
    l2=0
    CALL getngg(l2)
    lopeigr=MAX(l1,l2,2*ngg1*ngg2*ngg3)
    tag='MAX(2*NSTATE**2,2*NGG1*NGG2*NGG3)'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_opeigr
  ! ==================================================================
  SUBROUTINE getngg(l2)
    ! ==--------------------------------------------------------------==
    ! input
    INTEGER                                  :: l2

    INTEGER                                  :: i, j, k, kmax, kmin, nx1, &
                                                nx2, nx3
    REAL(real_8)                             :: g2, t1, t2, t3

! ==--------------------------------------------------------------==

    nx1=0
    nx2=0
    nx3=0
    ! Maximum K point amplitude (Kmax) in Brillouin zone.
    gcutka=gkpwf(gvec_com%b1,gvec_com%b2,gvec_com%b3,tkpts%tkpnt.AND.tsphere)
    IF (gcutka.EQ.0._real_8) THEN
       gcutwmin=gvec_com%gcutw
       gcutwmax=gvec_com%gcutw
    ELSE
       ! For each k and each G, |k+G|^2 < GCUTW
       ! so we define GCUTWMIN and GCUTWMAX as:
       ! GCUTWMIN < |G|^2 < GCUTWMAX we have to apply a mask.
       gcutwmin = MAX(0._real_8,SQRT(gvec_com%gcutw)-SQRT(gcutka))**2
       gcutwmax = (SQRT(gvec_com%gcutw)+SQRT(gcutka))**2
    ENDIF
    ! ..construction of g-vectors
    ! !    IG=0
    ! mb...I=0
    !$omp parallel do private(J,K,KMIN,KMAX,T1,T2,T3,G2) &
    !$omp  reduction(max:NX2,NX3)
    DO j=0,spar%nr2s-1
       kmin=-spar%nr3s+1
       kmax=spar%nr3s-1
       IF (j.EQ.0) kmin=0
       DO k=kmin,kmax
          t1=REAL(j,kind=real_8)*gvec_com%b2(1)+REAL(k,kind=real_8)*gvec_com%b3(1)
          t2=REAL(j,kind=real_8)*gvec_com%b2(2)+REAL(k,kind=real_8)*gvec_com%b3(2)
          t3=REAL(j,kind=real_8)*gvec_com%b2(3)+REAL(k,kind=real_8)*gvec_com%b3(3)
          g2=t1*t1+t2*t2+t3*t3
          IF (g2.LT.gcutwmax) THEN
             ! !          IG=IG+1
             ! !          NX1=MAX(0,NX1)
             nx2=MAX(ABS(j),nx2)
             nx3=MAX(ABS(k),nx3)
          ENDIF
       ENDDO
    ENDDO
    ! mb...I<>0
    !$omp parallel do private(I,J,K,T1,T2,T3,G2) &
    !$omp  reduction(max:NX1,NX2,NX3)
    DO i=1,spar%nr1s-1
       DO j=-spar%nr2s+1,spar%nr2s-1
          DO k=-spar%nr3s+1,spar%nr3s-1
             t1=REAL(i,kind=real_8)*gvec_com%b1(1)+REAL(j,kind=real_8)*gvec_com%b2(1)+REAL(k,kind=real_8)*gvec_com%b3(1)
             t2=REAL(i,kind=real_8)*gvec_com%b1(2)+REAL(j,kind=real_8)*gvec_com%b2(2)+REAL(k,kind=real_8)*gvec_com%b3(2)
             t3=REAL(i,kind=real_8)*gvec_com%b1(3)+REAL(j,kind=real_8)*gvec_com%b2(3)+REAL(k,kind=real_8)*gvec_com%b3(3)
             g2=t1*t1+t2*t2+t3*t3
             IF (g2.LT.gcutwmax) THEN
                ! !            IG=IG+1
                nx1=MAX(ABS(i),nx1)
                nx2=MAX(ABS(j),nx2)
                nx3=MAX(ABS(k),nx3)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    ! mb
    nx1=2*(nx1+1)+1
    nx2=2*(nx2+1)+1
    nx3=2*(nx3+1)+1
    l2=nx1*nx2*nx3*2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE getngg
  ! ==================================================================

END MODULE opeigr_utils

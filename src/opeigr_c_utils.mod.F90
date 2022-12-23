MODULE opeigr_c_utils
  USE ddip,                            ONLY: lenbk,&
                                             ngg1,&
                                             ngg2,&
                                             ngg3,&
                                             ngwmax
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE opeigr_utils,                    ONLY: opapp
  USE parac,                           ONLY: parai,&
                                             paral
  USE reshaper,                        ONLY: reshape_inplace
  USE spin,                            ONLY: lspin2,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             nkpt,&
                                             parap,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: determ
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: opeigr_c

CONTAINS

  ! ==================================================================
  SUBROUTINE opeigr_c(c0,c2,sc0,nstate,mapful,mapcol,ddmat,&
       nop1,nop2,dd)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,*)
    COMPLEX(real_8), TARGET                  :: c2(nkpt%ngwk,*), &
                                                sc0(spar%ngwks,*)
    INTEGER                                  :: nstate, mapful(2,*), &
                                                mapcol(ngwmax,*)
    COMPLEX(real_8)                          :: ddmat(nstate,nstate)
    INTEGER                                  :: nop1, nop2
    COMPLEX(real_8)                          :: dd

    CHARACTER(*), PARAMETER                  :: procedureN = 'opeigr_c'

    COMPLEX(real_8)                          :: dds
    COMPLEX(real_8), ALLOCATABLE             :: aux(:)
    COMPLEX(real_8), POINTER                 :: c2f(:,:), c2s(:,:), sc0s(:,:)
    INTEGER                                  :: i, ierr, ig, igii, ii, ip, &
                                                isub, n1, nggp, nn

    IF (paral%io_parent)&
         WRITE(6,*)'NGWK=',nkpt%ngwk,'NSTATE=',nstate
    CALL tiset('    OPEIGR',isub)
    CALL reshape_inplace(c2, (/spar%ngwks, nstate/), c2f)
    CALL zeroing(ddmat)!,nstate*nstate)
    CALL reshape_inplace(c2, (/lenbk, parai%nproc/), c2s)
    CALL reshape_inplace(sc0, (/lenbk, parai%nproc/), sc0s)
    DO ip=1,parai%nproc
       nn=parap%nst12(ip-1,2)-parap%nst12(ip-1,1)+1
       n1=parap%nst12(ip-1,1)
       CALL dcopy(2*nkpt%ngwk*nn,c0(1,n1),1,c2s(1,ip),1)
    ENDDO
    CALL my_trans(c2s,sc0s,16*lenbk,1)
    nn=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
    DO ip=1,parai%nproc
       nggp=parap%sparm(3,ip-1)
#ifdef __VECTOR 
       !$omp parallel do private(IGII,II,IG)
       DO igii=1,nn*nggp
          ii=1+(igii-1)/nggp
          ig=1+MOD(igii-1,nggp)
          c2f(mapcol(ig,ip),ii)=sc0s(igii,ip)
       ENDDO
#else 
       igii=0
       DO ii=1,nn
          DO ig=1,nggp
             igii=igii+1
             c2f(mapcol(ig,ip),ii)=sc0s(igii,ip)
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
       !$omp parallel do private(IG)
       DO ig=1,spar%ngwks
          aux(mapful(1,ig))=c2f(ig,ii)
          aux(mapful(2,ig))=CONJG(c2f(ig,ii))
       ENDDO
       CALL opapp(aux,ngg1,ngg2,ngg3,nop1)
       CALL opapp(aux,ngg1,ngg2,ngg3,nop2)
       !$omp parallel do private(IG)
       DO ig=1,spar%ngwks
          sc0(ig,ii)=aux(mapful(1,ig))
       ENDDO
    ENDDO
    nn=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
    DO ip=1,parai%nproc
       nggp=parap%sparm(3,ip-1)
#ifdef __VECTOR 
       !$omp parallel do private(IGII,II,IG)
       DO igii=1,nn*nggp
          ii=1+(igii-1)/nggp
          ig=1+MOD(igii-1,nggp)
          c2s(igii,ip)=sc0(mapcol(ig,ip),ii)
       ENDDO
#else 
       igii=0
       DO ii=1,nn
          DO ig=1,nggp
             igii=igii+1
             c2s(igii,ip)=sc0(mapcol(ig,ip),ii)
          ENDDO
       ENDDO
#endif 
    ENDDO
    CALL my_trans(c2s,sc0s,16*lenbk,1)
    DO ip=1,parai%nproc
       nn=parap%nst12(ip-1,2)-parap%nst12(ip-1,1)+1
       n1=parap%nst12(ip-1,1)
       CALL dcopy(2*nkpt%ngwk*nn,sc0s(1,ip),1,c2(1,n1),1)
    ENDDO
    ! IF(cntl%tlsd) THEN
    ! CALL AZZERO(A,NSTATE*NSTATE)
    ! ..Alpha spin
    ! CALL ZGEMM('C','N',NSUP,NSUP,NGWK,ZONE,C1(1,1),NGWK,C2(1,1),
    ! &     NGWK,ZZERO,A(1,1),NSTATE)
    ! ..Beta spin
    ! CALL ZGEMM('C','N',NSDOWN,NSDOWN,NGWK,ZONE,C1(1,NSUP+1),
    ! &     NGWK,C2(1,NSUP+1),NGWK,ZZERO,A(NSUP+1,NSUP+1),NSTATE)
    ! ELSE
    ! CALL ZGEMM('C','N',NSTATE,NSTATE,NGWK,ZONE,C1(1,1),NGWK,
    ! &     C2(1,1),NGWK,ZZERO,A(1,1),NSTATE)
    ! ENDIF

    IF (cntl%tlsd) THEN
       CALL zgemm('C','N',spin_mod%nsup,spin_mod%nsup,nkpt%ngwk,CMPLX(1._real_8,0._real_8,kind=real_8),c0(1,1),&
            nkpt%ngwk,c2(1,1),nkpt%ngwk,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat(1,1),nstate)
       CALL zgemm('C','N',spin_mod%nsdown,spin_mod%nsdown,nkpt%ngwk,CMPLX(1._real_8,0._real_8,kind=real_8),c0(1,&
            spin_mod%nsup+1),nkpt%ngwk,c2(1,spin_mod%nsup+1),nkpt%ngwk,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat(spin_mod%nsup+&
            1,spin_mod%nsup+1),nstate)
    ELSE
       CALL zgemm('C','N',nstate,nstate,nkpt%ngwk,CMPLX(1._real_8,0._real_8,kind=real_8),c0,nkpt%ngwk,&
            c2,nkpt%ngwk,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat,nstate)
    ENDIF
    CALL mp_sum(ddmat,nstate*nstate,parai%allgrp)

    DO i=1,nstate
       IF (paral%io_parent)&
            WRITE(6,*)'I=',i,'DDMAT=',ddmat(i,i)
    ENDDO

    IF (nop2.EQ.0) THEN
       CALL dcopy(2*nstate*nstate,ddmat,1,aux,1)
       CALL determ(aux,nstate,nstate,dd)
       IF (lspin2%tlse) THEN
          CALL dcopy(2*nstate*nstate,ddmat,1,aux,1)
          CALL determ(aux,nstate,nstate-2,dds)
          dd=dd*dds
       ENDIF
    ELSE
       dd=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
    ENDIF
    DEALLOCATE(aux,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    CALL tihalt('    OPEIGR',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE opeigr_c
  ! ==================================================================

END MODULE opeigr_c_utils

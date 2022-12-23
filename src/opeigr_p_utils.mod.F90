MODULE opeigr_p_utils
  USE ddip,                            ONLY: lenbk,&
                                             ngg1,&
                                             ngg2,&
                                             ngg3,&
                                             ngwmax
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE opeigr_utils,                    ONLY: opapp
  USE parac,                           ONLY: parai
  USE reshaper,                        ONLY: reshape_inplace
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             parap,&
                                             spar
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: opeigr_p

CONTAINS

  ! ==================================================================
  SUBROUTINE opeigr_p(c0,c2,sc0,nstate,mapful,mapcol,ddmat,&
       nop1,nop2,nop3,dd,cwork)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(ncpw%ngw,*)
    COMPLEX(real_8), TARGET                  :: c2(ncpw%ngw,*), &
                                                sc0(spar%ngws,*)
    INTEGER                                  :: nstate, mapful(2,*), &
                                                mapcol(ngwmax,*)
    COMPLEX(real_8)                          :: ddmat(nstate,nstate)
    INTEGER                                  :: nop1, nop2, nop3
    COMPLEX(real_8)                          :: dd, cwork(ncpw%ngw,nstate,2)

    COMPLEX(real_8), ALLOCATABLE             :: aux(:)
    COMPLEX(real_8), POINTER                 :: c2f(:,:), c2s(:,:), sc0s(:,:)
    INTEGER                                  :: i, ig, igii, ii, ip, istat, &
                                                length, n1, nggp, nn

! ==--------------------------------------------------------------==

    length=0
    ALLOCATE(aux(ngg1*ngg2*ngg3),stat=istat)
    IF (istat/=0) CALL stopgm('OPEIGR_P','allocation error',& 
         __LINE__,__FILE__)

    CALL reshape_inplace(c2, (/spar%ngws, nstate/), c2f)
    CALL zeroing(ddmat)!,nstate*nstate)
    CALL reshape_inplace(c2, (/lenbk, parai%nproc/), c2s)
    CALL reshape_inplace(sc0, (/lenbk, parai%nproc/), sc0s)
    DO ip=1,parai%nproc
       nn=parap%nst12(ip-1,2)-parap%nst12(ip-1,1)+1
       n1=parap%nst12(ip-1,1)
       CALL dcopy(2*ncpw%ngw*nn,c0(1,n1),1,c2s(1,ip),1)
    ENDDO
    CALL my_trans(c2s,sc0s,16*lenbk,1)
    nn=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
    !$omp parallel do private(IP,NGGP,IGII,II,IG)
    DO ip=1,parai%nproc
       nggp=parap%sparm(3,ip-1)
       igii=0
       DO ii=1,nn
          DO ig=1,nggp
             igii=igii+1
             c2f(mapcol(ig,ip),ii)=sc0s(igii,ip)
          ENDDO
       ENDDO
    ENDDO
    DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
       ii=i-parap%nst12(parai%mepos,1)+1
       CALL zeroing(aux)!,ngg1*ngg2*ngg3)
       !$omp parallel do private(IG)
       DO ig=1,spar%ngws
          aux(mapful(1,ig))=c2f(ig,ii)
          aux(mapful(2,ig))=CONJG(c2f(ig,ii))
       ENDDO
       CALL opapp(aux,ngg1,ngg2,ngg3,nop1)
       CALL opapp(aux,ngg1,ngg2,ngg3,nop2)
       CALL opapp(aux,ngg1,ngg2,ngg3,nop3)

       !$omp parallel do private(IG)
       DO ig=1,spar%ngws
          sc0(ig,ii)=aux(mapful(1,ig))
       ENDDO

    ENDDO

    nn=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
    !$omp parallel do private(IP,NGGP,IGII,II,IG)
    DO ip=1,parai%nproc
       nggp=parap%sparm(3,ip-1)
       igii=0
       DO ii=1,nn
          DO ig=1,nggp

             igii=igii+1

             c2s(igii,ip)=sc0(mapcol(ig,ip),ii)
          ENDDO
       ENDDO
    ENDDO
    CALL my_trans(c2s,sc0s,16*lenbk,1)
    DO ip=1,parai%nproc
       nn=parap%nst12(ip-1,2)-parap%nst12(ip-1,1)+1
       n1=parap%nst12(ip-1,1)
       IF (nn.GT.0) THEN
          CALL dcopy(2*ncpw%ngw*nn,sc0s(1,ip),1,c2(1,n1),1)
          CALL zcopy(ncpw%ngw*nn,sc0s(1,ip),1,cwork(1,n1,1),1)
       ENDIF
    ENDDO

    IF (cntl%tlsd) THEN
       CALL zgemm('C','N',spin_mod%nsup,spin_mod%nsup,ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),c0,ncpw%ngw,c2,&
            ncpw%ngw,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat,nstate)
       CALL zgemm('C','N',spin_mod%nsdown,spin_mod%nsdown,ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),c0(1,&
            spin_mod%nsup+1),ncpw%ngw,c2(1,spin_mod%nsup+1),ncpw%ngw,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat(spin_mod%nsup+1,&
            spin_mod%nsup+1),nstate)
    ELSE
       CALL zgemm('C','N',nstate,nstate,ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),c0,ncpw%ngw,&
            c2,ncpw%ngw,CMPLX(0._real_8,0._real_8,kind=real_8),ddmat,nstate)
    ENDIF
    DO ip=1,parai%nproc
       nn=parap%nst12(ip-1,2)-parap%nst12(ip-1,1)+1
       n1=parap%nst12(ip-1,1)
       CALL dcopy(2*ncpw%ngw*nn,c0(1,n1),1,c2s(1,ip),1)
    ENDDO
    CALL my_trans(c2s,sc0s,16*lenbk,1)
    nn=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
    !$omp parallel do private(IP,NGGP,IGII,II,IG)
    DO ip=1,parai%nproc
       nggp=parap%sparm(3,ip-1)
       igii=0
       DO ii=1,nn
          DO ig=1,nggp

             igii=igii+1

             c2f(mapcol(ig,ip),ii)=sc0s(igii,ip)
          ENDDO
       ENDDO
    ENDDO

    DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
       ii=i-parap%nst12(parai%mepos,1)+1
       CALL zeroing(aux)!,ngg1*ngg2*ngg3)
       !$omp parallel do private(IG)
       DO ig=1,spar%ngws
          aux(mapful(1,ig))=c2f(ig,ii)
          aux(mapful(2,ig))=CONJG(c2f(ig,ii))
       ENDDO
       CALL opapp(aux,ngg1,ngg2,ngg3,nop1)
       CALL opapp(aux,ngg1,ngg2,ngg3,nop2)
       CALL opapp(aux,ngg1,ngg2,ngg3,nop3)

       !$omp parallel do private(IG)
       DO ig=1,spar%ngws
          sc0(ig,ii)=aux(mapful(2,ig))
       ENDDO
       sc0(1,ii)=CMPLX(0._real_8,0._real_8,kind=real_8)

    ENDDO
    nn=parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1
    !$omp parallel do private(IP,NGGP,IGII,II,IG)
    DO ip=1,parai%nproc
       nggp=parap%sparm(3,ip-1)
       igii=0
       DO ii=1,nn
          DO ig=1,nggp

             igii=igii+1

             c2s(igii,ip)=sc0(mapcol(ig,ip),ii)
          ENDDO
       ENDDO
    ENDDO
    CALL my_trans(c2s,sc0s,16*lenbk,1)
    DO ip=1,parai%nproc
       nn=parap%nst12(ip-1,2)-parap%nst12(ip-1,1)+1
       n1=parap%nst12(ip-1,1)
       IF (nn.GT.0) THEN
          CALL dcopy(2*ncpw%ngw*nn,sc0s(1,ip),1,c2(1,n1),1)
          CALL zcopy(ncpw%ngw*nn,sc0s(1,ip),1,cwork(1,n1,2),1)
       ENDIF
    ENDDO

    IF (cntl%tlsd) THEN
       CALL zgemm('T','N',spin_mod%nsup,spin_mod%nsup,ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),c0,ncpw%ngw,c2,&
            ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),ddmat,nstate)
       CALL zgemm('T','N',spin_mod%nsdown,spin_mod%nsdown,ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),c0(1,&
            spin_mod%nsup+1),ncpw%ngw,c2(1,spin_mod%nsup+1),ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),ddmat(spin_mod%nsup+1,&
            spin_mod%nsup+1),nstate)
    ELSE
       CALL zgemm('T','N',nstate,nstate,ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),c0,ncpw%ngw,&
            c2,ncpw%ngw,CMPLX(1._real_8,0._real_8,kind=real_8),ddmat,nstate)
    ENDIF
    CALL mp_sum(ddmat,nstate*nstate,parai%allgrp)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(aux,stat=istat)
    IF (istat/=0) CALL stopgm('OPEIGR_P','deallocation error',& 
         __LINE__,__FILE__)
    RETURN
  END SUBROUTINE opeigr_p

END MODULE opeigr_p_utils

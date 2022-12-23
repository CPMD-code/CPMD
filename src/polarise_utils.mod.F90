MODULE polarise_utils
  USE atimesmod,                       ONLY: atimes_eval,&
                                             ikind_atimes
  USE calc_pij_utils,                  ONLY: calc_pi,&
                                             calc_pij
  USE clinbcg_utils,                   ONLY: clinbcg,&
                                             give_scr_clinbcg
  USE cnst,                            ONLY: fpi,&
                                             uimag
  USE cppt,                            ONLY: &
       gk, hg, indz, indzs, inyh, nzh, nzhs, rhops, vps
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: rk,&
                                             wk
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_walltime
  USE parac,                           ONLY: parai,&
                                             paral
  USE pola,                            ONLY: &
       alphap, alphapav, itmax, itol, tol, tollocc, tzeff, zeff, zeffav
  USE projv_utils,                     ONLY: projv
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: fpar,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgthr
  USE wrener_utils,                    ONLY: wreigen
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: polarise
  PUBLIC :: polarproj
  PUBLIC :: give_scr_polarproj

CONTAINS

  ! ==================================================================
  SUBROUTINE polarise(c0,we,f,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! == ROUTINE TO COMPUTE THE POLARISABILITY (A. Alavi 1997)        ==
    ! == VIA DIRECT SUMMATION OVER EXCITED STATES 
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate, nkpoint
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate,nkpoint)
    REAL(real_8)                             :: we(nstate,nkpoint), &
                                                f(nstate,nkpoint)

    COMPLEX(real_8), PARAMETER               :: zzero = (0._real_8,0._real_8)

    COMPLEX(real_8)                          :: cpx, cpy, cpz
    INTEGER                                  :: i, ik, j
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: alphaxx, alphaxxa, alphayy, &
                                                alphayya, alphazz, alphazza, &
                                                deng, df, SUM(3)

    IF (ifirst.EQ.0) THEN
       ifirst=1
       alphaxxa=0._real_8
       alphayya=0._real_8
       alphazza=0._real_8
    ENDIF
    CALL zeroing(sum)!,3)
    DO ik = 1 , nkpoint
       DO i = 1 , nstate-1
          DO j = i+1,nstate
             df = f(i,ik) - f(j,ik)
             IF (df.LT.1.e-9_real_8) GOTO 999
             deng = we(j,ik)-we(i,ik)
             CALL calc_pij(c0(1,j,ik),c0(1,i,ik),cpx,cpy,cpz,ik)
             SUM(1) = SUM(1) + REAL(df*wk(ik)*cpx*CONJG(cpx)/deng**3)
             SUM(2) = SUM(2) + REAL(df*wk(ik)*cpy*CONJG(cpy)/deng**3)
             SUM(3) = SUM(3) + REAL(df*wk(ik)*cpz*CONJG(cpz)/deng**3)
999          CONTINUE
          ENDDO
       ENDDO
    ENDDO
    alphaxx=SUM(1)*parm%tpiba2
    alphayy=SUM(2)*parm%tpiba2
    alphazz=SUM(3)*parm%tpiba2
    alphaxxa=alphaxxa+alphaxx
    alphayya=alphayya+alphayy
    alphazza=alphazza+alphazz
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE polarise
  ! ==================================================================
  SUBROUTINE polarproj(c0,we,f,amu,&
       vpot,psi,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! == ROUTINE TO COMPUTE THE POLARISABILITY                        ==
    ! == BY USING PROJECTION METHOD                                   ==
    ! == ( P. Gianozzi, S de Gironcoli, P. Pavone and S. Baroni,      ==
    ! == Phys. Rev. B, 43, 7231, (1991)                               ==
    ! ==                                                              ==
    ! == VALENCE STATES ARE DEFINED AS THOSE FOR WHOM THE OCCUPATION  ==
    ! == NUMBERIS GREATER THAN TOLLOCC ( specified in pola.inc)       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: amu, vpot(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(fpar%nnr1)
    INTEGER                                  :: nstate, nkpoint
    COMPLEX(real_8)                          :: c0(nkpt%ngwk,nstate,nkpoint)
    REAL(real_8)                             :: we(nstate,nkpoint), &
                                                f(nstate,nkpoint)

    CHARACTER(*), PARAMETER                  :: procedureN = 'polarproj'
    COMPLEX(real_8), PARAMETER               :: zzero = (0._real_8,0._real_8)
    INTEGER, DIMENSION(6), PARAMETER         :: jwh = (/1,1,1,2,2,3/), &
                                                kwh = (/1,2,3,2,3,3/)

    CHARACTER(len=1), SAVE                   :: comp(3)
    COMPLEX(real_8)                          :: caux, ei123, zdotc
    COMPLEX(real_8), ALLOCATABLE             :: psi2(:), vt1pp(:)
    COMPLEX(real_8), ALLOCATABLE, SAVE       :: cb(:), cproj(:,:), cx(:,:)
    EXTERNAL                                 :: zdotc
    INTEGER                                  :: i, ia, ierr, ifft, ig, ik, &
                                                ir, is, isa, isub, j, jj, k, &
                                                kk, l, m
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: aux, auxx(3,3), err, tim, &
                                                tim1, tim2
    REAL(real_8), ALLOCATABLE, SAVE          :: vext(:,:)

#ifdef __NEC
    REAL(real_8)                             :: SUM(3,3)
    INTEGER                                  :: isa0

    INTEGER, ALLOCATABLE                     :: isv(:)
    INTEGER, ALLOCATABLE                     :: iav(:)
#else
    REAL(real_8)                             :: sum
#endif

    ! ==--------------------------------------------------------------==
    CALL tiset(' POLARPROJ',isub)
    CALL setfftn(0)
    IF (ifirst.EQ.0) THEN
       ifirst=1
       CALL zeroing(alphapav)!,9)
       ALLOCATE(cproj(nkpt%ngwk,3),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cb(nkpt%ngwk),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(cx(nkpt%ngwk,3),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       comp(1) = 'X'
       comp(2) = 'Y'
       comp(3) = 'Z'
       IF (tzeff) THEN
          ALLOCATE(zeff(3,3,maxsys%nax,ions1%nsp),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(zeffav(3,3,maxsys%nax,ions1%nsp),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(vt1pp(maxfft),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(vext(ncpw%nhg,ions1%nsp),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(psi2(fpar%nnr1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL zeroing(zeffav)!,9*maxsys%nax*ions1%nsp)
#ifdef __NEC
          isa0=0
          DO is=1,ions1%nsp
             isa0 = isa0 + ions0%na(is)
          ENDDO
          ALLOCATE(isv(isa0),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(iav(isa0),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          isa=0
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                isa=isa+1
                isv(isa) = is
                iav(isa) = ia
             ENDDO
          ENDDO
#endif
          ! Compute external potential in g-space. 
          ! VPS=short-range pseudpotential term in g-space
          ! RHOPS=e^-g^2, arising from the Gaussian smearing of nuclear charge. 
          DO is=1,ions1%nsp
             DO ig=2,ncpw%nhg
                vext(ig,is)=(vps(is,ig)+fpi*rhops(is,ig)/hg(ig)/parm%tpiba2)*&
                     parm%tpiba2
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    ifft=1
    IF (paral%io_parent)&
         WRITE(6,'(1X,64("<")/20X,A26/1X,64(">"))')&
         'CALCULATING POLARISABILITY'
    CALL wreigen(we,f,amu,nstate)
    CALL zeroing(alphap)!,9)
    DO ik = 1 , nkpoint
       IF (paral%io_parent)&
            WRITE(6,'(1X,64("="))')
       IF (paral%io_parent)&
            WRITE(6,'(A9,I4,3X,4(F10.7,3X))')&
            ' K POINT:',ik,rk(1,ik),rk(2,ik),rk(3,ik),wk(ik)
       DO i = 1 , nstate
          IF (f(i,ik).LT.tollocc) GOTO 999
          IF (paral%io_parent)&
               WRITE(6,'(A9,I4,4X,A11,F7.3)')&
               '   STATE:',i,'OCCUPATION:',f(i,ik)
          atimes_eval=we(i,ik)
          ikind_atimes=ik
          tim1=m_walltime()
          DO j=1,3
             ! Returns the j-th component of momentum operator in cproj
             CALL calc_pi(c0(1,i,ik),cproj(1,j),j,ik)
             CALL projv(cproj(1,j),c0(1,1,ik),f(1,ik),nstate)
             ! Multiply by -i 
             ! CALL ZSCAL(NGWK,-UIMAG,CPROJ(1,J),1)
             CALL dcopy(nkpt%ngwk*2,cproj(1,j),1,cb,1)
             CALL zeroing(cx(:,j))!,nkpt%ngwk)
             ! Solve A.x=b
             CALL clinbcg(cb,cx(1,j),c0(1,1,ik),vpot,psi,&
                  nstate,f(1,ik),itol,tol,itmax,err)
             ! Project out valence states 
             CALL projv(cx(1,j),c0(1,1,ik),f(1,ik),nstate)
             ! Store away 
             CALL dcopy(nkpt%ngwk*2,cx(1,j),1,cproj(1,j),1)
             CALL dcopy(nkpt%ngwk*2,cx(1,j),1,cb,1)
             CALL zeroing(cx(:,j))!,nkpt%ngwk)
             ! Solve Ax=b.   
             CALL clinbcg(cb,cx(1,j),c0(1,1,ik),vpot,psi,&
                  nstate,f(1,ik),itol,tol,itmax,err)
             CALL projv(cx(1,j),c0(1,1,ik),f(1,ik),nstate)
          ENDDO
          tim2=m_walltime()
          tim=(tim2-tim1)*0.001_real_8
          IF (tkpts%tkpnt) THEN
             DO jj=1,6
                j=jwh(jj)
                k=kwh(jj)
                aux=REAL(zdotc(nkpt%ngwk,cproj(1,j),1,cx(1,k),1))
                auxx(j,k)=2._real_8*aux*wk(ik)*f(i,ik)
                alphap(j,k)=alphap(j,k)+auxx(j,k)
             ENDDO
          ELSE
             DO jj=1,6
                j=jwh(jj)
                k=kwh(jj)
                aux=dotp(nkpt%ngwk,cproj(:,j),cx(:,k))
                auxx(j,k)=2._real_8*aux*wk(ik)*f(i,ik)
                alphap(j,k)=alphap(j,k)+auxx(j,k)
             ENDDO
          ENDIF
          IF (paral%io_parent)&
               WRITE(6,'(1X,3(A1,A1,F10.4,3X),3X,A6,F9.2)')&
               comp(jwh(1)),comp(kwh(1)),auxx(jwh(1),kwh(1)),&
               comp(jwh(4)),comp(kwh(4)),auxx(jwh(4),kwh(4)),&
               comp(jwh(6)),comp(kwh(6)),auxx(jwh(6),kwh(6)),&
               ' TIME:',tim
          ! ==----------------------------------------------------------==
          ! == Effective charge calculation                             ==
          ! ==----------------------------------------------------------==
          ! First transform c0 to real space and store in PSI 
          IF (tzeff) THEN
             CALL zeroing(psi)!,maxfft)
#ifdef __NEC
             !CDIR NODEP
#else
             !$omp parallel do private(IG)
#endif
             DO ig = 1 , ncpw%ngw
                psi(nzhs(ig) ) = c0(ig,i,ik)
                psi(indzs(ig)) = c0(ig+ncpw%ngw,i,ik)
             ENDDO
             IF (geq0) psi(nzhs(1)) = c0(1,i,ik)
             IF (ifft.EQ.1) THEN
                CALL  invfftn(psi,.TRUE.,parai%allgrp)
             ELSE
                CALL stopgm("THIS","SG_INVFFT NOT AVAILABLE ANYMORE",& 
                     __LINE__,__FILE__)
             ENDIF
             CALL dcopy(2*fpar%nnr1,psi,1,psi2,1)
             ! Next compute dV/dRI in g-space, transform to real space
             ! and form dV/dRI|C0>. 
             isa=0
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   isa=isa+1
                   DO j=1,3
                      ! Put dV/dRIj in g-space in VT1PP 
                      CALL zeroing(vt1pp)!,maxfft)
#ifdef __NEC
                      !CDIR NODEP
#else
                      !$omp parallel do private(IG,EI123,CAUX)
#endif
                      DO ig=2,ncpw%nhg
                         ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                              ei3(isa,inyh(3,ig))
                         caux=vext(ig,is)*CMPLX(0._real_8,gk(j,ig),kind=real_8)*&
                              parm%tpiba*ei123
                         vt1pp(nzh(ig))=caux
                         vt1pp(indz(ig))=CONJG(caux)
                      ENDDO
                      ! Transform to real space 
                      CALL  invfftn(vt1pp,.FALSE.,parai%allgrp)
                      ! Multiply by psi in real space  
                      !$omp parallel do private(IR)
                      DO ir=1,fpar%nnr1
                         vt1pp(ir)=vt1pp(ir)*psi2(ir)
                      ENDDO
                      ! FWFFT VT1PP 
                      IF (ifft.EQ.1) THEN
                         CALL  fwfftn(vt1pp,.TRUE.,parai%allgrp)
                      ELSE
                         CALL stopgm("THIS","SG_FWFFT NOT AVAILABLE ANYMORE",& 
                              __LINE__,__FILE__)
                      ENDIF
                      ! Gather VT1PP into CB  
                      IF (tkpts%tkpnt) THEN
                         !$omp parallel do private(IG)
                         DO ig=1,ncpw%ngw
                            cb(ig)=vt1pp(nzhs(ig))
                            cb(ig+ncpw%ngw)=vt1pp(indzs(ig))
                         ENDDO
                         IF (geq0) cb(ncpw%ngw+1)=zzero
                      ELSE
                         CALL zgthr(ncpw%ngw,vt1pp,cb,nzhs)
                      ENDIF
                      CALL zscal(nkpt%ngwk,-uimag,cb,1)
                      IF (tkpts%tkpnt) THEN
                         DO jj=1,3
                            aux=4._real_8*REAL(zdotc(nkpt%ngwk,cx(1,jj),1,cb,1))*&
                                 f(i,ik)*wk(ik)
                            zeff(jj,j,ia,is)=zeff(jj,j,ia,is)+aux
                         ENDDO
                      ELSE
                         DO jj=1,3
                            aux=4._real_8*dotp(ncpw%ngw,cx(:,jj),cb)*&
                                 f(i,ik)*wk(ik)
                            zeff(jj,j,ia,is)=zeff(jj,j,ia,is)+aux
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
999       CONTINUE
       ENDDO
    ENDDO
    DO kk=1,6
       alphapav(jwh(kk),kwh(kk))=alphapav(jwh(kk),kwh(kk))+&
            alphap(jwh(kk),kwh(kk))
    ENDDO
    alphap(2,1)=alphap(1,2)
    alphap(3,1)=alphap(1,3)
    alphap(3,2)=alphap(2,3)
    alphapav(2,1)=alphapav(1,2)
    alphapav(3,1)=alphapav(1,3)
    alphapav(3,2)=alphapav(2,3)
    IF (tzeff) THEN
       ! Apply acoustic sum rule.
#ifdef __NEC
       CALL zeroing(sum)!,3*3)
       DO isa = 1,isa0
          is= isv(isa)
          ia= iav(isa)
          DO l=1,3
             DO m=1,3
                SUM(l,m)=SUM(l,m)+zeff(l,m,ia,is)
             ENDDO
          ENDDO
       ENDDO
       DO l=1,3
          DO m=1,3
             SUM(l,m)=SUM(l,m)/REAL(ions1%nat,kind=real_8)
          ENDDO
       ENDDO
       !CDIR NODEP
       DO isa = 1,isa0
          is= isv(isa)
          ia= iav(isa)
          DO l=1,3
             DO m=1,3
                zeff(l,m,ia,is)=zeff(l,m,ia,is)-SUM(l,m)
             ENDDO
          ENDDO
       ENDDO
#else
       DO l=1,3
          DO m=1,3
             sum=0._real_8
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   sum=sum+zeff(l,m,ia,is)
                ENDDO
             ENDDO
             sum=sum/REAL(ions1%nat,kind=real_8)
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   zeff(l,m,ia,is)=zeff(l,m,ia,is)-sum
                ENDDO
             ENDDO
          ENDDO
       ENDDO
#endif
       CALL daxpy(9*maxsys%nax*ions1%nsp,1._real_8,zeff,1,zeffav,1)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt(' POLARPROJ',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE polarproj
  ! ==================================================================
  SUBROUTINE give_scr_polarproj(lpolarproj,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lpolarproj
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

! ==--------------------------------------------------------------==

    CALL give_scr_clinbcg(lpolarproj,tag,nstate)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_polarproj
  ! ==================================================================

END MODULE polarise_utils

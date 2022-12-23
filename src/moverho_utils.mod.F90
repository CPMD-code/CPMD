MODULE moverho_utils
  USE atrho_utils,                     ONLY: atdens
  USE atwf,                            ONLY: dmovrmix,&
                                             tmovr
  USE coor,                            ONLY: tau0
  USE cppt,                            ONLY: hg,&
                                             indz,&
                                             inyh,&
                                             nzh
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE fitpack_utils,                   ONLY: curv2
  USE geq0mod,                         ONLY: geq0
  USE gvec,                            ONLY: gvec_com
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: parai,&
                                             paral
  USE prmem_utils,                     ONLY: prmem
  USE qspl,                            ONLY: ggnh,&
                                             nsplpo
  USE ropt,                            ONLY: ropt_mod
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             maxsys,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: moverho
  PUBLIC :: give_scr_moverho

CONTAINS

  ! ==================================================================
  SUBROUTINE moverho(rhoe,psi)
    ! ==--------------------------------------------------------------==
    ! ==  EXPANDS THE DENSITY IN AN ATOM CENTERED BASIS OF STO        ==
    ! ==  MOVES THE EXPANDED DENSITY ACCORDING TO THE MOVEMENT        ==
    ! ==  OF THE ATOMS                                                ==
    ! ==  ROUTINE ASSUMES THAT SETBASIS HAS BEEN CALLED ONCE          ==
    ! ==  AND PHFAC HAS BEEN CALLED FOR THE TAUP                      ==
    ! ==  CONFIGURATION (AND STORED IN EI1(IA,IG) ETC).               ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: rhoe(fpar%nnr1,*)
    COMPLEX(real_8)                          :: psi(maxfftn)

    CHARACTER(*), PARAMETER                  :: procedureN = 'moverho'
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) ,&
                                                zzero = (0._real_8,0._real_8)

    COMPLEX(real_8)                          :: ctem1, ctem2, ctem3, ctep1, &
                                                ctep2, ctep3, ei10, ei20, &
                                                ei30, svtmpm, svtmpp
    COMPLEX(real_8), ALLOCATABLE             :: eig1(:,:), rhog(:), tsfac(:)
    INTEGER                                  :: i, ia, ierr, ig, ir, is, isa, &
                                                isa0, isub, j, k, nh1, nh2, &
                                                nh3, nn, nnx
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: ar, rstmp, rsum, rsum1, sa1, &
                                                sa2, sa3, vol
    REAL(real_8), ALLOCATABLE, SAVE          :: datom(:,:,:)

! ==--------------------------------------------------------------==

    CALL tiset('   MOVERHO',isub)
    CALL setfftn(0)
    ! ==--------------------------------------------------------------==
    ! TRANSFORM THE DENSITY TO G SPACE
    nnx=fpar%nnr1*clsd%nlsd
    ! First call, allocate DATOM.
    IF (ifirst.EQ.0) THEN
       ifirst=1
       ALLOCATE(datom(nsplpo,2,ions1%nsp),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (paral%parent) CALL prmem('   MOVERHO')
       DO is=1,ions1%nsp
          CALL atdens(is,datom(1,1,is))
       ENDDO
    ENDIF
    ! TODO align for BG
    nn = MAX(spar%nr1s,spar%nr2s,spar%nr3s)
    ALLOCATE(eig1(3,2*nn-1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(rhog(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(tsfac(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    vol=1._real_8/SQRT(parm%omega)
    nh1=spar%nr1s/2+1
    nh2=spar%nr2s/2+1
    nh3=spar%nr3s/2+1
    IF (.NOT.cntl%tlsd) THEN
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          psi(ir)=CMPLX(rhoe(ir,1),0._real_8,kind=real_8)
       ENDDO
       CALL  fwfftn(psi,.FALSE.,parai%allgrp)
       !$omp parallel do private(IG)
       DO ig=1,ncpw%nhg
          rhog(ig)=psi(nzh(ig))
       ENDDO
       isa0=0
       DO is=1,ions1%nsp
          CALL zeroing(tsfac)!,nhg)
          DO ia=1,ions0%na(is)
             ! Generate structure factor for TAU0 config.
             sa1=(gvec_com%b1(1)*tau0(1,ia,is)+gvec_com%b1(2)*tau0(2,ia,is)+&
                  gvec_com%b1(3)*tau0(3,ia,is))*parm%tpiba
             sa2=(gvec_com%b2(1)*tau0(1,ia,is)+gvec_com%b2(2)*tau0(2,ia,is)+&
                  gvec_com%b2(3)*tau0(3,ia,is))*parm%tpiba
             sa3=(gvec_com%b3(1)*tau0(1,ia,is)+gvec_com%b3(2)*tau0(2,ia,is)+&
                  gvec_com%b3(3)*tau0(3,ia,is))*parm%tpiba
             eig1(1,1)=CMPLX(1._real_8,0._real_8,kind=real_8)
             eig1(2,1)=CMPLX(1._real_8,0._real_8,kind=real_8)
             eig1(3,1)=CMPLX(1._real_8,0._real_8,kind=real_8)
             ctep1=CMPLX(COS(sa1),-SIN(sa1),kind=real_8)
             ctep2=CMPLX(COS(sa2),-SIN(sa2),kind=real_8)
             ctep3=CMPLX(COS(sa3),-SIN(sa3),kind=real_8)
             ctem1=CONJG(ctep1)
             ctem2=CONJG(ctep2)
             ctem3=CONJG(ctep3)

             svtmpp=ctep1
             svtmpm=ctem1
             DO i=2,spar%nr1s
                ei1(isa,i)=svtmpp
                svtmpp=svtmpp*ctep1
                ei1(isa,spar%nr1s+i-1)=svtmpm
                svtmpm=svtmpm*ctem1
             ENDDO

             svtmpp=ctep2
             svtmpm=ctem2
             DO j=2,spar%nr2s
                ei2(isa,j)=svtmpp
                svtmpp=svtmpp*ctep2
                ei2(isa,spar%nr2s+j-1)=svtmpm
                svtmpm=svtmpm*ctem2
             ENDDO

             svtmpp=ctep3
             svtmpm=ctem3
             DO k=2,spar%nr3s
                ei3(isa,k)=svtmpp
                svtmpp=svtmpp*ctep3
                ei3(isa,spar%nr3s+k-1)=svtmpm
                svtmpm=svtmpm*ctem3
             ENDDO

             ei10=CONJG(eig1(1,nh1))
             ei20=CONJG(eig1(2,nh2))
             ei30=CONJG(eig1(3,nh3))
             !$omp parallel do private(I)
#ifdef __SR8000
             !poption parallel, tlocal(I)
#endif 
             DO i=1,2*spar%nr1s-1
                eig1(1,i)=eig1(1,i)*ei10
             ENDDO
             !$omp parallel do private(J)
#ifdef __SR8000
             !poption parallel, tlocal(J)
#endif 
             DO j=1,2*spar%nr2s-1
                eig1(2,j)=eig1(2,j)*ei20
             ENDDO
             !$omp parallel do private(K)
#ifdef __SR8000
             !poption parallel, tlocal(K)
#endif 
             DO k=1,2*spar%nr3s-1
                eig1(3,k)=eig1(3,k)*ei30
             ENDDO
             isa=isa0+ia
             !$omp parallel do private(IG)
#ifdef __SR8000
             !poption parallel, tlocal(IG)
#endif 
             DO ig=1,ncpw%nhg
                tsfac(ig)=tsfac(ig)+&
                     ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))-&
                     eig1(1,inyh(1,ig))*eig1(2,inyh(2,ig))*&
                     eig1(3,inyh(3,ig))
             ENDDO
          ENDDO             ! End loop over IA
          isa0=isa0+ions0%na(is)
          ! Construct difference in atomic densities arising
          ! from motion of ions
          !$omp parallel do private(IG,AR)
          DO ig=1,ncpw%nhg
             ar=curv2(hg(ig),nsplpo,ggnh(1),datom(1,1,is),datom(1,2,is),&
                  0._real_8)
             rhog(ig)=rhog(ig) + dmovrmix*ar*vol*tsfac(ig)
          ENDDO
       ENDDO                 ! End loop over IS
       ! FFT to real space 
       CALL zeroing(psi)!,maxfft)
       !$omp parallel do private(IG)
       !CDIR NODEP
       DO ig=2,ncpw%nhg
          psi(nzh(ig))=rhog(ig)
          psi(indz(ig))=CONJG(rhog(ig))
       ENDDO
       IF (geq0) psi(nzh(1))=rhog(1)
       CALL  invfftn(psi,.FALSE.,parai%allgrp)
       rsum1=0._real_8
       !$omp parallel do private(IR,RSTMP) reduction(+:RSUM1)
       DO ir=1,fpar%nnr1
          rstmp=REAL(psi(ir))
          rhoe(ir,1)=rstmp
          rsum1 = rsum1 + rstmp
       ENDDO
       rsum1=rsum1*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
       rsum = 0.0_real_8
       IF (geq0) rsum = REAL(rhog(1)*parm%omega)
    ELSE
    ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(eig1,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(rhog,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(tsfac,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('   MOVERHO',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE moverho
  ! ==================================================================
  SUBROUTINE give_scr_moverho(lmoverho,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lmoverho
    CHARACTER(len=30)                        :: tag

    IF (ropt_mod%modens.OR.tmovr) THEN
       lmoverho=2*ncpw%nhg+       & ! RHOG
            2*ncpw%nhg+             & ! TSFAC
            2*3*(2*MAX(spar%nr1s,spar%nr2s,spar%nr3s)-1) ! EIG1
       lmoverho=MAX(lmoverho,2*maxsys%mmaxx)+10
       tag='4*NHG+6*MAX(NR1S,NR2S,NR3S)'
    ELSE
       lmoverho=0
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_moverho
  ! ==================================================================

END MODULE moverho_utils

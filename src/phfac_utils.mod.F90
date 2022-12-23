#if defined(__SR11000)
!option OPT(O(ss))
#endif

MODULE phfac_utils
  USE cppt,                            ONLY: inyh
  USE error_handling,                  ONLY: stopgm
  USE gvec,                            ONLY: gvec_com
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: eigkr,&
                                             eikr,&
                                             rk
  USE kpts,                            ONLY: tkpts
  USE parac,                           ONLY: paral
  USE prmem_utils,                     ONLY: prmem
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigr,&
                                             eigrb,&
                                             natx
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             kpbeg,&
                                             ncpw,&
                                             nkpbl,&
                                             nkpt,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: phfac
  PUBLIC :: calc_eigkr

CONTAINS

  ! ==================================================================
  SUBROUTINE phfac(tau0)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'phfac'

    COMPLEX(real_8)                          :: ctem1, ctem2, ctem3, ctep1, &
                                                ctep2, ctep3, ei10, ei20, &
                                                ei30, svtmpm, svtmpp, zsum
    INTEGER                                  :: i, ia, ierr, ig, ik, ikk, &
                                                ikpt, is, isa, isub, j, k, &
                                                nh1, nh2, nh3
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: ar1, ar2, ar3, sum, sum1, &
                                                sum2, sum3

! ==--------------------------------------------------------------==

    IF (spar%nr1s.LT.3) THEN
       CALL stopgm('PHFAC',' PHFAC: NR1 TOO SMALL ',& 
            __LINE__,__FILE__)
    ENDIF
    IF (spar%nr2s.LT.3) THEN
       CALL stopgm('PHFAC',' PHFAC: NR2 TOO SMALL ',& 
            __LINE__,__FILE__)
    ENDIF
    IF (spar%nr3s.LT.3) THEN
       CALL stopgm('PHFAC',' PHFAC: NR3 TOO SMALL ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tiset('     PHFAC',isub)
    IF (ifirst.EQ.0) THEN
       IF (cntl%bigmem) THEN
          ALLOCATE(eigrb(ncpw%nhg,ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ALLOCATE(eigr(ncpw%ngw,ions1%nat,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! FIXME deallocate missing
       IF (tkpts%tkpnt) THEN
          ALLOCATE(eikr(nkpt%nkpts,ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)

          ALLOCATE(eigkr(nkpt%ngwk,ions1%nat,nkpt%nkpnt),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ALLOCATE(ei1(natx,(2*spar%nr1s-1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ei2(natx,(2*spar%nr2s-1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ei3(natx,(2*spar%nr3s-1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ifirst  = 1
       IF (paral%parent) THEN
          CALL prmem('     PHFAC')
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    nh1=spar%nr1s/2+1
    nh2=spar%nr2s/2+1
    nh3=spar%nr3s/2+1
    !$omp parallel do private(ISA,IA,IS,SUM1,SUM2,SUM3,AR1,AR2,AR3) &
    !$omp private(CTEP1,CTEP2,CTEP3,CTEM1,CTEM2,CTEM3) &
    !$omp private(SVTMPP,SVTMPM,I,J,K,EI10,EI20,EI30) &
    !$omp shared(NH1,NH2,NH3)
    DO isa=1,ions1%nat
       ia=iatpt(1,isa)
       is=iatpt(2,isa)
       sum1=gvec_com%b1(1)*tau0(1,ia,is)+gvec_com%b1(2)*tau0(2,ia,is)+gvec_com%b1(3)*tau0(3,ia,is)
       sum2=gvec_com%b2(1)*tau0(1,ia,is)+gvec_com%b2(2)*tau0(2,ia,is)+gvec_com%b2(3)*tau0(3,ia,is)
       sum3=gvec_com%b3(1)*tau0(1,ia,is)+gvec_com%b3(2)*tau0(2,ia,is)+gvec_com%b3(3)*tau0(3,ia,is)
       ar1=parm%tpiba*sum1
       ar2=parm%tpiba*sum2
       ar3=parm%tpiba*sum3
       ei1(isa,1)=CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
       ei2(isa,1)=CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
       ei3(isa,1)=CMPLX(1.0_real_8,0.0_real_8,kind=real_8)
       ctep1=CMPLX(COS(ar1),-SIN(ar1),kind=real_8)
       ctep2=CMPLX(COS(ar2),-SIN(ar2),kind=real_8)
       ctep3=CMPLX(COS(ar3),-SIN(ar3),kind=real_8)
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

       ei10=1._real_8/ei1(isa,nh1)
       ei20=1._real_8/ei2(isa,nh2)
       ei30=1._real_8/ei3(isa,nh3)
#ifdef __SR8000
       !poption parallel, tlocal(I)
#endif 
       DO i=1,2*spar%nr1s-1
          ei1(isa,i)=ei1(isa,i)*ei10
       ENDDO
#ifdef __SR8000
       !poption parallel, tlocal(J)
#endif 
       DO j=1,2*spar%nr2s-1
          ei2(isa,j)=ei2(isa,j)*ei20
       ENDDO
#ifdef __SR8000
       !poption parallel, tlocal(K)
#endif 
       DO k=1,2*spar%nr3s-1
          ei3(isa,k)=ei3(isa,k)*ei30
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
#ifdef __SR8000 
    !poption parallel 
#endif 
#ifdef _vpp_
    !OCL NOALIAS
#endif
    !$omp parallel do private (ISA,IG)
    DO isa=1,ions1%nat
       DO ig=1,ncpw%ngw
          eigr(ig,isa,1)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
               ei3(isa,inyh(3,ig))
       ENDDO
    ENDDO
    IF (cntl%bigmem) THEN
#ifdef __SR8000 
       !poption parallel 
#endif 
#ifdef _vpp_
       !OCL NOALIAS
#endif
       !$omp parallel do private (ISA,IG)
       DO isa=1,ions1%nat
          DO ig=1,ncpw%nhg
             eigrb(ig,isa)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*ei3(&
                  isa,inyh(3,ig))
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (tkpts%tkpnt) THEN
       DO ikpt=1,nkpt%nblkp
1         DO ik=1,nkpbl(ikpt)
             ikk=kpbeg(ikpt)+ik
             !$omp parallel do private(ISA,IA,IS,SUM1,SUM2,SUM3,SUM,ZSUM,IG)
             DO isa=1,ions1%nat
                ia=iatpt(1,isa)
                is=iatpt(2,isa)
                sum1=gvec_com%b1(1)*tau0(1,ia,is)+gvec_com%b1(2)*tau0(2,ia,is)+gvec_com%b1(3)*tau0(3,&
                     ia,is)
                sum2=gvec_com%b2(1)*tau0(1,ia,is)+gvec_com%b2(2)*tau0(2,ia,is)+gvec_com%b2(3)*tau0(3,&
                     ia,is)
                sum3=gvec_com%b3(1)*tau0(1,ia,is)+gvec_com%b3(2)*tau0(2,ia,is)+gvec_com%b3(3)*tau0(3,&
                     ia,is)
                sum=rk(1,ikk)*sum1+rk(2,ikk)*sum2+rk(3,ikk)*sum3
                zsum=CMPLX(COS(sum),SIN(sum),kind=real_8)
                eikr(ikk,isa)=zsum
                DO ig=1,ncpw%ngw
                   eigkr(ig,isa,ik)=eigr(ig,isa,1)*zsum
                   eigkr(ig+ncpw%ngw,isa,ik)=CONJG(eigr(ig,isa,1))*zsum
                ENDDO
             ENDDO
          ENDDO
          IF (tkpts%tkblock) THEN
             CALL wkpt_swap(eigkr,1,ikpt,'EIGKR')
          ENDIF
       ENDDO
    ELSE
       ! kpt..not to duplicate code.
       ! NB: it is impossible to achieve this using Fortran 90, because
       ! EIGKR has rank 3, and EIGR has rank 2. Even if 
       ! EIGKR has dimensions (:,:,1) it cannot point to (:,:)
       ! To refactor this we will have to either make EIGR rank 3 with last dim = 1,
       ! or allocate memory for EIGKR and copy...
       ! 
       ! semantically this should be the following
       ! EIGKR => EIGR
       ! in principle, EIGKR should be just removed (?), and EIGR should have rank 3
       ! where last dim is NKPNT

       ! TODO: fix this!
       eigkr => eigr
    ENDIF
    CALL tihalt('     PHFAC',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE phfac
  ! ==================================================================
  SUBROUTINE calc_eigkr(ikpt)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES EIGKR (USED WITH THE OPTION TKBLOCK)              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ikpt

    COMPLEX(real_8)                          :: zsum
    INTEGER                                  :: ia, ig, ik, ikind, is, isa, &
                                                isub

    CALL tiset('CALC_EIGKR',isub)
    DO ikind=1,nkpbl(ikpt)
       ik=kpbeg(ikpt)+ikind
       isa=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             isa=isa+1
             zsum=eikr(ik,isa)
             DO ig=1,ncpw%ngw
                eigkr(ig,isa,ikind)=eigr(ig,isa,1)*zsum
                eigkr(ig+ncpw%ngw,isa,ikind)=CONJG(eigr(ig,isa,1))*zsum
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL tihalt('CALC_EIGKR',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calc_eigkr
  ! ==================================================================

END MODULE phfac_utils

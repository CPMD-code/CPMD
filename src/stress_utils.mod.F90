#include "cpmd_global.h"

MODULE stress_utils
  USE cnst,                            ONLY: pi
  USE cppt,                            ONLY: gk,&
                                             hg
  USE cvan,                            ONLY: dvan
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kdp,                             ONLY: bkdp,&
                                             pkdp,&
                                             wkdp,&
                                             xkdp
  USE kdp_stress_kin_utils,            ONLY: kdp_stress_kin
  USE kdpc,                            ONLY: nkdp,&
                                             tkdp
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: hgkm,&
                                             hgkp,&
                                             rk,&
                                             wk
  USE kpts,                            ONLY: tkpts
  USE metr,                            ONLY: metr_com
  USE nlcc,                            ONLY: corel
  USE nlps,                            ONLY: imagp,&
                                             ndfnl,&
                                             nghtol,&
                                             nlps_com,&
                                             wsg
  USE nlsm1_s_utils,                   ONLY: nlsm1_s
  USE parac,                           ONLY: parai,&
                                             paral
  USE pbc_utils,                       ONLY: pbc
  USE prcp,                            ONLY: prcp_com
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE putbet_utils,                    ONLY: putbet
  USE ragg,                            ONLY: raggio
  USE rnlsm_utils,                     ONLY: give_scr_rnlsm,&
                                             rnlsm
  USE ropt,                            ONLY: iteropt
  USE sfac,                            ONLY: fnl
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE special_functions,               ONLY: cp_erfc
  USE spin,                            ONLY: clsd
  USE str2,                            ONLY: becs,&
                                             betaa,&
                                             dqg,&
                                             drhovg,&
                                             gagk,&
                                             gagkm,&
                                             gagkp
  USE strs,                            ONLY: &
       alpha, beta, decc, degc, dehc, deht, dekin, delta, denl, depst, desr, &
       dexc, dlam
  USE system,                          ONLY: cntl,&
                                             kpbeg,&
                                             maxsys,&
                                             ncpw,&
                                             nkpbl,&
                                             nkpt,&
                                             parap,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: unitmx
  USE vdwcmod,                         ONLY: vdwl,&
                                             vdwr
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: stress
  PUBLIC :: give_scr_stress

CONTAINS

  ! ==================================================================
  SUBROUTINE stress(c0,tau0,f,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE STRESS TENSOR FROM THE DFT TOTAL ENERGY                 ==
    ! ==  (THIS IS ESSENTIALLY A COPY OF PAOLO FOCHERS CODE)          ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: tau0(:,:,:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate,nkpt%nkpts)
    COMPLEX(real_8) :: c0(nkpt%ngwk,nstate,nkpt%nkpnt)

    CHARACTER(*), PARAMETER                  :: procedureN = 'stress'
    REAL(real_8), PARAMETER                  :: argmax = 20._real_8 , &
                                                deltakin = 1.e-10_real_8 

    COMPLEX(real_8), ALLOCATABLE             :: auxc(:), ctmp(:), qg(:)
    INTEGER :: i, ierr, ig, ii, ikind, ikk, ikpt, il_auxc, il_ctmp, il_gam, &
      il_qg, infm, is, isa0, ishft, isub, iv, ix, iy, iz, j, jv, k, kbeg, &
      kend, ki, kinc, kj, kk, l, l2, lbecs, li, lj, m, nkpoint
    INTEGER, SAVE                            :: istart = 0
    LOGICAL                                  :: tzero
    REAL(real_8) :: addesr, addpre, arg, arg1, arg2, erre2, esrtzero, fxx, &
      rckj, rckk, repand, rlm, rxlm(3), sgc(6), sum, topia, weight, xlm, &
      xskin, ylm, zlm, zv2, xlm_, ylm_, zlm_
    REAL(real_8), ALLOCATABLE                :: gam(:)
    REAL(real_8), EXTERNAL                   :: ddot

#if defined(__VECTOR)
    REAL(real_8)                             :: t1,t2,t3,t4,t5,t6
#endif 
    ! ==--------------------------------------------------------------==
    CALL tiset('    STRESS',isub)
    IF (cntl%tfdist) CALL stopgm('STRESS','TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (corel%tinlc) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' THIS NEEDS FIRST SOME PROGRAMING ! '
       CALL stopgm('STRESS',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! 07/2008 AK
    ! Stress tensor calculation is broken when using vanderbilt ultra-soft
    ! pseudopotentials and probably will never be fixed. use Q-E instead.
    ! 08/2017 MB-TI - The above statement is wrong. Fixed on 28/08/2017. 
  ! IF (pslo_com%tivan) THEN
  !    CALL stopgm('STRESS',&
  !         'STRESS TENSOR IS BROKEN FOR VANDERBILT USPP',& 
  !         __LINE__,__FILE__)
  ! ENDIF
    ! ==--------------------------------------------------------------==
    ! SCR partition
    il_auxc=2*nkpt%ngwk*maxsys%nax         ! NLSM1_S
    il_gam =imagp*maxsys%nax*MAX(parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1,nstate) ! NLSM1_S
    IF (pslo_com%tivan) THEN
       il_qg   =2*ncpw%nhg         ! DRHOV
       il_ctmp =2*ncpw%nhg         ! DRHOV
    ELSE
       il_qg   =0
       il_ctmp =0
    ENDIF
    ALLOCATE(auxc(nkpt%ngwk * maxsys%nax),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    call zeroing( auxc )
    ALLOCATE(gam(il_gam),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    call zeroing( gam )
    ALLOCATE(qg(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ctmp(ncpw%nhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! Initialize
    CALL unitmx(delta,3)
    CALL zeroing(degc)!,6)
    CALL zeroing(dlam)!,6)
    IF (.NOT.corel%tinlc) CALL zeroing(decc)!,6)
    ! Initialize memory            
    IF (istart.EQ.0) THEN
       istart=1
       alpha(1)=1
       alpha(2)=2
       alpha(3)=3
       alpha(4)=2
       alpha(5)=3
       alpha(6)=3
       beta(1) =1
       beta(2) =1
       beta(3) =1
       beta(4) =2
       beta(5) =2
       beta(6) =3

       ! NSTATE=NST12(MEPOS,2)-NST12(MEPOS,1)+1  !=NDFN
       ! This works also for parallel
       lbecs=imagp*ions1%nat*maxsys%nhxs*6*(parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,1)+1)*nkpt%nkpnt
       IF (lbecs.LE.0) lbecs=1
       ALLOCATE(becs(imagp,ions1%nat,maxsys%nhxs,6,ndfnl,nkpt%nkpnt),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(betaa(nkpt%ngwk,maxsys%nhxs,6,maxsys%nsx,nkpt%nkpnt),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(gagk(ncpw%nhg,6),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       IF (tkpts%tkpnt) THEN
          ALLOCATE(gagkp(ncpw%ngw,6,nkpt%nkpnt),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(gagkm(ncpw%ngw,6,nkpt%nkpnt),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF (pslo_com%tivan) THEN
          ALLOCATE(drhovg(ncpw%nhg,6*clsd%nlsd),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       IF (paral%parent) CALL prmem('    STRESS')
    ENDIF
    ! Compute GAGK = G_A * G_K
#if defined(__VECTOR)
    !OCL NOALIAS
#ifdef __SR8000
    !poption parallel, tlocal(IG)
#endif 
    !$omp parallel do private(IG)
    DO ig=1,ncpw%nhg
       gagk(ig,1) = gk(alpha(1),ig)*gk(beta(1),ig)*parm%tpiba2
       gagk(ig,2) = gk(alpha(2),ig)*gk(beta(2),ig)*parm%tpiba2
       gagk(ig,3) = gk(alpha(3),ig)*gk(beta(3),ig)*parm%tpiba2
       gagk(ig,4) = gk(alpha(4),ig)*gk(beta(4),ig)*parm%tpiba2
       gagk(ig,5) = gk(alpha(5),ig)*gk(beta(5),ig)*parm%tpiba2
       gagk(ig,6) = gk(alpha(6),ig)*gk(beta(6),ig)*parm%tpiba2
    ENDDO
#else 
    !$omp parallel do private(KK,IG) schedule(static)
    DO ig=1,ncpw%nhg
       DO kk=1,6
          gagk(ig,kk) = gk(alpha(kk),ig)*gk(beta(kk),ig)*parm%tpiba2
       ENDDO
    ENDDO
#endif 
    ! Initialize the stress arrays
    CALL zeroing(dekin)!,6)
    dehc=0.0_real_8
    CALL zeroing(deht )!,6)
    CALL zeroing(dexc )!,6)
    CALL zeroing(desr )!,6)
    depst=0.0_real_8
    CALL zeroing(denl)!, 6)
    CALL zeroing(degc)!, 6)
    CALL zeroing(dlam)!, 6)
    IF (.NOT.corel%tinlc) CALL zeroing(decc)!,6)
    IF (.NOT.vdwl%vdwc)  CALL zeroing(vdwr%devdw)!,6)
    IF (.NOT.vdwl%vdwd)  CALL zeroing(vdwr%devdw)!,6)
    ! ==--------------------------------------------------------------==
    ! == Loop over k points                                           ==
    ! ==--------------------------------------------------------------==
    CALL inq_swap(kbeg,kend,kinc)
    DO ikpt=kbeg,kend,kinc
       nkpoint=nkpbl(ikpt)
       IF (tkpts%tkblock) THEN
          CALL rkpt_swap(c0,nstate,ikpt,&
               'HGKP HGKM MASKGW EIGKR TWNL C0')
       ENDIF
       IF (tkpts%tkpnt) THEN
          ! Initialisation of GAGKP and GAGKM
#if defined(__VECTOR) || defined(__PRIMEHPC)
          DO ikind=1,nkpoint
             ikk=kpbeg(ikpt)+ikind
             !OCL NOALIAS
#ifdef __SR8000
             !poption parallel, tlocal(KK,IG)
#endif 
             !$omp parallel do private(KK,IG) __COLLAPSE2
             DO kk=1,6
                DO ig=1,ncpw%ngw
                   gagkp(ig,kk,ikind) =&
                        (rk(alpha(kk),ikk) + gk(alpha(kk),ig)) *&
                        (rk(beta (kk),ikk) + gk(beta (kk),ig)) * parm%tpiba2
                   gagkm(ig,kk,ikind) =&
                        (rk(alpha(kk),ikk) - gk(alpha(kk),ig)) *&
                        (rk(beta (kk),ikk) - gk(beta (kk),ig)) * parm%tpiba2
                ENDDO
             ENDDO
          ENDDO
#else
          !OCL NOALIAS
          !$omp  parallel do private(IG,KK,IKIND,IKK) &
          !$omp  schedule(static)
          DO ig=1,ncpw%ngw
             DO ikind=1,nkpoint
                ikk=kpbeg(ikpt)+ikind
                DO kk=1,6
                   gagkp(ig,kk,ikind) =&
                        (rk(alpha(kk),ikk) + gk(alpha(kk),ig)) *&
                        (rk(beta (kk),ikk) + gk(beta (kk),ig)) * parm%tpiba2
                   gagkm(ig,kk,ikind) =&
                        (rk(alpha(kk),ikk) - gk(alpha(kk),ig)) *&
                        (rk(beta (kk),ikk) - gk(beta (kk),ig)) * parm%tpiba2
                ENDDO
             ENDDO
          ENDDO
#endif
       ENDIF
       ! Initialize BETAA(deriv. of TWNL)
       CALL zeroing(betaa)!(1,1,1,1,1),nkpt%ngwk*maxsys%nhxs*maxsys%nsx*6*nkpoint)
       CALL putbet(ikpt)
       DO ikind=1,nkpoint
          ! Initialize FNL
          CALL rnlsm(c0(:,:,ikind),nstate,&
               1,ikind,.FALSE.)
          ! Initialize BECS(deriv. of FNL), DRHOVG
          ! (deriv. of Vanderbilt aug. charge)
          CALL nlsm1_s(c0(1,1,ikind),auxc,gam,nstate,ikind)
          isa0=0
          DO is=1,ions1%nsp
             DO kk=1,6
                IF (delta(alpha(kk),beta(kk)).GT.0.5_real_8) THEN
                   DO iv=1,nlps_com%ngh(is)
                      DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                         ii=i-parap%nst12(parai%mepos,1)+1
                         CALL daxpy(imagp*ions0%na(is),-0.5_real_8,&
                              fnl (1,isa0+1,iv,i,ikind),1,&
                              becs(1,isa0+1,iv,kk,ii,ikind),1)
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
       IF (pslo_com%tivan) CALL drhov(nstate,qg,ctmp)
       ! ==------------------------------------------------------------==
       ! ==  KINETIC ENERGY CONTRIBUTION                               ==
       ! ==------------------------------------------------------------==
       ! 
       ! k.p kinetic part o the stress (jk)
       IF (tkdp) THEN
          CALL kdp_stress_kin(nstate,c0,gagk,wkdp,xkdp,pkdp,bkdp,nkdp)
          GOTO 1111
       ENDIF
       ! 
       DO ikind=1,nkpoint
          ikk=kpbeg(ikpt)+ikind
          DO i=1,nstate
             IF (f(i,ikk).NE.0._real_8) THEN
                CALL zeroing(sgc)!,6)
                !OCL NOALIAS
#ifdef __SR8000
                !poption parallel, tlocal(IG)
#endif
#if defined(__VECTOR)
                !$omp parallel do private(IG)
#else
                !$omp parallel do private(IG) schedule(static)
#endif
                DO ig=1,nkpt%ngwk
                   dqg(ig,1)=REAL(c0(ig,i,ikind))*REAL(c0(ig,i,ikind))+&
                        AIMAG(c0(ig,i,ikind))*AIMAG(c0(ig,i,ikind))
                ENDDO
                IF (prcp_com%akin.GT.deltakin) THEN
                   xskin=1._real_8/prcp_com%gskin
                   topia=2._real_8*prcp_com%gakin/SQRT(pi)
                   IF (tkpts%tkpnt) THEN
                      !OCL NOALIAS
#ifdef __SR8000
                      !poption parallel, tlocal(IG,ARG)
#endif
#if defined(__VECTOR)
                      !$omp parallel do private(IG,ARG) shared(TOPIA,XSKIN)
                      DO ig=1,ncpw%ngw
                         arg=(hgkp(ig,ikind)-prcp_com%gckin)*xskin
                         dqg(ig,1)=dqg(ig,1)*&
                              (1._real_8+topia*xskin*EXP(-arg*arg))
                      ENDDO
                      !OCL NOALIAS
                      !OCL NOVREC
#ifdef __SR8000
                      !poption parallel, tlocal(IG,ARG)
#endif
                      !$omp parallel do private(IG,ARG) shared(TOPIA,XSKIN)
                      DO ig=1,ncpw%ngw
                         arg=(hgkm(ig,ikind)-prcp_com%gckin)*xskin
                         dqg(ncpw%ngw+ig,1)=dqg(ncpw%ngw+ig,1)&
                              *(1._real_8+topia*xskin*EXP(-arg*arg))
                      ENDDO
#else
                      !$omp parallel do private(IG,ARG1,ARG2) shared(TOPIA,XSKIN) &
                      !$omp  schedule(static)
                      DO ig=1,ncpw%ngw
                         arg1=(hgkp(ig,ikind)-prcp_com%gckin)*xskin
                         arg2=(hgkm(ig,ikind)-prcp_com%gckin)*xskin

                         dqg(ig,1)=dqg(ig,1)*&
                              (1._real_8+topia*xskin*EXP(-arg1*arg1))
                         dqg(ncpw%ngw+ig,1)=dqg(ncpw%ngw+ig,1)&
                              *(1._real_8+topia*xskin*EXP(-arg2*arg2))
                      ENDDO
#endif
                   ELSE
                      !OCL NOALIAS
#ifdef __SR8000
                      !poption parallel, tlocal(IG,ARG)
#endif 
#if defined(__VECTOR)
                      !$omp parallel do private(IG,ARG) shared(TOPIA,XSKIN)
#else
                      !$omp parallel do private(IG,ARG) shared(TOPIA,XSKIN) schedule(static)
#endif
                      DO ig=1,ncpw%ngw
                         arg=(hg(ig)-prcp_com%gckin)*xskin
                         dqg(ig,1)=dqg(ig,1)*(1._real_8+topia*xskin*EXP(-arg*&
                              arg))
                      ENDDO
                   ENDIF
                ENDIF
                IF (tkpts%tkpnt) THEN
#ifdef __SR8000
                   !poption parallel, tlocal(KK,IG)
#endif 
#if defined(__VECTOR)
                   !$omp parallel do private(KK,IG)
                   DO kk=1,6
                      DO ig=1,ncpw%ngw
                         sgc(kk)=sgc(kk)+REAL(dqg(ig,1)*gagkp(ig,kk,ikind))&
                              +REAL(dqg(ncpw%ngw+ig,1)*GAGKM(IG,KK,IKIND))
                      ENDDO
                   ENDDO
#else
                   !$omp parallel do private(KK,IG) reduction(+:SGC) schedule(static)
                   DO ig=1,ncpw%ngw
                      DO kk=1,6
                         sgc(kk)=sgc(kk)+REAL(dqg(ig,1)*gagkp(ig,kk,ikind))+&
                              REAL(dqg(ncpw%ngw+ig,1)*gagkm(ig,kk,ikind))
                      ENDDO
                   ENDDO
#endif
                ELSE
#if defined(__VECTOR)
                   t1=0.0_real_8
                   t2=0.0_real_8
                   t3=0.0_real_8
                   t4=0.0_real_8
                   t5=0.0_real_8
                   t6=0.0_real_8
                   !$omp parallel do private(IG) &
                   !$omp  reduction(+:T1,T2,T3,T4,T5,T6)
                   DO ig=1,ncpw%ngw
                      t1=t1+REAL(dqg(ig,1)*gagk(ig,1))
                      t2=t2+REAL(dqg(ig,1)*gagk(ig,2))
                      t3=t3+REAL(dqg(ig,1)*gagk(ig,3))
                      t4=t4+REAL(dqg(ig,1)*gagk(ig,4))
                      t5=t5+REAL(dqg(ig,1)*gagk(ig,5))
                      t6=t6+REAL(dqg(ig,1)*gagk(ig,6))
                   ENDDO
                   sgc(1)=sgc(1)+t1
                   sgc(2)=sgc(2)+t2
                   sgc(3)=sgc(3)+t3
                   sgc(4)=sgc(4)+t4
                   sgc(5)=sgc(5)+t5
                   sgc(6)=sgc(6)+t6
#else 
                   !$omp parallel do private(KK,IG) reduction(+:sgc)
                   DO ig=1,ncpw%ngw
                      DO kk=1,6
                         sgc(kk)=sgc(kk)+REAL(dqg(ig,1)*gagk(ig,kk))
                      ENDDO
                   ENDDO
#endif 
                ENDIF
                weight=wk(ikk)*f(i,ikk)
                CALL daxpy(6,weight,sgc(1),1,dekin(1),1)
             ENDIF
          ENDDO
       ENDDO
1111   CONTINUE
       ! ==------------------------------------------------------------==
       ! ==  NON-LOCAL CONTRIBUTION dENL/dH                            ==
       ! ==------------------------------------------------------------==
       DO ikind=1,nkpoint
          ikk=kpbeg(ikpt)+ikind
          isa0=0
          DO is=1,ions1%nsp
             IF (pslo_com%tvan(is)) THEN
                ! VANDERBILT PSEUDOPOTENTIAL
                DO iv=1,nlps_com%ngh(is)
                   DO jv=1,nlps_com%ngh(is)
                      DO kk=1,6
                         sum=0.0_real_8
                         DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                            ii=i-parap%nst12(parai%mepos,1)+1
                            IF (f(i,ikk).NE.0._real_8) THEN
                               weight=wk(ikk)*f(i,ikk)
                               sum=sum+weight*ddot(imagp*ions0%na(is),&
                                    becs(1,isa0+1,iv,kk,ii,ikind),1,&
                                    fnl (1,isa0+1,jv,i,ikind),1)
                            ENDIF
                         ENDDO
                         denl(kk) = denl(kk) + 2.0_real_8*sum*dvan(jv,iv,is)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF (sgpp1%tsgp(is)) THEN
                DO iv=1,nlps_com%ngh(is)
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   DO jv=1,nlps_com%ngh(is)
                      l2=nghtol(jv,is)+1
                      lj=sgpp2%lpval(jv,is)
                      IF (l.EQ.l2.AND.li.EQ.lj) THEN
                         kj=sgpp2%lfval(jv,is)
                         DO kk=1,6
                            sum=0.0_real_8
                            DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                               ii=i-parap%nst12(parai%mepos,1)+1
                               IF (f(i,ikk).NE.0._real_8) THEN
                                  weight=wk(ikk)*f(i,ikk)
                                  sum=sum+weight*ddot(imagp*ions0%na(is),&
                                       becs(1,isa0+1,iv,kk,ii,ikind),1,&
                                       fnl (1,isa0+1,jv,i,ikind),1)
                               ENDIF
                            ENDDO
                            denl(kk) = denl(kk) + 2._real_8*sum*sgpp2%hlsg(ki,kj,l,is)
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ELSE
                ! BHS AND RELATED
                DO iv=1,nlps_com%ngh(is)
                   DO kk=1,6
                      sum=0.0_real_8
                      DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                         ii=i-parap%nst12(parai%mepos,1)+1
                         IF (f(i,ikk).NE.0._real_8) THEN
                            weight=wk(ikk)*f(i,ikk)
                            sum=sum+weight*ddot(imagp*ions0%na(is),&
                                 becs(1,isa0+1,iv,kk,ii,ikind),1,&
                                 fnl (1,isa0+1,iv,i,ikind),1)
                         ENDIF
                      ENDDO
                      denl(kk) = denl(kk) + 2.0_real_8*sum*wsg(is,iv)
                   ENDDO
                ENDDO
             ENDIF
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    IF (tkpts%tkpnt) THEN
       CALL dscal(6,-1.0_real_8,dekin(1),1)
    ELSE IF (.NOT.tkdp) THEN
       CALL dscal(6,-2.0_real_8,dekin(1),1)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==  ESR CONTRIBUTION                                            ==
    ! ==--------------------------------------------------------------==
    IF (paral%parent) THEN
       DO k=1,ions1%nsp
          rckk = raggio(k)*raggio(k)
          DO j=k,ions1%nsp
             zv2 = ions0%zv(k)*ions0%zv(j)
             rckj = SQRT(rckk + raggio(j)*raggio(j))
             DO l=1,ions0%na(k)
                infm = 1
                IF (k.EQ.j)  infm=l
                DO m=infm,ions0%na(j)
                   IF (l.EQ.m.AND.k.EQ.j) THEN
                      xlm = 0.0_real_8
                      ylm = 0.0_real_8
                      zlm = 0.0_real_8
                      tzero=.TRUE.
                   ELSE
                      tzero=.FALSE.
                      xlm_ = tau0(1,l,k) - tau0(1,m,j)
                      ylm_ = tau0(2,l,k) - tau0(2,m,j)
                      zlm_ = tau0(3,l,k) - tau0(3,m,j)
                      CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
                   ENDIF
                   DO ix=-iteropt%iesr,iteropt%iesr
                      DO iy=-iteropt%iesr,iteropt%iesr
                         DO iz=-iteropt%iesr,iteropt%iesr
                            ishft=ix*ix+iy*iy+iz*iz
                            IF (.NOT.(tzero.AND.ishft.EQ.0)) THEN
                               rxlm(1)=xlm+ix*metr_com%ht(1,1)+iy*metr_com%ht(2,1)+iz*metr_com%ht(3,1)
                               rxlm(2)=ylm+ix*metr_com%ht(1,2)+iy*metr_com%ht(2,2)+iz*metr_com%ht(3,2)
                               rxlm(3)=zlm+ix*metr_com%ht(1,3)+iy*metr_com%ht(2,3)+iz*metr_com%ht(3,3)
                               erre2 = rxlm(1)**2 + rxlm(2)**2 + rxlm(3)**2
                               rlm = SQRT(erre2)
                               IF (tzero) THEN
                                  esrtzero=0.5_real_8
                               ELSE
                                  esrtzero=1._real_8
                               ENDIF
                               arg = rlm/rckj
                               IF (arg.LE.argmax) THEN! ADDESR,ADDPRE /= 0
                                  addesr = zv2*cp_erfc(arg)/rlm
                                  addpre=(2._real_8*zv2/SQRT(pi))&
                                       *EXP(-arg*arg)/rckj
                                  repand = esrtzero*(addesr+addpre)/erre2
                                  DO kk=1,6
                                     fxx = repand*rxlm(alpha(kk))*rxlm(beta(kk))
                                     desr(kk) = desr(kk) - fxx
                                  ENDDO
                               ENDIF
                            ENDIF
                         ENDDO! IZ
                      ENDDO! IY
                   ENDDO! IX
                ENDDO! M
             ENDDO    ! L
          ENDDO        ! J
       ENDDO            ! K
    ENDIF
    ! ==--------------------------------------------------------------==
    ! TODO maybe we can deallocate earlier?
    DEALLOCATE(auxc,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(gam,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(qg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ctmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tihalt('    STRESS',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE stress
  ! ==================================================================
  SUBROUTINE give_scr_stress(lstress,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lstress
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lrnlsm

! Variables
! ==--------------------------------------------------------------==
! NLSM1_S

    lstress=2*ncpw%ngw*maxsys%nax+&
         maxsys%nax*(parap%nst12(parai%mepos,2)-parap%nst12(parai%mepos,2)+1)
    ! DRHOV
    IF (pslo_com%tivan) lstress=MAX(lstress,4*ncpw%nhg)
    IF (tkpts%tkblock) THEN
       CALL give_scr_rnlsm(lrnlsm,tag,crge%n,.FALSE.)
       lstress=MAX(lstress,lrnlsm)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_stress
  ! ==================================================================

END MODULE stress_utils

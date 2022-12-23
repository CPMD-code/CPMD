MODULE ehpsi_utils
  USE cppt,                            ONLY: hg,&
                                             indzs,&
                                             nzhs
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: ngrm
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fint,                            ONLY: afnl,&
                                             alm,&
                                             biln,&
                                             cnl,&
                                             emk2,&
                                             fint1,&
                                             fint4,&
                                             fint5
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpclean_utils,                   ONLY: c_clean
  USE kpnt,                            ONLY: hgkm,&
                                             hgkp
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_prod,&
                                             mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlm,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai
  USE reshaper,                        ONLY: type_cast
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE spin,                            ONLY: clsd
  USE system,                          ONLY: &
       cntr, fpar, group, ncpw, nhx, nkpbl, nkpt, parap, parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  !  PUBLIC :: ehpsi
  PUBLIC :: ehpsi_c
  PUBLIC :: epot
  PUBLIC :: eg2
  PUBLIC :: calc_emk2
  PUBLIC :: calc_biln
  PUBLIC :: set_b2l
  PUBLIC :: change_b2l
  PUBLIC :: change_betap
  PUBLIC :: evpsi
  !public :: rcalc_biln

CONTAINS

  ! ==================================================================
  SUBROUTINE ehpsi_c(c0,c2,exp_vpot,psi,nstate,ikind,maxfft)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ACTION OF APPLYING e^{(-beta/P)H}  TO A WAVEFUNCTION C0,==
    ! ==  RETURNING THE ANSWER IN C2.                                 ==
    ! ==--------------------------------------------------------------==
    REAL(real_8), TARGET                     :: exp_vpot(fpar%nnr1)
    COMPLEX(real_8), TARGET                  :: psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(nkpt%ngwk,nstate), &
                                                c0(nkpt%ngwk,nstate)
    INTEGER                                  :: ikind, maxfft

    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) 

    INTEGER                                  :: i, ia, ib, ibb, ifft, ig, &
                                                isub, ix, ix1, ixp, lead, &
                                                nl2, nnrx, nsta
    REAL(real_8), POINTER                    :: exp_vpotx(:), psix(:)

! ==--------------------------------------------------------------==
! == INITIALIZE ELEC-FORCE ARRAYS                                 ==
! ==--------------------------------------------------------------==

    CALL tiset('   EHPSI_C',isub)
    ! ==--------------------------------------------------------------==
    IF (group%nogrp.EQ.1) THEN
       lead  = fpar%kr1s*parai%ngrays
       ifft  = 1
       nnrx=fpar%nnr1
    ELSE
       lead  = fpar%kr1s*ngrm
       ifft=2
       nnrx=fpar%krx*fpar%kr2s*fpar%kr3s
    ENDIF
    CALL type_cast(psi, SIZE(psi), psix)
    exp_vpotx => exp_vpot
    IF (group%nogrp.GT.1) THEN
       ! Summation of density within orbital split
       nl2=parap%nlink(group%nolist(group%nogrp))
       ix1=parap%nrxpl(parai%mepos,1)-parap%nrxpl(nl2,1)
       ! EXP(-betaP*0) = 1!!!
       !$omp parallel do private(I)
       DO i=1,nnrx
          psix(i)=1._real_8
       ENDDO
       DO ixp=1,parm%nr1
          ix=ix1+ixp
          CALL dcopy(fpar%kr2s*fpar%kr3s,exp_vpot(ixp),fpar%kr1,psix(ix),fpar%krx)
       ENDDO
       ! WARNING: COMBINE with MULTIPLICATION operation
       CALL mp_prod(psix,exp_vpotx,nnrx,group%meogrp)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == APPLY THE KINETIC ENERGY PROPAGATOR TO C0, RETURNING THE     ==
    ! == ANSWER IN PSI. ARRAY EMK2 CONTAINS e^{-beta g^2/4P}          ==
    ! ==--------------------------------------------------------------==
    !$omp parallel do private(I,IG)
    DO i=1,nstate
       DO ig = 1,nkpt%ngwk
          c2(ig,i) = emk2(ig,ikind)*c0(ig,i)
       ENDDO
    ENDDO
    IF (nlm.GT.0) THEN
       CALL ovlap2_c(nkpt%ngwk,nlm,nstate,afnl,cnl(1,1),c2(1,1))
       CALL mp_sum(afnl,2*nlm*nstate,parai%allgrp)
       CALL calc_biln(nstate,ikind)
    ENDIF
    ! ==--------------------------------------------------------------==
    DO ia=1,nstate,group%nogrp
       nsta=MIN((nstate-ia+1),group%nogrp)
       CALL zeroing(psi)!,maxfft*group%nogrp)
       DO ib=1,nsta
          i=ia+(ib-1)
          ibb=(ib-1)*lead
          !$omp parallel do private(IG)
          DO ig = 1 , ncpw%ngw
             psi(nzhs(ig)+ibb) = c2(ig,i)
             ! mb          ENDDO          ! useless to repeat
             ! mb          DO IG=1,NGW    ! twice the loop
             psi(indzs(ig)+ibb) = c2(ig+ncpw%ngw,i)
          ENDDO
          IF (geq0) psi(nzhs(1)+ibb) = c2(1,i)
       ENDDO
       ! ==------------------------------------------------------------==
       ! == INVERSE FFT PSI TO OBTAIN PSI ON REAL SPACE MESH           ==
       ! ==------------------------------------------------------------==
       IF (ifft.EQ.1) THEN
          CALL  invfftn(psi,.TRUE.,parai%allgrp)
       ELSE
          CALL stopgm("THIS","SG_INVFFT NOT AVAILABLE ANYMORE",& 
               __LINE__,__FILE__)
       ENDIF
       ! ==------------------------------------------------------------==
       ! == APPLY THE EXPONENTIATED POTENTIAL TO PSI                   ==
       ! == APPLY e^(-b/P V) IN REAL SPACE, AND RETURN ANSWER IN PSI.  ==
       ! == ARRAY POT CONTAINS e^(-b/P V(r)).                          ==
       ! ==------------------------------------------------------------==
       CALL evpsi(nnrx,psi,exp_vpotx)
       ! ==------------------------------------------------------------==
       ! == TRANSFORM BACK TO FOURIER SPACE                            ==
       ! ==------------------------------------------------------------==
       IF (ifft.EQ.1) THEN
          CALL  fwfftn(psi,.TRUE.,parai%allgrp)
       ELSE
          CALL stopgm("THIS","SG_FWFFT NOT AVAILABLE ANYMORE",& 
               __LINE__,__FILE__)
       ENDIF
       ! ==------------------------------------------------------------==
       DO ib=1,nsta
          i=ia+(ib-1)
          ibb=(ib-1)*lead
          !$omp parallel do private(IG)
          DO ig = 1 , ncpw%ngw
             c2(ig,i) = psi(nzhs(ig)+ibb)
             ! mb          ENDDO            ! useless to repeat
             ! mb          DO IG = 1 , NGW  ! twice the loop
             c2(ig+ncpw%ngw,i) = psi(indzs(ig)+ibb)
          ENDDO
          IF (geq0) c2(ncpw%ngw+1,i) = CMPLX(0._real_8,0._real_8,kind=real_8)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    ! == ADD NON-LOCAL PP PART                                        ==
    ! ==--------------------------------------------------------------==
    IF (nlm.GT.0) THEN
       CALL zgemm('N','N',nkpt%ngwk,nstate,nlm,zone,cnl(1,1),nkpt%ngwk,&
            biln(1,1,1),nlm,zone,c2(1,1),nkpt%ngwk)
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == APPLY THE KINETIC ENERGY PROPAGATOR TO C0, RETURNING THE     ==
    ! == ANSWER IN C2                                                 ==
    ! ==--------------------------------------------------------------==
    DO i=1,nstate
       !$omp parallel do private(IG)
       DO ig = 1,nkpt%ngwk
          c2(ig,i) = c2(ig,i)*emk2(ig,ikind)
       ENDDO
       IF (geq0) c2(ncpw%ngw+1,i) = CMPLX(0._real_8,0._real_8,kind=real_8)
    ENDDO
    CALL c_clean(c2,nstate,ikind)
    ! Reestablish initial conditions for the potential
    IF (group%nogrp.GT.1) THEN
       nl2=parap%nlink(group%nolist(group%nogrp))
       ix1=parap%nrxpl(parai%mepos,1)-parap%nrxpl(nl2,1)
       CALL dcopy(nnrx,exp_vpotx(1),1,psix(1),1)
       CALL zeroing(exp_vpotx)!,nnrx)
       DO ixp=1,parm%nr1
          ix=ix1+ixp
          CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,exp_vpot(ixp),fpar%kr1)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tihalt('   EHPSI_C',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ehpsi_c
  ! ==================================================================
  SUBROUTINE epot(vpot)
    ! ==--------------------------------------------------------------==
    ! ==  RETURNS THE EXPONENTIATED POTENTIAL IN VPOT                 ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vpot(fpar%nnr1,clsd%nlsd)

    INTEGER                                  :: ii, ir

    !$omp parallel do private(II,IR)
    DO ii = 1 , clsd%nlsd
       DO ir = 1 , fpar%nnr1
          vpot(ir,ii)=EXP(-fint1%betap*vpot(ir,ii))
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE epot
  ! ==================================================================
  SUBROUTINE eg2
    ! ==--------------------------------------------------------------==
    ! ==  RETURNS THE EXPONENTIATED KINETIC ENERGY IN G-SPACE         ==
    ! ==  IN ARRAY EMK2                                               ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'eg2'

    INTEGER                                  :: ierr, ig, ikpt, kbeg, kend, &
                                                kinc
    REAL(real_8)                             :: rs1

! ==--------------------------------------------------------------==
! TODO free EMK2

    ALLOCATE(emk2(nkpt%ngwk,nkpt%nkpnt),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (tkpts%tkpnt) THEN
       CALL inq_swap(kbeg,kend,kinc)
       DO ikpt=kbeg,kend,kinc
          IF (tkpts%tkblock) CALL rkpt_swap(rs1,1,ikpt,'HGKP HGKM')
          CALL calc_emk2(ikpt)
          IF (tkpts%tkblock) CALL wkpt_swap(emk2,1,ikpt,'EMK2')
       ENDDO
    ELSE
       !$omp parallel do private(IG)
       DO ig = 1, ncpw%ngw
          emk2(ig,1) = EXP(-fint1%betap*hg(ig)*parm%tpiba2*0.25_real_8)
       ENDDO
       IF (geq0) emk2(1,1) = 1._real_8
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE eg2
  ! ==================================================================
  SUBROUTINE calc_emk2(ikpt)
    ! ==--------------------------------------------------------------==
    ! == CALCULATE EMK2 FOR A SET OF K POINTS (SEE EG2)               ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ikpt

    INTEGER                                  :: ig, ik

    DO ik=1,nkpbl(ikpt)
       !$omp parallel do private(IG)
       DO ig = 1, ncpw%ngw
          emk2(ig,ik)     = EXP(-fint1%betap*hgkp(ig,ik)*parm%tpiba2*0.25_real_8)
          emk2(ncpw%ngw+ig,ik) = EXP(-fint1%betap*hgkm(ig,ik)*parm%tpiba2*0.25_real_8)
       ENDDO
       IF (geq0) emk2(ncpw%ngw+1,ik) = 0._real_8
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calc_emk2
  ! ==================================================================
  SUBROUTINE calc_biln(nxx,ikind)
    ! ==--------------------------------------------------------------==
    ! == COMPUTES BILN                                                ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nxx, ikind

    COMPLEX(real_8), PARAMETER               :: zzero = (0._real_8,0._real_8) 

    COMPLEX(real_8)                          :: zbetap
    INTEGER                                  :: i, ia, im, inlm, is, isub, &
                                                iv, j, jv, ki, kj, l, l2, li, &
                                                lj
    REAL(real_8)                             :: bfnl(2,nhx)

! ==--------------------------------------------------------------==

    CALL tiset(' CALC_BILN',isub)
    DO i=1,nxx
       inlm=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             IF (sgpp1%tsgp(is)) THEN
                CALL zeroing(bfnl(:,1:nlps_com%ngh(is)))!,2*nlps_com%ngh(is))
                DO iv=1,nlps_com%ngh(is)
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   IF (tkpts%tkpnt) THEN
                      DO jv=1,nlps_com%ngh(is)
                         l2=nghtol(jv,is)+1
                         lj=sgpp2%lpval(jv,is)
                         IF (l.EQ.l2.AND.li.EQ.lj) THEN
                            kj=sgpp2%lfval(jv,is)
                            bfnl(1,iv)=bfnl(1,iv)+&
                                 sgpp2%hlsg(ki,kj,l,is)*afnl(1,inlm+jv,i)
                            bfnl(2,iv)=bfnl(2,iv)+&
                                 sgpp2%hlsg(ki,kj,l,is)*afnl(2,inlm+jv,i)
                         ENDIF
                      ENDDO
                   ELSE
                      DO jv=1,nlps_com%ngh(is)
                         l2=nghtol(jv,is)+1
                         lj=sgpp2%lpval(jv,is)
                         IF (l.EQ.l2.AND.li.EQ.lj) THEN
                            kj=sgpp2%lfval(jv,is)
                            bfnl(1,iv)=bfnl(1,iv)+&
                                 sgpp2%hlsg(ki,kj,l,is)*afnl(1,inlm+jv,i)
                         ENDIF
                      ENDDO
                   ENDIF
                ENDDO
                IF (tkpts%tkpnt) THEN
                   CALL dcopy(imagp*nlps_com%ngh(is),bfnl(1,1),1,afnl(1,inlm+1,i),1)
                ELSE
                   CALL dcopy(nlps_com%ngh(is),bfnl(1,1),2,afnl(1,inlm+1,i),1)
                ENDIF
                inlm=inlm+nlps_com%ngh(is)
             ELSE
                DO iv=1,nlps_com%ngh(is)
                   inlm=inlm+1
                   afnl(1,inlm,i)=afnl(1,inlm,i)*wsg(is,iv)
                ENDDO
                IF (tkpts%tkpnt) THEN
                   inlm=inlm-nlps_com%ngh(is)
                   DO iv=1,nlps_com%ngh(is)
                      inlm=inlm+1
                      afnl(2,inlm,i)=afnl(2,inlm,i)*wsg(is,iv)
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    IF (tkpts%tkpnt) THEN
       zbetap=CMPLX(0.5_real_8*fint1%betap,0._real_8,kind=real_8)
       CALL zgemm('N','N',nlm,nxx,nlm,zbetap,alm(1,1,1,ikind),&
            nlm,afnl,nlm,zzero,biln,nlm)
    ELSE
       CALL dgemm('N','N',nlm,nxx,nlm,0.5_real_8*fint1%betap,alm(1,1,1,ikind),&
            nlm,afnl,nlm,0._real_8 ,biln,nlm)
    ENDIF
    DO i=1,nxx
       inlm=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             IF (sgpp1%tsgp(is)) THEN
                CALL zeroing(bfnl(:,1:nlps_com%ngh(is)))!,2*nlps_com%ngh(is))
                DO iv=1,nlps_com%ngh(is)
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   IF (tkpts%tkpnt) THEN
                      DO jv=1,nlps_com%ngh(is)
                         l2=nghtol(jv,is)+1
                         lj=sgpp2%lpval(jv,is)
                         IF (l.EQ.l2.AND.li.EQ.lj) THEN
                            kj=sgpp2%lfval(jv,is)
                            bfnl(1,iv)=bfnl(1,iv)+&
                                 sgpp2%hlsg(ki,kj,l,is)*biln(1,inlm+jv,i)
                            bfnl(2,iv)=bfnl(2,iv)+&
                                 sgpp2%hlsg(ki,kj,l,is)*biln(2,inlm+jv,i)
                         ENDIF
                      ENDDO
                   ELSE
                      DO jv=1,nlps_com%ngh(is)
                         l2=nghtol(jv,is)+1
                         lj=sgpp2%lpval(jv,is)
                         IF (l.EQ.l2.AND.li.EQ.lj) THEN
                            kj=sgpp2%lfval(jv,is)
                            bfnl(1,iv)=bfnl(1,iv)+&
                                 sgpp2%hlsg(ki,kj,l,is)*biln(1,inlm+jv,i)
                         ENDIF
                      ENDDO
                   ENDIF
                ENDDO
                IF (tkpts%tkpnt) THEN
                   CALL dcopy(imagp*nlps_com%ngh(is),bfnl(1,1),1,biln(1,inlm+1,i),1)
                ELSE
                   CALL dcopy(nlps_com%ngh(is),bfnl(1,1),2,biln(1,inlm+1,i),1)
                ENDIF
                inlm=inlm+nlps_com%ngh(is)
             ELSE
                DO iv=1,nlps_com%ngh(is)
                   inlm=inlm+1
                   biln(1,inlm,i)=biln(1,inlm,i)*wsg(is,iv)
                ENDDO
                IF (tkpts%tkpnt) THEN
                   inlm=inlm-nlps_com%ngh(is)
                   DO iv=1,nlps_com%ngh(is)
                      inlm=inlm+1
                      biln(2,inlm,i)=biln(2,inlm,i)*wsg(is,iv)
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       DO j=1,nlm
          DO im=1,imagp
             biln(im,j,i)=fint1%betap*(biln(im,j,i)-afnl(im,j,i))
          ENDDO
       ENDDO
    ENDDO
    CALL tihalt(' CALC_BILN',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calc_biln
  ! ==================================================================
  SUBROUTINE set_b2l()
    ! ==--------------------------------------------------------------==
    ! Variables
    ! ==--------------------------------------------------------------==
    IF (fint4%ntabtrot.LE.1) THEN
       RETURN
    ENDIF
    fint4%itrot=1
    cntr%b2limit=fint4%b2trot(fint4%itrot)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE set_b2l
  ! ==================================================================
  SUBROUTINE change_b2l(drhomax,trefine)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: drhomax
    LOGICAL                                  :: trefine

    INTEGER                                  :: i
    INTEGER, SAVE                            :: ifirst = 0

! ==--------------------------------------------------------------==

    IF (fint4%ntabtrot.LE.1) THEN
       RETURN
    ENDIF
    IF (ifirst.EQ.0) THEN
       fint4%itrot=1
       ifirst=1
    ENDIF
    cntr%b2limit=fint4%b2trot(fint4%itrot)
    IF (.NOT.trefine) THEN
       fint4%itrot=MIN(fint4%itrot+1,fint4%ntabtrot)
    ENDIF
    DO i=fint4%itrot,fint4%ntabtrot
       IF (drhomax.LT.fint4%denstrot(i)) THEN
          cntr%b2limit=fint4%b2trot(i)
          fint4%itrot=i
       ENDIF
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE change_b2l
  ! ==================================================================
  SUBROUTINE change_betap(drhomax,trefine)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: drhomax
    LOGICAL                                  :: trefine

    REAL(real_8), PARAMETER                  :: betaelmax = 1.e33_real_8 

    INTEGER                                  :: i
    INTEGER, SAVE                            :: ifirst = 0

! ==--------------------------------------------------------------==

    IF (fint5%ntabbetap.LE.1) RETURN
    IF (ifirst.EQ.0) THEN
       fint5%ibetap=1
       ifirst=1
    ENDIF
    fint1%betap=fint5%tabbetap(fint5%ibetap)
    IF (.NOT.trefine) THEN
       fint5%ibetap=MIN(fint5%ibetap+1,fint5%ntabbetap)
    ENDIF
    DO i=fint5%ibetap,fint5%ntabbetap
       IF (drhomax.LT.fint5%densbetap(i)) THEN
          fint1%betap=fint5%tabbetap(i)
          fint5%ibetap=i
       ENDIF
    ENDDO
    IF (fint1%betael.LT.betaelmax) fint1%p = fint1%betael/fint1%betap
    ! Exponentiated kinetic energy in G-space.
    CALL eg2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE change_betap
  ! ==================================================================
  SUBROUTINE evpsi(nnr1,psi,exp_vpot)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nnr1
    COMPLEX(real_8)                          :: psi(nnr1)
    REAL(real_8)                             :: exp_vpot(nnr1)

    INTEGER                                  :: ii, ir, isub

    CALL tiset('     EVPSI',isub)
    ii=MOD(nnr1,4)
    !$omp parallel private(IR) default(shared)
    !$omp do
    DO ir=1,ii
       psi(ir) = psi(ir)*exp_vpot(ir)
    ENDDO
    !$omp enddo
    !$omp do
    DO ir=1+ii,nnr1,4
       psi(ir) = psi(ir)*exp_vpot(ir)
       psi(ir+1) = psi(ir+1)*exp_vpot(ir+1)
       psi(ir+2) = psi(ir+2)*exp_vpot(ir+2)
       psi(ir+3) = psi(ir+3)*exp_vpot(ir+3)
    ENDDO
    !$omp enddo
    !$omp end parallel 
    CALL tihalt('     EVPSI',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE evpsi
  ! ==================================================================
  SUBROUTINE rcalc_biln(nxx)
    ! ==--------------------------------------------------------------==
    ! == COMPUTES BILN FOR REAL SPACE CODE WITH KPOINTS               ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nxx

    INTEGER                                  :: i, ia, inlm, is, isub, iv, j, &
                                                jv, ki, kj, l, l2, li, lj
    REAL(real_8)                             :: bfnl(nhx)

    CALL tiset('RCALC_BILN',isub)
    DO i=1,nxx
       inlm=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             IF (sgpp1%tsgp(is)) THEN
                CALL zeroing(bfnl(1:nlps_com%ngh(is)))!,nlps_com%ngh(is))
                DO iv=1,nlps_com%ngh(is)
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   DO jv=1,nlps_com%ngh(is)
                      l2=nghtol(jv,is)+1
                      lj=sgpp2%lpval(jv,is)
                      IF (l.EQ.l2.AND.li.EQ.lj) THEN
                         kj=sgpp2%lfval(jv,is)
                         bfnl(iv)=bfnl(iv)+&
                              sgpp2%hlsg(ki,kj,l,is)*afnl(1,inlm+jv,i)
                      ENDIF
                   ENDDO
                ENDDO
                CALL dcopy(nlps_com%ngh(is),bfnl(1),1,afnl(1,inlm+1,i),1)
                inlm=inlm+nlps_com%ngh(is)
             ELSE
                DO iv=1,nlps_com%ngh(is)
                   inlm=inlm+1
                   afnl(1,inlm,i)=afnl(1,inlm,i)*wsg(is,iv)
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    CALL dgemm('N','N',nlm,nxx,nlm,0.5_real_8*fint1%betap,alm,&
         nlm,afnl,nlm,0._real_8 ,biln,nlm)
    DO i=1,nxx
       inlm=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             IF (sgpp1%tsgp(is)) THEN
                CALL zeroing(bfnl(1:nlps_com%ngh(is)))!,nlps_com%ngh(is))
                DO iv=1,nlps_com%ngh(is)
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   DO jv=1,nlps_com%ngh(is)
                      l2=nghtol(jv,is)+1
                      lj=sgpp2%lpval(jv,is)
                      IF (l.EQ.l2.AND.li.EQ.lj) THEN
                         kj=sgpp2%lfval(jv,is)
                         bfnl(iv)=bfnl(iv)+&
                              sgpp2%hlsg(ki,kj,l,is)*biln(1,inlm+jv,i)
                      ENDIF
                   ENDDO
                ENDDO
                CALL dcopy(nlps_com%ngh(is),bfnl(1),1,biln(1,inlm+1,i),1)
                inlm=inlm+nlps_com%ngh(is)
             ELSE
                DO iv=1,nlps_com%ngh(is)
                   inlm=inlm+1
                   biln(1,inlm,i)=biln(1,inlm,i)*wsg(is,iv)
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       DO j=1,nlm
          biln(1,j,i)=fint1%betap*(biln(1,j,i)-afnl(1,j,i))
       ENDDO
    ENDDO
    CALL tihalt('RCALC_BILN',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rcalc_biln
  ! ==================================================================

END MODULE ehpsi_utils



SUBROUTINE ehpsi(c0,c2,exp_vpot,psi,nstate)
  ! ==--------------------------------------------------------------==
  ! ==                        COMPUTES                              ==
  ! ==  THE ACTION OF APPLYING e^{(-beta/P)H}  TO A WAVEFUNCTION C0,==
  ! ==  RETURNING THE ANSWER IN C2.                                 ==
  ! ==--------------------------------------------------------------==
  USE zeroing_utils,                   ONLY: zeroing
  USE cnst , ONLY:uimag
  USE cppt , ONLY:indzs,nzhs
  USE ehpsi_utils,                     ONLY: calc_biln,&
       evpsi
  USE error_handling,                  ONLY: stopgm
  USE fft , ONLY:ngrm
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: fwfftn,&
       invfftn
  USE fint , ONLY:afnl,biln,cnl,emk2
  USE geq0mod , ONLY:geq0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_prod,&
       mp_sum
  USE nlps , ONLY:nlm
  USE parac,                           ONLY: parai
  USE reshaper , ONLY:reshape_inplace
  USE system , ONLY:group,fpar,ncpw,parm, parap
  USE timer,                           ONLY: tihalt,&
       tiset
  IMPLICIT NONE
  REAL(real_8), TARGET                       :: exp_vpot(fpar%nnr1)
  COMPLEX(real_8)                            :: psi(maxfftn)
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: c2(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

  CHARACTER(*), PARAMETER                    :: procedureN = 'ehpsi'

  COMPLEX(real_8)                            :: fm, fp
  INTEGER                                    :: i, ia, ib, ibb, ierr, ifft, &
                                                ig, is1, is2, isub, ix, ix1, &
                                                ixp, lead, msglen, njump, &
                                                nl2, nnrx, nsta
  REAL(real_8), ALLOCATABLE, DIMENSION(:)    :: psix
  REAL(real_8), POINTER                      :: emk2_1d(:), exp_vpotx(:)

!!  use ovlap_utils, only : ovlap2
! Variables
! EMK2_1D is temporary alias for array EMK2
! for setting locally 1D view in this subroutine

  CALL reshape_inplace(emk2, (/ncpw%ngw/), emk2_1d)
  ! ==--------------------------------------------------------------==
  ! == INITIALIZE ELEC-FORCE ARRAYS                                 ==
  ! ==--------------------------------------------------------------==
  CALL tiset('     EHPSI',isub)
  ! ==--------------------------------------------------------------==
  IF (group%nogrp.EQ.1) THEN
     lead  = fpar%kr1s*parai%ngrays
     ifft  = 1
     nnrx=fpar%nnr1
  ELSE
     lead  = fpar%kr1s*ngrm
     ifft=2
     nnrx=fpar%krx*fpar%kr2s*fpar%kr3s
  ENDIF
  exp_vpotx => exp_vpot
  IF (group%nogrp.GT.1) THEN
     ! Summation of density within orbital split
     nl2=parap%nlink(group%nolist(group%nogrp))
     ix1=parap%nrxpl(parai%mepos,1)-parap%nrxpl(nl2,1)
     msglen= 8 * nnrx
     ! EXP(-betaP*0) = 1!!!
     ALLOCATE(psix(nnrx),stat=ierr)
     IF (ierr.NE.0) CALL stopgm( 'ehpsi.F90',&
          'Allocation problem' ,& 
          __LINE__,__FILE__)
     !$omp parallel do private(I)
#ifdef _vpp_
     !OCL NOVREC
#endif
     DO i=1,nnrx
        psix(i)=1._real_8
     ENDDO
     DO ixp=1,parm%nr1
        ix=ix1+ixp
        CALL dcopy(fpar%kr2s*fpar%kr3s,exp_vpot(ixp),fpar%kr1,psix(ix),fpar%krx)
     ENDDO
     ! WARNING: COMBINE with MULTIPLICATION operation
     CALL mp_prod(psix,exp_vpotx,nnrx,group%meogrp)
     DEALLOCATE(psix,stat=ierr)
     IF (ierr.NE.0) CALL stopgm( 'ehpsi.F90',&
          'DeAllocation problem' ,& 
          __LINE__,__FILE__)
  ENDIF
  ! ==-----------------------------------------------------------==
  ! == APPLY THE KINETIC ENERGY PROPAGATOR TO C0, RETURNING THE  ==
  ! == ANSWER IN PSI. ARRAY EMK2 CONTAINS e^{-beta g^2/4P}       ==
  ! ==-----------------------------------------------------------==
  !$omp parallel do private(I,IG)
  DO i=1,nstate
     DO ig = 1 , ncpw%ngw
        c2(ig,i) = emk2_1d(ig)*c0(ig,i)
     ENDDO
  ENDDO
  IF (nlm.GT.0) THEN
     CALL ovlap2(ncpw%ngw,nlm,nstate,afnl,cnl,c2(1,1),.TRUE.)
     CALL mp_sum(afnl,nstate*nlm,parai%allgrp)
     CALL calc_biln(nstate,1)
  ENDIF
  ! ==--------------------------------------------------------------==
  ! == WARNING: 
  ! == HERE WE APPLY EH PSI TWO AT A TIME, MAKING USE OF THE FACT   ==
  ! == THE WAVEFUNCTIONS ARE REAL.                                  ==
  ! ==--------------------------------------------------------------==
  njump=2*group%nogrp
  DO ia = 1 , nstate , njump
     nsta=MIN((nstate-ia+2)/2,group%nogrp)
     CALL zeroing(psi)!,maxfft*group%nogrp)
     ! ==------------------------------------------------------------==
     ! == Store the wavefunctions in array PSI                       ==
     ! == in the following way:                                      ==
     ! == Since we use the Gamma point only, the wavefunctions in    ==
     ! == real space are real. This fact is used to save time by     ==
     ! == using complex arrays where the wavefunctions of            ==
     ! == 2 different states (i and i+1) are combined, and           ==
     ! == performing complex Fourier Transforms.                     ==
     ! == Warning: If you use K-points other than Gamma, or zone     ==
     ! == boundary, this symmetry is broken and the trick is not     ==
     ! == valid anymore. First step is to build the array PSI.       ==
     ! ==------------------------------------------------------------==
     DO ib=1,nsta
        i=ia+2*(ib-1)
        ibb=(ib-1)*lead
        is1 = i
        is2 = i+1
        IF (is2.GT.nstate) THEN
           ! One state
           !$omp parallel do private(IG)
           DO ig=1,ncpw%ngw
              psi(nzhs(ig)+ibb)=c2(ig,is1)
              ! mb            ENDDO         ! useless to repeat
              ! mb            DO IG=1,NGW   ! twice the loop
              psi(indzs(ig)+ibb)=CONJG(c2(ig,is1))
           ENDDO
           IF (geq0) psi(nzhs(1)+ibb)=c2(1,is1)
        ELSE
           ! Two states
           !$omp parallel do private(IG)
           DO ig=1,ncpw%ngw
              psi(nzhs(ig)+ibb)=c2(ig,is1)+uimag*c2(ig,is2)
              ! mb            ENDDO        ! useless to repeat
              ! mb            DO IG=1,NGW  ! twice the loop
              psi(indzs(ig)+ibb)=CONJG(c2(ig,is1))+&
                   uimag*CONJG(c2(ig,is2))
           ENDDO
           IF (geq0) psi(nzhs(1)+ibb)=c2(1,is1)+uimag*CONJG(c2(1,is2))
        ENDIF
     ENDDO
     ! ==------------------------------------------------------------==
     ! == INVERSE FFT PSI TO OBTAIN PSI ON REAL SPACE MESH           ==
     ! ==------------------------------------------------------------==
     IF (ifft.EQ.1) THEN
        CALL  invfftn(psi,.TRUE.,parai%allgrp)
     ELSE
        CALL stopgm("THIS","SG_INVFFT NOT AVAILABLE ANYMORE",& 
             __LINE__,__FILE__)
     ENDIF
     ! ==------------------------------------------------------------==
     ! == APPLY THE EXPONENTIATED POTENTIAL TO PSI                   ==
     ! == APPLY e^(-b/P V) IN REAL SPACE, AND RETURN ANSWER IN PSI.  ==
     ! == ARRAY POT CONTAINS e^(-b/P V(r)).                          ==
     ! ==------------------------------------------------------------==
     CALL evpsi(nnrx,psi,exp_vpotx)
     ! ==------------------------------------------------------------==
     ! == TRANSFORM BACK TO RECIPROCAL SPACE                         ==
     ! ==------------------------------------------------------------==
     IF (ifft.EQ.1) THEN
        CALL  fwfftn(psi,.TRUE.,parai%allgrp)
     ELSE
        CALL stopgm("THIS","SG_FWFFT NOT AVAILABLE ANYMORE",& 
             __LINE__,__FILE__)
     ENDIF
     ! ==------------------------------------------------------------==
     ! == Decode the combination of wavefunctions (now multiplied by ==
     ! == the potential), back to those of states i and i+1.         ==
     ! ==------------------------------------------------------------==
     DO ib=1,nsta
        i=ia+2*(ib-1)
        ibb=(ib-1)*lead
        is1 = i
        is2 = i+1
        IF (is2.GT.nstate) THEN
           !$omp parallel do private(IG)
           DO ig = 1 , ncpw%ngw
              c2(ig,is1) = psi(nzhs(ig)+ibb)
           ENDDO
        ELSE
           !$omp parallel do private(IG,FP,FM)
           DO ig = 1 , ncpw%ngw
              fp=(psi(nzhs(ig)+ibb)+psi(indzs(ig)+ibb))*0.5_real_8
              fm=(psi(nzhs(ig)+ibb)-psi(indzs(ig)+ibb))*0.5_real_8
              c2(ig,is1) = CMPLX(REAL(fp), AIMAG(fm),kind=real_8)
              c2(ig,is2) = CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  ! ==------------------------------------------------------------==
  ! == ADD NON-LOCAL PP PART                                      ==
  ! ==------------------------------------------------------------==
  IF (nlm.GT.0) THEN
     CALL dgemm('N','N',2*ncpw%ngw,nstate,nlm,1._real_8,cnl,2*ncpw%ngw,&
          biln(1,1,1),nlm,1._real_8,c2(1,1),2*ncpw%ngw)
  ENDIF
  ! ==--------------------------------------------------------------==
  ! == APPLY THE KINETIC ENERGY PROPAGATOR TO C0, RETURNING THE     ==
  ! == ANSWER IN C2                                                 ==
  ! ==--------------------------------------------------------------==
  DO i=1,nstate
     !$omp parallel do private(IG)
     DO ig = 1 , ncpw%ngw
        c2(ig,i)  = emk2_1d(ig)*c2(ig,i)
     ENDDO
     IF (geq0) THEN
        c2(1,i) = CMPLX(REAL(c2(1,i)),0._real_8,kind=real_8)
     ENDIF
  ENDDO
  ! Reestablish initial conditions for the potential
  IF (group%nogrp.GT.1) THEN
     nl2=parap%nlink(group%nolist(group%nogrp))
     ix1=parap%nrxpl(parai%mepos,1)-parap%nrxpl(nl2,1)
     CALL dcopy(nnrx,exp_vpotx(1),1,psix(1),1)
     CALL zeroing(exp_vpotx)!,nnrx)
     DO ixp=1,parm%nr1
        ix=ix1+ixp
        CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,exp_vpot(ixp),fpar%kr1)
     ENDDO
  ENDIF
  ! ==--------------------------------------------------------------==
  CALL tihalt('     EHPSI',isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE ehpsi
! ==================================================================

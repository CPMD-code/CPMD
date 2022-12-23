MODULE rho1ofr_utils
  USE cnst,                            ONLY: uimag
  USE cppt,                            ONLY: indzs,&
                                             nzhs
  USE dotp_utils,                      ONLY: dotp
  USE ener,                            ONLY: chrg
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: ngrm
  USE fftmain_utils,                   ONLY: invfftn
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE reshaper,                        ONLY: reshape_inplace,&
                                             type_cast
  USE spin,                            ONLY: clsd,&
                                             lspin2,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             group,&
                                             ncpw,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rho1ofr
  PUBLIC :: rhosofr
  PUBLIC :: rhoabofr

CONTAINS

  ! ==================================================================
  SUBROUTINE rho1ofr(c0,c1,focc,rhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == THE FIRST ORDER RESPONSE DENSITY RHO1                        ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: focc(:)
    REAL(real_8), TARGET                     :: rhoe(:,:)
    COMPLEX(real_8)                          :: psi(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    REAL(real_8), PARAMETER :: o3 = 0.33333333333333333_real_8

    INTEGER                                  :: i, ib, ibb, id, ifft, ig, ir, &
                                                ispin1, isub, ix, ix1, ixp, &
                                                l, lead, leadx, nl2, nnrx, &
                                                nsta
    REAL(real_8)                             :: coef3, r1, r2, ral, rbe, rsp, &
                                                rsumg, rsumr, rto
    REAL(real_8), POINTER                    :: psix(:), rhoeg(:)

    CALL tiset('   RHO1OFR',isub)
    ! ==--------------------------------------------------------------==
    IF (group%nogrp.EQ.1) THEN
       lead  = fpar%kr1s*parai%ngrays
       leadx = fpar%nnr1
       ifft=1
       nnrx=fpar%nnr1
    ELSE
       lead  = fpar%kr1s*ngrm
       leadx = fpar%krx*fpar%kr2s*fpar%kr3s
       ifft=2
       nnrx=fpar%krx*fpar%kr2s*fpar%kr3s
    ENDIF
    CALL TYPE_CAST(PSI, SIZE(PSI), PSIX) ! ! PSI IS C(8), PSIX IS REAL(8)
    !CALL RESHAPE_INPLACE(RHOE, (/NNR1 * clsd%nlsd/), RHOEG) ! ! RHOEG(:) => RHOE(:,:)
    CALL RESHAPE_INPLACE(RHOE, (/SIZE(rhoe)/), RHOEG) ! ! RHOEG(:) => RHOE(:,:)
    ! Initialize
    CALL zeroing(rhoeg)!,clsd%nlsd*leadx)
    ! Loop over the electronic states
    DO id=1,nstate,group%nogrp
       nsta=MIN(nstate-id+1,group%nogrp)
       CALL zeroing(psi)!,maxfft*group%nogrp)
       DO ib=1,nsta
          i=id+(ib-1)
          ibb=(ib-1)*lead
          !CDIR NODEP
          DO ig=1,ncpw%ngw
             psi(nzhs(ig)+ibb)=c0(ig,i)+uimag*c1(ig,i)
             psi(indzs(ig)+ibb)=CONJG(c0(ig,i))+uimag*CONJG(c1(ig,i))
          ENDDO
          psi(nzhs(1)+ibb)=c0(1,i)+uimag*c1(1,i)
       ENDDO
       ! ==--------------------------------------------------------------==
       ! ==  Fourier transform the wave functions to real space.         ==
       ! ==--------------------------------------------------------------==
       IF (ifft.EQ.1) THEN
          CALL  invfftn(psi,.TRUE.,parai%allgrp)
       ELSE
          CALL stopgm("THIS","SG_INVFFT NOT AVAILABLE ANYMORE",& 
               __LINE__,__FILE__)
       ENDIF
       IF (group%nogrp.GT.1) THEN
          DO ib=1,nsta
             IF (group%nolist(ib).EQ.parai%me) THEN
                i=id+(ib-1)
             ENDIF
          ENDDO
       ENDIF
       coef3=2.0_real_8*focc(i)/parm%omega
       IF (cntl%tlsd) THEN
          ispin1=0
          IF (i.GT.spin_mod%nsup) ispin1=leadx
          DO ir=1,nnrx
             r1=REAL(psi(ir))
             r2=AIMAG(psi(ir))
             rhoeg(ir+ispin1)=rhoeg(ir+ispin1)+coef3*r1*r2
          ENDDO
       ELSEIF (lspin2%tlse) THEN
          DO ir=1,nnrx
             r1=REAL(psi(ir))
             r2=AIMAG(psi(ir))
             rhoeg(ir)=rhoeg(ir)+coef3*r1*r2
          ENDDO
          IF (i.EQ.clsd%ialpha) THEN
             DO l=1,nnrx
                r1=REAL(psi(ir))
                r2=AIMAG(psi(ir))
                rhoeg(l+2*leadx)=rhoeg(l+2*leadx)+coef3*r1*r2
             ENDDO
          ENDIF
          IF (i.EQ.clsd%ibeta) THEN
             DO l=1,nnrx
                r1=REAL(psi(ir))
                r2=AIMAG(psi(ir))
                rhoeg(l+leadx)=rhoeg(l+leadx)+coef3*r1*r2
             ENDDO
          ENDIF
       ELSE
          DO ir=1,nnrx
             r1=REAL(psi(ir))
             r2=AIMAG(psi(ir))
             rhoeg(ir)=rhoeg(ir)+coef3*r1*r2
          ENDDO
       ENDIF
    ENDDO                     ! End loop over the electronic states
    ! ==--------------------------------------------------------------==
    IF (group%nogrp.GT.1) THEN
       ! Summation of density within orbital split
       nl2=parap%nlink(group%nolist(group%nogrp))
       ix1=parap%nrxpl(parai%mepos,1)-parap%nrxpl(nl2,1)
       IF (cntl%tlsd) THEN
          CALL mp_sum(rhoeg,psix,nnrx,group%meogrp)
          CALL zeroing(rhoe(:,1))
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,rhoe(ixp,1),fpar%kr1)
          ENDDO
          CALL mp_sum(rhoeg(leadx+1:),psix,nnrx,group%meogrp)
          CALL zeroing(rhoe(:,2))
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,rhoe(ixp,2),fpar%kr1)
          ENDDO
       ELSEIF (lspin2%tlse) THEN
          CALL mp_sum(rhoeg,psix,nnrx,group%meogrp)
          CALL zeroing(rhoe(:,1))
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,rhoe(ixp,1),fpar%kr1)
          ENDDO
          CALL mp_sum(rhoeg(leadx+1:),psix,nnrx,group%meogrp)
          CALL zeroing(rhoe(:,2))
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,rhoe(ixp,2),fpar%kr1)
          ENDDO
          CALL mp_sum(rhoeg(2*leadx+1:),psix,nnrx,group%meogrp)
          CALL zeroing(rhoe(:,3))
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,rhoe(ixp,3),fpar%kr1)
          ENDDO
       ELSE
          CALL mp_sum(rhoeg,psix,nnrx,group%meogrp)
          CALL zeroing(rhoe(:,1))
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,rhoe(ixp,1),fpar%kr1)
          ENDDO
       ENDIF
    ENDIF
    ! CONTRIBUTION OF THE VANDERBILT PP TO RHOE
    IF (pslo_com%tivan) THEN
       CALL stopgm('RHO1OFR','VDB not implemented',& 
            __LINE__,__FILE__)
    ENDIF
    ! SYMMETRIZE DENSITY IF POINT GROUP SPECIFIED 
    ! (NEED FOR SPECIAL K-POINTS).
    IF (cntl%tsymrho) THEN
       CALL stopgm('RHO1OFR','SYMRHO not implemented',& 
            __LINE__,__FILE__)
    ENDIF
    ! ALPHA+BETA DENSITY IN RHOE(*,1), BETA DENSITY IN RHOE(*,2)
    chrg%csums=0._real_8
    IF (cntl%tlsd) THEN
       DO ir=1,fpar%nnr1
          rhoe(ir,1) = rhoe(ir,1) + rhoe(ir,2)
       ENDDO
       ! ALPHA+BETA DENSITY IN RHOE(*,1), BETA DENSITY IN RHOE(*,2) ; M STATE
       ! ALPHA+BETA DENSITY IN RHOE(*,3), BETA DENSITY IN RHOE(*,4) ; T STATE
    ELSEIF (lspin2%tlse) THEN
       IF (lspin2%tlsets) THEN
          DO ir=1,fpar%nnr1
             rto=rhoe(ir,1)
             rbe=rhoe(ir,2)
             ral=rhoe(ir,3)
             rsp=0.5_real_8*(rto-ral-rbe)
             rhoe(ir,2)=rsp+rbe+o3*ral
             rhoe(ir,3)=rto
             rhoe(ir,4)=rsp+rbe+2._real_8*o3*ral
          ENDDO
       ELSE
          DO ir=1,fpar%nnr1
             rto=rhoe(ir,1)
             rbe=rhoe(ir,2)
             ral=rhoe(ir,3)
             rsp=0.5_real_8*(rto-ral-rbe)
             rhoe(ir,2)=rsp+rbe
             rhoe(ir,3)=rto
             rhoe(ir,4)=rsp+ral+rbe
          ENDDO
       ENDIF
    ENDIF
    ! HERE TO CHECK THE INTEGRAL OF THE CHARGE DENSITY
    ! FIRST IN R SPACE
    rsumr=0._real_8
    DO ir=1,fpar%nnr1
       rsumr=rsumr+rhoe(ir,1)
    ENDDO
    rsumr=rsumr*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    CALL mp_sum(rsumr,parai%allgrp)
    ! NOW IN G SPACE
    rsumg=0._real_8
    DO i=1,nstate
       rsumg=rsumg+dotp(ncpw%ngw,c0(:,i),c1(:,i))
    ENDDO
    CALL mp_sum(rsumg,parai%allgrp)
    IF (ABS(rsumr).GT.1.e-8_real_8.AND.paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,1PE16.8)')&
            " WARNING! First order density integral (R space) ",rsumR
    ENDIF
    IF (ABS(rsumg).GT.1.e-8_real_8.AND.paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,1PE16.8)')&
            " WARNING! First order density integral (G space) ",rsumG
    ENDIF
    CALL tihalt('   RHO1OFR',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rho1ofr
  ! ==================================================================
  SUBROUTINE rhosofr(c0,rhoe,psi)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NORMALIZED ELECTRON DENSITY RHOE IN REAL SPACE          ==
    ! ==  FOR A SINGLE STATE C0                                       ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:)
    REAL(real_8)                             :: rhoe(:)
    COMPLEX(real_8)                          :: psi(:)

    REAL(real_8), PARAMETER                  :: delta = 1.e-8_real_8

    INTEGER                                  :: ig, ir, isub
    REAL(real_8)                             :: om, r1, rsum
    REAL(real_8), EXTERNAL                   :: dasum

    CALL tiset('   RHOSOFR',isub)
    om=1._real_8/parm%omega
    CALL zeroing(rhoe)!,nnr1)
    CALL zeroing(psi)!,maxfft)
    !CDIR NODEP
    DO ig=1,ncpw%ngw
       psi(nzhs(ig))=c0(ig)
       psi(indzs(ig))=CONJG(c0(ig))
    ENDDO
    IF (geq0) psi(nzhs(1))=c0(1)
    CALL  invfftn(psi,.TRUE.,parai%allgrp)
    DO ir=1,fpar%nnr1
       r1=REAL(psi(ir))
       rhoe(ir)=rhoe(ir)+r1*r1
    ENDDO
    CALL dscal(fpar%nnr1,om,rhoe,1)
    rsum=dasum(fpar%nnr1,rhoe,1)*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    CALL mp_sum(rsum,parai%allgrp)
    IF (paral%parent.AND.ABS(rsum-1._real_8).GT.delta) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,T45,F20.12)')&
            'ORBITAL DENSITY IN REAL SPACE:',rsum
       CALL stopgm('RHOSOFR','ORBITAL DENSITY SUM IS NOT 1',& 
            __LINE__,__FILE__)
    ENDIF
    CALL tihalt('   RHOSOFR',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rhosofr
  ! ==================================================================
  SUBROUTINE rhoabofr(nstate,ca,cb,rhoe,psi)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NORMALIZED ELECTRON DENSITY RHOE IN REAL SPACE          ==
    ! ==  FOR A PAIR OF WAVEFUNCTIONS CA AND CB                       ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: ca(:,:), cb(:,:)
    REAL(real_8)                             :: rhoe(:)
    COMPLEX(real_8)                          :: psi(:)

    INTEGER                                  :: ig, ir, is, isub
    REAL(real_8)                             :: om, r1, r2

    CALL tiset('  RHOABOFR',isub)
    om=1._real_8/parm%omega
    DO is=1,nstate
       CALL zeroing(psi)!,maxfft)
       !CDIR NODEP
       DO ig=1,ncpw%ngw
          psi(nzhs(ig))=ca(ig,is)+uimag*cb(ig,is)
          psi(indzs(ig))=CONJG(ca(ig,is))+uimag*CONJG(cb(ig,is))
       ENDDO
       IF (geq0) psi(nzhs(1))=ca(1,is)+uimag*cb(1,is)
       CALL invfftn(psi,.TRUE.,parai%allgrp)
       DO ir=1,fpar%nnr1
          r1=REAL(psi(ir))
          r2=AIMAG(psi(ir))
          rhoe(ir)=rhoe(ir)+r1*r2*om
       ENDDO
    ENDDO
    CALL tihalt('  RHOABOFR',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rhoabofr
  ! ==================================================================

END MODULE rho1ofr_utils

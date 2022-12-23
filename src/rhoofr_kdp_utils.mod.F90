MODULE rhoofr_kdp_utils
  USE augchg_utils,                    ONLY: augchg
  USE cnst,                            ONLY: uimag
  USE cppt,                            ONLY: indzs,&
                                             nzhs
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: chrg
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: ngrm
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: invfftn
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE moverho_utils,                   ONLY: give_scr_moverho,&
                                             moverho
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE reshaper,                        ONLY: reshape_inplace,&
                                             type_cast
  USE rhov_utils,                      ONLY: give_scr_rhov,&
                                             rhov
  USE ropt,                            ONLY: ropt_mod
  USE sfac,                            ONLY: fnl
  USE spin,                            ONLY: clsd,&
                                             lspin2,&
                                             spin_mod
  USE symm,                            ONLY: symmi
  USE symtrz_utils,                    ONLY: give_scr_symrho,&
                                             symrho
  USE system,                          ONLY: cntl,&
                                             dual00,&
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

  PUBLIC :: rhoofr_kdp
  PUBLIC :: give_scr_rhoofr_kdp

CONTAINS

  ! ==================================================================
  SUBROUTINE rhoofr_kdp(c0,ckdp,rhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NORMALIZED ELECTRON DENSITY RHOE IN REAL SPACE          ==
    ! ==  FOR THE K.DOT.P (OPTIMAL BASIS SET) METHOD                  ==
    ! ==  CO IS INPUT ORBITALS                                        ==
    ! ==  CKDP IS INPUT TRANSFORMED ORBITALS                          ==
    ! ==--------------------------------------------------------------==
    ! == WARNING: ALL WAVEFUNCTIONS C0 HAVE TO BE ORTHOGONAL          ==
    ! ==--------------------------------------------------------------==
    REAL(real_8), TARGET                     :: rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(maxfftn)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: ckdp(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rhoofr_kdp'
    REAL(real_8), PARAMETER                  :: delta = 1.e-6_real_8, &
                                                deltakin = 1.e-10_real_8

    COMPLEX(real_8)                          :: czr1, czr2
    INTEGER                                  :: i, ia, iat, ib, ibb, id, &
                                                ierr, ifft, ig, is, is1, &
                                                ispin1, isub, ix, ix1, ixp, &
                                                l, lead, leadx, msglen, nl2, &
                                                nnrx, nsta
    LOGICAL                                  :: tfcal
    REAL(real_8)                             :: coef3, r1, r2, rsum, rsum1, &
                                                rsum1abs, rsumv, tt
    REAL(real_8), ALLOCATABLE                :: aux(:)
    REAL(real_8), POINTER                    :: psix(:), rhoeg(:)

    CALL tiset('RHOOFR_KDP',isub)
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
    CALL RESHAPE_INPLACE(RHOE, (/fpar%nnr1 * clsd%nlsd/), RHOEG) ! ! RHOEG(:) => RHOE(:,:)
    ! Initialize
    CALL zeroing(rhoeg)!,clsd%nlsd*leadx)
    ! Accumulate the charge and kinetic energy
    rsum=0._real_8
    DO i=1,nstate
       rsum=rsum+dotp(ncpw%ngw,c0(:,i),ckdp(:,i))
    ENDDO
    ! Loop over the electronic states
    DO id=1,nstate,group%nogrp
       nsta=MIN(nstate-id+1,group%nogrp)
       tfcal=.TRUE.
       IF (tfcal) THEN
          CALL zeroing(psi)!,maxfft*group%nogrp)
#ifdef __SR8000
          !poption parallel
#endif
          DO ib=1,nsta
             i=id+(ib-1)
             ibb=(ib-1)*lead
             is1=i
#ifdef _vpp_
             !OCL NOVREC(PSI)
#endif
             !$omp parallel do private(IG,CZR1,CZR2)
             !CDIR NODEP
             DO ig=1,ncpw%ngw
                czr1=c0(ig,is1)
                czr2=ckdp(ig,is1)
                psi(nzhs(ig)+ibb)=czr1+uimag*czr2
                psi(indzs(ig)+ibb)=CONJG(czr1)+uimag*CONJG(czr2)
             ENDDO
             IF (geq0) THEN
                psi(nzhs(1)+ibb)=c0(1,is1)+uimag*ckdp(1,is1)
             ENDIF
          ENDDO
          ! ==--------------------------------------------------------------==
          ! ==  Fourier transform the wave functions to real space.         ==
          ! ==  In the array PSI was used also the fact that the wave       ==
          ! ==  functions at Gamma are real, to form a complex array (PSI)  ==
          ! ==  with the wave functions corresponding to two different      ==
          ! ==  states (i and i+1) as the real and imaginary part. This     ==
          ! ==  allows to call the FFT routine 1/2 of the times and save    ==
          ! ==  time.                                                       ==
          ! ==--------------------------------------------------------------==
          IF (ifft.EQ.1) THEN
             CALL  invfftn(psi,.TRUE.,parai%allgrp)
          ELSE
             CALL stopgm("THIS","SG_INVFFT NOT AVAILABLE ANYMORE",& 
                  __LINE__,__FILE__)
          ENDIF
          ! Compute the charge density from the wave functions
          ! in real space 
          IF (group%nogrp.GT.1) THEN
             DO ib=1,nsta
                IF (group%nolist(ib).EQ.parai%me) THEN
                   i=id+(ib-1)
                   is1=i
                ENDIF
             ENDDO
          ELSE
             is1=id
          ENDIF
          coef3=1._real_8/parm%omega
          IF (cntl%tlsd) THEN
             IF (is1.GT.spin_mod%nsup) THEN
                ispin1=leadx
             ELSE
                ispin1=0
             ENDIF
#ifdef __SR8000
             !poption parallel
#endif
#ifdef _vpp_
             !OCL NOVREC
#endif
             !$omp parallel do private(L,R1,R2,TT)
             DO l=ispin1+1,ispin1+nnrx
                r1=REAL(psi(l))
                r2=AIMAG(psi(l))
                tt=coef3*r1*r2
                rhoeg(l)=rhoeg(l)+tt
             ENDDO
          ELSEIF (lspin2%tlse) THEN
             CALL stopgm('RHOOFR_KDP',' NOT IMPLEMENTED TLSE',& 
                  __LINE__,__FILE__)
          ELSE
#ifdef __SR8000
             !poption parallel
#endif
#ifdef _vpp_
             !OCL NOVREC
#endif
             !$omp parallel do private(L,R1,R2,TT)
             DO l=1,nnrx
                r1=REAL(psi(l))
                r2=AIMAG(psi(l))
                tt=coef3*r1*r2
                rhoeg(l)=rhoeg(l)+tt
             ENDDO
          ENDIF
       ENDIF                 ! endif TFCAL
    ENDDO                     ! End loop over the electronic states
    ! ==--------------------------------------------------------------==
    IF (group%nogrp.GT.1) THEN
       ! Summation of density within orbital split
       nl2=parap%nlink(group%nolist(group%nogrp))
       ix1=parap%nrxpl(parai%mepos,1)-parap%nrxpl(nl2,1)
       msglen= 8 * nnrx
       IF (cntl%tlsd) THEN
          CALL mp_sum(rhoeg,psix,nnrx,group%meogrp)
          CALL zeroing(rhoe(:,1))!,nnr1)
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,rhoe(ixp,1),fpar%kr1)
          ENDDO
          CALL mp_sum(rhoeg(leadx+1:),psix,fpar%nnr1,group%meogrp)
          CALL zeroing(rhoe(:,2))!,nnr1)
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,rhoe(ixp,2),fpar%kr1)
          ENDDO
       ELSEIF (lspin2%tlse) THEN
          CALL stopgm('RHOOFR_KDP',' NOT IMPLEMENTED TLSE',& 
               __LINE__,__FILE__)
       ELSE
          CALL mp_sum(rhoeg,psix,fpar%nnr1,group%meogrp)
          CALL zeroing(rhoe(:,1))!,nnr1)
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,rhoe(ixp,1),fpar%kr1)
          ENDDO
       ENDIF
    ENDIF
    ! MOVE DENSITY ACCORDING TO MOVEMENT OF ATOMS
    IF (ropt_mod%modens) CALL moverho(rhoe,psi)
    ! CONTRIBUTION OF THE VANDERBILT PP TO RHOE
    IF (pslo_com%tivan) THEN
       IF (cntl%tlsd) THEN
          ! ALPHA SPIN
          CALL rhov(nstate,1,spin_mod%nsup,rsumv,psi)
          rsum=rsum+parm%omega*rsumv
          !$omp parallel do private(I)
          DO i=1,fpar%nnr1
             rhoe(i,1)=rhoe(i,1)+REAL(psi(i))
          ENDDO
          ! BETA SPIN
          CALL rhov(nstate,spin_mod%nsup+1,nstate,rsumv,psi)
          rsum=rsum+parm%omega*rsumv
          !$omp parallel do private(I)
          DO i=1,fpar%nnr1
             rhoe(i,2)=rhoe(i,2)+REAL(psi(i))
          ENDDO
       ELSE
          CALL rhov(nstate,1,nstate,rsumv,psi)
          rsum=rsum+parm%omega*rsumv
          !$omp parallel do private(I)
          DO i=1,fpar%nnr1
             rhoe(i,1)=rhoe(i,1)+REAL(psi(i))
          ENDDO
       ENDIF
       ! Vanderbilt Charges
       IF (paral%parent) THEN
          ALLOCATE(aux(ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          CALL zeroing(aux)!,ions1%nat)
          CALL augchg(fnl,crge%f,aux,ions1%nat)
          iat=0
          DO is=1,ions1%nsp
             chrg%vdbchg(is)=0.0_real_8
             DO ia=1,ions0%na(is)
                iat=iat+1
                chrg%vdbchg(is)=chrg%vdbchg(is)+aux(iat)
             ENDDO
             chrg%vdbchg(is)=chrg%vdbchg(is)/ions0%na(is)
          ENDDO
          DEALLOCATE(aux,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    ! ALPHA+BETA DENSITY IN RHOE(*,1), BETA DENSITY IN RHOE(*,2)
    chrg%csums=0._real_8
    chrg%csumsabs=0._real_8
    IF (cntl%tlsd) THEN
       rsum1=0._real_8
       rsum1abs=0._real_8
       DO i=1,fpar%nnr1
          rsum1 = rsum1 + rhoe(i,1) - rhoe(i,2)
          rsum1abs = rsum1abs + ABS(rhoe(i,1) - rhoe(i,2))
          rhoe(i,1) = rhoe(i,1) + rhoe(i,2)
       ENDDO
       chrg%csums=rsum1*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
       chrg%csumsabs=rsum1abs*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ELSEIF (lspin2%tlse) THEN
       CALL stopgm('RHOOFR_KDP',' NOT IMPLEMENTED TLSE',& 
            __LINE__,__FILE__)
    ENDIF
    ! SYMMETRIZE DENSITY IF POINT GROUP SPECIFIED
    ! (NEED FOR SPECIAL K-POINTS).
    IF (cntl%tsymrho) THEN
       IF (cntl%tlsd) THEN
          CALL symrho(rhoe(:,1),psi)
          CALL symrho(rhoe(:,2),psi)
       ELSE
          CALL symrho(rhoe(:,1),psi)
       ENDIF
    ENDIF
    ! HERE TO CHECK THE INTEGRAL OF THE CHARGE DENSITY
    ! RSUM1=DASUM(NNR1,RHOE(1,1),1)
    ! --> with VDB PP RHOE might be negative in some points
    rsum1=0._real_8
    !$omp parallel do private(I) reduction(+:RSUM1)
    DO i=1,fpar%nnr1
       rsum1=rsum1+rhoe(i,1)
    ENDDO
    rsum1=rsum1*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    chrg%csumg=rsum
    chrg%csumr=rsum1
    CALL mp_sum(chrg%csumg,parai%allgrp)
    CALL mp_sum(chrg%csumr,parai%allgrp)
    CALL mp_sum(chrg%csums,parai%allgrp)
    CALL mp_sum(chrg%csumsabs,parai%allgrp)

    IF (paral%parent.AND.ABS(chrg%csumr-chrg%csumg).GT.delta) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,T45,F20.12)') 'IN FOURIER SPACE:', chrg%csumg
       IF (paral%io_parent)&
            WRITE(6,'(A,T45,F20.12)') 'IN REAL SPACE:', chrg%csumr
       IF ((symmi%indpg.NE.0.AND.dual00%cdual.LT.4._real_8).AND.paral%io_parent)&
            WRITE(6,*) 'YOUR DUAL NUMBER ',dual00%cdual,&
            ' COULD BE TOO SMALL WITH DENSITY SYMMETRISATION'
       CALL stopgm('RHOOFR','TOTAL DENSITY SUMS ARE NOT EQUAL',& 
            __LINE__,__FILE__)
    ENDIF
    CALL tihalt('RHOOFR_KDP',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rhoofr_kdp
  ! ==================================================================
  SUBROUTINE give_scr_rhoofr_kdp(lrhoofr,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrhoofr
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lmoverho, lrhov, lsymrho

    lmoverho=0
    IF (pslo_com%tivan.AND.paral%parent) THEN ! RHOOFR and RHOOFR_C
       CALL give_scr_rhov(lrhov,tag)
       lrhoofr=MAX(ions1%nat,lrhov)! AUGCHG
    ELSE
       lrhoofr=0
    ENDIF
    CALL give_scr_symrho(lsymrho,tag)
    IF (ropt_mod%modens) CALL give_scr_moverho(lmoverho,tag)
    lrhoofr=MAX(lrhoofr,lsymrho,lmoverho)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rhoofr_kdp
  ! ==================================================================

END MODULE rhoofr_kdp_utils

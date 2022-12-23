MODULE difrho_utils
  USE cnst,                            ONLY: uimag
  USE cppt,                            ONLY: indzs,&
                                             nzhs
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: ngrm
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: invfftn
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE reshaper,                        ONLY: reshape_inplace
  USE rhov_utils,                      ONLY: give_scr_rhov,&
                                             rhov
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             group,&
                                             ncpw,&
                                             parap,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: difrho
  PUBLIC :: give_scr_difrho

CONTAINS

  ! ==================================================================
  SUBROUTINE difrho(c0,rhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==            ELECTRON DENSITY RHOE IN REAL SPACE               ==
    ! ==--------------------------------------------------------------==
    REAL(real_8), TARGET                     :: rhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8), TARGET                  :: psi(maxfftn)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw,nstate)

    COMPLEX(real_8), POINTER                 :: psire(:)
    INTEGER                                  :: i, ib, ibb, id, ierr, ifft, &
                                                ig, is1, is2, ispin1, ispin2, &
                                                isub, ix, ix1, ixp, l, lead, &
                                                leadx, msglen, nl2, nnrx, nsta
    REAL(real_8)                             :: coef3, coef4, r1, r2, rsumv, &
                                                tff
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: psix
    REAL(real_8), EXTERNAL                   :: dasum
    REAL(real_8), POINTER                    :: rhoeg(:)

    CALL tiset('    DIFRHO',isub)
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
    psire => psi
    CALL reshape_inplace(rhoe, (/fpar%nnr1 * clsd%nlsd/), rhoeg)

    ! Initialize
    CALL zeroing(rhoeg)!,clsd%nlsd*leadx)
    ! Loop over the electronic states
    DO id=1,nstate,2*group%nogrp
       nsta=MIN((nstate-id+2)/2,group%nogrp)
       CALL zeroing(psi)!,maxfft*group%nogrp)
       tff=dasum(2*nsta,crge%f(id,1),1)
       IF (tff.GT.1.e-5_real_8) THEN
          DO ib=1,nsta
             i=id+2*(ib-1)
             ibb=(ib-1)*lead
             is1=i
             is2=i+1
             IF (is2.GT.nstate) THEN
                !CDIR NODEP
                DO ig=1,ncpw%ngw
                   psi(nzhs(ig)+ibb)=c0(ig,is1)
                   psi(indzs(ig)+ibb)=CONJG(c0(ig,is1))
                ENDDO
             ELSE
                !CDIR NODEP
                DO ig=1,ncpw%ngw
                   psi(nzhs(ig)+ibb)=c0(ig,is1)+uimag*c0(ig,is2)
                   psi(indzs(ig)+ibb)=CONJG(c0(ig,is1))+&
                        uimag*CONJG(c0(ig,is2))
                ENDDO
             ENDIF
             IF (geq0.AND.(is2.GT.nstate)) THEN
                psi(nzhs(1)+ibb)=c0(1,is1)
             ELSEIF (geq0) THEN
                psi(nzhs(1)+ibb)=c0(1,is1)+uimag*c0(1,is2)
             ENDIF
          ENDDO
          ! Compute the charge density from the wave functions
          ! in real space 
          IF (group%nogrp.GT.1) THEN
             is2=0
             DO ib=1,nsta
                IF (group%nolist(ib).EQ.parai%me) THEN
                   i=id+2*(ib-1)
                   is1=i
                   is2=i+1
                ENDIF
             ENDDO
          ELSE
             is1=id
             is2=id+1
          ENDIF
          coef3=crge%f(is1,1)/parm%omega
          IF (is2.GT.nstate) THEN
             coef4=0.0_real_8
          ELSE
             coef4=crge%f(is2,1)/parm%omega
          ENDIF
          IF (is2.EQ.0) THEN
             coef3=0.0_real_8
             coef4=0.0_real_8
          ENDIF
          IF (ifft.EQ.1) THEN
             CALL  invfftn(psi,.TRUE.,parai%allgrp)
          ELSE
             CALL stopgm("THIS","SG_INVFFT NOT AVAILABLE ANYMORE",& 
                  __LINE__,__FILE__)
          ENDIF
          IF (cntl%tlsd) THEN
             ispin1=0
             ispin2=0
             IF (is1.GT.spin_mod%nsup) ispin1=leadx
             IF (is2.GT.spin_mod%nsup) ispin2=leadx
             !$omp parallel do private(L,R1,R2) shared(COEF3,COEF4)
             DO l=1,nnrx
                r1=REAL(psire(l))
                r2=AIMAG(psire(l))
                rhoeg(l+ispin1)=rhoeg(l+ispin1)+coef3*r1*r1
                rhoeg(l+ispin2)=rhoeg(l+ispin2)+coef4*r2*r2
             ENDDO
          ELSE
             !$omp parallel do private(L,R1,R2) shared(COEF3,COEF4)
             DO l=1,nnrx
                r1=REAL(psire(l))
                r2=AIMAG(psire(l))
                rhoeg(l)=rhoeg(l)+coef3*r1*r1+coef4*r2*r2
             ENDDO
          ENDIF
       ENDIF
       ! End loop over the electronic states
    ENDDO
    IF (group%nogrp.GT.1) THEN
       ! Summation of density within orbital split
       nl2=parap%nlink(group%nolist(group%nogrp))
       ix1=parap%nrxpl(parai%mepos,1)-parap%nrxpl(nl2,1)
       msglen= 8 * nnrx
       ALLOCATE(psix(nnrx),stat=ierr)
       IF (ierr.NE.0) CALL stopgm( 'difrho.F90',&
            'Allocation problem' ,& 
            __LINE__,__FILE__)
       IF (cntl%tlsd) THEN
          CALL mp_sum(rhoeg,psix,nnrx,group%meogrp)
          CALL zeroing(rhoe(:,1))!,nnr1)
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,rhoe(ixp,1),fpar%kr1)
          ENDDO
          CALL mp_sum(rhoeg(leadx+1:),psix,nnrx,group%meogrp)
          CALL zeroing(rhoe(:,2))!,nnr1)
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,rhoe(ixp,2),fpar%kr1)
          ENDDO
       ELSE
          CALL mp_sum(rhoeg,psix,nnrx,group%meogrp)
          CALL zeroing(rhoe(:,1))!,nnr1)
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,rhoe(ixp,1),fpar%kr1)
          ENDDO
       ENDIF
       DEALLOCATE(psix,stat=ierr)
       IF (ierr.NE.0) CALL stopgm( 'difrho.F90',&
            'DeAllocation problem' ,& 
            __LINE__,__FILE__)
    ENDIF
    ! CONTRIBUTION OF THE VANDERBILT PP TO RHOE
    IF (pslo_com%tivan) THEN
       IF (cntl%tlsd) THEN
          ! ALPHA SPIN
          CALL rhov(nstate,1,spin_mod%nsup,rsumv,psi)
          !$omp parallel do private(I)
          DO i=1,fpar%nnr1
             rhoe(i,1)=rhoe(i,1)+REAL(psire(i))
          ENDDO
          ! BETA SPIN
          CALL rhov(nstate,spin_mod%nsup+1,nstate,rsumv,psi)
          !$omp parallel do private(I)
          DO i=1,fpar%nnr1
             rhoe(i,2)=rhoe(i,2)+REAL(psire(i))
          ENDDO
       ELSE
          CALL rhov(nstate,1,nstate,rsumv,psi)
          !$omp parallel do private(I)
          DO i=1,fpar%nnr1
             rhoe(i,1)=rhoe(i,1)+REAL(psire(i))
          ENDDO
       ENDIF
    ENDIF
    ! ALPHA+BETA DENSITY IN RHOE(*,1), BETA DENSITY IN RHOE(*,2)
    IF (cntl%tlsd) THEN
       !$omp parallel do private(I)
       DO i=1,fpar%nnr1
          rhoe(i,1) = rhoe(i,1) + rhoe(i,2)
       ENDDO
    ENDIF
    CALL tihalt('    DIFRHO',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE difrho
  ! ==================================================================
  SUBROUTINE give_scr_difrho(ldifrho,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ldifrho
    CHARACTER(len=30)                        :: tag

! ==--------------------------------------------------------------==

    IF (pslo_com%tivan) THEN
       CALL give_scr_rhov(ldifrho,tag)
    ELSE
       ldifrho=0
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_difrho
  ! ==================================================================

END MODULE difrho_utils

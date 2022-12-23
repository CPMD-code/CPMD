MODULE rhoofr_p_utils
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
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE moverho_utils,                   ONLY: moverho
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE perturbation_p_utils,            ONLY: restrain_ngw_zero
  USE pslo,                            ONLY: pslo_com
  USE reshaper,                        ONLY: reshape_inplace,&
                                             type_cast
  USE response_pmod,                   ONLY: dmbi,&
                                             response1
  USE rhov_utils,                      ONLY: give_scr_rhov
  USE ropt,                            ONLY: ropt_mod
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE symtrz_utils,                    ONLY: give_scr_symrho,&
                                             symrho
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

  PUBLIC :: rhoofr_p
  PUBLIC :: give_scr_rhoofr_p

CONTAINS

  ! ==================================================================
  SUBROUTINE rhoofr_p(c0,c1,drhoe,psi,nstate)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  the first order electron density drhoe in real space.       ==
    ! ==  drhoe =  <psi_i1|r><r|psi_i0> + cc                          ==
    ! ==           - (<psi_i1|psi_j0>+<psi_i0|psi_j1>) *              ==
    ! ==               * <psi_j0|r><r|psi_i0>                         ==
    ! ==                                                              ==
    ! ==  the kinetic energy ekin. It is done in reciprocal space     ==
    ! ==  where the associated operators are diagonal.                ==
    ! ==--------------------------------------------------------------==
    REAL(real_8), TARGET                     :: drhoe(fpar%nnr1,clsd%nlsd)
    COMPLEX(real_8)                          :: psi(maxfftn)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)

    INTEGER                                  :: i, ib, ibb, id, ifft, ig, &
                                                is1, ispin1, isub, ix, ix1, &
                                                ixp, l, lead, leadx, msglen, &
                                                nl2, nnrx, nsta
    LOGICAL                                  :: tfcal
    REAL(real_8)                             :: coef3, r1, r2, rsum, rsum1, &
                                                rsum1abs, tt
    REAL(real_8), POINTER                    :: psix(:), rhoeg(:)

! Variables
! ==--------------------------------------------------------------==
! When calculating the magnetic response, the perturbation operator
! and thus the first order wavefunction are PURELY IMAGINARY. 
! Therefore, the first order density is zero by construction.

    IF (response1%tnmr .OR. response1%tepr) RETURN
    IF (dmbi%cutoff_restr) CALL restrain_ngw_zero(c1,dmbi%ngw_zero,ncpw%ngw,nstate)
    ! ==--------------------------------------------------------------==
    CALL tiset('    rhoofr_p',isub)
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
    CALL type_cast(psi, SIZE(psi), psix) ! ! psi is c(8), psix is real(8)
    CALL reshape_inplace(drhoe, (/fpar%nnr1 * clsd%nlsd/), rhoeg) ! ! rhoeg(:) => drhoe(:,:)
    CALL zeroing(rhoeg)!,clsd%nlsd*leadx)

    ! accumulate the charge and kinetic energy
    rsum=0._real_8
    DO i=1,nstate
       IF (crge%f(i,1).NE.0._real_8) THEN
          rsum=rsum+2._real_8*crge%f(i,1)*dotp(ncpw%ngw,c0(:,i),c1(:,i))
       ENDIF
    ENDDO

    ! loop over the electronic states
    DO id=1,nstate,group%nogrp
       nsta=MIN((nstate-id+2)/2,group%nogrp)
       tfcal=.FALSE.
       DO i=id,MIN(id+2*nsta-1,nstate)
          tfcal=tfcal.OR.(crge%f(i,1).NE.0._real_8)
       ENDDO
       IF (tfcal) THEN
          CALL zeroing(psi)!,maxfft*group%nogrp)
          DO ib=1,nsta
             i=id+2*(ib-1)
             ibb=(ib-1)*lead
             is1=i
#ifdef _vpp_
             !OCL NOVREC
#endif
             !CDIR NODEP
             !$omp parallel do private(IG)
             DO ig=1,ncpw%ngw
                psi(nzhs(ig)+ibb)=c0(ig,is1)+uimag*c1(ig,is1)
                psi(indzs(ig)+ibb)=CONJG(c0(ig,is1))+&
                     uimag*CONJG(c1(ig,is1))
             ENDDO
             IF (geq0) THEN
                psi(nzhs(1)+ibb)=c0(1,is1)+uimag*c1(1,is1)
             ENDIF
          ENDDO
          ! ==--------------------------------------------------------------==
          ! ==  Fourier transform the wave functions to real space.         ==
          ! ==  ASSUMING THAT phi_0 AND phi_1 ARE REAL!!
          ! ==--------------------------------------------------------------==
          IF (ifft.EQ.1) THEN
             CALL  invfftn(psi,.TRUE.,parai%allgrp)
          ELSE
             CALL stopgm("THIS","SG_INVFFT NOT AVAILABLE ANYMORE",& 
                  __LINE__,__FILE__)
          ENDIF
          ! compute the charge density from the wave functions
          ! in real space 
          IF (group%nogrp.GT.1) THEN
             DO ib=1,nsta
                IF (group%nolist(ib).EQ.parai%me) THEN
                   i=id+2*(ib-1)
                   is1=i
                ENDIF
             ENDDO
          ELSE
             is1=id
          ENDIF
          coef3=crge%f(is1,1)/parm%omega
          IF (cntl%tlsd) THEN
             ispin1=0
             IF (is1.GT.spin_mod%nsup) ispin1=leadx
#ifdef _vpp_
             !OCL NOVREC
#endif
             !$omp parallel do private(l,r1,r2)
             DO  l=1,nnrx
                r1=REAL(psi(l))
                r2=AIMAG(psi(l))
                rhoeg(l+ispin1)=rhoeg(l+ispin1)+coef3*r1*r2*2._real_8
             ENDDO
          ELSE
#ifdef _vpp_
             !OCL NOVREC
#endif
             !$omp parallel do private(l,r1,r2,tt)
             DO l=1,nnrx
                r1=REAL(psi(l))
                r2=AIMAG(psi(l))
                tt=coef3*r1*r2*2._real_8
                rhoeg(l)=rhoeg(l)+tt
             ENDDO
          ENDIF
       ENDIF              ! tfcal
    ENDDO                     ! end loop over the electronic states
    ! ==--------------------------------------------------------------==
    IF (group%nogrp.GT.1) THEN
       ! summation of density within orbital split
       nl2=parap%nlink(group%nolist(group%nogrp))
       ix1=parap%nrxpl(parai%mepos,1)-parap%nrxpl(nl2,1)
       msglen= 8 * nnrx
       IF (cntl%tlsd) THEN
          CALL mp_sum(rhoeg,psix,nnrx,group%meogrp)
          CALL zeroing(drhoe(:,1))
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,drhoe(ixp,1),fpar%kr1)
          ENDDO
          CALL mp_sum(rhoeg(leadx+1:),psix,nnrx,group%meogrp)
          CALL zeroing(drhoe(:,2))
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,drhoe(ixp,2),fpar%kr1)
          ENDDO
       ELSE
          CALL mp_sum(rhoeg,psix,nnrx,group%meogrp)
          CALL zeroing(drhoe(:,1))!,nnr1)
          DO ixp=1,parm%nr1
             ix=ix1+ixp
             CALL dcopy(fpar%kr2s*fpar%kr3s,psix(ix),fpar%krx,drhoe(ixp,1),fpar%kr1)
          ENDDO
       ENDIF
    ENDIF

    ! move density according to movement of atoms
    IF (ropt_mod%modens)  THEN
       CALL stopgm('rhoofr_p','MOVE DENSITY not supported.',& 
            __LINE__,__FILE__)
       CALL moverho(drhoe,psi)
    ENDIF
    ! contribution of the vanderbilt pp to rhoe
    IF (pslo_com%tivan) THEN
       CALL stopgm('rhoofr_p','raman vanderbilt not implemented',& 
            __LINE__,__FILE__)
    ENDIF
    ! symmetrize density if point group specified 
    ! (need for special k-points.
    IF (cntl%tlsd) THEN
       CALL stopgm('rhoofr_p','raman: LSD not implemented',& 
            __LINE__,__FILE__)
    ELSE
       CALL symrho(drhoe(:,1),psi)
    ENDIF
    ! alpha+beta density in rhoe(*,1), beta density in rhoe(*,2)
    chrg%csums=0._real_8
    chrg%csumsabs=0._real_8
    IF (cntl%tlsd) THEN
       rsum1=0._real_8
       rsum1abs=0._real_8
       DO i=1,fpar%nnr1
          rsum1 = rsum1 + drhoe(i,1) - drhoe(i,2)
          rsum1abs = rsum1abs + ABS(drhoe(i,1)-drhoe(i,2))
          drhoe(i,1) = drhoe(i,1) + drhoe(i,2)
       ENDDO
       chrg%csums=rsum1*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
       chrg%csumsabs=rsum1abs*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    ENDIF
    ! here to check the integral of the charge density
    rsum1=0._real_8
    !$omp parallel do private(i) reduction(+:rsum1)
    DO i=1,fpar%nnr1
       rsum1=rsum1+drhoe(i,1)
    ENDDO
    rsum1=rsum1*parm%omega/REAL(spar%nr1s*spar%nr2s*spar%nr3s,kind=real_8)
    chrg%csumg=rsum
    chrg%csumr=rsum1
    CALL mp_sum(chrg%csumg,parai%allgrp)
    CALL mp_sum(chrg%csumr,parai%allgrp)
    CALL mp_sum(chrg%csums,parai%allgrp)
    CALL mp_sum(chrg%csumsabs,parai%allgrp)


    IF (paral%parent.AND.ABS(chrg%csumr-chrg%csumg).GT.1.e-6_real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,'(a,t45,f20.12)') 'in fourier space:', chrg%csumg
       IF (paral%io_parent)&
            WRITE(6,'(a,t45,f20.12)') 'in real space:', chrg%csumr
       CALL stopgm('rhoofr_p','total density sums are not equal',& 
            __LINE__,__FILE__)
    ENDIF
    IF (paral%parent.AND.ABS(chrg%csumr).GT.1.e-6_real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,'(a,t45,f20.12)') 'in real space:', chrg%csumr
       ! call stopgm('rhoofr_p','density in real space not zero')
    ENDIF
    IF (paral%parent.AND.ABS(chrg%csumg).GT.1.e-6_real_8) THEN
       IF (paral%io_parent)&
            WRITE(6,'(a,t45,f20.12)') 'in fourier space:', chrg%csumg
       ! call stopgm('rhoofr_p','density in g space not zero')
    ENDIF

    CALL tihalt('    rhoofr_p',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rhoofr_p
  ! ==================================================================
  SUBROUTINE give_scr_rhoofr_p(lrhoofr,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lrhoofr
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: lrhov, lsymrho_f

! variables
! ==--------------------------------------------------------------==

    IF (pslo_com%tivan.AND.paral%parent) THEN ! rhoofr and rhoofr_c
       CALL give_scr_rhov(lrhov,tag)
       lrhoofr=MAX(ions1%nat,lrhov)! augchg
    ELSE
       lrhoofr=0
    ENDIF
    CALL give_scr_symrho(lsymrho_f,tag)
    lrhoofr=MAX(lrhoofr,lsymrho_f)
    IF (ropt_mod%modens) lrhoofr = MAX(lrhoofr,fpar%nnr1*2) ! moverho
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_rhoofr_p
  ! ==================================================================

END MODULE rhoofr_p_utils

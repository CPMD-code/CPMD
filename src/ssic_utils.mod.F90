MODULE ssic_utils
  USE cppt,                            ONLY: indz,&
                                             nzh,&
                                             scg
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE system,                          ONLY: cntr,&
                                             fpar,&
                                             ncpw,&
                                             parm
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ssic

CONTAINS

  ! ==================================================================
  SUBROUTINE ssic(ehsic,rhoe,v)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == SELF INTERACTION CORRECTION (SIC) TO THE HARTREE ENERGY AND  ==
    ! == RELATED POTENTIAL CONTRIBUTIONS TO THE FORCES                ==
    ! == default: asic=0.2_real_8 bsic=0.0 after PCCP 7, 1363-1367 (2005)  ==
    ! == First version: Tsukuba, 01 September 2006                    ==
    ! == Revision     : Tsukuba, 13 October   2006                    ==
    ! ==--------------------------------------------------------------==
    ! == Input:                                                       ==
    ! ==   RHOE  Density in R-space (input)                           ==
    ! ==   V     Potential in G-space (input)                         ==
    ! ==   REMINDER - SPINDEN(i,j,k)=RHOE(i,j,k,1)-2*RHOE(i,j,k,2)    ==
    ! ==              RHOE(i,j,k,1) = rho_alpha + rho_beta            ==
    ! ==              RHOE(i,j,k,2) = rho_beta                        ==
    ! == Output:                                                      ==
    ! ==   EHSIC SIC correction to the e-e Hartree term               ==
    ! ==   V     V+VSIC potential in R-space (spin-up)                ==
    ! ==         V-VSIC potential in R-space (spin-down)              ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: ehsic, rhoe(:,:)
    COMPLEX(real_8)                          :: v(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'ssic'

    COMPLEX(real_8)                          :: ehsicg, rhet, rhets, &
                                                vsict(ncpw%nhg)
    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:)                           :: sping, vsic
    INTEGER                                  :: ierr, ig, ig1, ir
    REAL(real_8)                             :: ehsicc(2)

!(nnr1,clsd%nlsd)
!(maxfft,clsd%nlsd)

    CALL setfftn(0)
    ALLOCATE(sping(maxfft),vsic(maxfft),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! == Transform SPINDEN(IR) into SPINDEN(IG)                       ==
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       sping(ir)=CMPLX(rhoe(ir,1)-2.0_real_8*rhoe(ir,2),0._real_8,kind=real_8)
    ENDDO
    CALL  fwfftn(sping,.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    IF (geq0) THEN
       ig1=2
       rhet=sping(nzh(1))
       rhets=scg(1)*rhet
       ehsicg=0.5_real_8*rhets*rhet
       vsict(1)=rhets
    ELSE
       ig1=1
       ehsicg=(0.0_real_8,0.0_real_8)
    ENDIF
    !$omp  parallel do private(IG,RHET,RHETS) &
    !$omp  reduction(+:EHSICG)
#ifdef __SR11000
    !poption parallel
    !poption tlocal(IG,RHET,RHETS)
    !poption psum(EHSICG)
#endif 
    DO ig=ig1,ncpw%nhg
       ! ==------------------------------------------------------------==
       ! == In VSICT we store temporarily the potential in G-space     ==
       ! ==------------------------------------------------------------==
       rhet=sping(nzh(ig))
       rhets=scg(ig)*rhet
       ! Hartree e-e SIC potential
       vsict(ig)=rhets
       ! Hartree e-e SIC energy
       ehsicg=ehsicg+rhets*CONJG(rhet)
    ENDDO
    ehsicc(1)=REAL(ehsicg)
    ehsicc(2)=AIMAG(ehsicg)
    CALL mp_sum(ehsicc,2,parai%allgrp)
    ehsic=-cntr%asic*ehsicc(1)*parm%omega
    ! ==------------------------------------------------------------==
    ! == Go back in R-space and add/subtract VSIC to V(i,j,k,spin)  ==
    ! ==------------------------------------------------------------==
    CALL zeroing(vsic)!,maxfft)
    !CDIR NODEP
    DO ig=1,ncpw%nhg
       vsic(indz(ig))=CONJG(vsict(ig))
       vsic(nzh(ig)) =vsict(ig)
    ENDDO
    IF (geq0) THEN
       vsic(nzh(1))=vsict(1)
    ENDIF
    CALL  invfftn(vsic,.FALSE.,parai%allgrp)
    ! mb      IF(bsic.eq.0._real_8) THEN
    !$omp parallel do private(IR)
    DO ir=1,fpar%nnr1
       v(ir,1)=v(ir,1)-cntr%asic*vsic(ir)! alpha spin
       v(ir,2)=v(ir,2)+cntr%asic*vsic(ir)! beta  spin
    ENDDO
    ! mb      ELSE
    ! mb   put here the +/-(asic*VSIC(IR)+bsic*VXC(IR)) stuff
    ! mb      ENDIF
    ! ==--------------------------------------------------------------==
    DEALLOCATE(sping,vsic,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ssic
  ! ==================================================================

END MODULE ssic_utils
